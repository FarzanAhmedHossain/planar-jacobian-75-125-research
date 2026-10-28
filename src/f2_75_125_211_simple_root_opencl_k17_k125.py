#!/usr/bin/env python3
"""Intel OpenCL K17--K125 local valuation DP at simple roots of h.

The exact K9--K16 screen leaves precisely profiles with
ord(U3)>=2 and ord(U5/t)>=1.  This program sends that surviving finite grid
through the same exact min-plus partition recurrence on the Intel Iris Xe,
then checks every output against an independent vectorized NumPy recurrence.

This is local order-profile pruning only.  A unique pole is credited only
when its unique minimizing monomial uses U2,...,U5 and lies below every
possible first-integral resonance correction.
"""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pyopencl as cl


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_k11_k16_zero_simple_root_screen_20261024.json"
OUT = HERE / "f2_75_125_211_simple_root_opencl_k17_k125_20261024.json"
FIRST, LAST, MAX_N = 17, 125, 125
N_BANDS = LAST - FIRST + 1
GRID_COUNT = 22_176
INFINITY_ORDER = 10_000
INF = 1 << 29


KERNEL = r"""
__kernel void profile_dp(
    __global const int *profiles, const int nprofiles,
    const int infinity_order, const int inf,
    __global int *main_min, __global int *main_ways,
    __global uchar *main_small, __global int *corr_min,
    __global int *corr_argmin)
{
    const int p = get_global_id(0);
    if (p >= nprofiles) return;
    int dp[126], ways[126];
    uchar small[126];
    for (int s=0; s<=125; ++s) { dp[s]=inf; ways[s]=0; small[s]=0; }
    dp[0]=0; ways[0]=1; small[0]=1;
    for (int i=2; i<=125; ++i) {
        int oi=0;
        if (i<=5) {
            oi=profiles[p*5+i-1];
            if (oi>=infinity_order) continue;
        }
        const int cost=oi-3;
        for (int s=i; s<=125; ++s) {
            const int prev=s-i;
            if (ways[prev]==0) continue;
            const int candidate=dp[prev]+cost;
            if (candidate<dp[s]) {
                dp[s]=candidate; ways[s]=ways[prev];
                small[s]=(uchar)((i<=5) && small[prev]);
            } else if (candidate==dp[s]) {
                ways[s]=min(2,ways[s]+ways[prev]); small[s]=0;
            }
        }
    }
    for (int b=0; b<109; ++b) {
        const int n=17+b, out=p*109+b;
        main_min[out]=5+dp[n]; main_ways[out]=ways[n];
        main_small[out]=small[n];
        int best=inf, best_r=-1;
        for (int r=1; r<=n/5; ++r) {
            const int total=n-5*r;
            int value;
            if (r==2) {
                // F^1 contributes the single slice F_total; total>=7 here,
                // whose conservative simple-root order is 0.
                value=0;
            } else if (r==5) {
                // F^0 has only its constant coefficient at total=0.
                if (total!=0) continue;
                value=0;
            } else {
                if (dp[total]==inf) continue;
                value=(5-r)+dp[total];
            }
            if (value<best) { best=value; best_r=r; }
        }
        corr_min[out]=best; corr_argmin[out]=best_r;
    }
}
"""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def profile_grid() -> np.ndarray:
    q_values = list(range(2, 12)) + [INFINITY_ORDER]
    r_values = list(range(2, 12)) + [INFINITY_ORDER]
    s_values = list(range(1, 12)) + [INFINITY_ORDER]
    v_values = list(range(1, 13)) + [INFINITY_ORDER]
    rows = [
        [INFINITY_ORDER, q, r, s, v]
        for q in q_values
        for r in r_values
        for s in s_values
        for v in v_values
    ]
    profiles = np.asarray(rows, dtype=np.int32)
    if profiles.shape != (18_876, 5):
        raise AssertionError(profiles.shape)
    return profiles


def cpu_dp(profiles: np.ndarray) -> tuple[np.ndarray, ...]:
    count = len(profiles)
    dp = np.full((count, MAX_N + 1), INF, dtype=np.int32)
    ways = np.zeros((count, MAX_N + 1), dtype=np.uint8)
    small = np.zeros((count, MAX_N + 1), dtype=np.bool_)
    dp[:, 0] = 0
    ways[:, 0] = 1
    small[:, 0] = True
    for part in range(2, MAX_N + 1):
        if part <= 5:
            orders = profiles[:, part - 1]
            allowed = orders < INFINITY_ORDER
            costs = orders.astype(np.int32) - 3
        else:
            allowed = np.ones(count, dtype=np.bool_)
            costs = np.full(count, -3, dtype=np.int32)
        for total in range(part, MAX_N + 1):
            previous_ways = ways[:, total - part]
            valid = allowed & (previous_ways != 0)
            candidate = dp[:, total - part] + costs
            improve = valid & (candidate < dp[:, total])
            tie = valid & (candidate == dp[:, total])
            dp[improve, total] = candidate[improve]
            ways[improve, total] = previous_ways[improve]
            small[improve, total] = (part <= 5) & small[improve, total - part]
            ways[tie, total] = np.minimum(
                2,
                ways[tie, total].astype(np.int16)
                + previous_ways[tie].astype(np.int16),
            ).astype(np.uint8)
            small[tie, total] = False

    shape = (count, N_BANDS)
    main_min = np.empty(shape, dtype=np.int32)
    main_ways = np.empty(shape, dtype=np.int32)
    main_small = np.empty(shape, dtype=np.uint8)
    corr_min = np.full(shape, INF, dtype=np.int32)
    corr_argmin = np.full(shape, -1, dtype=np.int32)
    for b, n in enumerate(range(FIRST, LAST + 1)):
        main_min[:, b] = 5 + dp[:, n]
        main_ways[:, b] = ways[:, n]
        main_small[:, b] = small[:, n]
        for resonance in range(1, n // 5 + 1):
            total = n - 5 * resonance
            if resonance == 2:
                value = np.zeros(count, dtype=np.int32)
                valid = np.ones(count, dtype=np.bool_)
            elif resonance == 5:
                if total != 0:
                    continue
                value = np.zeros(count, dtype=np.int32)
                valid = np.ones(count, dtype=np.bool_)
            else:
                value = (5 - resonance) + dp[:, total]
                valid = dp[:, total] != INF
            improve = valid & (value < corr_min[:, b])
            corr_min[improve, b] = value[improve]
            corr_argmin[improve, b] = resonance
    return main_min, main_ways, main_small, corr_min, corr_argmin


def partitions(n: int) -> int:
    ways = [0] * (n + 1)
    ways[0] = 1
    for part in range(1, n + 1):
        for total in range(part, n + 1):
            ways[total] += ways[total - part]
    return ways[n]


def main() -> None:
    started = time.monotonic()
    prior = json.loads(SOURCE.read_text(encoding="utf-8"))
    if prior.get("profiles_surviving_K9_K16") != 18_876:
        raise AssertionError("K9--K16 simple-root survivor count changed")
    if prior.get("profiles_first_excluded_K9_K10") != 3_300:
        raise AssertionError("K9--K10 simple-root exclusion count changed")
    if not prior.get("survivor_characterization_asserted_profile_by_profile"):
        raise AssertionError("K9--K16 profile characterization was not checked")
    profiles = profile_grid()
    nprofiles = len(profiles)

    devices = [
        device
        for platform in cl.get_platforms()
        for device in platform.get_devices()
        if device.type & cl.device_type.GPU
    ]
    if not devices:
        raise RuntimeError("no OpenCL GPU is available")
    device = next((item for item in devices if "Iris" in item.name), devices[0])
    context = cl.Context([device])
    queue = cl.CommandQueue(
        context, properties=cl.command_queue_properties.PROFILING_ENABLE
    )
    program = cl.Program(context, KERNEL).build()

    shape = (nprofiles, N_BANDS)
    gm = np.empty(shape, np.int32)
    gw = np.empty(shape, np.int32)
    gs = np.empty(shape, np.uint8)
    gc = np.empty(shape, np.int32)
    gr = np.empty(shape, np.int32)
    pbuf = cl.Buffer(
        context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
        hostbuf=profiles,
    )
    outputs = (gm, gw, gs, gc, gr)
    buffers = [
        cl.Buffer(context, cl.mem_flags.WRITE_ONLY, array.nbytes)
        for array in outputs
    ]
    global_size = ((nprofiles + 127) // 128) * 128
    event = program.profile_dp(
        queue, (global_size,), (128,), pbuf,
        np.int32(nprofiles), np.int32(INFINITY_ORDER), np.int32(INF), *buffers,
    )
    event.wait()
    for array, buffer in zip(outputs, buffers):
        cl.enqueue_copy(queue, array, buffer)
    queue.finish()
    for buffer in buffers:
        buffer.release()
    pbuf.release()

    cpu = cpu_dp(profiles)
    for label, gpu_array, cpu_array in zip(
        ("main minimum", "tie count", "small-support flag",
         "correction minimum", "correction argmin"), outputs, cpu,
    ):
        if not np.array_equal(gpu_array, cpu_array):
            mismatch = np.argwhere(gpu_array != cpu_array)[0].tolist()
            raise AssertionError(f"GPU/CPU mismatch in {label} at {mismatch}")

    flags = (gm < 0) & (gw == 1) & (gs != 0) & (gm < gc)
    active = np.ones(nprofiles, dtype=np.bool_)
    counts: dict[str, int] = {}
    records: dict[str, dict[str, int]] = {}
    examples: dict[str, list[dict[str, object]]] = {}
    new_excluded = 0
    for b, band in enumerate(range(FIRST, LAST + 1)):
        entering = active.copy()
        selected = np.flatnonzero(entering & flags[:, b])
        counts[str(band)] = int(len(selected))
        new_excluded += len(selected)
        band_examples = []
        for index in selected[:8].tolist():
            row = profiles[index].tolist()
            band_examples.append({
                "orders_U1_to_U5": ["infinity" if x >= INFINITY_ORDER else x for x in row],
                "main_order": int(gm[index, b]),
                "correction_order": int(gc[index, b]),
                "correction_r": int(gr[index, b]),
            })
        examples[str(band)] = band_examples
        records[str(band)] = {
            "profiles_entering_band": int(np.count_nonzero(entering)),
            "profiles_excluded": int(len(selected)),
            "profiles_remaining": int(np.count_nonzero(entering) - len(selected)),
            "profiles_with_unique_main_minimizer": int(np.count_nonzero(
                (gm[entering, b] < INF) & (gw[entering, b] == 1)
            )),
            "profiles_with_unique_minimizer_supported_on_U2_U5": int(np.count_nonzero(
                (gm[entering, b] < INF) & (gw[entering, b] == 1) & (gs[entering, b] != 0)
            )),
        }
        active[selected] = False

    survivors = [
        ["infinity" if x >= INFINITY_ORDER else int(x) for x in profiles[i].tolist()]
        for i in np.flatnonzero(active).tolist()
    ]
    prior_excluded = int(prior["profiles_first_excluded_K9_K10"])
    cumulative_excluded = prior_excluded + new_excluded
    full_grid = GRID_COUNT
    result = {
        "scope": "exact local simple-root valuation profile screen for F2 2+1+1, U1=0; not whole-branch exclusion",
        "root_multiplicity": 1,
        "field": "integer dynamic programming for exact local orders and partition minima",
        "input_receipts": {
            SOURCE.name: sha256(SOURCE),
        },
        "script_sha256": sha256(Path(__file__)),
        "gpu": {
            "platform": device.platform.name,
            "device": device.name,
            "driver": device.driver_version,
        },
        "gpu_cpu_parity": True,
        "gpu_cpu_profile_results_compared": int(nprofiles * N_BANDS),
        "gpu_kernel_seconds": float(event.profile.end - event.profile.start) / 1e9,
        "wall_seconds": round(time.monotonic() - started, 3),
        "finite_local_profile_grid_count": full_grid,
        "profiles_entering_K17": nprofiles,
        "orders_fixed_by_exact_K9_K16_screen": {
            "U1": "zero",
            "ord_U2": "2..11 or infinity",
            "ord_U3": "2..11 or infinity",
            "ord_U4": "1..11 or infinity",
            "ord_U5_over_t": "1..12 or infinity",
        },
        "band_range": [FIRST, LAST],
        "first_new_exclusions_by_band": counts,
        "new_profiles_excluded_K17_K125": int(new_excluded),
        "profiles_excluded_K9_K16": prior_excluded,
        "cumulative_profiles_excluded_K9_K125": cumulative_excluded,
        "cumulative_percent_of_finite_simple_root_profile_grid": round(
            100 * cumulative_excluded / full_grid, 6
        ),
        "profiles_surviving_K9_K125": len(survivors),
        "band_records": records,
        "survivor_profiles": survivors,
        "examples_K17_K125": examples,
        "whole_2+1+1_branch_excluded_percent": 0.0,
        "status": "finite local-profile pruning only; no whole-branch exclusion",
    }
    temporary = OUT.with_suffix(OUT.suffix + ".tmp")
    temporary.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    temporary.replace(OUT)
    print(json.dumps({k: result[k] for k in (
        "gpu", "gpu_cpu_parity", "gpu_cpu_profile_results_compared",
        "gpu_kernel_seconds", "wall_seconds", "profiles_entering_K17",
        "new_profiles_excluded_K17_K125", "cumulative_profiles_excluded_K9_K125",
        "cumulative_percent_of_finite_simple_root_profile_grid",
        "profiles_surviving_K9_K125", "whole_2+1+1_branch_excluded_percent",
    )}, indent=2))
    print(f"report={OUT}")


if __name__ == "__main__":
    main()
