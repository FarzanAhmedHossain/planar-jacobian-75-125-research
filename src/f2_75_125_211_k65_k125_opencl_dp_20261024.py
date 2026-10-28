"""Exact compressed K65-K125 continuation of the F2 2+1+1 profile screen.

The same unbounded-partition valuation recurrence used through K64 applies at
each later band: parts 1..5 use the recorded local orders, and parts above 5
have screen order zero. The OpenCL integer DP is checked entrywise against an
independent NumPy integer DP. This is finite local-profile pruning only; it is
not a whole-branch exclusion.
"""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pyopencl as cl

from f2_75_125_k17_nonzero_first_slice import INFINITY


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_k55_k64_opencl_dp_20261024.json"
OUT = HERE / "f2_75_125_211_k65_k125_opencl_dp_20261024.json"
FIRST, LAST, MAX_N = 65, 125, 125
N_BANDS = LAST - FIRST + 1
GRID_COUNT = 51_870
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
    for (int i=1; i<=125; ++i) {
        int cost=-6;
        if (i<=5) {
            const int oi=profiles[p*5+i-1];
            if (oi>=infinity_order) continue;
            cost=oi-6;
        }
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
    for (int b=0; b<61; ++b) {
        const int n=65+b, out=p*61+b;
        main_min[out]=10+dp[n]; main_ways[out]=ways[n];
        main_small[out]=small[n];
        int best=inf, best_r=-1;
        for (int r=1; r<=n/5; ++r) {
            const int total=n-5*r;
            if (r==5) continue;
            int value;
            if (r==2) value=0;
            else {
                if (dp[total]==inf) continue;
                value=2*(5-r)+dp[total];
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


def cpu_dp(profiles: np.ndarray) -> tuple[np.ndarray, ...]:
    count = len(profiles)
    dp = np.full((count, MAX_N + 1), INF, dtype=np.int32)
    ways = np.zeros((count, MAX_N + 1), dtype=np.uint8)
    small = np.zeros((count, MAX_N + 1), dtype=np.bool_)
    dp[:, 0] = 0
    ways[:, 0] = 1
    small[:, 0] = True
    for part in range(1, MAX_N + 1):
        if part <= 5:
            orders = profiles[:, part - 1]
            allowed = orders < INFINITY
            costs = orders.astype(np.int32) - 6
        else:
            allowed = np.ones(count, dtype=np.bool_)
            costs = np.full(count, -6, dtype=np.int32)
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
        main_min[:, b] = 10 + dp[:, n]
        main_ways[:, b] = ways[:, n]
        main_small[:, b] = small[:, n]
        for r in range(1, n // 5 + 1):
            total = n - 5 * r
            if r == 5:
                continue
            if r == 2:
                values = np.zeros(count, dtype=np.int32)
                valid = np.ones(count, dtype=np.bool_)
            else:
                values = 2 * (5 - r) + dp[:, total]
                valid = dp[:, total] != INF
            improve = valid & (values < corr_min[:, b])
            corr_min[improve, b] = values[improve]
            corr_argmin[improve, b] = r
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
    if prior.get("last_complete_band") != FIRST - 1:
        raise AssertionError("source receipt is not complete through K64")
    raw = prior["survivor_profiles"]
    if len(raw) != prior.get("survivor_profiles_after_last_complete_band"):
        raise AssertionError("K64 survivor count mismatch")
    if prior.get("finite_profile_grid_count") != GRID_COUNT:
        raise AssertionError("unexpected finite profile grid")
    profiles = np.asarray(
        [[INFINITY if value == "infinity" else int(value) for value in row]
         for row in raw],
        dtype=np.int32,
    )
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
    event = program.profile_dp(
        queue, (((nprofiles + 127) // 128) * 128,), (128,), pbuf,
        np.int32(nprofiles), np.int32(INFINITY), np.int32(INF), *buffers,
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

    flags = (gm < 0) & (gw == 1) & (gs != 0) & (gc != INF) & (gm < gc)
    active = np.ones(nprofiles, dtype=np.bool_)
    old_excluded = int(prior["cumulative_profiles_excluded_K9_K64"])
    cumulative_new = 0
    counts = dict(prior["first_excluding_band_counts"])
    by_u2 = dict(prior["first_excluding_band_by_ord_U2"])
    term_counts = dict(prior["term_counts_by_band"])
    records: dict[str, dict[str, object]] = {}
    examples: dict[str, list[dict[str, object]]] = {}
    for b, band in enumerate(range(FIRST, LAST + 1)):
        entering = active.copy()
        selected = np.flatnonzero(entering & flags[:, b])
        counts[str(band)] = int(len(selected))
        cumulative_new += int(len(selected))
        histogram: dict[str, int] = {}
        band_examples = []
        for index in selected.tolist():
            u2 = "infinity" if raw[index][1] == "infinity" else str(raw[index][1])
            histogram[u2] = histogram.get(u2, 0) + 1
            if len(band_examples) < 8:
                band_examples.append({
                    "orders_U1_to_U5": raw[index],
                    "main_order": int(gm[index, b]),
                    "correction_order": int(gc[index, b]),
                    "correction_r": int(gr[index, b]),
                })
        by_u2[str(band)] = histogram
        examples[str(band)] = band_examples
        term_counts[str(band)] = {
            "main_F_5_3": partitions(band),
            "corrections_by_r": {
                str(r): (1 if r == 2 else 0 if r == 5 else partitions(band - 5 * r))
                for r in range(1, band // 5 + 1)
            },
        }
        records[str(band)] = {
            "profiles_entering_band": int(np.count_nonzero(entering)),
            "profiles_excluded": int(len(selected)),
            "profiles_remaining": int(np.count_nonzero(entering) - len(selected)),
            "main_min_order_range": [
                int(gm[entering, b].min()), int(gm[entering, b].max())
            ],
            "unique_main_minimizers": int(np.count_nonzero(gw[entering, b] == 1)),
            "unique_small_support_minimizers": int(np.count_nonzero(
                (gw[entering, b] == 1) & (gs[entering, b] != 0)
            )),
            "profiles_with_finite_correction": int(
                np.count_nonzero(gc[entering, b] != INF)
            ),
        }
        active[selected] = False

    survivors = [raw[index] for index in np.flatnonzero(active).tolist()]
    total = old_excluded + cumulative_new
    result = {
        "scope": "exact compressed K65-K125 finite local valuation screen for F2 2+1+1 profiles",
        "field": "integer dynamic programming; exact profile orders and partition minima",
        "source_receipts": {SOURCE.name: sha256(SOURCE)},
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
        "finite_profile_grid_count": GRID_COUNT,
        "profiles_entering_K65": nprofiles,
        "last_complete_band": LAST,
        "first_excluding_band_counts": counts,
        "first_excluding_band_by_ord_U2": by_u2,
        "term_counts_by_band": term_counts,
        "band_records": records,
        "new_profiles_excluded_K65_K125": cumulative_new,
        "cumulative_profiles_excluded_K9_K125": total,
        "cumulative_percent_of_finite_profile_grid": round(100 * total / GRID_COUNT, 6),
        "survivor_profiles_after_last_complete_band": len(survivors),
        "survivor_profiles": survivors,
        "examples_K65_K125": examples,
        "status": "finite local-profile pruning only; whole 2+1+1 branch exclusion remains unproved",
    }
    temporary = OUT.with_suffix(OUT.suffix + ".tmp")
    temporary.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    temporary.replace(OUT)
    print(json.dumps({
        "bands_with_new_exclusions": {
            str(band): counts[str(band)]
            for band in range(FIRST, LAST + 1)
            if counts[str(band)]
        },
        "newly_excluded": cumulative_new,
        "cumulative_excluded": total,
        "finite_grid_percent": result["cumulative_percent_of_finite_profile_grid"],
        "survivors": len(survivors),
        "gpu": device.name,
        "gpu_cpu_parity": True,
        "receipt": OUT.name,
    }, separators=(",", ":")), flush=True)


if __name__ == "__main__":
    main()
