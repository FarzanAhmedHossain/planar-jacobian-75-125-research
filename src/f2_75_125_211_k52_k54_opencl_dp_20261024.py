"""Compressed exact K52-K54 profile screen for F2 2+1+1.

For a partition term with k parts, its order is
    2*base - 6*k + sum_i multiplicity_i * ord(U_i).
After the screen's convention ord(U_i)=0 for i>5, this is an unbounded
coin-change minimum with cost ord(U_i)-6 for part i. This script computes
that minimum, its capped tie count, and whether the unique minimizer uses
only U1,...,U5. The OpenCL result is checked against a separate NumPy DP.

Scope is only the finite local valuation grid. It is not a whole-branch
exclusion or an attachment proof for the degree pair.
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
SOURCE = HERE / "f2_75_125_211_k50_k54_nonzero_first_slice_20261023.json"
OUT = HERE / "f2_75_125_211_k52_k54_opencl_dp_20261024.json"
FIRST_BAND = 52
LAST_BAND = 54
MAX_N = LAST_BAND
GRID_COUNT = 51_870
INF = 1 << 29


KERNEL = r"""
__kernel void profile_dp(
    __global const int *profiles,
    const int nprofiles,
    const int infinity_order,
    const int inf,
    __global int *main_min,
    __global int *main_ways,
    __global uchar *main_small,
    __global int *corr_min,
    __global int *corr_argmin)
{
    const int p = get_global_id(0);
    if (p >= nprofiles) return;

    int dp[55];
    int ways[55];
    uchar small[55];
    for (int s = 0; s <= 54; ++s) {
        dp[s] = inf;
        ways[s] = 0;
        small[s] = 0;
    }
    dp[0] = 0;
    ways[0] = 1;
    small[0] = 1;

    // Coin order is part size. Ascending s permits repeated use of a part;
    // the outer coin loop counts unordered partitions exactly once.
    for (int i = 1; i <= 54; ++i) {
        int cost = -6;
        if (i <= 5) {
            const int oi = profiles[p * 5 + i - 1];
            if (oi >= infinity_order) continue;
            cost = oi - 6;
        }
        for (int s = i; s <= 54; ++s) {
            const int prev = s - i;
            if (ways[prev] == 0) continue;
            const int candidate = dp[prev] + cost;
            if (candidate < dp[s]) {
                dp[s] = candidate;
                ways[s] = ways[prev];
                small[s] = (uchar)((i <= 5) && small[prev]);
            } else if (candidate == dp[s]) {
                ways[s] = min(2, ways[s] + ways[prev]);
                small[s] = 0;
            }
        }
    }

    for (int b = 0; b < 3; ++b) {
        const int n = 52 + b;
        const int out = p * 3 + b;
        main_min[out] = 10 + dp[n];
        main_ways[out] = ways[n];
        main_small[out] = small[n];

        int best = inf;
        int best_r = -1;
        for (int r = 1; r <= n / 5; ++r) {
            const int total = n - 5 * r;
            if (r == 5) continue; // alpha=0; total>0, coefficient vanishes.
            int value;
            if (r == 2) {
                // alpha=1 leaves only the one-part term; total>5 here, so
                // its variable has screen order zero and its order is 0.
                value = 0;
            } else {
                if (dp[total] == inf) continue;
                value = 2 * (5 - r) + dp[total];
            }
            if (value < best) {
                best = value;
                best_r = r;
            }
        }
        corr_min[out] = best;
        corr_argmin[out] = best_r;
    }
}
"""


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def cpu_profile_dp(profiles: np.ndarray) -> tuple[np.ndarray, ...]:
    """Vectorized exact-integer unbounded-partition DP, independent of OpenCL."""
    count = len(profiles)
    dp = np.full((count, MAX_N + 1), INF, dtype=np.int32)
    ways = np.zeros((count, MAX_N + 1), dtype=np.uint8)
    small = np.zeros((count, MAX_N + 1), dtype=np.bool_)
    dp[:, 0] = 0
    ways[:, 0] = 1
    small[:, 0] = True

    for i in range(1, MAX_N + 1):
        if i <= 5:
            orders = profiles[:, i - 1]
            allowed = orders < INFINITY
            costs = orders.astype(np.int32) - 6
        else:
            allowed = np.ones(count, dtype=np.bool_)
            costs = np.full(count, -6, dtype=np.int32)
        for s in range(i, MAX_N + 1):
            prev_ways = ways[:, s - i]
            valid = allowed & (prev_ways != 0)
            candidate = dp[:, s - i] + costs
            improve = valid & (candidate < dp[:, s])
            tie = valid & (candidate == dp[:, s])
            if np.any(improve):
                dp[improve, s] = candidate[improve]
                ways[improve, s] = prev_ways[improve]
                small[improve, s] = (i <= 5) & small[improve, s - i]
            if np.any(tie):
                ways[tie, s] = np.minimum(
                    2, ways[tie, s].astype(np.int16) + prev_ways[tie].astype(np.int16)
                ).astype(np.uint8)
                small[tie, s] = False

    bands = np.arange(FIRST_BAND, LAST_BAND + 1, dtype=np.int32)
    main_min = np.empty((count, 3), dtype=np.int32)
    main_ways = np.empty((count, 3), dtype=np.int32)
    main_small = np.empty((count, 3), dtype=np.uint8)
    corr_min = np.full((count, 3), INF, dtype=np.int32)
    corr_argmin = np.full((count, 3), -1, dtype=np.int32)
    for b, n in enumerate(bands.tolist()):
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
                valid = dp[:, total] != INF
                values = 2 * (5 - r) + dp[:, total]
            improve = valid & (values < corr_min[:, b])
            corr_min[improve, b] = values[improve]
            corr_argmin[improve, b] = r
    return main_min, main_ways, main_small, corr_min, corr_argmin


def count_partitions(n: int) -> int:
    ways = [0] * (n + 1)
    ways[0] = 1
    for part in range(1, n + 1):
        for total in range(part, n + 1):
            ways[total] += ways[total - part]
    return ways[n]


def recover_small_partition(profile: np.ndarray, total: int) -> list[int]:
    """Recover the unique minimizing partition restricted to parts 1..5."""
    old = [INF] * (total + 1)
    old[0] = 0
    layers = [old]
    for part in range(1, 6):
        prev = layers[-1]
        curr = prev.copy()
        order = int(profile[part - 1])
        if order < INFINITY:
            cost = order - 6
            for amount in range(part, total + 1):
                if curr[amount - part] == INF:
                    continue
                candidate = curr[amount - part] + cost
                if candidate < curr[amount]:
                    curr[amount] = candidate
        layers.append(curr)
    if layers[5][total] == INF:
        return []
    parts: list[int] = []
    amount = total
    part = 5
    while amount > 0 and part > 0:
        current = layers[part][amount]
        without = layers[part - 1][amount]
        if current == without:
            part -= 1
        else:
            parts.append(part)
            amount -= part
    if amount:
        raise AssertionError("failed to recover supported unique minimizer")
    return sorted(parts)


def main() -> None:
    started = time.monotonic()
    prior = json.loads(SOURCE.read_text(encoding="utf-8"))
    if prior.get("last_complete_band") != FIRST_BAND - 1:
        raise AssertionError("input receipt is not complete through K51")
    raw_profiles = prior["survivor_profiles"]
    if len(raw_profiles) != prior.get("survivor_profiles_after_last_complete_band"):
        raise AssertionError("K51 survivor count does not match its list")
    if prior.get("finite_profile_grid_count") != GRID_COUNT:
        raise AssertionError("unexpected finite profile grid")
    profiles = np.asarray([
        [INFINITY if value == "infinity" else int(value) for value in raw]
        for raw in raw_profiles
    ], dtype=np.int32)

    platforms = cl.get_platforms()
    devices = [d for platform in platforms for d in platform.get_devices()
               if d.type & cl.device_type.GPU]
    if not devices:
        raise RuntimeError("no OpenCL GPU is available")
    device = next((d for d in devices if "Iris" in d.name), devices[0])
    context = cl.Context([device])
    queue = cl.CommandQueue(
        context, properties=cl.command_queue_properties.PROFILING_ENABLE
    )
    program = cl.Program(context, KERNEL).build()

    n = len(profiles)
    shape = (n, 3)
    gm = np.empty(shape, dtype=np.int32)
    gw = np.empty(shape, dtype=np.int32)
    gs = np.empty(shape, dtype=np.uint8)
    gc = np.empty(shape, dtype=np.int32)
    gr = np.empty(shape, dtype=np.int32)
    mf = cl.mem_flags
    pbuf = cl.Buffer(context, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=profiles)
    buffers = [
        cl.Buffer(context, mf.WRITE_ONLY, obj.nbytes)
        for obj in (gm, gw, gs, gc, gr)
    ]
    event = program.profile_dp(
        queue, (((n + 127) // 128) * 128,), (128,), pbuf,
        np.int32(n), np.int32(INFINITY), np.int32(INF), *buffers,
    )
    event.wait()
    for obj, buffer in zip((gm, gw, gs, gc, gr), buffers):
        cl.enqueue_copy(queue, obj, buffer)
    queue.finish()
    for buffer in buffers:
        buffer.release()
    pbuf.release()

    cpu_result = cpu_profile_dp(profiles)
    gpu_result = (gm, gw, gs, gc, gr)
    for label, gpu, cpu in zip(
        ("main minimum", "tie count", "small-support flag", "correction minimum", "correction argmin"),
        gpu_result,
        cpu_result,
    ):
        if not np.array_equal(gpu, cpu):
            mismatch = np.argwhere(gpu != cpu)[0].tolist()
            raise AssertionError(f"OpenCL/NumPy mismatch in {label} at {mismatch}")

    excluded = (gm < 0) & (gw == 1) & (gs != 0) & (gm < gc) & (gc != INF)
    # Verify the profile count and classifications are deterministic and
    # retain complete counts, not only a percentage.
    source_hash = sha256(SOURCE)
    script_hash = sha256(Path(__file__))
    counts_by_band = dict(prior["first_excluding_band_counts"])
    by_u2 = dict(prior["first_excluding_band_by_ord_U2"])
    old_excluded = int(prior["cumulative_profiles_excluded_K9_K54"])
    cumulative_new = 0
    term_counts = dict(prior["term_counts_by_band"])
    examples: dict[str, list[dict[str, object]]] = {}
    band_records: dict[str, dict[str, object]] = {}

    active = np.ones(n, dtype=np.bool_)
    for b, band in enumerate(range(FIRST_BAND, LAST_BAND + 1)):
        entering = active.copy()
        selected = np.flatnonzero(entering & excluded[:, b])
        cumulative_new += int(len(selected))
        counts_by_band[str(band)] = int(len(selected))
        band_u2: dict[str, int] = {}
        band_examples = []
        for idx in selected.tolist():
            raw = raw_profiles[idx]
            u2 = "infinity" if raw[1] == "infinity" else str(raw[1])
            band_u2[u2] = band_u2.get(u2, 0) + 1
            if len(band_examples) < 8:
                unique_partition = recover_small_partition(profiles[idx], band)
                band_examples.append({
                    "orders_U1_to_U5": raw,
                    "main_min_order": int(gm[idx, b]),
                    "main_unique_small_partition": unique_partition,
                    "strictly_lower_than_corrections": int(gc[idx, b]),
                    "lowest_correction_r": int(gr[idx, b]),
                })
        by_u2[str(band)] = band_u2
        examples[str(band)] = band_examples
        term_counts[str(band)] = {
            "main_F_5_3": count_partitions(band),
            "corrections_by_r": {
                str(r): (1 if r == 2 else 0 if r == 5 else count_partitions(band - 5 * r))
                for r in range(1, band // 5 + 1)
            },
        }
        band_records[str(band)] = {
            "profiles_entering_band": int(np.count_nonzero(entering)),
            "profiles_excluded": int(len(selected)),
            "profiles_remaining": int(np.count_nonzero(entering) - len(selected)),
            "main_min_order_range": [int(gm[entering, b].min()), int(gm[entering, b].max())],
            "unique_main_minimizers": int(np.count_nonzero(gw[entering, b] == 1)),
            "unique_small_support_minimizers": int(np.count_nonzero((gw[entering, b] == 1) & (gs[entering, b] != 0))),
            "profiles_with_finite_correction": int(np.count_nonzero(gc[entering, b] != INF)),
        }
        active[selected] = False

    sequential_survivors = [raw_profiles[i] for i in np.flatnonzero(active).tolist()]

    total_excluded = old_excluded + cumulative_new
    result = {
        "scope": "exact compressed K52-K54 finite local valuation screen for F2 2+1+1 nonzero-U1 profiles",
        "field": "integer dynamic programming; exact profile orders and exact combinatorial partition minima",
        "source_receipts": {SOURCE.name: source_hash},
        "script_sha256": script_hash,
        "gpu": {"platform": device.platform.name, "device": device.name, "driver": device.driver_version},
        "gpu_cpu_parity": True,
        "gpu_cpu_profile_results_compared": int(n * 3),
        "gpu_kernel_seconds": float(event.profile.end - event.profile.start) / 1e9,
        "wall_seconds": round(time.monotonic() - started, 3),
        "finite_profile_grid_count": GRID_COUNT,
        "profiles_entering_K52": int(n),
        "last_complete_band": LAST_BAND,
        "first_excluding_band_counts": counts_by_band,
        "first_excluding_band_by_ord_U2": by_u2,
        "term_counts_by_band": term_counts,
        "band_records": band_records,
        "new_profiles_excluded_K52_K54": int(cumulative_new),
        "cumulative_profiles_excluded_K9_K54": int(total_excluded),
        "cumulative_percent_of_finite_profile_grid": round(100 * total_excluded / GRID_COUNT, 6),
        "survivor_profiles_after_last_complete_band": int(len(sequential_survivors)),
        "survivor_profiles": sequential_survivors,
        "examples_K52_K54": examples,
        "status": "finite local-profile pruning only; whole 2+1+1 branch exclusion remains unproved",
    }
    temp = OUT.with_suffix(OUT.suffix + ".tmp")
    temp.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    temp.replace(OUT)
    print(json.dumps({
        "bands": {str(b): counts_by_band[str(b)] for b in range(FIRST_BAND, LAST_BAND + 1)},
        "new_profiles_excluded": result["new_profiles_excluded_K52_K54"],
        "cumulative_excluded": total_excluded,
        "finite_grid_percent": result["cumulative_percent_of_finite_profile_grid"],
        "survivors": len(sequential_survivors),
        "gpu": device.name,
        "gpu_cpu_parity": True,
        "receipt": OUT.name,
    }), flush=True)


if __name__ == "__main__":
    main()
