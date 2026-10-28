"""Extend the exact U1=0 double-root valuation screen from K24 to K125.

This is a local order-profile screen for the F2 2+1+1 carrier.  It is not a
proof excluding the double-root branch or the full 2+1+1 branch.  The
min-plus recurrence follows the exact template valuation

    2 * (a - 3 * sum(m_i)) + sum(m_i * ord(U_i)),

with a=5 for the main term and a=5-r for resonance r.  For i>=6 the
unknown order is relaxed to zero, so this remains a conservative lower bound.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEPS = HERE / "f2_75_125_flint_deps"
if DEPS.is_dir():
    sys.path.insert(0, str(DEPS))

import numpy as np
import pyopencl as cl

from f2_75_125_k17_nonzero_first_slice import templates


SOURCE = HERE / "f2_75_125_k15_k24_zero_u1_double_root.json"
OUT = HERE / "f2_75_125_211_double_root_opencl_k25_k125.json"
FIRST, LAST, MAX_N = 25, 125, 125
N_BANDS = LAST - FIRST + 1
INFINITY_ORDER = 10_000
INF = 1 << 29
EXPECTED_PROFILES = 12_128
DIRECT_CHECK_LAST = 30


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
    for (int part=2; part<=125; ++part) {
        int cost;
        if (part<=5) {
            const int oi=profiles[p*4+part-2];
            if (oi>=infinity_order) continue;
            cost=oi-6;
        } else {
            // Unknown higher U_i orders are relaxed to zero.
            cost=-6;
        }
        for (int total=part; total<=125; ++total) {
            const int prev=total-part;
            if (ways[prev]==0) continue;
            const int candidate=dp[prev]+cost;
            if (candidate<dp[total]) {
                dp[total]=candidate; ways[total]=ways[prev];
                small[total]=(uchar)((part<=5) && small[prev]);
            } else if (candidate==dp[total]) {
                ways[total]=min(2,ways[total]+ways[prev]);
                small[total]=0;
            }
        }
    }
    for (int b=0; b<101; ++b) {
        const int n=25+b, out=p*101+b;
        main_min[out]=10+dp[n]; main_ways[out]=ways[n];
        main_small[out]=small[n];
        int best=inf, best_r=-1;
        for (int r=1; r<=n/5; ++r) {
            const int total=n-5*r;
            int value;
            if (r==2) {
                // alpha=1 leaves the single, untracked U_total term;
                // its exact template base order is zero.
                value=0;
            } else if (r==5) {
                // alpha=0 contributes only the constant at total=0.
                if (total!=0) continue;
                value=0;
            } else {
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


def load_profiles() -> tuple[dict, np.ndarray]:
    prior = json.loads(SOURCE.read_text(encoding="utf-8"))
    rows = prior.get("survivor_profiles")
    if not isinstance(rows, list) or len(rows) != EXPECTED_PROFILES:
        raise AssertionError("K15--K24 double-root survivor count changed")
    profiles = np.asarray(rows, dtype=np.int32)
    if profiles.shape != (EXPECTED_PROFILES, 4):
        raise AssertionError(f"unexpected profile shape: {profiles.shape}")
    return prior, profiles


def cpu_dp(profiles: np.ndarray) -> tuple[np.ndarray, ...]:
    count = len(profiles)
    dp = np.full((count, MAX_N + 1), INF, dtype=np.int32)
    ways = np.zeros((count, MAX_N + 1), dtype=np.uint8)
    small = np.zeros((count, MAX_N + 1), dtype=np.uint8)
    dp[:, 0] = 0
    ways[:, 0] = 1
    small[:, 0] = 1
    for part in range(2, MAX_N + 1):
        if part <= 5:
            orders = profiles[:, part - 2]
            allowed = orders < INFINITY_ORDER
            costs = orders - 6
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
            small[tie, total] = 0

    shape = (count, N_BANDS)
    main_min = np.empty(shape, dtype=np.int32)
    main_ways = np.empty(shape, dtype=np.int32)
    main_small = np.empty(shape, dtype=np.uint8)
    corr_min = np.full(shape, INF, dtype=np.int32)
    corr_argmin = np.full(shape, -1, dtype=np.int32)
    for band_index, n in enumerate(range(FIRST, LAST + 1)):
        main_min[:, band_index] = 10 + dp[:, n]
        main_ways[:, band_index] = ways[:, n]
        main_small[:, band_index] = small[:, n]
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
                value = 2 * (5 - resonance) + dp[:, total]
                valid = dp[:, total] != INF
            improve = valid & (value < corr_min[:, band_index])
            corr_min[improve, band_index] = value[improve]
            corr_argmin[improve, band_index] = resonance
    return main_min, main_ways, main_small, corr_min, corr_argmin


def gpu_dp(profiles: np.ndarray) -> tuple[tuple[np.ndarray, ...], str]:
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
    queue = cl.CommandQueue(context)
    program = cl.Program(context, KERNEL).build()
    shape = (len(profiles), N_BANDS)
    outputs = (
        np.empty(shape, np.int32), np.empty(shape, np.int32),
        np.empty(shape, np.uint8), np.empty(shape, np.int32),
        np.empty(shape, np.int32),
    )
    pbuf = cl.Buffer(
        context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
        hostbuf=profiles,
    )
    buffers = [
        cl.Buffer(context, cl.mem_flags.WRITE_ONLY, array.nbytes)
        for array in outputs
    ]
    global_size = ((len(profiles) + 127) // 128) * 128
    event = program.profile_dp(
        queue, (global_size,), (128,), pbuf,
        np.int32(len(profiles)), np.int32(INFINITY_ORDER), np.int32(INF),
        *buffers,
    )
    event.wait()
    for array, buffer in zip(outputs, buffers):
        cl.enqueue_copy(queue, array, buffer)
    queue.finish()
    for buffer in buffers:
        buffer.release()
    pbuf.release()
    return outputs, device.name


def direct_minimum(terms, weights: tuple[int, ...], profile: np.ndarray):
    best = None
    ways = 0
    best_small = False
    for multiplicities, coefficient, base_order in terms:
        if not coefficient:
            continue
        order = base_order
        supported = True
        finite = True
        for weight, multiplicity in zip(weights, multiplicities):
            if not multiplicity:
                continue
            if weight <= 5:
                observed_order = int(profile[weight - 2])
                if observed_order >= INFINITY_ORDER:
                    finite = False
                    break
                order += multiplicity * observed_order
            else:
                supported = False
        if not finite:
            continue
        if best is None or order < best:
            best, ways, best_small = order, 1, supported
        elif order == best:
            ways = min(2, ways + 1)
            best_small = False
    return INF if best is None else best, ways, best_small


def validate_recurrence(profiles: np.ndarray, arrays: tuple[np.ndarray, ...]) -> int:
    """Compare selected profile/band cells with direct exact templates."""
    selected = sorted(set(
        list(range(0, min(24, len(profiles)), 3))
        + list(range(257, len(profiles), 997))
        + [len(profiles) - 1]
    ))
    main_min, main_ways, main_small, corr_min, _corr_argmin = arrays
    checked = 0
    for n in range(FIRST, DIRECT_CHECK_LAST + 1):
        main_weights = tuple(range(2, n + 1))
        main_terms = templates(n, Fraction(5, 3), main_weights, 5)
        corrections = []
        for resonance in range(1, n // 5 + 1):
            total = n - 5 * resonance
            weights = tuple(range(2, total + 1))
            terms = templates(
                total, Fraction(5 - resonance, 3), weights, 5 - resonance
            )
            corrections.append((weights, terms))
        band_index = n - FIRST
        for profile_index in selected:
            profile = profiles[profile_index]
            direct_main = direct_minimum(main_terms, main_weights, profile)
            observed_main = (
                int(main_min[profile_index, band_index]),
                int(main_ways[profile_index, band_index]),
                bool(main_small[profile_index, band_index]),
            )
            if direct_main != observed_main:
                raise AssertionError(
                    f"direct main-template mismatch at profile {profile_index}, K={n}: "
                    f"{direct_main} != {observed_main}"
                )
            direct_correction = min(
                (direct_minimum(terms, weights, profile)[0]
                 for weights, terms in corrections),
                default=INF,
            )
            observed_correction = int(corr_min[profile_index, band_index])
            if direct_correction != observed_correction:
                raise AssertionError(
                    f"direct resonance-template mismatch at profile {profile_index}, "
                    f"K={n}: {direct_correction} != {observed_correction}"
                )
            checked += 1
    return checked


def main() -> None:
    started = time.monotonic()
    prior, profiles = load_profiles()
    gpu, device_name = gpu_dp(profiles)
    cpu = cpu_dp(profiles)
    labels = (
        "main minimum", "tie count", "small-support flag",
        "correction minimum", "correction argmin",
    )
    for label, gpu_array, cpu_array in zip(labels, gpu, cpu):
        if not np.array_equal(gpu_array, cpu_array):
            mismatch = np.argwhere(gpu_array != cpu_array)[0].tolist()
            raise AssertionError(f"GPU/CPU mismatch in {label} at {mismatch}")
    directly_checked = validate_recurrence(profiles, gpu)

    main_min, main_ways, main_small, corr_min, corr_argmin = gpu
    flags = (main_min < 0) & (main_ways == 1) & (main_small != 0) & (main_min < corr_min)
    active = np.ones(len(profiles), dtype=np.bool_)
    counts: dict[str, int] = {}
    records: dict[str, dict[str, int]] = {}
    examples: dict[str, list[dict[str, object]]] = {}
    new_excluded = 0
    for band_index, band in enumerate(range(FIRST, LAST + 1)):
        entering = active.copy()
        selected = np.flatnonzero(entering & flags[:, band_index])
        counts[str(band)] = int(len(selected))
        new_excluded += len(selected)
        examples[str(band)] = [
            {
                "orders_U2_to_U5": [
                    "infinity" if int(value) >= INFINITY_ORDER else int(value)
                    for value in profiles[index].tolist()
                ],
                "main_order": int(main_min[index, band_index]),
                "correction_order": int(corr_min[index, band_index]),
                "correction_resonance": int(corr_argmin[index, band_index]),
            }
            for index in selected[:8].tolist()
        ]
        records[str(band)] = {
            "profiles_entering_band": int(np.count_nonzero(entering)),
            "profiles_excluded": int(len(selected)),
            "profiles_remaining": int(np.count_nonzero(entering) - len(selected)),
        }
        active[selected] = False

    survivors = [
        profiles[index].tolist() for index in np.flatnonzero(active).tolist()
    ]
    result = {
        "scope": "K25-K125 local min-plus valuation screen on U1=0 deeper-double-root profiles of the F2 2+1+1 quartic carrier",
        "status": "local order-profile pruning only; not a branch exclusion",
        "field": "integer valuation recurrence with exact rational-template coefficients; GPU/CPU recurrence parity",
        "valuation_model": "main order 10 + min sum((ord(U_i)-6)*m_i) for i=2..5 and -6*m_i for i>=6; resonance r has base 2*(5-r), except r=2 uses the single untracked U_(K-10) term with base 0 and r=5 contributes only at K=25",
        "input_file": SOURCE.name,
        "input_sha256": sha256(SOURCE),
        "input_scope": prior.get("scope"),
        "input_profiles_after_K24": len(profiles),
        "bands": [FIRST, LAST],
        "gpu_device": device_name,
        "cpu_gpu_parity": True,
        "direct_exact_template_cells_checked": directly_checked,
        "direct_exact_template_bands": [FIRST, DIRECT_CHECK_LAST],
        "profiles_newly_excluded_K25_K125": int(new_excluded),
        "profiles_remaining_after_K125": len(survivors),
        "first_exclusions_by_band": counts,
        "band_records": records,
        "survivor_profiles": survivors,
        "examples": examples,
        "elapsed_seconds": round(time.monotonic() - started, 3),
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    concise = {key: value for key, value in result.items() if key not in {"survivor_profiles"}}
    print(json.dumps(concise, indent=2), flush=True)
    print(f"receipt={OUT} sha256={sha256(OUT)}", flush=True)


if __name__ == "__main__":
    main()
