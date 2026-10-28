"""GPU-accelerated exact-integer K52-K54 finite-profile screen.

The rational binomial coefficients are generated exactly on the CPU. The
OpenCL kernel evaluates only their integer valuation orders and tie counts;
an independent batched NumPy calculation must match every GPU profile result.
This is a finite local-profile screen, not a whole-branch exclusion.
"""

from __future__ import annotations

import hashlib
import json
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import pyopencl as cl

from f2_75_125_k17_nonzero_first_slice import INFINITY, templates
from f2_75_125_k17_k24_nonzero_first_slice import packed_terms, term_orders


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_k50_k54_nonzero_first_slice_20261023.json"
OUT = HERE / "f2_75_125_211_k52_k54_opencl_screen_20261024.json"
GRID_COUNT = 51_870
BATCH = 2048


KERNEL = r"""
__kernel void screen_profiles(
    __global const int *profiles,
    __global const int *exponents,
    __global const int *base_orders,
    __global const uchar *main_flags,
    __global const uchar *main_support,
    const int nprofiles,
    const int nterms,
    const int infinity_order,
    __global int *excluded,
    __global int *winner,
    __global int *main_min_out,
    __global int *corr_min_out)
{
    const int p = get_global_id(0);
    if (p >= nprofiles) return;

    int best_main = 2147483647;
    int best_corr = 2147483647;
    int main_count = 0;
    int main_winner = -1;
    int winner_supported = 0;

    for (int j = 0; j < nterms; ++j) {
        int order = base_orders[j];
        int finite = 1;
        for (int k = 0; k < 5; ++k) {
            const int e = exponents[j * 5 + k];
            const int profile_order = profiles[p * 5 + k];
            if (profile_order >= infinity_order && e != 0) {
                finite = 0;
                break;
            }
            order += e * profile_order;
        }
        if (!finite) continue;

        if (main_flags[j]) {
            if (order < best_main) {
                best_main = order;
                main_count = 1;
                main_winner = j;
                winner_supported = main_support[j] != 0;
            } else if (order == best_main) {
                ++main_count;
            }
        } else if (order < best_corr) {
            best_corr = order;
        }
    }

    const int yes = (best_main < 0 && main_count == 1 && winner_supported
                     && best_corr != 2147483647 && best_main < best_corr);
    excluded[p] = yes;
    winner[p] = main_winner;
    main_min_out[p] = best_main;
    corr_min_out[p] = best_corr;
}
"""


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path: Path, value: dict) -> None:
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    temp.replace(path)


def terms_for_band(n: int):
    main = packed_terms(
        templates(n, Fraction(5, 3), tuple(range(1, n + 1)), 5)
    )
    corrections = []
    for r in range(1, n // 5 + 1):
        total = n - 5 * r
        alpha = Fraction(5 - r, 3)
        corrections.append((r, packed_terms(
            templates(total, alpha, tuple(range(1, total + 1)), 5 - r)
        )))
    return main, corrections


def concatenate_terms(main, corrections):
    exponents = [main[0]]
    base = [main[1]]
    flags = [np.ones(len(main[1]), dtype=np.uint8)]
    support = [main[2].astype(np.uint8, copy=False)]
    for _r, correction in corrections:
        exponents.append(correction[0])
        base.append(correction[1])
        flags.append(np.zeros(len(correction[1]), dtype=np.uint8))
        support.append(np.zeros(len(correction[1]), dtype=np.uint8))
    exp = np.ascontiguousarray(np.concatenate(exponents, axis=0), dtype=np.int32)
    bases = np.ascontiguousarray(np.concatenate(base), dtype=np.int32)
    main_flags = np.ascontiguousarray(np.concatenate(flags), dtype=np.uint8)
    main_support = np.ascontiguousarray(np.concatenate(support), dtype=np.uint8)
    if len(exp) != len(bases) or len(bases) != len(main_flags):
        raise AssertionError("packed OpenCL term arrays have inconsistent lengths")
    bound = int(np.max(np.abs(bases), initial=0))
    bound += sum(
        int(np.max(exp[:, j], initial=0)) * INFINITY for j in range(5)
    )
    if bound >= 2**31:
        raise OverflowError(f"integer order bound {bound} exceeds int32")
    return exp, bases, main_flags, main_support, bound


def cpu_batched_check(exponents, base, main_flags, support, profiles):
    """Independent NumPy order evaluation for the same finite profile batch."""
    nprofiles = len(profiles)
    excluded = np.zeros(nprofiles, dtype=np.bool_)
    winners = np.full(nprofiles, -1, dtype=np.int32)
    main_minima = np.full(nprofiles, np.iinfo(np.int32).max, dtype=np.int32)
    corr_minima = np.full(nprofiles, np.iinfo(np.int32).max, dtype=np.int32)
    for start in range(0, nprofiles, BATCH):
        stop = min(start + BATCH, nprofiles)
        block = profiles[start:stop].astype(np.int32, copy=False)
        orders = np.broadcast_to(base[:, None], (len(base), len(block))).copy()
        invalid = np.zeros(orders.shape, dtype=np.bool_)
        for j in range(5):
            e = exponents[:, j, None]
            o = block[None, :, j]
            orders += e * o
            invalid |= (e != 0) & (o >= INFINITY)
        orders[invalid] = np.iinfo(np.int32).max

        main_values = orders.copy()
        main_values[main_flags == 0, :] = np.iinfo(np.int32).max
        main_min = main_values.min(axis=0)
        main_idx = main_values.argmin(axis=0)
        main_idx = np.where(main_min == np.iinfo(np.int32).max, -1, main_idx)
        main_ties = (main_values == main_min[None, :]).sum(axis=0)
        corr_values = orders.copy()
        corr_values[main_flags != 0, :] = np.iinfo(np.int32).max
        corr_min = corr_values.min(axis=0)
        ok = (
            (main_min < 0)
            & (main_ties == 1)
            & (corr_min != np.iinfo(np.int32).max)
            & (main_min < corr_min)
            & support[np.maximum(main_idx, 0)]
        )
        excluded[start:stop] = ok
        winners[start:stop] = main_idx
        main_minima[start:stop] = main_min
        corr_minima[start:stop] = corr_min
    return excluded, winners, main_minima, corr_minima


def opencl_screen(program, context, queue, exponents, base, main_flags,
                  main_support, profiles):
    result = np.zeros(len(profiles), dtype=np.int32)
    winners = np.full(len(profiles), -1, dtype=np.int32)
    main_minima = np.full(len(profiles), np.iinfo(np.int32).max, dtype=np.int32)
    corr_minima = np.full(len(profiles), np.iinfo(np.int32).max, dtype=np.int32)

    mf = cl.mem_flags
    exp_buf = cl.Buffer(context, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=exponents)
    base_buf = cl.Buffer(context, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=base)
    flag_buf = cl.Buffer(context, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=main_flags)
    support_buf = cl.Buffer(context, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=main_support)
    for start in range(0, len(profiles), BATCH):
        stop = min(start + BATCH, len(profiles))
        block = np.ascontiguousarray(profiles[start:stop], dtype=np.int32)
        out = np.zeros(len(block), dtype=np.int32)
        win = np.full(len(block), -1, dtype=np.int32)
        mm = np.full(len(block), np.iinfo(np.int32).max, dtype=np.int32)
        cm = np.full(len(block), np.iinfo(np.int32).max, dtype=np.int32)
        pbuf = cl.Buffer(context, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=block)
        obuf = cl.Buffer(context, mf.WRITE_ONLY, out.nbytes)
        wbuf = cl.Buffer(context, mf.WRITE_ONLY, win.nbytes)
        mmbuf = cl.Buffer(context, mf.WRITE_ONLY, mm.nbytes)
        cmbuf = cl.Buffer(context, mf.WRITE_ONLY, cm.nbytes)
        global_size = ((len(block) + 127) // 128 * 128,)
        program.screen_profiles(
            queue, global_size, (128,), pbuf, exp_buf, base_buf, flag_buf,
            support_buf, np.int32(len(block)), np.int32(len(base)),
            np.int32(INFINITY), obuf, wbuf, mmbuf, cmbuf,
        )
        cl.enqueue_copy(queue, out, obuf)
        cl.enqueue_copy(queue, win, wbuf)
        cl.enqueue_copy(queue, mm, mmbuf)
        cl.enqueue_copy(queue, cm, cmbuf).wait()
        result[start:stop] = out
        winners[start:stop] = win
        main_minima[start:stop] = mm
        corr_minima[start:stop] = cm
        pbuf.release()
        obuf.release()
        wbuf.release()
        mmbuf.release()
        cmbuf.release()
    exp_buf.release()
    base_buf.release()
    flag_buf.release()
    support_buf.release()
    return result.astype(np.bool_), winners, main_minima, corr_minima


def main() -> None:
    started = time.monotonic()
    prior = json.loads(SOURCE.read_text(encoding="utf-8"))
    if prior.get("last_complete_band") != 51:
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
    devices = [device for platform in platforms for device in platform.get_devices()
               if device.type & cl.device_type.GPU]
    if not devices:
        raise RuntimeError("no OpenCL GPU is available")
    device = next((d for d in devices if "Iris" in d.name), devices[0])
    context = cl.Context([device])
    queue = cl.CommandQueue(context)
    program = cl.Program(context, KERNEL).build()

    counts = {n: 0 for n in range(52, 55)}
    by_u2 = {n: {} for n in counts}
    examples = {n: [] for n in counts}
    term_counts = dict(prior["term_counts_by_band"])
    first_counts = dict(prior["first_excluding_band_counts"])
    first_by_u2 = dict(prior["first_excluding_band_by_ord_U2"])
    cumulative_before = int(prior["cumulative_profiles_excluded_K9_K54"])
    initial_profiles = len(profiles)
    source_hash = sha256(SOURCE)
    script_hash = sha256(Path(__file__))
    crosscheck_profiles = 0
    gpu_seconds = 0.0

    for n in range(52, 55):
        main_terms, corrections = terms_for_band(n)
        term_counts[str(n)] = {
            "main_F_5_3": len(main_terms[3]),
            "corrections_by_r": {str(r): len(packed[3]) for r, packed in corrections},
        }
        exponents, bases, flags, support, bound = concatenate_terms(main_terms, corrections)
        band_start = time.monotonic()
        gpu_result, gpu_winners, gpu_main, gpu_corr = opencl_screen(
            program, context, queue, exponents, bases, flags, support, profiles
        )
        gpu_seconds += time.monotonic() - band_start
        cpu_result, cpu_winners, cpu_main, cpu_corr = cpu_batched_check(
            exponents, bases, flags, support, profiles
        )
        if not (
            np.array_equal(gpu_result, cpu_result)
            and np.array_equal(gpu_winners, cpu_winners)
            and np.array_equal(gpu_main, cpu_main)
            and np.array_equal(gpu_corr, cpu_corr)
        ):
            raise AssertionError(f"GPU/CPU exact order mismatch at K{n}")
        crosscheck_profiles += len(profiles)

        winners_now = np.flatnonzero(gpu_result)
        survivors = []
        winner_set = set(winners_now.tolist())
        for index, raw in enumerate(raw_profiles):
            if index not in winner_set:
                survivors.append(raw)
                continue
            counts[n] += 1
            key = "infinity" if raw[1] == "infinity" else str(raw[1])
            by_u2[n][key] = by_u2[n].get(key, 0) + 1
            if len(examples[n]) < 8:
                profile = tuple(INFINITY if value == "infinity" else int(value)
                                for value in raw)
                main_order, ties, _finite, _main_support, main_raw = term_orders(
                    main_terms, profile
                )
                if main_order is None or len(ties) != 1:
                    raise AssertionError("GPU winner fails exact main-term replay")
                main_index = int(ties[0])
                multiplicities, coefficient, _base = main_raw[main_index]
                corr_orders = {}
                for r, correction in corrections:
                    order, _tied, _mask, _sup, _raw = term_orders(correction, profile)
                    corr_orders[str(r)] = order
                examples[n].append({
                    "orders_U1_to_U5": raw,
                    "main_pole_order": -int(main_order),
                    "main_term_multiplicities": list(multiplicities),
                    "main_coefficient": str(coefficient),
                    "correction_orders_by_r": corr_orders,
                })

        profiles = np.asarray([
            [INFINITY if value == "infinity" else int(value) for value in raw]
            for raw in survivors
        ], dtype=np.int32).reshape((-1, 5))
        raw_profiles = survivors
        first_counts[str(n)] = counts[n]
        first_by_u2[str(n)] = by_u2[n]
        cumulative = cumulative_before + sum(counts.values())
        result = {
            "scope": "K52-K54 exact rational local valuation screen at the distinguished double root of the F2 2+1+1 nonzero-U1 branch",
            "field": "exact Q generalized-binomial coefficients and integer order arithmetic",
            "source_receipts": {SOURCE.name: source_hash},
            "script_sha256": script_hash,
            "finite_profile_grid_count": GRID_COUNT,
            "profiles_entering_K52": initial_profiles,
            "last_complete_band": n,
            "first_excluding_band_counts": first_counts,
            "first_excluding_band_by_ord_U2": first_by_u2,
            "term_counts_by_band": term_counts,
            "new_profiles_excluded_K52_K54": sum(counts.values()),
            "cumulative_profiles_excluded_K9_K54": cumulative,
            "cumulative_percent_of_finite_profile_grid": round(100 * cumulative / GRID_COUNT, 6),
            "survivor_profiles_after_last_complete_band": len(raw_profiles),
            "survivor_profiles": raw_profiles,
            "examples_K52_K54": {str(j): examples[j] for j in examples},
            "computation": {
                "device": device.name,
                "device_global_memory_bytes": int(device.global_mem_size),
                "device_max_allocation_bytes": int(device.max_mem_alloc_size),
                "batch_profiles": BATCH,
                "int32_order_bound": bound,
                "gpu_cpu_all_profile_parity": True,
                "profile_band_evaluations_crosschecked": crosscheck_profiles,
                "gpu_kernel_seconds": round(gpu_seconds, 6),
                "elapsed_seconds": round(time.monotonic() - started, 6),
            },
            "status": "finite local-profile pruning only; whole 2+1+1 branch exclusion remains unproved",
        }
        atomic_json(OUT, result)
        print(json.dumps({
            "band": n,
            "newly_excluded": counts[n],
            "survivors": len(raw_profiles),
            "cumulative_excluded": cumulative,
            "finite_grid_percent": result["cumulative_percent_of_finite_profile_grid"],
            "gpu_seconds": result["computation"]["gpu_kernel_seconds"],
            "device": device.name,
            "receipt": OUT.name,
        }), flush=True)


if __name__ == "__main__":
    main()
