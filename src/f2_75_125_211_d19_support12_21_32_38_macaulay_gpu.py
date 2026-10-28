#!/usr/bin/env python3
"""GPU Macaulay row-span screen for the d19 2+1+1, U1=0 subsystem.

All membership results are over GF(251) and are reconnaissance only.  Any
useful pure-power membership must be reconstructed and verified over QQ.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import time
from functools import lru_cache
from pathlib import Path
import sys

import numpy as np
import pyopencl as cl
import sympy as sp
from sympy.polys.domains import GF


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import f2_75_125_211_d19_support12_21_32_38_combined_groebner as builder

PRIME = 251
OUTPUT_STEM = HERE / "f2_75_125_211_d19_support12_21_32_38_macaulay_gpu_w"


class _CapturedSystem(Exception):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


@lru_cache(None)
def monomials_of_weight(weights: tuple[int, ...], target: int) -> tuple[tuple[int, ...], ...]:
    @lru_cache(None)
    def rec(index: int, remaining: int) -> tuple[tuple[int, ...], ...]:
        if index == len(weights):
            return ((),) if remaining == 0 else ()
        out = []
        for exponent in range(remaining // weights[index] + 1):
            out.extend((exponent,) + tail
                       for tail in rec(index + 1, remaining - exponent * weights[index]))
        return tuple(out)
    return rec(0, target)


def capture_compressed_system():
    captured = {}
    original = builder.ring_groebner

    def stop_after_construction(inputs, ring, **kwargs):
        captured["inputs"] = list(inputs)
        captured["ring"] = ring
        raise _CapturedSystem

    builder.QQ = GF(PRIME)
    builder.ring_groebner = stop_after_construction
    try:
        builder.main()
    except _CapturedSystem:
        pass
    finally:
        builder.ring_groebner = original
    if len(captured.get("inputs", ())) != 18:
        raise RuntimeError("failed to capture the checked 18-generator system")
    return captured["inputs"], captured["ring"]


KERNEL_SOURCE = r"""
__kernel void gather_columns(__global const uchar *matrix, __global uchar *out,
                             int nrows, int ncols, int first_col, int width) {
    int offset = (int)get_global_id(0);
    int row = (int)get_global_id(1);
    if (offset < width && row < nrows)
        out[offset * nrows + row] = matrix[row * ncols + first_col + offset];
}
__kernel void swap_and_normalize(__global uchar *matrix, int ncols,
                                 int rank, int pivot, int col,
                                 int prime, int inverse) {
    int c = (int)get_global_id(0);
    if (c >= ncols) return;
    if (pivot != rank) {
        uchar tmp = matrix[rank * ncols + c];
        matrix[rank * ncols + c] = matrix[pivot * ncols + c];
        matrix[pivot * ncols + c] = tmp;
    }
    if (c >= col) {
        int value = (int)matrix[rank * ncols + c];
        matrix[rank * ncols + c] = (uchar)((value * inverse) % prime);
    }
}
__kernel void eliminate_below(__global uchar *matrix, int nrows, int ncols,
                              int rank, int col, int prime) {
    int col_offset = (int)get_global_id(0);
    int row_offset = (int)get_global_id(1);
    int width = ncols - col - 1;
    int active_rows = nrows - rank - 1;
    if (col_offset >= width || row_offset >= active_rows) return;
    int row = rank + 1 + row_offset;
    int c = col + 1 + col_offset;
    int factor = (int)matrix[row * ncols + col];
    int value = (int)matrix[row * ncols + c]
              - factor * (int)matrix[rank * ncols + c];
    value %= prime;
    if (value < 0) value += prime;
    matrix[row * ncols + c] = (uchar)value;
}
__kernel void zero_pivot_column(__global uchar *matrix, int ncols,
                                int nrows, int rank, int col) {
    int offset = (int)get_global_id(0);
    int row = rank + 1 + offset;
    if (row < nrows) matrix[row * ncols + col] = (uchar)0;
}
__kernel void reduce_targets(__global const uchar *matrix,
                             __global uchar *targets, int ncols,
                             int rank, int col, int n_targets, int prime) {
    int col_offset = (int)get_global_id(0);
    int target = (int)get_global_id(1);
    int c = col + col_offset;
    if (target >= n_targets || c >= ncols) return;
    int base = target * ncols;
    int factor = (int)targets[base + col];
    int value = (int)targets[base + c]
              - factor * (int)matrix[rank * ncols + c];
    value %= prime;
    if (value < 0) value += prime;
    targets[base + c] = (uchar)value;
}
"""


def run(target_weight: int) -> dict:
    inputs, ring = capture_compressed_system()
    analysis_path = HERE / "f2_75_125_211_d19_support38_analysis_20261025.json"
    analysis38 = json.loads(analysis_path.read_text(encoding="utf-8"))
    names = tuple(analysis38["selected_variables"])
    weights = tuple(int(value) for value in analysis38["primitive_positive_integer_weights"])
    if len(names) != 28 or len(weights) != 28 or len(ring.gens) != len(names):
        raise AssertionError("captured ring/weight metadata mismatch")

    generator_terms = []
    generator_degrees = []
    for poly in inputs:
        terms = [(tuple(monomial), int(coefficient))
                 for monomial, coefficient in poly.terms()]
        term_degrees = {sum(w * e for w, e in zip(weights, monomial))
                        for monomial, _ in terms}
        if len(term_degrees) != 1:
            raise AssertionError("captured generator is not weight homogeneous")
        generator_terms.append(terms)
        generator_degrees.append(term_degrees.pop())

    target_monomials = monomials_of_weight(weights, target_weight)
    if not target_monomials:
        raise ValueError(f"no monomials of target weight {target_weight}")
    target_index = {monomial: index for index, monomial in enumerate(target_monomials)}
    multipliers = [monomials_of_weight(weights, target_weight - degree)
                   if degree <= target_weight else ()
                   for degree in generator_degrees]
    nrows = sum(map(len, multipliers))
    ncols = len(target_monomials)
    if not nrows:
        raise ValueError("target weight is below all generator degrees")

    target_specs = []
    target_vectors = np.zeros((len(names), ncols), dtype=np.uint8)
    for variable_index, (name, weight) in enumerate(zip(names, weights)):
        if target_weight % weight == 0:
            exponent = [0] * len(names)
            exponent[variable_index] = target_weight // weight
            target_specs.append((f"{name}^{target_weight // weight}", tuple(exponent)))
            target_vectors[len(target_specs) - 1, target_index[tuple(exponent)]] = 1
    target_vectors = np.ascontiguousarray(target_vectors[:len(target_specs)])

    matrix = np.zeros((nrows, ncols), dtype=np.uint8)
    row = 0
    for terms, multiplier_list in zip(generator_terms, multipliers):
        for multiplier in multiplier_list:
            for exponent, coefficient in terms:
                shifted = tuple(a + b for a, b in zip(exponent, multiplier))
                matrix[row, target_index[shifted]] = coefficient
            row += 1
    if row != nrows:
        raise AssertionError("Macaulay row census mismatch")
    matrix_hash = hashlib.sha256(matrix.tobytes(order="C")).hexdigest()
    matrix_bytes = int(matrix.nbytes)
    print(json.dumps({"stage": "matrix-built", "target_weight": target_weight,
                      "shape": [nrows, ncols], "bytes": matrix_bytes,
                      "matrix_sha256": matrix_hash}), flush=True)

    platforms = cl.get_platforms()
    candidates = [(platform, device) for platform in platforms
                  for device in platform.get_devices()
                  if device.type & cl.device_type.GPU]
    if not candidates:
        raise RuntimeError("no OpenCL GPU device is available")
    platform, device = next(((p, d) for p, d in candidates if "Intel" in d.name), candidates[0])
    if matrix_bytes > int(device.max_mem_alloc_size) or matrix_bytes > int(device.global_mem_size):
        raise MemoryError("Macaulay matrix exceeds the selected GPU's allocation limits")

    context = cl.Context([device])
    queue = cl.CommandQueue(context, device)
    program = cl.Program(context, KERNEL_SOURCE).build()
    gather_columns = cl.Kernel(program, "gather_columns")
    swap_and_normalize = cl.Kernel(program, "swap_and_normalize")
    eliminate_below = cl.Kernel(program, "eliminate_below")
    zero_pivot_column = cl.Kernel(program, "zero_pivot_column")
    reduce_targets = cl.Kernel(program, "reduce_targets")
    flags = cl.mem_flags
    # Iris Xe is an integrated GPU.  USE_HOST_PTR lets its OpenCL runtime use
    # the contiguous host allocation directly when supported, avoiding a
    # second full-size copy of the dense Macaulay matrix.
    matrix_buffer = cl.Buffer(context, flags.READ_WRITE | flags.USE_HOST_PTR,
                              hostbuf=matrix)
    target_buffer = cl.Buffer(context, flags.READ_WRITE | flags.COPY_HOST_PTR,
                              hostbuf=target_vectors)
    # Keep `matrix` alive for the entire lifetime of the USE_HOST_PTR buffer.
    del target_vectors
    gc.collect()

    batch_size = min(32, ncols)
    pivot_batch = np.empty((batch_size, nrows), dtype=np.uint8)
    pivot_buffer = cl.Buffer(context, flags.WRITE_ONLY, size=batch_size * nrows)
    pivot_columns = []
    rank = 0
    scan_col = 0
    started = time.monotonic()
    while rank < nrows and scan_col < ncols:
        width = min(batch_size, ncols - scan_col)
        gather_columns(queue, (width, nrows), None, matrix_buffer,
                       pivot_buffer,
                       np.int32(nrows), np.int32(ncols), np.int32(scan_col), np.int32(width))
        cl.enqueue_copy(queue, pivot_batch[:width], pivot_buffer).wait()
        active_columns = np.flatnonzero(np.any(pivot_batch[:width, rank:] != 0, axis=1))
        if not len(active_columns):
            scan_col += width
            continue
        offset = int(active_columns[0])
        col = scan_col + offset
        nonzero = np.flatnonzero(pivot_batch[offset, rank:] != 0)
        pivot = rank + int(nonzero[0])
        pivot_value = int(pivot_batch[offset, pivot])
        pivot_columns.append(col)
        scan_col = col + 1
        inverse = pow(pivot_value, PRIME - 2, PRIME)
        swap_and_normalize(queue, (ncols,), None, matrix_buffer,
                           np.int32(ncols), np.int32(rank), np.int32(pivot),
                           np.int32(col), np.int32(PRIME), np.int32(inverse))
        active_rows = nrows - rank - 1
        width_after = ncols - col - 1
        if active_rows and width_after:
            eliminate_below(queue, (width_after, active_rows), None,
                             matrix_buffer, np.int32(nrows), np.int32(ncols),
                             np.int32(rank), np.int32(col), np.int32(PRIME))
        if active_rows:
            zero_pivot_column(queue, (active_rows,), None, matrix_buffer,
                              np.int32(ncols), np.int32(nrows),
                              np.int32(rank), np.int32(col))
        target_width = ncols - col
        if target_width and target_specs:
            reduce_targets(queue, (target_width, len(target_specs)), None,
                            matrix_buffer, target_buffer, np.int32(ncols),
                            np.int32(rank), np.int32(col),
                            np.int32(len(target_specs)), np.int32(PRIME))
        rank += 1
        if rank % 100 == 0:
            print(f"rank_progress={rank} columns_scanned={col + 1}", flush=True)

    residuals = np.empty((len(target_specs), ncols), dtype=np.uint8)
    cl.enqueue_copy(queue, residuals, target_buffer).wait()
    matrix_buffer.release()
    del matrix_buffer, matrix
    gc.collect()
    elapsed = time.monotonic() - started
    membership = {
        name: {"in_row_span_mod_prime": bool(not np.any(residuals[i])),
               "residual_term_count": int(np.count_nonzero(residuals[i]))}
        for i, (name, _) in enumerate(target_specs)
    }

    source_files = [builder.FAMILY21, builder.ANALYSIS21, builder.FAMILY32,
                    builder.ANALYSIS32, builder.FAMILY38, builder.ANALYSIS38,
                    builder.SUPPORT32_REDUNDANCY, builder.SUPPORT38_LOW_REDUCTION,
                    builder.SUPPORT12_GB, builder.SUPPORT12_INPUT, builder.SUPPORT12_REUSE]
    result = {
        "scope": "GF(251) Macaulay row-span reconnaissance for the d19 support-12/21/32/38 necessary subsystem; not a QQ membership proof or branch exclusion",
        "branch": "2+1+1, U1=0",
        "prime": PRIME,
        "target_weight": target_weight,
        "variables": list(names),
        "weights": list(weights),
        "generator_count": len(inputs),
        "generator_weight_degrees": generator_degrees,
        "macaulay_shape": [nrows, ncols],
        "macaulay_matrix_sha256": matrix_hash,
        "macaulay_matrix_bytes": matrix_bytes,
        "rank": rank,
        "pivot_columns": pivot_columns,
        "pure_power_membership": membership,
        "gpu": {"platform": platform.name, "device": device.name},
        "gpu_elimination_seconds": elapsed,
        "origin_common_zero": True,
        "input_sha256": {str(path.name): sha256_file(path) for path in source_files},
        "warning": "Modular row-span membership is reconnaissance only. Exact QQ witnesses and a proof that the forced origin violates the branch nonzero conditions are needed before claiming exclusion.",
    }
    output = Path(f"{OUTPUT_STEM}{target_weight}_20261028.json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "complete", "target_weight": target_weight,
                      "rank": rank, "shape": [nrows, ncols],
                      "pure_power_membership": membership,
                      "gpu": result["gpu"], "seconds": elapsed,
                      "report": str(output)}, indent=2), flush=True)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-weight", type=int, default=20)
    args = parser.parse_args()
    run(args.target_weight)


if __name__ == "__main__":
    main()
