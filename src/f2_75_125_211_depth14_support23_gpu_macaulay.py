#!/usr/bin/env python3
"""Intel OpenCL Macaulay reconnaissance for the 35 support-23 equations.

The modular row-space and monomial-membership results are reconnaissance only.
Any proposed characteristic-zero membership must be replayed exactly over Q.
"""

from __future__ import annotations

import hashlib
import json
import time
import argparse
import gc
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import ijson
import numpy as np
import pyopencl as cl


SOURCE = Path("research/f2_75_125_211_depth14_global_compatibility_reduction_20261021.json")
OUTPUT_STEM = "research/f2_75_125_211_depth14_support23_gpu_macaulay"
PRIME = 32749
SELECTED_ROWS = set(range(86, 121))
VARIABLES = tuple([*range(5), *range(6, 13), *range(14, 25)])
WEIGHT_BY_VARIABLE = {
    0: 1, 1: 2,
    **{i: 3 for i in range(2, 5)},
    **{i: 4 for i in range(6, 13)},
    **{i: 5 for i in range(14, 20)},
    **{i: 6 for i in range(20, 25)},
}
WEIGHTS = tuple(WEIGHT_BY_VARIABLE[i] for i in VARIABLES)


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


def mod_fraction(value: Fraction) -> int:
    denominator = value.denominator % PRIME
    if denominator == 0:
        raise ValueError("the modular prime divides a coefficient denominator")
    return (value.numerator % PRIME) * pow(denominator, -1, PRIME) % PRIME


KERNEL_SOURCE = r"""
__kernel void gather_column(__global const int *a, __global int *v,
                            int nrows, int ncols, int col) {
    int r = (int)get_global_id(0);
    if (r < nrows) v[r] = a[r * ncols + col];
}
__kernel void swap_and_normalize(__global int *a, int ncols, int rank,
                                 int pivot, int col, int prime, int inverse) {
    int c = (int)get_global_id(0);
    if (c >= ncols) return;
    if (pivot != rank) {
        int tmp = a[rank * ncols + c];
        a[rank * ncols + c] = a[pivot * ncols + c];
        a[pivot * ncols + c] = tmp;
    }
    if (c >= col) a[rank * ncols + c] = (a[rank * ncols + c] * inverse) % prime;
}
__kernel void eliminate_below(__global int *a, int nrows, int ncols,
                              int rank, int col, int prime) {
    int width = ncols - col - 1;
    int active_rows = nrows - rank - 1;
    int i = (int)get_global_id(0);
    if (i >= active_rows * width) return;
    int r = rank + 1 + i / width;
    int c = col + 1 + i % width;
    int factor = a[r * ncols + col];
    int value = a[r * ncols + c] - factor * a[rank * ncols + c];
    value %= prime;
    if (value < 0) value += prime;
    a[r * ncols + c] = value;
}
__kernel void zero_pivot_column(__global int *a, int nrows, int ncols,
                                int rank, int col) {
    int offset = (int)get_global_id(0);
    int r = rank + 1 + offset;
    if (r < nrows) a[r * ncols + col] = 0;
}
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-weight", type=int, default=18)
    args = parser.parse_args()
    target_weight = args.target_weight
    output = Path(f"{OUTPUT_STEM}_w{target_weight}.json")
    source_hash = sha256_file(SOURCE)
    symbols = [f"p{i}" for i in VARIABLES]
    target_monomials = monomials_of_weight(WEIGHTS, target_weight)
    target_index = {monomial: j for j, monomial in enumerate(target_monomials)}
    multiplier_monomials = monomials_of_weight(WEIGHTS, target_weight - 14)
    if len(VARIABLES) != 23 or target_weight < 14:
        raise AssertionError("unexpected weighted monomial census")
    expected_census = {18: (15614, 13), 20: (34819, 41), 21: (51226, 71)}
    if target_weight in expected_census and (
        len(target_monomials), len(multiplier_monomials)
    ) != expected_census[target_weight]:
        raise AssertionError("weighted monomial census changed")

    equations = []
    tags = []
    with SOURCE.open("rb") as stream:
        for equation in ijson.items(stream, "global_compatibility_equations.item"):
            if equation["constraint_row_zero_based"] not in SELECTED_ROWS:
                continue
            if equation["support_variables"] != symbols:
                raise AssertionError(f"support mismatch at {equation['constraint_tag']}")
            terms = []
            for sparse, coefficient in equation["terms"]:
                exponent = [0] * len(VARIABLES)
                for variable, power in sparse:
                    if variable not in VARIABLES:
                        raise AssertionError(f"unexpected variable index {variable}")
                    exponent[VARIABLES.index(variable)] = power
                if sum(w * e for w, e in zip(WEIGHTS, exponent)) != 14:
                    raise AssertionError(f"wrong weighted degree at {equation['constraint_tag']}")
                terms.append((tuple(exponent), mod_fraction(Fraction(coefficient))))
            equations.append(terms)
            tags.append(equation["constraint_tag"])

    if len(equations) != 35:
        raise AssertionError(f"expected 35 equations, found {len(equations)}")

    matrix = np.zeros((len(equations) * len(multiplier_monomials), len(target_monomials)),
                      dtype=np.int32)
    row = 0
    for equation in equations:
        for multiplier in multiplier_monomials:
            for exponent, coefficient in equation:
                shifted = tuple(a + b for a, b in zip(exponent, multiplier))
                matrix[row, target_index[shifted]] = coefficient
            row += 1
    matrix_hash = hashlib.sha256(matrix.tobytes(order="C")).hexdigest()
    print(json.dumps({
        "stage": "matrix-built", "source_sha256": source_hash,
        "equations": len(equations), "terms": sum(map(len, equations)),
        "shape": list(matrix.shape), "matrix_sha256": matrix_hash,
    }), flush=True)

    platforms = cl.get_platforms()
    candidates = [(platform, device) for platform in platforms for device in platform.get_devices()
                  if device.type & cl.device_type.GPU]
    if not candidates:
        raise RuntimeError("no OpenCL GPU device is available")
    platform, device = next(((p, d) for p, d in candidates if "Intel" in d.name), candidates[0])
    context = cl.Context([device])
    queue = cl.CommandQueue(context, device)
    program = cl.Program(context, KERNEL_SOURCE).build()
    gather_column = cl.Kernel(program, "gather_column")
    swap_and_normalize = cl.Kernel(program, "swap_and_normalize")
    eliminate_below = cl.Kernel(program, "eliminate_below")
    zero_pivot_column = cl.Kernel(program, "zero_pivot_column")
    flags = cl.mem_flags
    matrix_buffer = cl.Buffer(context, flags.READ_WRITE | flags.COPY_HOST_PTR, hostbuf=matrix)
    nrows, ncols = matrix.shape
    del matrix
    gc.collect()
    column_buffer = cl.Buffer(context, flags.WRITE_ONLY, size=nrows * np.dtype(np.int32).itemsize)
    column = np.empty(nrows, dtype=np.int32)
    row_ids = np.arange(nrows, dtype=np.int32)
    pivot_columns = []
    selected_rows = []
    rank = 0
    started = time.monotonic()
    for col in range(ncols):
        if rank == nrows:
            break
        gather_column(queue, (nrows,), None, matrix_buffer, column_buffer,
                      np.int32(nrows), np.int32(ncols), np.int32(col))
        cl.enqueue_copy(queue, column, column_buffer).wait()
        nonzero = np.flatnonzero(column[rank:])
        if not len(nonzero):
            continue
        pivot = rank + int(nonzero[0])
        pivot_value = int(column[pivot])
        selected_rows.append(int(row_ids[pivot]))
        if pivot != rank:
            row_ids[rank], row_ids[pivot] = row_ids[pivot], row_ids[rank]
        pivot_columns.append(col)
        inverse = pow(pivot_value, PRIME - 2, PRIME)
        swap_and_normalize(queue, (ncols,), None, matrix_buffer,
                           np.int32(ncols), np.int32(rank), np.int32(pivot),
                           np.int32(col), np.int32(PRIME), np.int32(inverse))
        active_rows = nrows - rank - 1
        width = ncols - col - 1
        if active_rows and width:
            eliminate_below(queue, (active_rows * width,), None, matrix_buffer,
                            np.int32(nrows), np.int32(ncols), np.int32(rank),
                            np.int32(col), np.int32(PRIME))
        if active_rows:
            zero_pivot_column(queue, (active_rows,), None, matrix_buffer,
                              np.int32(nrows), np.int32(ncols),
                              np.int32(rank), np.int32(col))
        rank += 1
        if rank % 50 == 0:
            print(f"rank_progress={rank} columns_scanned={col + 1}", flush=True)

    reduced = np.empty((rank, ncols), dtype=np.int32)
    cl.enqueue_copy(queue, reduced, matrix_buffer).wait()
    elapsed = time.monotonic() - started

    pure_power_membership = {}
    for variable_pos, variable in enumerate(VARIABLES):
        weight = WEIGHTS[variable_pos]
        if target_weight % weight:
            continue
        power = target_weight // weight
        exponent = [0] * len(VARIABLES)
        exponent[variable_pos] = power
        target = np.zeros(ncols, dtype=np.int32)
        target[target_index[tuple(exponent)]] = 1
        for pivot_row, pivot_col in enumerate(pivot_columns):
            factor = int(target[pivot_col])
            if factor:
                target[pivot_col:] = (target[pivot_col:].astype(np.int64)
                                      - factor * reduced[pivot_row, pivot_col:].astype(np.int64)) % PRIME
        residual_support = np.flatnonzero(target)
        pure_power_membership[f"p{variable}^{power}"] = {
            "in_row_span_mod_prime": len(residual_support) == 0,
            "residual_term_count": int(len(residual_support)),
        }

    result = {
        "source": str(SOURCE), "source_sha256": source_hash,
        "scope": f"35 support-23 compatibility equations, weighted Macaulay degree {target_weight}; modular reconnaissance only",
        "branch": "2+1+1", "field": f"GF({PRIME})", "prime": PRIME,
        "support_variables": symbols, "weights": list(WEIGHTS),
        "selected_rows_zero_based": sorted(SELECTED_ROWS), "selected_tags": tags,
        "input_term_count": sum(map(len, equations)),
        "macaulay_shape": [int(nrows), int(ncols)],
        "multiplier_monomial_count": len(multiplier_monomials),
        "target_weight": target_weight, "rank": rank,
        "pivot_columns": pivot_columns, "selected_original_rows": selected_rows,
        "matrix_sha256": matrix_hash, "gpu": {"platform": platform.name, "device": device.name},
        "gpu_elimination_seconds": elapsed,
        "pure_power_membership": pure_power_membership,
        "warning": "A modular membership result is not a characteristic-zero ideal membership or branch exclusion. Exact rational certificates and normalization nonzero conditions are still required.",
    }
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"rank": rank, "shape": [nrows, ncols],
                      "pure_power_membership": pure_power_membership,
                      "gpu": result["gpu"], "seconds": elapsed}, indent=2), flush=True)
    print(f"receipt={output} sha256={sha256_file(output)}", flush=True)


if __name__ == "__main__":
    main()
