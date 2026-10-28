#!/usr/bin/env python3
"""Intel OpenCL modular Macaulay rank for the full active depth-13 2+1+1 system.

This is a replayable rank reconnaissance calculation, not a branch-exclusion
certificate. A nonzero modular minor certifies a rational lower bound on the
Macaulay rank; failure to reach full rank has no exclusion consequence.
"""

from __future__ import annotations

import hashlib
import json
import time
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np
import pyopencl as cl

from f2_75_125_211_degree_structure import streamed_array


SOURCE = Path("research/f2_75_125_211_p43zero_depth13_constraints.jsonfrag")
TARGET_WEIGHT = 14
PRIME = 32749
VARIABLES = tuple(range(36)) + tuple(range(37, 43))
WEIGHT_BY_VARIABLE = {
    **{0: 1, 1: 2},
    **{i: 3 for i in range(2, 5)},
    **{i: 4 for i in range(5, 13)},
    **{i: 5 for i in range(13, 20)},
    **{i: 6 for i in range(20, 25)},
    **{i: 7 for i in range(25, 30)},
    **{i: 8 for i in range(30, 36)},
    **{i: 9 for i in range(37, 43)},
}
WEIGHTS = tuple(WEIGHT_BY_VARIABLE[i] for i in VARIABLES)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def monomials_of_weight(weights: tuple[int, ...], target: int) -> list[tuple[int, ...]]:
    @lru_cache(None)
    def rec(index: int, remaining: int) -> tuple[tuple[int, ...], ...]:
        if index == len(weights):
            return ((),) if remaining == 0 else ()
        out: list[tuple[int, ...]] = []
        for exponent in range(remaining // weights[index] + 1):
            for tail in rec(index + 1, remaining - exponent * weights[index]):
                out.append((exponent,) + tail)
        return tuple(out)

    return list(rec(0, target))


@lru_cache(8192)
def inverse_mod(denominator: int) -> int:
    denominator %= PRIME
    if denominator == 0:
        raise ValueError(f"prime {PRIME} divides a source coefficient denominator")
    return pow(denominator, -1, PRIME)


def coefficient_mod(value: Fraction) -> int:
    return (value.numerator % PRIME) * inverse_mod(value.denominator) % PRIME


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
    if (c >= col)
        a[rank * ncols + c] = (a[rank * ncols + c] * inverse) % prime;
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
    if not is_prime(PRIME):
        raise AssertionError("chosen modular characteristic is not prime")
    if len(VARIABLES) != 42 or len(WEIGHT_BY_VARIABLE) != 42:
        raise AssertionError("unexpected active-variable list")
    source_hash = sha256_file(SOURCE)
    variable_index = {variable: i for i, variable in enumerate(VARIABLES)}
    expected_counts = {8: 112, 9: 112, 10: 115, 11: 115, 12: 149, 13: 151, 14: 0}
    monomials = monomials_of_weight(WEIGHTS, TARGET_WEIGHT)
    monomial_index = {monomial: i for i, monomial in enumerate(monomials)}
    multipliers = {
        degree: monomials_of_weight(WEIGHTS, TARGET_WEIGHT - degree)
        for degree in range(8, TARGET_WEIGHT + 1)
    }
    expected_rows = sum(expected_counts[d] * len(multipliers[d]) for d in multipliers)
    matrix = np.zeros((expected_rows, len(monomials)), dtype=np.int32)
    equation_counts: Counter[int] = Counter()
    row = 0
    term_count = 0
    for tag, saved in streamed_array(SOURCE):
        terms: defaultdict[tuple[int, ...], Fraction] = defaultdict(Fraction)
        for sparse, values in saved:
            if len(values) != 1:
                raise ValueError(f"non-rational coefficient encoding in {tag}")
            exponent_map = {i - 1: e for i, e in sparse if 1 <= i <= 75 and e}
            if not exponent_map.keys() <= variable_index.keys():
                raise ValueError(f"parameter outside the 42-variable active set in {tag}")
            exponent = tuple(exponent_map.get(i, 0) for i in VARIABLES)
            terms[exponent] += Fraction(values[0])
            term_count += 1
        terms = defaultdict(Fraction, {e: c for e, c in terms.items() if c})
        if not terms:
            continue
        degrees = {sum(w * e for w, e in zip(WEIGHTS, exponent)) for exponent in terms}
        if len(degrees) != 1:
            raise AssertionError(f"nonhomogeneous source equation {tag}: {sorted(degrees)}")
        degree = next(iter(degrees))
        if degree not in expected_counts:
            raise AssertionError(f"unexpected source weighted degree {degree}: {tag}")
        equation_counts[degree] += 1
        if degree <= TARGET_WEIGHT:
            modular_terms = [(e, coefficient_mod(c)) for e, c in terms.items()]
            for multiplier in multipliers[degree]:
                for exponent, coefficient in modular_terms:
                    shifted = tuple(a + b for a, b in zip(exponent, multiplier))
                    matrix[row, monomial_index[shifted]] = coefficient
                row += 1

    actual_counts = {d: equation_counts[d] for d in expected_counts}
    if actual_counts != expected_counts:
        raise AssertionError(f"weighted generator counts changed: {actual_counts}")
    if row != expected_rows:
        raise AssertionError("Macaulay matrix row construction mismatch")
    matrix_hash = hashlib.sha256(matrix.tobytes(order="C")).hexdigest()
    print(json.dumps({
        "stage": "matrix-built",
        "source_sha256": source_hash,
        "equations": actual_counts,
        "source_terms_scanned": term_count,
        "shape": list(matrix.shape),
        "monomials_by_weight": len(monomials),
        "modular_matrix_sha256": matrix_hash,
    }), flush=True)

    platforms = cl.get_platforms()
    candidates = [(p, d) for p in platforms for d in p.get_devices()
                  if d.type & cl.device_type.GPU]
    if not candidates:
        raise RuntimeError("no OpenCL GPU device is available")
    platform, device = next(
        ((p, d) for p, d in candidates if "Intel" in d.name), candidates[0]
    )
    context = cl.Context([device])
    queue = cl.CommandQueue(context, device)
    program = cl.Program(context, KERNEL_SOURCE).build()
    gather_column = cl.Kernel(program, "gather_column")
    swap_and_normalize = cl.Kernel(program, "swap_and_normalize")
    eliminate_below = cl.Kernel(program, "eliminate_below")
    zero_pivot_column = cl.Kernel(program, "zero_pivot_column")
    flags = cl.mem_flags
    matrix_buffer = cl.Buffer(context, flags.READ_WRITE | flags.COPY_HOST_PTR,
                              hostbuf=matrix)
    matrix_shape = matrix.shape
    del matrix
    column_buffer = cl.Buffer(context, flags.WRITE_ONLY,
                               size=matrix_shape[0] * np.dtype(np.int32).itemsize)
    column = np.empty(matrix_shape[0], dtype=np.int32)
    row_ids = np.arange(matrix_shape[0], dtype=np.int32)
    selected_rows: list[int] = []
    pivot_columns: list[int] = []
    pivot_values: list[int] = []
    swaps = 0
    rank = 0
    started = time.monotonic()
    nrows, ncols = matrix_shape
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
        if pivot_value == 0:
            raise AssertionError("GPU pivot search selected a zero")
        selected_rows.append(int(row_ids[pivot]))
        if pivot != rank:
            swaps += 1
            row_ids[rank], row_ids[pivot] = row_ids[pivot], row_ids[rank]
        pivot_columns.append(col)
        pivot_values.append(pivot_value)
        inverse = pow(pivot_value, PRIME - 2, PRIME)
        swap_and_normalize(
            queue, (ncols,), None, matrix_buffer,
            np.int32(ncols), np.int32(rank), np.int32(pivot), np.int32(col),
            np.int32(PRIME), np.int32(inverse),
        )
        active_rows = nrows - rank - 1
        width = ncols - col - 1
        if active_rows and width:
            eliminate_below(
                queue, (active_rows * width,), None, matrix_buffer,
                np.int32(nrows), np.int32(ncols), np.int32(rank),
                np.int32(col), np.int32(PRIME),
            )
        if active_rows:
            zero_pivot_column(
                queue, (active_rows,), None, matrix_buffer,
                np.int32(nrows), np.int32(ncols), np.int32(rank),
                np.int32(col),
            )
        rank += 1
        if rank % 100 == 0:
            queue.finish()
            print(json.dumps({"stage": "gpu-elimination", "rank": rank,
                              "column": col, "seconds": round(time.monotonic()-started, 2)}),
                  flush=True)
    queue.finish()

    selected_columns = pivot_columns[:]
    minor_determinant = None if not rank else 1
    if rank:
        # At each pivot, the unnormalized pivot is the next Schur-complement
        # factor for the minor on the selected source rows and pivot columns.
        for pivot_value in pivot_values:
            minor_determinant = (minor_determinant * pivot_value) % PRIME
        if minor_determinant == 0:
            raise AssertionError("selected modular pivot minor has zero determinant")

    result = {
        "source": str(SOURCE),
        "source_sha256": source_hash,
        "branch": "2+1+1",
        "scope": "weight-14 Macaulay matrix for the 42 active variables in the detached depth-13 necessary subsystem only; not a whole-branch exclusion",
        "coefficient_field": "Q rank lower bound certified by a nonzero minor modulo a prime",
        "prime": PRIME,
        "device": {"platform": platform.name, "vendor": platform.vendor,
                   "name": device.name, "driver": device.driver_version,
                   "global_memory_bytes": int(device.global_mem_size)},
        "positive_integer_weights": list(WEIGHTS),
        "variables": [f"p{i}" for i in VARIABLES],
        "selected_generators": actual_counts,
        "source_terms_scanned": term_count,
        "target_weight": TARGET_WEIGHT,
        "monomial_count": len(monomials),
        "matrix_shape": list(matrix_shape),
        "modular_matrix_sha256": matrix_hash,
        "rank_mod_prime": rank,
        "independent_minor_rows": selected_rows,
        "independent_minor_columns": selected_columns,
        "independent_minor_determinant_mod_prime": minor_determinant,
        "row_swap_count": swaps,
        "elapsed_gpu_seconds": round(time.monotonic() - started, 3),
        "interpretation": "The nonzero minor proves the rational Macaulay rank is at least the reported modular rank. Only a full column rank would show that every weight-14 form belongs to this necessary-equation subideal; even that would not by itself exclude the full 2+1+1 branch.",
    }
    receipt = Path("research/f2_75_125_211_gpu_weight14_macaulay_20261021.json")
    receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
