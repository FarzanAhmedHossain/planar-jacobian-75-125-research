#!/usr/bin/env python3
"""One-byte Intel OpenCL Macaulay reconnaissance for support-23 equations.

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
BASIS_CERTIFICATE = Path("research/f2_75_125_211_support23_exact_rank4_certificate.json")
OUTPUT_STEM = "research/f2_75_125_211_depth14_support23_gpu_macaulay_active_rows"
PRIME = 251
SELECTED_ROWS = {86, 87, 88, 89}
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
inline uchar read_entry(__global const uchar *a0, __global const uchar *a1,
                        int split, int ncols, int row, int col) {
    return row < split ? a0[row * ncols + col]
                       : a1[(row - split) * ncols + col];
}
inline void write_entry(__global uchar *a0, __global uchar *a1,
                        int split, int ncols, int row, int col, uchar value) {
    if (row < split) a0[row * ncols + col] = value;
    else a1[(row - split) * ncols + col] = value;
}
__kernel void gather_column(__global const uchar *a0, __global const uchar *a1,
                            __global int *v, int split, int nrows,
                            int ncols, int col) {
    int r = (int)get_global_id(0);
    if (r < nrows) v[r] = (int)read_entry(a0, a1, split, ncols, r, col);
}
__kernel void gather_columns(__global const uchar *a0, __global const uchar *a1,
                             __global uchar *v, int split, int nrows,
                             int ncols, int first_col, int width) {
    int offset = (int)get_global_id(0);
    int row = (int)get_global_id(1);
    if (offset < width && row < nrows) {
        v[offset * nrows + row] = read_entry(a0, a1, split, ncols,
                                             row, first_col + offset);
    }
}
__kernel void swap_and_normalize(__global uchar *a0, __global uchar *a1,
                                 int split, int ncols, int rank, int pivot,
                                 int col, int prime, int inverse) {
    int c = (int)get_global_id(0);
    if (c >= ncols) return;
    if (pivot != rank) {
        uchar tmp = read_entry(a0, a1, split, ncols, rank, c);
        uchar val = read_entry(a0, a1, split, ncols, pivot, c);
        write_entry(a0, a1, split, ncols, rank, c, val);
        write_entry(a0, a1, split, ncols, pivot, c, tmp);
    }
    if (c >= col) {
        int val = (int)read_entry(a0, a1, split, ncols, rank, c);
        write_entry(a0, a1, split, ncols, rank, c,
                    (uchar)((val * inverse) % prime));
    }
}
__kernel void eliminate_below(__global uchar *a0, __global uchar *a1,
                              __global const int *active_row_ids,
                              int nactive, int split, int ncols,
                              int rank, int col, int prime) {
    int active_offset = (int)get_global_id(1);
    int col_offset = (int)get_global_id(0);
    int width = ncols - col - 1;
    if (active_offset >= nactive || col_offset >= width) return;
    int r = active_row_ids[active_offset];
    int c = col + 1 + col_offset;
    int factor = (int)read_entry(a0, a1, split, ncols, r, col);
    int value = (int)read_entry(a0, a1, split, ncols, r, c)
              - factor * (int)read_entry(a0, a1, split, ncols, rank, c);
    value %= prime;
    if (value < 0) value += prime;
    write_entry(a0, a1, split, ncols, r, c, (uchar)value);
}
__kernel void zero_pivot_column(__global uchar *a0, __global uchar *a1,
                                int split, int nrows, int ncols,
                                int rank, int col) {
    int offset = (int)get_global_id(0);
    int r = rank + 1 + offset;
    if (r < nrows) write_entry(a0, a1, split, ncols, r, col, (uchar)0);
}
__kernel void reduce_targets(__global const uchar *a0, __global const uchar *a1,
                             __global uchar *targets, int split, int ncols,
                             int rank, int col, int n_targets) {
    int col_offset = (int)get_global_id(0);
    int target = (int)get_global_id(1);
    int c = col + col_offset;
    if (target >= n_targets || c >= ncols) return;
    int base = target * ncols;
    int factor = (int)targets[base + col];
    int value = (int)targets[base + c]
              - factor * (int)read_entry(a0, a1, split, ncols, rank, c);
    value %= 251;
    if (value < 0) value += 251;
    targets[base + c] = (uchar)value;
}
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-weight", type=int, default=24)
    args = parser.parse_args()
    target_weight = args.target_weight
    output = Path(f"{OUTPUT_STEM}_w{target_weight}.json")
    source_hash = sha256_file(SOURCE)
    basis_certificate = json.loads(BASIS_CERTIFICATE.read_text(encoding="utf-8"))
    basis_certificate_hash = sha256_file(BASIS_CERTIFICATE)
    if basis_certificate["source_sha256"] != source_hash:
        raise AssertionError("rank-four certificate belongs to a different source")
    certified_rows = {int(row["row"]) for row in basis_certificate["selected_source_rows"]}
    if certified_rows != SELECTED_ROWS or basis_certificate["exact_rank"] != 4:
        raise AssertionError("unexpected exact rank-four basis rows")
    if not basis_certificate["all_source_rows_exactly_replayed_mod_251"]:
        raise AssertionError("the selected basis is not certified equivalent modulo 251")
    for relation in basis_certificate["row_relations"]:
        for coefficient in relation["coefficients"]:
            if Fraction(coefficient).denominator % PRIME == 0:
                raise AssertionError("a span relation denominator is not invertible modulo 251")
    symbols = [f"p{i}" for i in VARIABLES]
    target_monomials = monomials_of_weight(WEIGHTS, target_weight)
    target_index = {monomial: j for j, monomial in enumerate(target_monomials)}
    target_specs = []
    for variable_pos, variable in enumerate(VARIABLES):
        weight = WEIGHTS[variable_pos]
        if target_weight % weight == 0:
            power = target_weight // weight
            exponent = [0] * len(VARIABLES)
            exponent[variable_pos] = power
            target_specs.append((f"p{variable}^{power}", tuple(exponent)))
    multiplier_monomials = monomials_of_weight(WEIGHTS, target_weight - 14)
    if len(VARIABLES) != 23 or target_weight < 14:
        raise AssertionError("unexpected weighted monomial census")
    expected_census = {
        18: (15614, 13), 20: (34819, 41), 21: (51226, 71),
        22: (74986, 136), 23: (108560, 233), 24: (156376, 396),
    }
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

    if len(equations) != 4:
        raise AssertionError(f"expected the four certified basis equations, found {len(equations)}")

    nrows = len(equations) * len(multiplier_monomials)
    ncols = len(target_monomials)
    target_vectors = np.zeros((len(target_specs), ncols), dtype=np.uint8)
    for target_pos, (_, exponent) in enumerate(target_specs):
        target_vectors[target_pos, target_index[exponent]] = 1
    upload_rows = 16
    split_rows = (nrows // 2 // upload_rows) * upload_rows
    if split_rows <= 0 or split_rows >= nrows:
        raise AssertionError("invalid segmented-buffer split")
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
    gather_columns = cl.Kernel(program, "gather_columns")
    swap_and_normalize = cl.Kernel(program, "swap_and_normalize")
    eliminate_below = cl.Kernel(program, "eliminate_below")
    zero_pivot_column = cl.Kernel(program, "zero_pivot_column")
    reduce_targets = cl.Kernel(program, "reduce_targets")
    flags = cl.mem_flags
    bytes_per_entry = np.dtype(np.uint8).itemsize
    buffer0_bytes = split_rows * ncols * bytes_per_entry
    buffer1_bytes = (nrows - split_rows) * ncols * bytes_per_entry
    if max(buffer0_bytes, buffer1_bytes) > int(device.max_mem_alloc_size):
        raise MemoryError("a segmented matrix buffer exceeds CL_DEVICE_MAX_MEM_ALLOC_SIZE")
    if buffer0_bytes + buffer1_bytes > int(device.global_mem_size):
        raise MemoryError("the segmented matrix exceeds CL_DEVICE_GLOBAL_MEM_SIZE")
    matrix_buffer0 = cl.Buffer(context, flags.READ_WRITE, size=buffer0_bytes)
    matrix_buffer1 = cl.Buffer(context, flags.READ_WRITE, size=buffer1_bytes)
    target_buffer = cl.Buffer(context, flags.READ_WRITE, size=target_vectors.nbytes)
    active_row_ids_buffer = cl.Buffer(
        context, flags.READ_ONLY, size=nrows * np.dtype(np.int32).itemsize
    )
    pivot_batch_size = 64
    pivot_batch_buffer = cl.Buffer(
        context, flags.WRITE_ONLY, size=pivot_batch_size * nrows
    )
    pivot_batch = np.empty((pivot_batch_size, nrows), dtype=np.uint8)
    cl.enqueue_copy(queue, target_buffer, target_vectors).wait()
    matrix_digest = hashlib.sha256()

    def upload_rows_at(start_row: int, rows: np.ndarray) -> None:
        rows = np.ascontiguousarray(rows, dtype=np.uint8)
        matrix_digest.update(rows.tobytes(order="C"))
        if start_row < split_rows:
            first_count = min(len(rows), split_rows - start_row)
            if first_count:
                cl.enqueue_copy(
                    queue, matrix_buffer0, rows[:first_count],
                    dst_offset=start_row * ncols * bytes_per_entry,
                ).wait()
        else:
            first_count = 0
        if first_count < len(rows):
            second_start = max(start_row, split_rows) - split_rows
            cl.enqueue_copy(
                queue, matrix_buffer1, rows[first_count:],
                dst_offset=second_start * ncols * bytes_per_entry,
            ).wait()

    batch = np.zeros((upload_rows, ncols), dtype=np.uint8)
    row = 0
    for equation in equations:
        for multiplier in multiplier_monomials:
            for exponent, coefficient in equation:
                shifted = tuple(a + b for a, b in zip(exponent, multiplier))
                batch[row % upload_rows, target_index[shifted]] = coefficient
            row += 1
            if row % upload_rows == 0:
                upload_rows_at(row - upload_rows, batch)
                batch.fill(0)
    if row % upload_rows:
        final_count = row % upload_rows
        upload_rows_at(row - final_count, batch[:final_count])
    if row != nrows:
        raise AssertionError(f"matrix builder made {row} rows, expected {nrows}")
    matrix_hash = matrix_digest.hexdigest()
    print(json.dumps({
        "stage": "matrix-built", "source_sha256": source_hash,
        "basis_certificate_sha256": basis_certificate_hash,
        "equations": len(equations), "terms": sum(map(len, equations)),
        "shape": [nrows, ncols], "matrix_sha256": matrix_hash,
        "storage_dtype": "uint8", "buffer_split_rows": split_rows,
        "buffer_bytes": [buffer0_bytes, buffer1_bytes],
    }), flush=True)

    row_ids = np.arange(nrows, dtype=np.int32)
    pivot_columns = []
    selected_rows = []
    rank = 0
    scan_col = 0
    started = time.monotonic()
    while rank < nrows and scan_col < ncols:
        width = min(pivot_batch_size, ncols - scan_col)
        gather_columns(queue, (width, nrows), None,
                       matrix_buffer0, matrix_buffer1, pivot_batch_buffer,
                       np.int32(split_rows), np.int32(nrows),
                       np.int32(ncols), np.int32(scan_col), np.int32(width))
        cl.enqueue_copy(queue, pivot_batch[:width], pivot_batch_buffer).wait()
        active_columns = np.flatnonzero(
            np.any(pivot_batch[:width, rank:] != 0, axis=1)
        )
        if not len(active_columns):
            scan_col += width
            continue
        offset = int(active_columns[0])
        col = scan_col + offset
        column = pivot_batch[offset]
        nonzero = np.flatnonzero(column[rank:])
        pivot = rank + int(nonzero[0])
        pivot_value = int(column[pivot])
        active_row_ids = np.flatnonzero(column[rank + 1:] != 0).astype(np.int32)
        active_row_ids += np.int32(rank + 1)
        if pivot != rank:
            # The selected pivot row is swapped into `rank`; its old physical
            # row slot receives the former zero row at `rank`, so it is no
            # longer an elimination target.
            active_row_ids = active_row_ids[active_row_ids != pivot]
        selected_rows.append(int(row_ids[pivot]))
        if pivot != rank:
            row_ids[rank], row_ids[pivot] = row_ids[pivot], row_ids[rank]
        pivot_columns.append(col)
        scan_col = col + 1
        inverse = pow(pivot_value, PRIME - 2, PRIME)
        swap_and_normalize(queue, (ncols,), None,
                           matrix_buffer0, matrix_buffer1,
                           np.int32(split_rows), np.int32(ncols),
                           np.int32(rank), np.int32(pivot), np.int32(col),
                           np.int32(PRIME), np.int32(inverse))
        active_rows = nrows - rank - 1
        width = ncols - col - 1
        if len(active_row_ids) and width:
            cl.enqueue_copy(queue, active_row_ids_buffer, active_row_ids).wait()
            eliminate_below(queue, (width, len(active_row_ids)), None,
                            matrix_buffer0, matrix_buffer1, active_row_ids_buffer,
                            np.int32(len(active_row_ids)), np.int32(split_rows),
                            np.int32(ncols), np.int32(rank), np.int32(col),
                            np.int32(PRIME))
        if active_rows:
            zero_pivot_column(queue, (active_rows,), None,
                              matrix_buffer0, matrix_buffer1,
                              np.int32(split_rows), np.int32(nrows),
                              np.int32(ncols), np.int32(rank), np.int32(col))
        target_width = ncols - col
        if target_width and target_specs:
            reduce_targets(queue, (target_width, len(target_specs)), None,
                           matrix_buffer0, matrix_buffer1, target_buffer,
                           np.int32(split_rows), np.int32(ncols),
                           np.int32(rank), np.int32(col),
                           np.int32(len(target_specs)))
        rank += 1
        if rank % 25 == 0:
            print(f"rank_progress={rank} columns_scanned={col + 1}", flush=True)

    residuals = np.empty_like(target_vectors)
    cl.enqueue_copy(queue, residuals, target_buffer).wait()
    elapsed = time.monotonic() - started

    pure_power_membership = {}
    for target_pos, (target_name, _) in enumerate(target_specs):
        residual_support = np.flatnonzero(residuals[target_pos])
        pure_power_membership[target_name] = {
            "in_row_span_mod_prime": len(residual_support) == 0,
            "residual_term_count": int(len(residual_support)),
        }

    result = {
        "source": str(SOURCE), "source_sha256": source_hash,
        "scope": f"four exact basis equations for the 35 support-23 compatibility equations, weighted Macaulay degree {target_weight}; modular reconnaissance only",
        "branch": "2+1+1", "field": f"GF({PRIME})", "prime": PRIME,
        "basis_certificate": str(BASIS_CERTIFICATE),
        "basis_certificate_sha256": basis_certificate_hash,
        "support_variables": symbols, "weights": list(WEIGHTS),
        "selected_rows_zero_based": sorted(SELECTED_ROWS), "selected_tags": tags,
        "input_term_count": sum(map(len, equations)),
        "macaulay_shape": [int(nrows), int(ncols)],
        "multiplier_monomial_count": len(multiplier_monomials),
        "target_weight": target_weight, "rank": rank,
        "pivot_columns": pivot_columns, "selected_original_rows": selected_rows,
        "matrix_sha256": matrix_hash, "gpu": {"platform": platform.name, "device": device.name},
        "gpu_elimination_seconds": elapsed,
        "matrix_storage_dtype": "uint8",
        "elimination_indexing": "coalesced 2D column offsets over only nonzero pivot rows, with streaming target residuals",
        "buffer_split_rows": split_rows,
        "matrix_buffer_bytes": [buffer0_bytes, buffer1_bytes],
        "pure_power_membership": pure_power_membership,
        "warning": "The four generators are exactly row-span equivalent to the 35 selected equations, but a modular membership result is not a characteristic-zero ideal membership or branch exclusion. Exact rational certificates and normalization nonzero conditions are still required.",
    }
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"rank": rank, "shape": [nrows, ncols],
                      "pure_power_membership": pure_power_membership,
                      "gpu": result["gpu"], "seconds": elapsed}, indent=2), flush=True)
    print(f"receipt={output} sha256={sha256_file(output)}", flush=True)


if __name__ == "__main__":
    main()
