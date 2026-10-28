#!/usr/bin/env python3
"""Independently reconstruct and verify the GPU Macaulay minor modulo p."""

from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

from f2_75_125_211_degree_structure import streamed_array


SOURCE = Path("research/f2_75_125_211_p43zero_depth13_constraints.jsonfrag")
GPU_RECEIPT = Path("research/f2_75_125_211_gpu_weight14_macaulay_20261021.json")
OUT = Path("research/f2_75_125_211_gpu_weight14_minor_check_20261021.json")
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
TARGET_WEIGHT = 14


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def monomials_of_weight(weights: tuple[int, ...], target: int):
    output: list[tuple[int, ...]] = []
    current = [0] * len(weights)

    def visit(index: int, remaining: int) -> None:
        if index == len(weights):
            if remaining == 0:
                output.append(tuple(current))
            return
        for exponent in range(remaining // weights[index] + 1):
            current[index] = exponent
            visit(index + 1, remaining - exponent * weights[index])
        current[index] = 0

    visit(0, target)
    return output


def determinant_mod(matrix: list[list[int]], prime: int) -> int:
    a = [row[:] for row in matrix]
    n = len(a)
    if any(len(row) != n for row in a):
        raise ValueError("minor is not square")
    determinant = 1
    for col in range(n):
        pivot = next((row for row in range(col, n) if a[row][col] % prime), None)
        if pivot is None:
            return 0
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            determinant = -determinant
        value = a[col][col] % prime
        determinant = determinant * value % prime
        inverse = pow(value, -1, prime)
        for row in range(col + 1, n):
            factor = a[row][col] * inverse % prime
            if factor:
                for j in range(col + 1, n):
                    a[row][j] = (a[row][j] - factor * a[col][j]) % prime
            a[row][col] = 0
    return determinant % prime


def main() -> None:
    source_hash = sha256_file(SOURCE)
    gpu = json.loads(GPU_RECEIPT.read_text(encoding="utf-8"))
    if source_hash != gpu["source_sha256"]:
        raise ValueError("GPU receipt and source equation hash do not match")
    prime = int(gpu["prime"])
    selected_rows = [int(i) for i in gpu["independent_minor_rows"]]
    selected_columns = [int(i) for i in gpu["independent_minor_columns"]]
    expected_det = int(gpu["independent_minor_determinant_mod_prime"])
    rank = int(gpu["rank_mod_prime"])
    if len(selected_rows) != rank or len(selected_columns) != rank:
        raise ValueError("minor dimensions do not match the reported modular rank")

    variable_index = {variable: i for i, variable in enumerate(VARIABLES)}
    monomials = monomials_of_weight(WEIGHTS, TARGET_WEIGHT)
    monomial_index = {monomial: i for i, monomial in enumerate(monomials)}
    multipliers = {
        degree: monomials_of_weight(WEIGHTS, TARGET_WEIGHT - degree)
        for degree in range(8, TARGET_WEIGHT + 1)
    }
    row_position = {source_row: position for position, source_row in enumerate(selected_rows)}
    column_position = {source_col: position for position, source_col in enumerate(selected_columns)}
    minor = [[0] * rank for _ in range(rank)]
    equation_counts: Counter[int] = Counter()
    source_terms = 0
    macaulay_row = 0

    for tag, saved in streamed_array(SOURCE):
        terms: defaultdict[tuple[int, ...], Fraction] = defaultdict(Fraction)
        for sparse, values in saved:
            if len(values) != 1:
                raise ValueError(f"non-rational coefficient encoding in {tag}")
            exponent_map = {i - 1: e for i, e in sparse if 1 <= i <= 75 and e}
            if not exponent_map.keys() <= variable_index.keys():
                raise ValueError(f"parameter outside the active set in {tag}")
            exponent = tuple(exponent_map.get(i, 0) for i in VARIABLES)
            terms[exponent] += Fraction(values[0])
            source_terms += 1
        terms = defaultdict(Fraction, {e: c for e, c in terms.items() if c})
        if not terms:
            continue
        degrees = {sum(w * e for w, e in zip(WEIGHTS, exponent)) for exponent in terms}
        if len(degrees) != 1:
            raise AssertionError(f"nonhomogeneous source equation {tag}")
        degree = next(iter(degrees))
        equation_counts[degree] += 1
        if degree > TARGET_WEIGHT:
            continue

        modular_terms: list[tuple[tuple[int, ...], int]] = []
        for exponent, coefficient in terms.items():
            denominator = coefficient.denominator % prime
            if denominator == 0:
                raise ValueError(f"source denominator divisible by {prime}")
            value = coefficient.numerator % prime * pow(denominator, -1, prime) % prime
            modular_terms.append((exponent, value))

        for multiplier in multipliers[degree]:
            position = row_position.get(macaulay_row)
            if position is not None:
                for exponent, coefficient in modular_terms:
                    shifted = tuple(a + b for a, b in zip(exponent, multiplier))
                    source_column = monomial_index[shifted]
                    target_column = column_position.get(source_column)
                    if target_column is not None:
                        minor[position][target_column] = coefficient
            macaulay_row += 1

    expected_counts = {8: 112, 9: 112, 10: 115, 11: 115, 12: 149, 13: 151, 14: 0}
    actual_counts = {degree: equation_counts[degree] for degree in expected_counts}
    if actual_counts != expected_counts or macaulay_row != int(gpu["matrix_shape"][0]):
        raise AssertionError("reconstructed equation or Macaulay row counts differ")

    computed_det = determinant_mod(minor, prime)
    result = {
        "source_sha256": source_hash,
        "gpu_receipt_sha256": sha256_file(GPU_RECEIPT),
        "prime": prime,
        "selected_minor_shape": [rank, rank],
        "selected_macaulay_rows": selected_rows,
        "selected_monomial_columns": selected_columns,
        "equation_counts": actual_counts,
        "macaulay_rows_reconstructed": macaulay_row,
        "source_terms_scanned": source_terms,
        "gpu_reported_determinant_mod_prime": expected_det,
        "independently_reconstructed_determinant_mod_prime": computed_det,
        "matches": computed_det == expected_det,
        "consequence": "The independently replayed nonzero minor certifies rational Macaulay rank at least the reported value for this truncated necessary system only.",
    }
    if computed_det == 0 or computed_det != expected_det:
        raise AssertionError("independent minor determinant check failed")
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
