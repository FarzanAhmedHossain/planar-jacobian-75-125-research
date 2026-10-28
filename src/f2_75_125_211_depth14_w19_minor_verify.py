"""Independently verify the selected modular Macaulay minor at weight 19.

The OpenCL run records 88 pivot rows and 88 pivot columns for a 770 by 23,365
Macaulay matrix over GF(32749).  This checker reconstructs only that 88 by 88
minor from the hashed rational source and verifies its determinant directly.
A nonzero determinant modulo a prime not dividing any input denominator gives
an exact characteristic-zero lower bound of 88 on the rational matrix rank.
It does not turn the modular row-space screen into a branch exclusion.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import ijson


SOURCE = Path(
    "research/f2_75_125_211_depth14_global_compatibility_reduction_20261021.json"
)
GPU_RECEIPT = Path(
    "research/f2_75_125_211_depth14_support23_gpu_macaulay_w19.json"
)
OUTPUT = Path(
    "research/f2_75_125_211_depth14_w19_minor_verification.json"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


@lru_cache(None)
def monomials_of_weight(weights: tuple[int, ...], target: int):
    @lru_cache(None)
    def rec(index: int, remaining: int):
        if index == len(weights):
            return ((),) if remaining == 0 else ()
        out = []
        for exponent in range(remaining // weights[index] + 1):
            out.extend(
                (exponent,) + tail
                for tail in rec(index + 1, remaining - exponent * weights[index])
            )
        return tuple(out)

    return rec(0, target)


def reduce_fraction(value: str | int, prime: int) -> int:
    rational = Fraction(value)
    denominator = rational.denominator % prime
    if denominator == 0:
        raise ValueError("prime divides an input coefficient denominator")
    return (rational.numerator % prime) * pow(denominator, -1, prime) % prime


def determinant_mod(matrix: list[list[int]], prime: int) -> int:
    n = len(matrix)
    if any(len(row) != n for row in matrix):
        raise ValueError("minor is not square")
    determinant = 1
    for col in range(n):
        pivot = next((row for row in range(col, n) if matrix[row][col]), None)
        if pivot is None:
            return 0
        if pivot != col:
            matrix[col], matrix[pivot] = matrix[pivot], matrix[col]
            determinant = (-determinant) % prime
        value = matrix[col][col] % prime
        determinant = determinant * value % prime
        inverse = pow(value, prime - 2, prime)
        for row in range(col + 1, n):
            entry = matrix[row][col] % prime
            if not entry:
                continue
            factor = entry * inverse % prime
            for j in range(col + 1, n):
                matrix[row][j] = (
                    matrix[row][j] - factor * matrix[col][j]
                ) % prime
            matrix[row][col] = 0
    return determinant


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite {OUTPUT}")
    receipt = json.loads(GPU_RECEIPT.read_text(encoding="utf-8"))
    source_hash = sha256(SOURCE)
    gpu_receipt_hash = sha256(GPU_RECEIPT)
    if source_hash != receipt["source_sha256"]:
        raise AssertionError("source hash differs from the OpenCL receipt")
    if receipt["branch"] != "2+1+1" or receipt["target_weight"] != 19:
        raise AssertionError("unexpected branch or Macaulay weight")

    prime = int(receipt["prime"])
    weights = tuple(map(int, receipt["weights"]))
    variable_ids = tuple(int(name[1:]) for name in receipt["support_variables"])
    target_weight = int(receipt["target_weight"])
    target_monomials = monomials_of_weight(weights, target_weight)
    multiplier_monomials = monomials_of_weight(weights, target_weight - 14)
    if [len(receipt["selected_rows_zero_based"]) * len(multiplier_monomials),
        len(target_monomials)] != receipt["macaulay_shape"]:
        raise AssertionError("Macaulay shape does not reconstruct")

    selected_sources = tuple(map(int, receipt["selected_rows_zero_based"]))
    selected_set = set(selected_sources)
    equations: dict[int, list[tuple[tuple[int, ...], int]]] = {}
    positions = {variable: index for index, variable in enumerate(variable_ids)}
    with SOURCE.open("rb") as stream:
        for equation in ijson.items(stream, "global_compatibility_equations.item"):
            source_index = int(equation["constraint_row_zero_based"])
            if source_index not in selected_set:
                continue
            if equation["support_variables"] != receipt["support_variables"]:
                raise AssertionError("selected equation support changed")
            terms = []
            for sparse, coefficient in equation["terms"]:
                exponent = [0] * len(variable_ids)
                for variable, power in sparse:
                    exponent[positions[int(variable)]] = int(power)
                if sum(w * e for w, e in zip(weights, exponent)) != 14:
                    raise AssertionError("selected equation is not weight 14")
                terms.append((tuple(exponent), reduce_fraction(coefficient, prime)))
            equations[source_index] = terms

    if tuple(sorted(equations)) != selected_sources:
        raise AssertionError("not all selected equations were reconstructed")

    target_index = {monomial: index for index, monomial in enumerate(target_monomials)}
    multiplier_count = len(multiplier_monomials)
    row_ids = list(map(int, receipt["selected_original_rows"]))
    pivot_columns = list(map(int, receipt["pivot_columns"]))
    if len(row_ids) != 88 or len(pivot_columns) != 88:
        raise AssertionError("unexpected pivot count")
    if len(set(row_ids)) != 88 or len(set(pivot_columns)) != 88:
        raise AssertionError("pivot rows or columns are not distinct")

    minor = [[0] * 88 for _ in range(88)]
    for minor_row, original_row in enumerate(row_ids):
        equation_position, multiplier_position = divmod(
            original_row, multiplier_count
        )
        source_index = selected_sources[equation_position]
        multiplier = multiplier_monomials[multiplier_position]
        for exponent, coefficient in equations[source_index]:
            shifted = tuple(a + b for a, b in zip(exponent, multiplier))
            column = target_index[shifted]
            minor_column = next(
                (j for j, pivot_column in enumerate(pivot_columns)
                 if pivot_column == column),
                None,
            )
            if minor_column is not None:
                minor[minor_row][minor_column] = coefficient

    determinant = determinant_mod(minor, prime)
    if determinant == 0:
        raise AssertionError("the selected GPU minor vanishes on replay")

    result = {
        "source": str(SOURCE),
        "source_sha256": source_hash,
        "gpu_receipt": str(GPU_RECEIPT),
        "gpu_receipt_sha256": gpu_receipt_hash,
        "branch": "2+1+1",
        "coefficient_field": "Q, reduced to GF(32749)",
        "macaulay_weight": target_weight,
        "minor_shape": [88, 88],
        "selected_original_rows": row_ids,
        "selected_columns": pivot_columns,
        "determinant_mod_prime": determinant,
        "consequence": (
            "The selected rational Macaulay matrix has rank at least 88 in "
            "characteristic zero. This is a rank certificate for the selected "
            "35-equation slice only, not a branch exclusion."
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
