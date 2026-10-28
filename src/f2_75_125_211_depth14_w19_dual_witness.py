"""Try to lift the weight-19 modular nonmembership screen to Q exactly.

The script reconstructs the 770 by 23,365 Macaulay matrix only on the 88
recorded pivot columns and the p0^19 target column.  It solves for a relation
expressing the target column through the pivot columns and checks that relation
in every row over Q.  If it succeeds, the resulting kernel vector has target
coordinate 1 and is an exact dual witness that p0^19 is not in the rational
row span of this selected 35-equation Macaulay slice.

This is only a finite-degree linear certificate, not an exclusion of the
2+1+1 branch.
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
RECEIPT = Path(
    "research/f2_75_125_211_depth14_support23_gpu_macaulay_w19.json"
)
OUTPUT = Path(
    "research/f2_75_125_211_depth14_w19_dual_witness.json"
)
EXPECTED_SOURCE_SHA256 = (
    "372ba14b06d04297b25e098ccbdfc03c3d03dffd6ae70f9060655be608a0689b"
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


def solve_square_fraction(
    matrix: list[list[Fraction]], rhs: list[Fraction]
) -> list[Fraction]:
    n = len(matrix)
    if len(rhs) != n or any(len(row) != n for row in matrix):
        raise ValueError("expected a square system")
    augmented = [row[:] + [value] for row, value in zip(matrix, rhs)]
    for col in range(n):
        pivot = next((row for row in range(col, n) if augmented[row][col]), None)
        if pivot is None:
            raise ArithmeticError(f"singular rational basis at column {col}")
        if pivot != col:
            augmented[col], augmented[pivot] = augmented[pivot], augmented[col]
        pivot_value = augmented[col][col]
        for j in range(col, n + 1):
            augmented[col][j] /= pivot_value
        for row in range(n):
            if row == col:
                continue
            factor = augmented[row][col]
            if factor:
                for j in range(col, n + 1):
                    augmented[row][j] -= factor * augmented[col][j]
    return [augmented[row][n] for row in range(n)]


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite {OUTPUT}")
    source_hash = sha256(SOURCE)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise AssertionError("global compatibility source hash changed")
    receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
    if receipt["source_sha256"] != source_hash:
        raise AssertionError("GPU receipt has a different source hash")
    if receipt["branch"] != "2+1+1" or receipt["target_weight"] != 19:
        raise AssertionError("unexpected branch or target weight")

    prime = int(receipt["prime"])
    variable_names = receipt["support_variables"]
    variable_ids = tuple(int(name[1:]) for name in variable_names)
    weights = tuple(map(int, receipt["weights"]))
    target_weight = int(receipt["target_weight"])
    target_monomials = monomials_of_weight(weights, target_weight)
    target_index = {monomial: index for index, monomial in enumerate(target_monomials)}
    target_exponent = [0] * len(variable_ids)
    target_exponent[variable_ids.index(0)] = target_weight
    target_col = target_index[tuple(target_exponent)]

    selected_sources = tuple(map(int, receipt["selected_rows_zero_based"]))
    selected_set = set(selected_sources)
    source_positions = {variable: pos for pos, variable in enumerate(variable_ids)}
    equations: dict[int, list[tuple[tuple[int, ...], Fraction]]] = {}
    with SOURCE.open("rb") as stream:
        for equation in ijson.items(stream, "global_compatibility_equations.item"):
            source_row = int(equation["constraint_row_zero_based"])
            if source_row not in selected_set:
                continue
            if equation["support_variables"] != variable_names:
                raise AssertionError("selected support changed")
            terms = []
            for sparse, coefficient in equation["terms"]:
                exponent = [0] * len(variable_ids)
                for variable, power in sparse:
                    exponent[source_positions[int(variable)]] = int(power)
                if sum(w * e for w, e in zip(weights, exponent)) != 14:
                    raise AssertionError("selected equation is not weight 14")
                terms.append((tuple(exponent), Fraction(coefficient)))
            equations[source_row] = terms
    if tuple(sorted(equations)) != selected_sources:
        raise AssertionError("selected equation census does not match receipt")

    multiplier_monomials = monomials_of_weight(weights, target_weight - 14)
    pivot_columns = list(map(int, receipt["pivot_columns"]))
    pivot_rows = list(map(int, receipt["selected_original_rows"]))
    if len(pivot_columns) != 88 or len(pivot_rows) != 88:
        raise AssertionError("expected the recorded 88 pivots")
    if target_col in pivot_columns:
        raise AssertionError("target monomial unexpectedly is a pivot column")
    if len(set(pivot_columns)) != 88 or len(set(pivot_rows)) != 88:
        raise AssertionError("pivot rows or columns are repeated")

    needed_columns = set(pivot_columns)
    needed_columns.add(target_col)
    matrix_rows: list[dict[int, Fraction]] = []
    for source_row in selected_sources:
        terms = equations[source_row]
        for multiplier in multiplier_monomials:
            restricted: dict[int, Fraction] = {}
            for exponent, coefficient in terms:
                shifted = tuple(a + b for a, b in zip(exponent, multiplier))
                column = target_index[shifted]
                if column in needed_columns:
                    restricted[column] = restricted.get(column, Fraction(0)) + coefficient
            matrix_rows.append(restricted)
    if len(matrix_rows) != receipt["macaulay_shape"][0]:
        raise AssertionError("Macaulay row count changed")
    if len(pivot_rows) and max(pivot_rows) >= len(matrix_rows):
        raise AssertionError("pivot row is out of range")

    basis_matrix = [
        [matrix_rows[row].get(column, Fraction(0)) for column in pivot_columns]
        for row in pivot_rows
    ]
    rhs = [-matrix_rows[row].get(target_col, Fraction(0)) for row in pivot_rows]
    coefficients = solve_square_fraction(basis_matrix, rhs)

    failures = []
    for row_index, row in enumerate(matrix_rows):
        residual = row.get(target_col, Fraction(0))
        residual += sum(
            coefficient * row.get(column, Fraction(0))
            for coefficient, column in zip(coefficients, pivot_columns)
        )
        if residual:
            failures.append((row_index, str(residual)))
            if len(failures) >= 10:
                break
    if failures:
        raise ArithmeticError(
            "the modular column relation does not lift to Q; first failures: "
            + repr(failures)
        )

    witness = {
        "source": str(SOURCE),
        "source_sha256": source_hash,
        "gpu_receipt": str(RECEIPT),
        "gpu_receipt_sha256": sha256(RECEIPT),
        "branch": "2+1+1",
        "field": "Q",
        "macaulay_weight": target_weight,
        "selected_equation_rows_zero_based": list(selected_sources),
        "selected_equation_count": len(selected_sources),
        "macaulay_shape": receipt["macaulay_shape"],
        "target_monomial": "p0^19",
        "target_column_zero_based": target_col,
        "kernel_vector_support_columns": [*pivot_columns, target_col],
        "kernel_vector_coefficients": [
            {"column": col, "coefficient": str(value)}
            for col, value in zip(pivot_columns, coefficients)
        ]
        + [{"column": target_col, "coefficient": "1"}],
        "verified_zero_residual_rows": len(matrix_rows),
        "consequence": (
            "An exact rational Macaulay-kernel vector has p0^19 coordinate 1. "
            "Therefore p0^19 is not in the Q-row-span of this selected "
            "35-equation, weight-19 Macaulay slice. This is not an ideal or "
            "branch exclusion."
        ),
    }
    OUTPUT.write_text(json.dumps(witness, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "EXACT_DUAL_WITNESS_VERIFIED",
        "target_monomial": witness["target_monomial"],
        "target_column_zero_based": target_col,
        "macaulay_shape": witness["macaulay_shape"],
        "exactly_verified_rows": len(matrix_rows),
        "witness_support_size": len(pivot_columns) + 1,
        "output": str(OUTPUT),
        "output_sha256": sha256(OUTPUT),
    }, indent=2))


if __name__ == "__main__":
    main()
