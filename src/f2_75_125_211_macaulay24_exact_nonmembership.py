#!/usr/bin/env python3
"""Build modular-minor witnesses for the weight-24 pure-power nonmemberships."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from f2_75_125_211_macaulay22_23_exact_nonmembership import (
    RANK4,
    SOURCE,
    load_rows,
    replay_rank_four_relations,
    sha256_file,
    monomials_of_weight,
)


HERE = Path(__file__).resolve().parent
GPU_RECEIPT = HERE / "f2_75_125_211_depth14_support23_gpu_macaulay_chunked8_basis4_w24.json"
OUT = HERE / "f2_75_125_211_macaulay24_exact_nonmembership_certificate.json"
PRIME = 251
EXPECTED_BASIS_ROWS = [86, 87, 88, 89]
SUPPORT_POSITIONS = tuple([*range(5), *range(6, 13), *range(14, 25)])


def solve_many_and_det(matrix: np.ndarray, rhs: np.ndarray) -> tuple[int, np.ndarray]:
    """Solve matrix*x=rhs mod PRIME, returning det(matrix) and all solutions."""
    size, right_count = rhs.shape
    if matrix.shape != (size, size):
        raise AssertionError("unexpected pivot minor dimensions")
    work = np.concatenate((matrix.T.copy(), rhs.copy()), axis=1).astype(np.int64, copy=False)
    determinant = 1
    for column in range(size):
        candidates = np.flatnonzero(work[column:, column])
        if len(candidates) == 0:
            raise AssertionError("selected modular pivot minor is singular")
        pivot = column + int(candidates[0])
        if pivot != column:
            work[[column, pivot]] = work[[pivot, column]]
            determinant = (-determinant) % PRIME
        pivot_value = int(work[column, column])
        determinant = (determinant * pivot_value) % PRIME
        inverse = pow(pivot_value, -1, PRIME)
        work[column, column:] = (work[column, column:] * inverse) % PRIME
        if column + 1 < size:
            factors = work[column + 1:, column].copy()
            work[column + 1:, column + 1:] = (
                work[column + 1:, column + 1:]
                - factors[:, None] * work[column, column + 1:][None, :]
            ) % PRIME
            work[column + 1:, column] = 0

    solutions = np.zeros((size, right_count), dtype=np.int64)
    for row in range(size - 1, -1, -1):
        if row + 1 < size:
            tail = work[row, row + 1:size] @ solutions[row + 1:, :]
        else:
            tail = np.zeros(right_count, dtype=np.int64)
        solutions[row, :] = (work[row, size:] - tail) % PRIME
    return determinant, solutions


def main() -> None:
    receipt = json.loads(GPU_RECEIPT.read_text(encoding="utf-8"))
    source_hash = sha256_file(SOURCE)
    rank4 = json.loads(RANK4.read_text(encoding="utf-8"))
    rank4_hash = sha256_file(RANK4)
    if receipt["source_sha256"] != source_hash or rank4["source_sha256"] != source_hash:
        raise AssertionError("GPU receipt or exact span certificate has a different source")
    basis_rows = [int(item["row"]) for item in rank4["selected_source_rows"]]
    if basis_rows != EXPECTED_BASIS_ROWS:
        raise AssertionError("unexpected four-row exact-Q basis")
    if rank4["coefficient_field"] != "Q" or rank4["exact_rank"] != 4:
        raise AssertionError("missing exact-Q rank-four certificate")
    if not rank4.get("all_source_rows_exactly_replayed_over_Q"):
        raise AssertionError("the 35 source equations are not certified over Q")
    if receipt["selected_rows_zero_based"] != basis_rows:
        raise AssertionError("GPU matrix did not use the certified four basis rows")
    if receipt["prime"] != PRIME or receipt["target_weight"] != 24:
        raise AssertionError("unexpected modular prime or Macaulay degree")

    data, weights, source_rows = load_rows()
    replay_rank_four_relations(source_rows, rank4)
    selected_equations = [source_rows[index] for index in basis_rows]
    support_weights = tuple(weights[index] for index in SUPPORT_POSITIONS)
    degree = 24
    target_monomials = monomials_of_weight(support_weights, degree)
    target_index = {monomial: index for index, monomial in enumerate(target_monomials)}
    multiplier_monomials = monomials_of_weight(support_weights, degree - 14)
    multiplier_count = len(multiplier_monomials)
    ncols = len(target_monomials)
    nrows = len(selected_equations) * multiplier_count
    if (nrows, ncols) != tuple(receipt["macaulay_shape"]):
        raise AssertionError("weight-24 Macaulay monomial census differs from the GPU receipt")
    if multiplier_count != 396 or receipt["rank"] != nrows:
        raise AssertionError("weight-24 modular rank does not meet the exact four-row bound")

    target_names = list(receipt["pure_power_membership"])
    if not target_names or any(receipt["pure_power_membership"][name]["in_row_span_mod_prime"]
                               for name in target_names):
        raise AssertionError("a target was reported in the modular row span")
    target_columns = []
    for name in target_names:
        variable, power_text = name.split("^")
        variable_position = int(variable[1:])
        support_position = SUPPORT_POSITIONS.index(variable_position)
        exponent = [0] * len(SUPPORT_POSITIONS)
        exponent[support_position] = int(power_text)
        target_columns.append(target_index[tuple(exponent)])

    selected_rows = [int(value) for value in receipt["selected_original_rows"]]
    pivot_columns = [int(value) for value in receipt["pivot_columns"]]
    if len(selected_rows) != nrows or len(pivot_columns) != nrows:
        raise AssertionError("the GPU receipt lacks a full pivot basis")
    pivot_position = {column: index for index, column in enumerate(pivot_columns)}
    if len(pivot_position) != nrows:
        raise AssertionError("pivot columns are not distinct")
    pivot_map = np.full(ncols, -1, dtype=np.int32)
    pivot_map[np.asarray(pivot_columns, dtype=np.int32)] = np.arange(nrows, dtype=np.int32)

    matrix = np.zeros((nrows, nrows), dtype=np.int64)
    sparse_rows = []
    for row_id in selected_rows:
        equation_position, multiplier_position = divmod(row_id, multiplier_count)
        equation = selected_equations[equation_position]
        multiplier = multiplier_monomials[multiplier_position]
        columns = np.empty(len(equation), dtype=np.int32)
        coefficients = np.empty(len(equation), dtype=np.int64)
        for term_position, (exponent, coefficient) in enumerate(equation.items()):
            projected = tuple(exponent[index] for index in SUPPORT_POSITIONS)
            if any(exponent[index] for index in range(len(exponent))
                   if index not in SUPPORT_POSITIONS):
                raise AssertionError("basis equation uses an unexpected variable")
            shifted = tuple(a + b for a, b in zip(projected, multiplier))
            columns[term_position] = target_index[shifted]
            denominator = coefficient.denominator % PRIME
            if denominator == 0:
                raise AssertionError("basis coefficient denominator is divisible by 251")
            coefficients[term_position] = (
                (coefficient.numerator % PRIME) * pow(denominator, -1, PRIME) % PRIME
            )
        mapped = pivot_map[columns]
        mask = mapped >= 0
        matrix[len(sparse_rows), mapped[mask]] = coefficients[mask]
        sparse_rows.append((columns, coefficients))

    rhs = np.zeros((nrows, len(target_names)), dtype=np.int64)
    for target_position, target_column in enumerate(target_columns):
        pivot_position_for_target = pivot_position.get(target_column)
        if pivot_position_for_target is not None:
            rhs[pivot_position_for_target, target_position] = 1
    pivot_determinant, combinations = solve_many_and_det(matrix, rhs)
    if pivot_determinant == 0:
        raise AssertionError("the selected pivot minor is zero")

    residuals = np.zeros((len(target_names), ncols), dtype=np.int64)
    for target_position, target_column in enumerate(target_columns):
        residuals[target_position, target_column] = 1
    for row_position, (columns, coefficients) in enumerate(sparse_rows):
        for target_position in range(len(target_names)):
            scalar = int(combinations[row_position, target_position])
            if scalar:
                residuals[target_position, columns] = (
                    residuals[target_position, columns] - scalar * coefficients
                ) % PRIME

    witnesses = {}
    for target_position, name in enumerate(target_names):
        support = np.flatnonzero(residuals[target_position])
        expected_support = int(receipt["pure_power_membership"][name]["residual_term_count"])
        if len(support) != expected_support or len(support) == 0:
            raise AssertionError(f"modular residual mismatch for {name}")
        extra_column = int(support[0])
        residual_value = int(residuals[target_position, extra_column])
        bordered_determinant = pivot_determinant * residual_value % PRIME
        if bordered_determinant == 0:
            raise AssertionError(f"zero bordered minor for {name}")
        witnesses[name] = {
            "pivot_minor_determinant_mod_prime": int(pivot_determinant),
            "extra_column_zero_based": extra_column,
            "target_residual_mod_prime": residual_value,
            "bordered_augmented_minor_determinant_mod_prime": int(bordered_determinant),
        }

    output = {
        "branch": "2+1+1",
        "scope": (
            "Exact-Q non-membership of all 17 weight-24 pure powers tested by the GPU receipt "
            "in the ideal generated by the four basis equations, equivalently the 35 support-23 "
            "equations. This does not exclude their zero set or the full 2+1+1 branch."
        ),
        "source": SOURCE.name,
        "source_sha256": source_hash,
        "exact_rank4_certificate": RANK4.name,
        "exact_rank4_certificate_sha256": rank4_hash,
        "gpu_receipt": GPU_RECEIPT.name,
        "gpu_receipt_sha256": sha256_file(GPU_RECEIPT),
        "gpu_receipt_basis_certificate_sha256_at_run": receipt["basis_certificate_sha256"],
        "current_exact_basis_rows": basis_rows,
        "exact_Q_relations_replayed": len(rank4["row_relations"]),
        "coefficient_field_for_conclusion": "Q",
        "prime": PRIME,
        "weighted_degree": degree,
        "macaulay_shape": [nrows, ncols],
        "multiplier_monomial_count": multiplier_count,
        "exact_Q_rank_upper_bound": 4 * multiplier_count,
        "pivot_minor_determinant_mod_prime": int(pivot_determinant),
        "exact_Q_rank": nrows,
        "targets": witnesses,
        "whole_branch_exclusion_percent": 0,
        "argument": (
            "The exact-Q relations show that the 35 generators span the four basis equations. "
            "Thus the weight-24 Macaulay rank over Q is at most 4*396=1584. The selected "
            "1584-by-1584 pivot minor and each bordered target minor are nonzero modulo 251; "
            "all denominators in those matrices are invertible modulo 251. These determinants "
            "prove the exact-Q rank is 1584 and each target raises the rank to 1585, hence "
            "each target is outside the ideal's weight-24 component."
        ),
    }
    OUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "scope": output["scope"],
        "exact_Q_rank": output["exact_Q_rank"],
        "pivot_minor_determinant_mod_prime": output["pivot_minor_determinant_mod_prime"],
        "target_count": len(witnesses),
        "targets": {
            name: item["bordered_augmented_minor_determinant_mod_prime"]
            for name, item in witnesses.items()
        },
        "receipt": str(OUT),
        "receipt_sha256": sha256_file(OUT),
    }, indent=2))


if __name__ == "__main__":
    main()
