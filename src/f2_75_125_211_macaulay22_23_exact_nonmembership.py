#!/usr/bin/env python3
"""Replay exact-Q consequences of the weight-22/23 support-23 Macaulay runs.

The GPU ranks are modular, but an exact-Q rank-four certificate supplies a
matching upper bound on the rational Macaulay ranks.  Together these imply
exact rational non-membership of the listed pure powers in the ideal generated
by the 35 support-23 compatibility polynomials.  This is a scoped ideal fact,
not a full F2 branch exclusion.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

import ijson
import numpy as np


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_depth14_global_compatibility_reduction_20261021.json"
RANK4 = HERE / "f2_75_125_211_support23_exact_rank4_certificate.json"
RECEIPTS = {
    22: HERE / "f2_75_125_211_depth14_support23_gpu_macaulay_chunked16_w22.json",
    23: HERE / "f2_75_125_211_depth14_support23_gpu_macaulay_chunked16_w23.json",
}
OUT = HERE / "f2_75_125_211_macaulay22_23_exact_nonmembership_certificate.json"
PRIME = 32749
ROWS = set(range(86, 121))
EXPECTED_TARGETS = {22: ("p0^22", "p1^11"), 23: ("p0^23",)}
EXPECTED_RANKS = {22: 544, 23: 932}
EXPECTED_MULTIPLIERS = {22: 136, 23: 233}
SUPPORT_POSITIONS = tuple([*range(5), *range(6, 13), *range(14, 25)])


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_rows():
    with SOURCE.open("rb") as stream:
        data = json.load(stream)
    weights = tuple(int(weight) for weight in data["old_parameter_weights"])
    if len(weights) != 42:
        raise AssertionError("expected 42 old parameters")
    selected = {}
    with SOURCE.open("rb") as stream:
        for row in ijson.items(stream, "global_compatibility_equations.item"):
            index = int(row["constraint_row_zero_based"])
            if index not in ROWS:
                continue
            polynomial = {}
            for sparse, coefficient in row["terms"]:
                exponent = [0] * len(weights)
                for variable, power in sparse:
                    exponent[int(variable)] = int(power)
                exponent = tuple(exponent)
                if sum(w * e for w, e in zip(weights, exponent)) != 14:
                    raise AssertionError(f"wrong weighted degree in source row {index}")
                q = Fraction(coefficient)
                if q.denominator % PRIME == 0:
                    raise AssertionError(f"bad reduction denominator in row {index}")
                polynomial[exponent] = q
            selected[index] = polynomial
    if set(selected) != ROWS:
        raise AssertionError("did not read all 35 support-23 rows")
    return data, weights, selected


def replay_rank_four_relations(rows, certificate):
    basis_rows = [int(item["row"]) for item in certificate["selected_source_rows"]]
    if basis_rows != [86, 87, 88, 89]:
        raise AssertionError("unexpected exact basis rows")
    relation_map = {int(item["row"]): item for item in certificate["row_relations"]}
    if set(relation_map) != ROWS:
        raise AssertionError("rank-four certificate does not list every source row")
    for row_index in sorted(ROWS):
        relation = relation_map[row_index]
        coefficients = [Fraction(c) for c in relation["coefficients"]]
        if len(coefficients) != 4:
            raise AssertionError(f"bad relation width at row {row_index}")
        replayed = {}
        for basis_index, coefficient in zip(basis_rows, coefficients):
            if not coefficient:
                continue
            for exponent, basis_coefficient in rows[basis_index].items():
                replayed[exponent] = replayed.get(exponent, Fraction(0)) + coefficient * basis_coefficient
        replayed = {exponent: coefficient for exponent, coefficient in replayed.items() if coefficient}
        if replayed != rows[row_index]:
            raise AssertionError(f"exact Q relation failed for source row {row_index}")
    return basis_rows


def monomials_of_weight(weights: tuple[int, ...], target: int):
    from functools import lru_cache

    @lru_cache(None)
    def rec(index: int, remaining: int):
        if index == len(weights):
            return ((),) if remaining == 0 else ()
        out = []
        for exponent in range(remaining // weights[index] + 1):
            out.extend((exponent,) + tail for tail in rec(
                index + 1, remaining - exponent * weights[index]
            ))
        return tuple(out)

    return rec(0, target)


def modular_solve_and_determinant(matrix: np.ndarray, rhs: np.ndarray):
    """Return det(matrix), x with matrix*x=rhs, using Gauss--Jordan mod p."""
    size = int(matrix.shape[0])
    work = np.concatenate((matrix.astype(np.int64).copy(), rhs.reshape(-1, 1)), axis=1)
    determinant = 1
    for column in range(size):
        candidates = np.flatnonzero(work[column:, column])
        if len(candidates) == 0:
            raise AssertionError("selected pivot minor is singular modulo the prime")
        pivot = column + int(candidates[0])
        if pivot != column:
            work[[column, pivot]] = work[[pivot, column]]
            determinant = (-determinant) % PRIME
        pivot_value = int(work[column, column])
        determinant = (determinant * pivot_value) % PRIME
        inverse = pow(pivot_value, -1, PRIME)
        work[column, column:] = (work[column, column:] * inverse) % PRIME
        factors = work[:, column].copy()
        factors[column] = 0
        active = np.flatnonzero(factors)
        for row in active:
            work[row, column:] = (
                work[row, column:] - factors[row] * work[column, column:]
            ) % PRIME
    return determinant, work[:, -1].copy()


def modular_minor_witnesses(degree, receipt, weights, source_rows):
    support_weights = tuple(weights[index] for index in SUPPORT_POSITIONS)
    target_monomials = monomials_of_weight(support_weights, degree)
    target_index = {monomial: index for index, monomial in enumerate(target_monomials)}
    multiplier_monomials = monomials_of_weight(support_weights, degree - 14)
    multiplier_count = len(multiplier_monomials)
    selected_rows = [int(value) for value in receipt["selected_original_rows"]]
    pivot_columns = [int(value) for value in receipt["pivot_columns"]]
    rank = len(selected_rows)
    if len(pivot_columns) != rank:
        raise AssertionError("row and column pivot counts disagree")
    pivot_column_index = {column: index for index, column in enumerate(pivot_columns)}
    if len(pivot_column_index) != rank:
        raise AssertionError("pivot columns are not distinct")

    selected_equation_rows = []
    selected_multipliers = []
    for row_id in selected_rows:
        equation_position, multiplier_position = divmod(row_id, multiplier_count)
        equation_index = 86 + equation_position
        selected_equation_rows.append(source_rows[equation_index])
        selected_multipliers.append(multiplier_monomials[multiplier_position])

    matrix = np.zeros((rank, rank), dtype=np.int64)
    row_term_caches = []
    for row_position, (equation, multiplier) in enumerate(
        zip(selected_equation_rows, selected_multipliers)
    ):
        shifted_terms = []
        for exponent, coefficient in equation.items():
            projected = tuple(exponent[index] for index in SUPPORT_POSITIONS)
            if any(exponent[index] for index in range(len(exponent))
                   if index not in SUPPORT_POSITIONS):
                raise AssertionError("support-23 source row uses an unexpected variable")
            shifted = tuple(a + b for a, b in zip(projected, multiplier))
            column = target_index[shifted]
            denominator = coefficient.denominator % PRIME
            if denominator == 0:
                raise AssertionError("selected matrix coefficient has bad denominator")
            value = (coefficient.numerator % PRIME) * pow(denominator, -1, PRIME) % PRIME
            shifted_terms.append((column, value))
            pivot_position = pivot_column_index.get(column)
            if pivot_position is not None:
                matrix[row_position, pivot_position] = value
        row_term_caches.append(shifted_terms)

    witnesses = {}
    for target in EXPECTED_TARGETS[degree]:
        variable, exponent_text = target.split("^")
        variable_index = int(variable[1:])
        exponent_value = int(exponent_text)
        projected_position = SUPPORT_POSITIONS.index(variable_index)
        target_exponent = [0] * len(SUPPORT_POSITIONS)
        target_exponent[projected_position] = exponent_value
        target_column = target_index[tuple(target_exponent)]

        rhs = np.fromiter(
            (1 if column == target_column else 0 for column in pivot_columns),
            dtype=np.int64,
            count=rank,
        )
        determinant, combination = modular_solve_and_determinant(matrix.T, rhs)
        residual = np.zeros(len(target_monomials), dtype=np.int64)
        residual[target_column] = 1
        for scalar, terms in zip(combination, row_term_caches):
            scalar = int(scalar)
            if scalar:
                for column, coefficient in terms:
                    residual[column] = (residual[column] - scalar * coefficient) % PRIME
        support = np.flatnonzero(residual)
        if len(support) != int(receipt["pure_power_membership"][target]["residual_term_count"]):
            raise AssertionError(f"residual support mismatch for {target}, weight {degree}")
        if len(support) == 0:
            raise AssertionError(f"target {target} lies in the modular row span")
        extra_column = int(support[0])
        residual_value = int(residual[extra_column])
        bordered_determinant = (determinant * residual_value) % PRIME
        if bordered_determinant == 0:
            raise AssertionError("bordered augmented minor is zero modulo the prime")
        witnesses[target] = {
            "target_column_zero_based": target_column,
            "pivot_minor_determinant_mod_prime": int(determinant),
            "extra_column_zero_based": extra_column,
            "target_residual_at_extra_column_mod_prime": residual_value,
            "bordered_augmented_minor_determinant_mod_prime": int(bordered_determinant),
            "independent_original_macaulay_rows": selected_rows,
            "independent_pivot_columns": pivot_columns,
        }
    return witnesses


def main() -> None:
    source_hash = sha256_file(SOURCE)
    rank4 = json.loads(RANK4.read_text(encoding="utf-8"))
    if rank4["source_sha256"] != source_hash:
        raise AssertionError("rank-four certificate is for a different source")
    if rank4["coefficient_field"] != "Q" or rank4["exact_rank"] != 4:
        raise AssertionError("missing exact-Q rank-four assertion")
    if not rank4.get("all_source_rows_exactly_replayed_over_Q"):
        raise AssertionError("rank-four certificate lacks exact-Q relation replay")

    data, weights, rows = load_rows()
    basis_rows = replay_rank_four_relations(rows, rank4)
    source_row_term_count = sum(len(polynomial) for polynomial in rows.values())
    if source_row_term_count != rank4["source_term_count"]:
        raise AssertionError("source term count does not match the exact rank-four certificate")

    conclusions = []
    for degree, receipt_path in RECEIPTS.items():
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if receipt["source_sha256"] != source_hash:
            raise AssertionError(f"weight-{degree} receipt is for a different source")
        if receipt["prime"] != PRIME or receipt["target_weight"] != degree:
            raise AssertionError(f"unexpected field or degree in weight-{degree} receipt")
        if set(receipt["selected_rows_zero_based"]) != ROWS:
            raise AssertionError(f"weight-{degree} receipt has a different row block")
        nrows, ncols = receipt["macaulay_shape"]
        multiplier_count = nrows // len(ROWS)
        if nrows % len(ROWS) or multiplier_count != EXPECTED_MULTIPLIERS[degree]:
            raise AssertionError(f"bad weight-{degree} Macaulay row census")
        modular_rank = int(receipt["rank"])
        exact_upper_bound = len(basis_rows) * multiplier_count
        if modular_rank != EXPECTED_RANKS[degree] or modular_rank != exact_upper_bound:
            raise AssertionError(f"weight-{degree} rank does not meet the exact-Q upper bound")

        targets = {}
        witnesses = modular_minor_witnesses(degree, receipt, weights, rows)
        for target in EXPECTED_TARGETS[degree]:
            status = receipt["pure_power_membership"][target]
            if status["in_row_span_mod_prime"] or status["residual_term_count"] < 1:
                raise AssertionError(f"weight-{degree} target {target} is not modularly excluded")
            targets[target] = {
                "modular_row_nonmembership": True,
                "augmented_modular_rank_lower_bound": modular_rank + 1,
                "exact_Q_row_nonmembership": True,
                "modular_minor_witness": witnesses[target],
                "argument": (
                    "Every source generator is an exact Q-linear combination of four basis rows, "
                    "so the rational Macaulay rank is at most 4 times the multiplier count. "
                    "The modular rank reaches that upper bound, while the target raises the "
                    "modular augmented rank by one. All source coefficient denominators are "
                    "invertible modulo the stated prime, so these nonzero modular minors also "
                    "prove the corresponding rational ranks."
                ),
            }
        conclusions.append({
            "weighted_degree": degree,
            "matrix_shape": [nrows, ncols],
            "multiplier_monomial_count": multiplier_count,
            "rank_upper_bound_over_Q_from_four_generators": exact_upper_bound,
            "rank_lower_bound_over_GF_from_nonzero_pivot_minor": modular_rank,
            "exact_Q_rank": exact_upper_bound,
            "prime": PRIME,
            "targets": targets,
            "modular_matrix_sha256": receipt["matrix_sha256"],
            "gpu": receipt["gpu"],
        })

    result = {
        "branch": "2+1+1",
        "scope": (
            "Exact-Q non-membership of p0^22, p1^11, and p0^23 in the ideal generated by "
            "the 35 support-23 depth-14 compatibility equations. This does not exclude the "
            "support-23 variety or the full 2+1+1 branch."
        ),
        "coefficient_field_for_conclusion": "Q",
        "source": SOURCE.name,
        "source_sha256": source_hash,
        "rank_four_certificate": RANK4.name,
        "rank_four_certificate_sha256": sha256_file(RANK4),
        "exact_Q_relations_replayed": len(ROWS),
        "exact_Q_basis_rows": basis_rows,
        "source_row_term_count": source_row_term_count,
        "all_source_denominators_invertible_mod_prime": True,
        "prime": PRIME,
        "weight_results": conclusions,
        "whole_branch_exclusion_percent": 0,
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    summary = {
        "scope": result["scope"],
        "exact_Q_relations_replayed": result["exact_Q_relations_replayed"],
        "weight_results": [
            {
                "weighted_degree": item["weighted_degree"],
                "exact_Q_rank": item["exact_Q_rank"],
                "targets": {
                    name: {
                        "pivot_minor_determinant_mod_prime": target["modular_minor_witness"][
                            "pivot_minor_determinant_mod_prime"
                        ],
                        "bordered_augmented_minor_determinant_mod_prime": target[
                            "modular_minor_witness"
                        ]["bordered_augmented_minor_determinant_mod_prime"],
                        "exact_Q_row_nonmembership": target["exact_Q_row_nonmembership"],
                    }
                    for name, target in item["targets"].items()
                },
            }
            for item in result["weight_results"]
        ],
        "whole_branch_exclusion_percent": result["whole_branch_exclusion_percent"],
    }
    print(json.dumps(summary, indent=2))
    print(f"receipt={OUT} sha256={sha256_file(OUT)}")


if __name__ == "__main__":
    main()
