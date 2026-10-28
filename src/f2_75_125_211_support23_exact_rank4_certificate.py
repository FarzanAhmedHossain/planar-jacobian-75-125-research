#!/usr/bin/env python3
"""Certify the exact four-dimensional span of the 35 support-23 equations.

The weight-24 GPU Macaulay job was built from 35 weighted-degree-14
compatibility equations.  This checker proves over Q that those equations
span only four polynomials, records an exact expression of every row in four
selected source rows, and checks that the same reduction is valid modulo 251.
That permits a much smaller but row-space-equivalent modular Macaulay replay.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

import ijson


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_depth14_global_compatibility_reduction_20261021.json"
OUT = HERE / "f2_75_125_211_support23_exact_rank4_certificate.json"
ROWS = set(range(86, 121))
VARIABLES = [*range(5), *range(6, 13), *range(14, 25)]
POSITION = {i: j for j, i in enumerate(VARIABLES)}
WEIGHT_BY_VARIABLE = {
    0: 1, 1: 2,
    2: 3, 3: 3, 4: 3,
    **{i: 4 for i in range(6, 13)},
    **{i: 5 for i in range(14, 20)},
    **{i: 6 for i in range(20, 25)},
}
PRIME = 251


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def mod_fraction(q: Fraction) -> int:
    denominator = q.denominator % PRIME
    if denominator == 0:
        raise ValueError("a relation denominator is not invertible modulo 251")
    return (q.numerator % PRIME) * pow(denominator, -1, PRIME) % PRIME


def solve_square(matrix: list[list[Fraction]], rhs: list[Fraction]) -> list[Fraction]:
    n = len(rhs)
    augmented = [list(matrix[i]) + [rhs[i]] for i in range(n)]
    for col in range(n):
        pivot = next((r for r in range(col, n) if augmented[r][col]), None)
        if pivot is None:
            raise AssertionError("four pivot monomials did not give an invertible matrix")
        augmented[col], augmented[pivot] = augmented[pivot], augmented[col]
        scale = augmented[col][col]
        augmented[col] = [x / scale for x in augmented[col]]
        for r in range(n):
            if r == col or not augmented[r][col]:
                continue
            scale = augmented[r][col]
            augmented[r] = [x - scale * y for x, y in zip(augmented[r], augmented[col])]
    return [augmented[i][-1] for i in range(n)]


def determinant_mod(matrix: list[list[Fraction]]) -> int:
    a = [[mod_fraction(x) for x in row] for row in matrix]
    n = len(a)
    determinant = 1
    for col in range(n):
        pivot = next((r for r in range(col, n) if a[r][col]), None)
        if pivot is None:
            return 0
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            determinant = -determinant
        value = a[col][col]
        determinant = determinant * value % PRIME
        inverse = pow(value, -1, PRIME)
        for r in range(col + 1, n):
            if not a[r][col]:
                continue
            factor = a[r][col] * inverse % PRIME
            for c in range(col + 1, n):
                a[r][c] = (a[r][c] - factor * a[col][c]) % PRIME
            a[r][col] = 0
    return determinant % PRIME


def mod_poly(poly: dict[tuple[int, ...], Fraction]) -> dict[tuple[int, ...], int]:
    result = {}
    for monomial, coefficient in poly.items():
        value = mod_fraction(coefficient)
        if value:
            result[monomial] = value
    return result


def main() -> None:
    source_rows: list[tuple[int, str, dict[tuple[int, ...], Fraction]]] = []
    with SOURCE.open("rb") as stream:
        for row in ijson.items(stream, "global_compatibility_equations.item"):
            if row["constraint_row_zero_based"] not in ROWS:
                continue
            if row["support_variables"] != [f"p{i}" for i in VARIABLES]:
                raise AssertionError(f"unexpected support at {row['constraint_tag']}")
            polynomial: dict[tuple[int, ...], Fraction] = {}
            for sparse, coefficient in row["terms"]:
                exponent = [0] * len(VARIABLES)
                for variable, power in sparse:
                    exponent[POSITION[variable]] = power
                monomial = tuple(exponent)
                value = Fraction(coefficient)
                if monomial in polynomial:
                    value += polynomial[monomial]
                if value:
                    polynomial[monomial] = value
                else:
                    polynomial.pop(monomial, None)
            source_rows.append((row["constraint_row_zero_based"], row["constraint_tag"], polynomial))

    if len(source_rows) != 35:
        raise AssertionError(f"expected 35 source equations, got {len(source_rows)}")

    axis_report = {}
    for variable in VARIABLES:
        weight = WEIGHT_BY_VARIABLE[variable]
        exponent = 14 // weight if 14 % weight == 0 else None
        nonzero_rows = []
        if exponent is not None:
            position = POSITION[variable]
            pure_monomial = tuple(exponent if j == position else 0 for j in range(len(VARIABLES)))
            nonzero_rows = [index for index, _tag, polynomial in source_rows
                            if polynomial.get(pure_monomial, Fraction(0))]
        axis_report[f"p{variable}"] = {
            "weight": weight,
            "tested_pure_power_exponent": exponent,
            "rows_with_nonzero_axis_restriction": nonzero_rows,
            "coordinate_axis_survives_this_35_equation_subideal": not nonzero_rows,
        }

    # Incremental sparse row reduction.  Every stored row is normalized and
    # supported strictly below its pivot in lexicographic monomial order.
    echelon: dict[tuple[int, ...], dict[tuple[int, ...], Fraction]] = {}
    selected: list[tuple[int, str, dict[tuple[int, ...], Fraction]]] = []
    pivots: list[tuple[int, ...]] = []
    for row in source_rows:
        reduced = row[2].copy()
        while reduced:
            pivot = max(reduced)
            if pivot not in echelon:
                scale = reduced[pivot]
                echelon[pivot] = {m: c / scale for m, c in reduced.items() if c / scale}
                selected.append(row)
                pivots.append(pivot)
                break
            factor = reduced[pivot]
            for monomial, coefficient in echelon[pivot].items():
                value = reduced.get(monomial, Fraction(0)) - factor * coefficient
                if value:
                    reduced[monomial] = value
                else:
                    reduced.pop(monomial, None)

    if len(selected) != 4 or len(pivots) != 4:
        raise AssertionError(f"expected rank four, found {len(selected)}")

    # Solve for each source row using the four certified pivot coordinates,
    # then replay the proposed identity on every monomial exactly over Q.
    coefficient_matrix = [
        [polynomial.get(pivot, Fraction(0)) for pivot in pivots]
        for _index, _tag, polynomial in selected
    ]
    coefficient_matrix_transpose = [
        [coefficient_matrix[r][c] for r in range(4)] for c in range(4)
    ]
    relations = []
    max_denominator_bits = 0
    all_relation_denominators_invertible = True
    for index, tag, polynomial in source_rows:
        rhs = [polynomial.get(pivot, Fraction(0)) for pivot in pivots]
        coefficients = solve_square(coefficient_matrix_transpose, rhs)
        replay: dict[tuple[int, ...], Fraction] = {}
        for scalar, (_basis_index, _basis_tag, basis_poly) in zip(coefficients, selected):
            if not scalar:
                continue
            for monomial, value in basis_poly.items():
                replay[monomial] = replay.get(monomial, Fraction(0)) + scalar * value
                if not replay[monomial]:
                    del replay[monomial]
        if replay != polynomial:
            raise AssertionError(f"exact span identity failed for source row {index}")
        for scalar in coefficients:
            max_denominator_bits = max(max_denominator_bits, scalar.denominator.bit_length())
            if scalar.denominator % PRIME == 0:
                all_relation_denominators_invertible = False
        relations.append({"row": index, "tag": tag, "coefficients": [str(c) for c in coefficients]})

    if not all_relation_denominators_invertible:
        raise AssertionError("the exact span identities do not all reduce modulo 251")

    modular_pivot_matrix_determinant = determinant_mod(coefficient_matrix)
    if modular_pivot_matrix_determinant == 0:
        raise AssertionError("the four selected source rows lose independence modulo 251")

    # Explicitly replay every relation modulo 251 as a second check.
    basis_mod = [mod_poly(row[2]) for row in selected]
    for source_index, row in enumerate(source_rows):
        coefficients = [Fraction(value) for value in relations[source_index]["coefficients"]]
        replay_mod: dict[tuple[int, ...], int] = {}
        for scalar, basis_poly in zip(coefficients, basis_mod):
            factor = mod_fraction(scalar)
            if not factor:
                continue
            for monomial, value in basis_poly.items():
                total = (replay_mod.get(monomial, 0) + factor * value) % PRIME
                if total:
                    replay_mod[monomial] = total
                else:
                    replay_mod.pop(monomial, None)
        if replay_mod != mod_poly(row[2]):
            raise AssertionError(f"mod-251 span identity failed for source row {row[0]}")

    result = {
        "source": SOURCE.name,
        "source_sha256": sha256_file(SOURCE),
        "branch": "2+1+1",
        "scope": "the 35 support-23 weighted-degree-14 global compatibility polynomials span exactly four polynomials over Q",
        "coefficient_field": "Q",
        "selected_source_rows": [
            {"row": index, "tag": tag} for index, tag, _polynomial in selected
        ],
        "pivot_monomials_by_support_order": [
            [[VARIABLES[j], int(power)] for j, power in enumerate(monomial) if power]
            for monomial in pivots
        ],
        "exact_rank": len(selected),
        "source_equation_count": len(source_rows),
        "source_term_count": sum(len(row[2]) for row in source_rows),
        "coordinate_axis_screen": axis_report,
        "coordinate_axes_excluded_by_this_subideal": [
            name for name, data in axis_report.items()
            if not data["coordinate_axis_survives_this_35_equation_subideal"]
        ],
        "coordinate_axes_surviving_this_subideal": [
            name for name, data in axis_report.items()
            if data["coordinate_axis_survives_this_35_equation_subideal"]
        ],
        "all_source_rows_exactly_replayed_over_Q": True,
        "row_relations": relations,
        "all_relation_denominators_invertible_mod_251": all_relation_denominators_invertible,
        "selected_rows_pivot_determinant_mod_251": modular_pivot_matrix_determinant,
        "all_source_rows_exactly_replayed_mod_251": True,
        "modular_equivalent_macaulay_shape_at_weight_24": [4 * 396, 156376],
        "original_macaulay_shape_at_weight_24": [35 * 396, 156376],
        "interpretation": (
            "At every higher weighted degree, the full 35-equation Macaulay row span equals the row span generated by the four selected source equations, because all 35 source equations are their exact Q-linear combinations and those coefficients reduce modulo 251. This reduces weight-24 construction from 13,860 to 1,584 rows; it is a computation reduction, not a branch exclusion."
        ),
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in result.items() if key != "row_relations"}, indent=2), flush=True)
    print(f"receipt={OUT} sha256={sha256_file(OUT)}", flush=True)


if __name__ == "__main__":
    main()
