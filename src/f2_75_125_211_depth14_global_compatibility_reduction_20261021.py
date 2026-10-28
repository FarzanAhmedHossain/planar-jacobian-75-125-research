#!/usr/bin/env python3
"""Globally eliminate row dependencies in the depth-14 affine system.

The 141 constraints are affine-linear in p44,...,p50.  This checker proves
that all 141 coefficient rows lie in the *constant* Q-span of rows 1 and 2,
not merely their rational-function-field span.  It then writes the exact 139
forcing compatibility equations in the 42 old parameters.  These equations
are necessary on every pivot stratum.  After they hold, the only remaining
new-parameter system is 2-by-7; the locus where its rank drops still needs
separate treatment.
"""

from __future__ import annotations

import hashlib
import json
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import ijson


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE / "f2_75_125_flint_deps"))
from flint import fmpq, fmpq_mpoly_ctx, __version__ as FLINT_VERSION  # noqa: E402

SOURCE = HERE / "f2_75_125_211_p43zero_depth14_work_flint_20261021.json"
OUT = HERE / "f2_75_125_211_depth14_global_compatibility_reduction_20261021.json"
OLD_NAMES = [*(f"p{i}" for i in range(36)), *(f"p{i}" for i in range(37, 43))]
NEW_NAMES = [f"p{i}" for i in range(44, 51)]
WEIGHT_BY_PARAMETER = {
    **{0: 1, 1: 2, 2: 3, 3: 3, 4: 3},
    **{i: 4 for i in range(5, 13)},
    **{i: 5 for i in range(13, 20)},
    **{i: 6 for i in range(20, 25)},
    **{i: 7 for i in range(25, 30)},
    **{i: 8 for i in range(30, 36)},
    **{i: 9 for i in range(37, 43)},
}
OLD_WEIGHTS = [WEIGHT_BY_PARAMETER[int(name[1:])] for name in OLD_NAMES]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def as_mpoly(ctx, terms: dict[tuple[int, ...], Fraction]):
    return ctx.from_dict(
        {
            exponent: fmpq(coefficient.numerator, coefficient.denominator)
            for exponent, coefficient in terms.items()
            if coefficient
        }
    )


def encode_poly(poly) -> list[list[object]]:
    result = []
    for exponents, coefficient in sorted(poly.to_dict().items()):
        sparse = [[i, int(power)] for i, power in enumerate(exponents) if power]
        result.append([sparse, str(coefficient)])
    return result


def weighted_degrees(poly) -> set[int]:
    return {
        sum(weight * power for weight, power in zip(OLD_WEIGHTS, exponents))
        for exponents in poly.to_dict()
    }


def coefficient_at(poly, exponents: tuple[int, ...]):
    return poly.to_dict().get(exponents, fmpq(0))


def main() -> None:
    with SOURCE.open("rb") as stream:
        generators = list(ijson.items(stream, "generators.item"))
    old_index = {generators.index(name): i for i, name in enumerate(OLD_NAMES)}
    new_index = {generators.index(name): i for i, name in enumerate(NEW_NAMES)}
    ctx = fmpq_mpoly_ctx.get(tuple(OLD_NAMES), "lex")

    raw_rows: list[list[dict[tuple[int, ...], Fraction]]] = []
    tags: list[str] = []
    max_new_degree = 0
    with SOURCE.open("rb") as stream:
        for tag, terms in ijson.items(stream, "constraints.item"):
            entries: list[dict[tuple[int, ...], Fraction]] = [
                {} for _ in range(len(NEW_NAMES) + 1)
            ]
            for sparse, values in terms:
                coefficient = Fraction(values[0])
                new_factors = [(index, exponent) for index, exponent in sparse if index in new_index]
                degree = sum(exponent for _, exponent in new_factors)
                max_new_degree = max(max_new_degree, degree)
                if degree > 1 or (new_factors and (len(new_factors) != 1 or new_factors[0][1] != 1)):
                    raise AssertionError(f"constraint {tag} is nonlinear in p44,...,p50")
                column = new_index[new_factors[0][0]] if new_factors else len(NEW_NAMES)
                exponents = [0] * len(OLD_NAMES)
                for index, exponent in sparse:
                    if index in old_index:
                        exponents[old_index[index]] += exponent
                    elif index in new_index:
                        continue
                    elif generators[index] == "t":
                        raise ValueError("unexpected active t variable in a constraint")
                    else:
                        raise ValueError(f"unclassified generator: {generators[index]}")
                monomial = tuple(exponents)
                entries[column][monomial] = entries[column].get(monomial, Fraction(0)) + coefficient
            raw_rows.append(entries)
            tags.append(tag)

    if len(raw_rows) != 141 or max_new_degree != 1:
        raise AssertionError("unexpected depth-14 constraint shape")
    matrix = [[as_mpoly(ctx, terms) for terms in row] for row in raw_rows]
    row0, row1 = 1, 2

    # View each polynomial coefficient as a coordinate in one large QQ-vector
    # space. Find two coordinates proving the pivot rows independent over Q.
    coordinates = []
    for column in range(len(NEW_NAMES)):
        d0 = matrix[row0][column].to_dict()
        d1 = matrix[row1][column].to_dict()
        for exponent in sorted(set(d0) | set(d1)):
            a0 = d0.get(exponent, fmpq(0))
            a1 = d1.get(exponent, fmpq(0))
            if a0 or a1:
                coordinates.append((column, exponent, a0, a1))

    scalar_pivot = None
    for i, first in enumerate(coordinates):
        for second in coordinates[i + 1 :]:
            determinant = first[2] * second[3] - second[2] * first[3]
            if determinant:
                scalar_pivot = (first, second, determinant)
                break
        if scalar_pivot:
            break
    if scalar_pivot is None:
        raise AssertionError("pivot rows are not independent over Q")
    first, second, scalar_det = scalar_pivot

    row_combinations = []
    residuals = []
    residual_support_histogram: Counter[int] = Counter()
    residual_digests: set[str] = set()
    for row in range(len(matrix)):
        x = coefficient_at(matrix[row][first[0]], first[1])
        y = coefficient_at(matrix[row][second[0]], second[1])
        alpha = (x * second[3] - y * first[3]) / scalar_det
        beta = (first[2] * y - second[2] * x) / scalar_det

        for column in range(len(NEW_NAMES)):
            remainder = matrix[row][column] - alpha * matrix[row0][column] - beta * matrix[row1][column]
            if remainder:
                raise AssertionError(f"row {row} is not in the constant Q-span at column {column}")
        row_combinations.append([str(alpha), str(beta)])

        residual = matrix[row][7] - alpha * matrix[row0][7] - beta * matrix[row1][7]
        if not residual:
            continue
        degrees = weighted_degrees(residual)
        if degrees != {14}:
            raise AssertionError(f"forcing residual row {row} has unexpected weight {degrees}")
        terms = encode_poly(residual)
        support = sorted({
            i for exponents, _coefficient in residual.to_dict().items()
            for i, power in enumerate(exponents) if power
        })
        residual_support_histogram[len(support)] += 1
        digest = sha256_bytes(str(residual).encode("utf-8"))
        if digest in residual_digests:
            raise AssertionError("duplicate forcing residual")
        residual_digests.add(digest)
        residuals.append({
            "constraint_row_zero_based": row,
            "constraint_tag": tags[row],
            "term_count": len(terms),
            "support_variables": [OLD_NAMES[i] for i in support],
            "sha256_of_flint_text": digest,
            "terms": terms,
        })

    if len(residuals) != 139 or len(residual_digests) != 139:
        raise AssertionError(f"expected 139 distinct compatibility equations, found {len(residuals)}")
    if row_combinations[row0] != ["1", "0"] or row_combinations[row1] != ["0", "1"]:
        raise AssertionError("pivot rows did not normalize to the identity basis")

    basis_rows = []
    for row in (row0, row1):
        basis_rows.append({
            "constraint_row_zero_based": row,
            "constraint_tag": tags[row],
            "new_parameter_coefficients": [encode_poly(matrix[row][j]) for j in range(7)],
            "forcing_polynomial": encode_poly(matrix[row][7]),
        })

    result = {
        "source": str(SOURCE.relative_to(ROOT)),
        "source_sha256": sha256_file(SOURCE),
        "branch": "2+1+1",
        "scope": "global exact row-space elimination of p44,...,p50; necessary for every pivot/rank stratum",
        "coefficient_field": "Q",
        "flint_version": FLINT_VERSION,
        "constraint_count": len(matrix),
        "old_parameters": OLD_NAMES,
        "old_parameter_weights": OLD_WEIGHTS,
        "new_parameters": NEW_NAMES,
        "max_new_parameter_degree": max_new_degree,
        "coefficient_row_rank_over_Q_rational_function_field": 2,
        "coefficient_row_rank_over_constant_Q_span": 2,
        "constant_span_coordinate_pivot": {
            "basis_rows_zero_based": [row0, row1],
            "basis_row_tags": [tags[row0], tags[row1]],
            "first_coordinate": {
                "column": NEW_NAMES[first[0]],
                "monomial": [[i, int(power)] for i, power in enumerate(first[1]) if power],
                "row_values": [str(first[2]), str(first[3])],
            },
            "second_coordinate": {
                "column": NEW_NAMES[second[0]],
                "monomial": [[i, int(power)] for i, power in enumerate(second[1]) if power],
                "row_values": [str(second[2]), str(second[3])],
            },
            "determinant": str(scalar_det),
        },
        "row_combination_coefficients_in_Q": row_combinations,
        "coefficient_basis_rows": basis_rows,
        "global_compatibility_equation_count": len(residuals),
        "global_compatibility_support_size_histogram": dict(sorted(residual_support_histogram.items())),
        "distinct_compatibility_equation_count": len(residual_digests),
        "global_compatibility_equations": residuals,
        "interpretation": (
            "The 141-by-7 coefficient matrix factors as C*U with C a constant 141-by-2 rational matrix. "
            "The saved 139 equations are exactly the conditions that the forcing column lie in im(C); "
            "they are globally necessary and sufficient for row consistency. On that locus the original "
            "affine system reduces to two equations in seven new parameters. Points where the resulting "
            "2-by-7 matrix U has rank below 2 still require a separate check."
        ),
    }
    OUT.write_text(json.dumps(result, separators=(",", ":")) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "source_sha256", "constraint_count", "coefficient_row_rank_over_Q_rational_function_field",
        "coefficient_row_rank_over_constant_Q_span", "constant_span_coordinate_pivot",
        "global_compatibility_equation_count", "global_compatibility_support_size_histogram",
        "distinct_compatibility_equation_count", "interpretation",
    )}, separators=(",", ":")), flush=True)
    print(f"reduction={OUT} bytes={OUT.stat().st_size} sha256={sha256_file(OUT)}", flush=True)


if __name__ == "__main__":
    main()
