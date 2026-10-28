#!/usr/bin/env python3
"""Construct exact compatibility residuals on the depth-14 pivot chart.

Input: the 141 exact depth-14 constraints, affine-linear in p44,...,p50.
Output: the pivot Delta and every nonzero quotient R_r from

    det(rows 1,2,r; columns p44,p45,forcing) = Delta * R_r.

The 139 R_r are exact necessary equations on the chart Delta != 0.  They are
not imposed on Delta=0; that pivot locus requires a separate analysis.  The
script also rechecks the exact rank-two structure of the coefficient matrix.
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
OUT = HERE / "f2_75_125_211_depth14_compatibility_atlas_20261021.json"
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


def determinant_pivot_forcing(matrix, row: int, row0: int = 1, row1: int = 2):
    return (
        matrix[row0][0] * (matrix[row1][1] * matrix[row][7] - matrix[row1][7] * matrix[row][1])
        - matrix[row0][1] * (matrix[row1][0] * matrix[row][7] - matrix[row1][7] * matrix[row][0])
        + matrix[row0][7] * (matrix[row1][0] * matrix[row][1] - matrix[row1][1] * matrix[row][0])
    )


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
    delta = matrix[row0][0] * matrix[row1][1] - matrix[row0][1] * matrix[row1][0]
    if not delta or weighted_degrees(delta) != {8}:
        raise AssertionError("pivot minor is zero or has an unexpected weight")

    # Prove every coefficient row lies in the span of the two pivot rows over
    # QQ(old parameters), using the nonzero 2x2 pivot and five exact minors per row.
    cross_minor_count = 0
    for row in range(len(matrix)):
        for column in range(2, len(NEW_NAMES)):
            minor = (
                matrix[row0][0]
                * (matrix[row1][1] * matrix[row][column] - matrix[row1][column] * matrix[row][1])
                - matrix[row0][1]
                * (matrix[row1][0] * matrix[row][column] - matrix[row1][column] * matrix[row][0])
                + matrix[row0][column]
                * (matrix[row1][0] * matrix[row][1] - matrix[row1][1] * matrix[row][0])
            )
            cross_minor_count += 1
            if minor:
                raise AssertionError(f"nonzero coefficient cross-minor at row={row}, col={column}")

    residuals = []
    zero_rows = []
    support_histogram: Counter[int] = Counter()
    distinct_hashes: set[str] = set()
    for row in range(len(matrix)):
        numerator = determinant_pivot_forcing(matrix, row)
        if not numerator:
            zero_rows.append(row)
            continue
        if numerator % delta:
            raise AssertionError(f"compatibility minor at row {row} is not divisible by the pivot")
        residual = numerator // delta
        if weighted_degrees(residual) != {14}:
            raise AssertionError(f"row {row} residual is not weighted-homogeneous of degree 14")
        polynomial_terms = encode_poly(residual)
        support = sorted({
            i for exponents, _coefficient in residual.to_dict().items()
            for i, power in enumerate(exponents) if power
        })
        support_histogram[len(support)] += 1
        canonical = str(residual).encode("utf-8")
        digest = sha256_bytes(canonical)
        if digest in distinct_hashes:
            raise AssertionError("unexpected duplicate compatibility residual")
        distinct_hashes.add(digest)
        residuals.append({
            "constraint_row_zero_based": row,
            "constraint_tag": tags[row],
            "new_coefficient_row_nonzero": any(matrix[row][column] for column in range(7)),
            "term_count": len(polynomial_terms),
            "support_variables": [OLD_NAMES[i] for i in support],
            "sha256_of_flint_text": digest,
            "terms": polynomial_terms,
        })

    if len(residuals) != 139 or zero_rows != [1, 2] or cross_minor_count != 705:
        raise AssertionError("unexpected compatibility atlas dimensions")

    delta_terms = encode_poly(delta)
    result = {
        "source": str(SOURCE.relative_to(ROOT)),
        "source_sha256": sha256_file(SOURCE),
        "branch": "2+1+1",
        "scope": "exact elimination of p44,...,p50 on the single pivot chart Delta!=0; Delta=0 is not covered",
        "coefficient_field": "Q",
        "flint_version": FLINT_VERSION,
        "constraint_count": len(matrix),
        "old_parameters": OLD_NAMES,
        "old_parameter_weights": OLD_WEIGHTS,
        "new_parameters": NEW_NAMES,
        "max_new_parameter_degree": max_new_degree,
        "coefficient_matrix_rank_over_Q_rational_function_field": 2,
        "cross_minor_identities_checked_exactly_zero": cross_minor_count,
        "pivot": {
            "rows_zero_based": [row0, row1],
            "row_tags": [tags[row0], tags[row1]],
            "columns": [NEW_NAMES[0], NEW_NAMES[1]],
            "weighted_degree": 8,
            "term_count": len(delta_terms),
            "sha256_of_flint_text": sha256_bytes(str(delta).encode("utf-8")),
            "terms": delta_terms,
        },
        "zero_compatibility_rows": zero_rows,
        "nonzero_compatibility_residual_count": len(residuals),
        "residual_support_size_histogram": dict(sorted(support_histogram.items())),
        "distinct_residual_count": len(distinct_hashes),
        "residuals": residuals,
        "interpretation": (
            "For Delta!=0, the 141-equation affine system in seven new parameters is solvable exactly when all saved residuals vanish. "
            "The residuals are exact necessary equations on that chart. The Delta=0 locus requires a separate pivot/rank-drop analysis."
        ),
    }
    OUT.write_text(json.dumps(result, separators=(",", ":")) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "source_sha256", "constraint_count", "coefficient_matrix_rank_over_Q_rational_function_field",
        "cross_minor_identities_checked_exactly_zero", "pivot", "zero_compatibility_rows",
        "nonzero_compatibility_residual_count", "residual_support_size_histogram",
        "distinct_residual_count", "interpretation",
    )}, separators=(",", ":")), flush=True)
    print(f"atlas={OUT} bytes={OUT.stat().st_size} sha256={sha256_file(OUT)}", flush=True)


if __name__ == "__main__":
    main()
