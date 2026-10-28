#!/usr/bin/env python3
"""Exact affine-rank reduction of the depth-14 2+1+1 constraints.

The 141 retained constraints are affine-linear in p44,...,p50.  This checker
reconstructs their 141-by-7 coefficient matrix A and forcing column b over
QQ[p0,...,p35,p37,...,p42].  It proves rank(A)=2 over the rational-function
field by an exact nonzero 2-by-2 pivot and 705 exact zero cross-minors.  It
also proves rank([A|b]) >= 3 by one nonzero 3-by-3 minor.  Thus the system is
inconsistent on a nonempty Zariski-open set of old parameters; the exceptional
hypersurface remains unresolved, so this is not a branch exclusion.

Requires ijson.  The repository's bundled FLINT package is located below.
"""

from __future__ import annotations

import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

import ijson


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE / "f2_75_125_flint_deps"))
from flint import fmpq, fmpq_mpoly_ctx, __version__ as FLINT_VERSION  # noqa: E402

SOURCE = HERE / "f2_75_125_211_p43zero_depth14_work_flint_20261021.json"
OUT = HERE / "f2_75_125_211_depth14_affine_rank_reduction_20261021.json"
OLD_NAMES = [*(f"p{i}" for i in range(36)), *(f"p{i}" for i in range(37, 43))]
NEW_NAMES = [f"p{i}" for i in range(44, 51)]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def as_mpoly(ctx, terms: dict[tuple[int, ...], Fraction]):
    return ctx.from_dict(
        {
            exponent: fmpq(coefficient.numerator, coefficient.denominator)
            for exponent, coefficient in terms.items()
            if coefficient
        }
    )


def main() -> None:
    source_digest = sha256_file(SOURCE)
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
                new_degree = sum(exponent for _, exponent in new_factors)
                max_new_degree = max(max_new_degree, new_degree)
                if new_degree > 1 or (new_factors and (len(new_factors) != 1 or new_factors[0][1] != 1)):
                    raise AssertionError(f"constraint {tag} is not affine-linear in p44,...,p50")
                column = new_index[new_factors[0][0]] if new_factors else len(NEW_NAMES)
                exponents = [0] * len(OLD_NAMES)
                for index, exponent in sparse:
                    if index in old_index:
                        exponents[old_index[index]] += exponent
                    elif index in new_index:
                        continue
                    elif generators[index] == "t":
                        raise ValueError("unexpected active t variable in a constraint monomial")
                    else:
                        raise ValueError(f"unclassified generator in constraints: {generators[index]}")
                monomial = tuple(exponents)
                entries[column][monomial] = entries[column].get(monomial, Fraction(0)) + coefficient
            raw_rows.append(entries)
            tags.append(tag)

    matrix = [[as_mpoly(ctx, terms) for terms in row] for row in raw_rows]
    if len(matrix) != 141 or max_new_degree != 1:
        raise AssertionError("unexpected depth-14 constraint shape")

    # Pivot on constraint rows 1 and 2 and new-parameter columns p44,p45.
    r0, r1 = 1, 2
    c0, c1 = 0, 1
    pivot = matrix[r0][c0] * matrix[r1][c1] - matrix[r0][c1] * matrix[r1][c0]
    if not pivot:
        raise AssertionError("the proposed 2-by-2 pivot minor vanished")

    # With this pivot invertible in the rational-function field, it suffices
    # to check each remaining column against each row using these five minors.
    cross_minor_count = 0
    for row in range(len(matrix)):
        for column in range(2, len(NEW_NAMES)):
            determinant = (
                matrix[r0][c0]
                * (matrix[r1][c1] * matrix[row][column] - matrix[r1][column] * matrix[row][c1])
                - matrix[r0][c1]
                * (matrix[r1][c0] * matrix[row][column] - matrix[r1][column] * matrix[row][c0])
                + matrix[r0][column]
                * (matrix[r1][c0] * matrix[row][c1] - matrix[r1][c1] * matrix[row][c0])
            )
            cross_minor_count += 1
            if determinant:
                raise AssertionError(f"nonzero coefficient cross-minor at row={row}, col={column}")

    # A nonzero 3-by-3 minor of [A|b], using rows 1,2,0 and columns p44,p45,b.
    r2 = 0
    augmented_minor = (
        matrix[r0][0] * (matrix[r1][1] * matrix[r2][7] - matrix[r1][7] * matrix[r2][1])
        - matrix[r0][1] * (matrix[r1][0] * matrix[r2][7] - matrix[r1][7] * matrix[r2][0])
        + matrix[r0][7] * (matrix[r1][0] * matrix[r2][1] - matrix[r1][1] * matrix[r2][0])
    )
    if not augmented_minor:
        raise AssertionError("the proposed augmented minor vanished")

    pivot_text = str(pivot).encode("utf-8")
    augmented_text = str(augmented_minor).encode("utf-8")
    result = {
        "source": str(SOURCE.relative_to(ROOT)),
        "source_sha256": source_digest,
        "branch": "2+1+1",
        "constraint_count": len(matrix),
        "old_parameters": OLD_NAMES,
        "new_parameters": NEW_NAMES,
        "max_new_parameter_degree": max_new_degree,
        "coefficient_matrix_rank_over_Q_rational_function_field": 2,
        "pivot_minor": {
            "rows_zero_based": [r0, r1],
            "row_tags": [tags[r0], tags[r1]],
            "columns": [NEW_NAMES[c0], NEW_NAMES[c1]],
            "term_count": len(pivot),
            "sha256_of_flint_text": sha256_bytes(pivot_text),
        },
        "cross_minor_identities_checked_exactly_zero": cross_minor_count,
        "augmented_matrix_rank_lower_bound_over_Q_rational_function_field": 3,
        "augmented_minor_third_row_is_pure_old": not any(matrix[0][j] for j in range(7)),
        "augmented_minor": {
            "rows_zero_based": [r0, r1, r2],
            "row_tags": [tags[r0], tags[r1], tags[r2]],
            "columns": [NEW_NAMES[0], NEW_NAMES[1], "forcing_column"],
            "term_count": len(augmented_minor),
            "sha256_of_flint_text": sha256_bytes(augmented_text),
        },
        "flint_version": FLINT_VERSION,
        "interpretation": (
            "The displayed augmented minor uses row 0, whose new-parameter coefficient row is zero, so it is the coefficient pivot times an already-present pure-old forcing equation. "
            "The local rank discrepancy is correct but adds no exclusion beyond that existing equation. "
            "The globally useful reduction is recorded in f2_75_125_211_depth14_global_compatibility_reduction_20261021.json."
        ),
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
