#!/usr/bin/env python3
"""Compute the exact Q-linear span of all 139 depth-14 compatibility rows."""

from __future__ import annotations

import hashlib
import json
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import ijson


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_depth14_global_compatibility_reduction_20261021.json"
OUT = HERE / "f2_75_125_211_global_compatibility_exact_span.json"
sys.path.insert(0, str(HERE / "f2_75_125_flint_deps"))
from flint import fmpq, fmpq_mpoly_ctx  # noqa: E402


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    with SOURCE.open("rb") as stream:
        names = list(ijson.items(stream, "old_parameters.item"))
    with SOURCE.open("rb") as stream:
        weights = list(ijson.items(stream, "old_parameter_weights.item"))
    if len(names) != 42 or len(weights) != 42:
        raise AssertionError("expected the 42 old parameters and their weights")
    ctx = fmpq_mpoly_ctx.get(tuple(names), "lex")

    echelon = {}
    selected = []
    support_histogram: Counter[int] = Counter()
    term_count = 0
    dependent_rows = []
    independent_rows = []
    with SOURCE.open("rb") as stream:
        for row in ijson.items(stream, "global_compatibility_equations.item"):
            terms = {}
            support = set()
            for sparse, coefficient in row["terms"]:
                exponent = [0] * len(names)
                for variable_position, power in sparse:
                    exponent[variable_position] = power
                    if power:
                        support.add(variable_position)
                q = Fraction(coefficient)
                if sum(w * e for w, e in zip(weights, exponent)) != 14:
                    raise AssertionError(f"wrong weighted degree at {row['constraint_tag']}")
                terms[tuple(exponent)] = fmpq(q.numerator, q.denominator)
            polynomial = ctx.from_dict(terms)
            support_histogram[len(support)] += 1
            term_count += len(terms)

            reduced = polynomial
            while not reduced.is_zero():
                coefficients = reduced.to_dict()
                pivot = max(coefficients)
                coefficient = coefficients[pivot]
                if pivot in echelon:
                    reduced -= echelon[pivot] * coefficient
                else:
                    echelon[pivot] = reduced * (fmpq(1) / coefficient)
                    selected.append({
                        "equation_index": row["constraint_row_zero_based"],
                        "tag": row["constraint_tag"],
                        "pivot_monomial": [
                            [names[j], int(power)]
                            for j, power in enumerate(pivot) if power
                        ],
                    })
                    independent_rows.append(row["constraint_row_zero_based"])
                    break
            else:
                dependent_rows.append(row["constraint_row_zero_based"])

    if len(independent_rows) + len(dependent_rows) != 139:
        raise AssertionError("not all 139 compatibility rows were processed")

    result = {
        "source": SOURCE.name,
        "source_sha256": sha256_file(SOURCE),
        "branch": "2+1+1",
        "scope": "exact Q-linear span of the 139 weighted-degree-14 global compatibility equations",
        "coefficient_field": "Q",
        "monomial_order_for_row_reduction": "lexicographic on the saved 42-parameter exponent tuples",
        "equation_count": 139,
        "term_count": term_count,
        "support_size_histogram": dict(sorted(support_histogram.items())),
        "exact_vector_space_rank": len(independent_rows),
        "independent_source_rows": independent_rows,
        "dependent_source_rows": dependent_rows,
        "independent_row_pivots": selected,
        "interpretation": (
            "This is Q-linear row-space reduction of polynomial generators, not Gröbner-basis reduction and not a branch exclusion. It may replace the 139 displayed equations by an exact smaller generating set for their Q-linear span; additional nonlinear ideal consequences and the full 2+1+1 conditions remain to be analyzed."
        ),
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    print(f"receipt={OUT} sha256={sha256_file(OUT)}", flush=True)


if __name__ == "__main__":
    main()
