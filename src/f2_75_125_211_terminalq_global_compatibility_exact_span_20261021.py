#!/usr/bin/env python3
"""Compute the exact Q-linear span of the 149 extended compatibility rows."""

from __future__ import annotations

import hashlib
import json
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import ijson


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "f2_75_125_211_terminalq_global_compatibility_reduction_20261021.json"
OUTPUT = ROOT / "f2_75_125_211_terminalq_global_compatibility_exact_span_20261021.json"
EXPECTED_SHA256 = "1171348e859a13cc47540b825e12a9de5eef858f5020ee9fa49bdaae6b5ce8c9"
WEIGHT_BY_PARAMETER = {
    **{0: 1, 1: 2, 2: 3, 3: 3, 4: 3},
    **{i: 4 for i in range(5, 13)},
    **{i: 5 for i in range(13, 20)},
    **{i: 6 for i in range(20, 25)},
    **{i: 7 for i in range(25, 30)},
    **{i: 8 for i in range(30, 36)},
    **{i: 9 for i in range(37, 43)},
}

sys.path.insert(0, str(ROOT / "f2_75_125_flint_deps"))
from flint import fmpq, fmpq_mpoly_ctx, __version__ as FLINT_VERSION  # noqa: E402


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    source_hash = sha256_file(SOURCE)
    if source_hash != EXPECTED_SHA256:
        raise SystemExit(f"source SHA mismatch: {source_hash}")
    with SOURCE.open("rb") as stream:
        names = list(ijson.items(stream, "old_parameters.item"))
    weights = [WEIGHT_BY_PARAMETER[int(name[1:])] for name in names]
    if len(names) != 42 or len(weights) != 42:
        raise AssertionError("expected 42 old parameters and weights")
    ctx = fmpq_mpoly_ctx.get(tuple(names), "lex")

    echelon = {}
    selected = []
    support_histogram: Counter[int] = Counter()
    term_count = 0
    dependent_rows = []
    independent_rows = []
    with SOURCE.open("rb") as stream:
        for row in ijson.items(stream, "compatibility_equations.item"):
            terms = {}
            support = set()
            for sparse, coefficient in row["forcing_polynomial"]:
                exponent = [0] * len(names)
                for variable_position, power in sparse:
                    exponent[int(variable_position)] = int(power)
                    if power:
                        support.add(int(variable_position))
                q = Fraction(coefficient)
                if sum(w * e for w, e in zip(weights, exponent)) != 14:
                    raise AssertionError(f"wrong weighted degree at {row['tag']}")
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
                        "equation_row": row["row"],
                        "tag": row["tag"],
                        "pivot_monomial": [
                            [names[j], int(power)]
                            for j, power in enumerate(pivot) if power
                        ],
                    })
                    independent_rows.append(row["row"])
                    break
            else:
                dependent_rows.append(row["row"])

    if len(independent_rows) + len(dependent_rows) != 149:
        raise AssertionError("not all 149 compatibility equations were processed")

    result = {
        "source": SOURCE.name,
        "source_sha256": source_hash,
        "branch": "2+1+1",
        "scope": "exact Q-linear span of the 149 weighted-degree-14 compatibility equations after adding distance-14 terminal-Q rows",
        "coefficient_field": "Q",
        "flint_version": FLINT_VERSION,
        "monomial_order_for_row_reduction": "lexicographic on the saved 42-parameter exponent tuples",
        "equation_count": 149,
        "term_count": term_count,
        "support_size_histogram": dict(sorted(support_histogram.items())),
        "exact_vector_space_rank": len(independent_rows),
        "independent_equation_rows": independent_rows,
        "dependent_equation_rows": dependent_rows,
        "independent_row_pivots": selected,
        "interpretation": "This is exact Q-linear row-space reduction of polynomial generators, not Gröbner-basis reduction and not a branch exclusion.",
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "source_sha256", "equation_count", "term_count", "support_size_histogram",
        "exact_vector_space_rank", "interpretation",
    )}, indent=2), flush=True)
    print(f"output={OUTPUT.name} sha256={sha256_file(OUTPUT)}", flush=True)


if __name__ == "__main__":
    main()
