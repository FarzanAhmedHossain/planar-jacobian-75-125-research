#!/usr/bin/env python3
"""Exact sparse scan for affine residual equations in an F2 lift checkpoint.

Only rows whose every monomial has parameter total degree at most one are
collected.  Their rational row space is reduced incrementally; no full JSON
checkpoint or multivariate polynomial ring is materialized.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

from f2_75_125_211_streaming_checkpoint_scan import JsonStream


def scan(source: Path, nparameters: int, source_sha256: str) -> dict:
    reader = JsonStream(source)
    basis: dict[int, dict] = {}
    degree_histogram: Counter[int] = Counter()
    row_tags: dict[int, str] = {}
    linear_rows = nonlinear_rows = auxiliary_rows = 0
    zero_rows = dependent_rows = 0
    contradiction = None

    try:
        metadata, _ = reader.top_level_constraints()
        if metadata.get("format") != "f2-sparse-ring-v1":
            raise ValueError("unexpected checkpoint format")
        if metadata.get("branch") != "2+1+1" or metadata.get("quadratic"):
            raise ValueError("expected rational 2+1+1 checkpoint")
        if int(metadata.get("nparameters", -1)) != nparameters:
            raise ValueError("parameter count disagrees with checkpoint")

        row_count = 0
        for row in reader.array_values():
            row_count += 1
            tag, terms = row
            row_tags[row_count] = tag
            coefficients: dict[int, Fraction] = {}
            constant = Fraction(0)
            maximum_degree = 0
            nonlinear = False
            auxiliary = False

            for sparse, values in terms:
                if len(values) != 1:
                    raise ValueError("expected rational coefficients")
                coefficient = Fraction(values[0])
                monomial_degree = sum(exponent for _, exponent in sparse)
                maximum_degree = max(maximum_degree, monomial_degree)
                if monomial_degree > 1:
                    nonlinear = True
                    continue
                if not sparse:
                    constant += coefficient
                    continue
                if len(sparse) != 1:
                    nonlinear = True
                    break
                generator, exponent = sparse[0]
                if exponent != 1:
                    nonlinear = True
                    break
                if not 1 <= generator <= nparameters:
                    auxiliary = True
                    continue
                variable = generator - 1
                coefficients[variable] = coefficients.get(variable, Fraction(0)) + coefficient

            degree_histogram[maximum_degree] += 1
            if nonlinear:
                nonlinear_rows += 1
                continue
            if auxiliary:
                auxiliary_rows += 1
                continue

            linear_rows += 1
            coefficients = {v: c for v, c in coefficients.items() if c}
            combination = {row_count: Fraction(1)}
            if not coefficients:
                if constant:
                    contradiction = {"tag": tag, "constant": str(constant)}
                    break
                zero_rows += 1
                continue

            for pivot in sorted(basis):
                factor = coefficients.get(pivot, Fraction(0))
                if not factor:
                    continue
                pivot_row = basis[pivot]
                for variable, value in pivot_row["coefficients"].items():
                    coefficients[variable] = coefficients.get(variable, Fraction(0)) - factor * value
                    if not coefficients[variable]:
                        del coefficients[variable]
                constant -= factor * pivot_row["constant"]
                for original_row, value in pivot_row["combination"].items():
                    combination[original_row] = combination.get(original_row, Fraction(0)) - factor * value
                    if not combination[original_row]:
                        del combination[original_row]

            if not coefficients:
                if constant:
                    contradiction = {
                        "row_index": row_count,
                        "tag": tag,
                        "constant": str(constant),
                        "combination": [
                            {"row_index": index, "tag": row_tags[index], "coefficient": str(value)}
                            for index, value in sorted(combination.items())
                        ],
                    }
                    break
                dependent_rows += 1
                continue

            pivot = min(coefficients)
            scale = coefficients[pivot]
            coefficients = {v: c / scale for v, c in coefficients.items()}
            constant /= scale
            combination = {index: value / scale for index, value in combination.items()}
            basis[pivot] = {
                "coefficients": coefficients,
                "constant": constant,
                "tag": tag,
                "combination": combination,
            }

        serialized_basis = []
        for pivot, row in sorted(basis.items()):
            serialized_basis.append({
                "pivot": f"p{pivot}",
                "tag": row["tag"],
                "constant": str(row["constant"]),
                "coefficients": {f"p{v}": str(c)
                                 for v, c in sorted(row["coefficients"].items())},
            })
        return {
            "source": str(source),
            "source_sha256": source_sha256,
            "branch": metadata["branch"],
            "nparameters": nparameters,
            "rows_scanned": row_count,
            "linear_rows": linear_rows,
            "nonlinear_rows": nonlinear_rows,
            "auxiliary_rows": auxiliary_rows,
            "degree_histogram": {str(k): v for k, v in sorted(degree_histogram.items())},
            "linear_rank": len(basis),
            "zero_rows": zero_rows,
            "dependent_rows": dependent_rows,
            "linear_contradiction": contradiction,
            "linear_basis": serialized_basis,
            "scope": "Exact affine consequences of the saved necessary system only; a nonempty linear solution space is not a branch exclusion.",
        }
    finally:
        reader.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("nparameters", type=int)
    parser.add_argument("source_sha256")
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    result = scan(args.source, args.nparameters, args.source_sha256)
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "source_sha256", "rows_scanned", "linear_rows", "nonlinear_rows",
        "auxiliary_rows", "degree_histogram", "linear_rank", "zero_rows",
        "dependent_rows", "linear_contradiction",
    )}, indent=2))
    print(f"report={args.report}")


if __name__ == "__main__":
    main()
