"""Memory-bounded exact profile for a rational sparse F2 checkpoint.

The checkpoint text is kept as one string, but its constraint array is decoded
one polynomial at a time.  The report is diagnostic: it does not claim that
the full branch or any Keller pair has been excluded.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter
from fractions import Fraction
from pathlib import Path

from f2_75_125_quadratic_affine_profile_stream import (
    iter_json_array,
    top_metadata_and_constraints,
)


def exact_rank(rows: list[list[Fraction]]) -> int:
    basis: dict[int, list[Fraction]] = {}
    for source in rows:
        row = list(source)
        for pivot in sorted(basis):
            factor = row[pivot]
            if factor:
                row = [a - factor * b for a, b in zip(row, basis[pivot])]
        pivot = next((i for i, value in enumerate(row) if value), None)
        if pivot is not None:
            scale = row[pivot]
            basis[pivot] = [value / scale for value in row]
    return len(basis)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def profile(source: Path) -> dict:
    text = source.read_text(encoding="utf-8")
    metadata, constraint_pos = top_metadata_and_constraints(text)
    if metadata.get("format") != "f2-sparse-ring-v1":
        raise ValueError("unexpected sparse checkpoint format")
    if metadata.get("quadratic"):
        raise ValueError("expected a rational coefficient checkpoint")

    marker = re.search(r',\s*"constraints"\s*:', text)
    if marker is None:
        raise ValueError("could not locate the constraint-array header")
    header = json.loads(text[: marker.start()] + "}")
    numeric_steps = [
        row["step"] for row in header.get("receipt", [])
        if isinstance(row.get("step"), int)
    ]
    depth = max(numeric_steps) if numeric_steps else None

    names = metadata["generators"]
    nparameters = int(metadata["nparameters"])
    parameter_indices = {
        names.index(f"p{i}"): i for i in range(nparameters)
    }

    affine_rows: list[list[Fraction]] = []
    constraint_count = term_count = zero_count = 0
    support_histogram: Counter[int] = Counter()
    degree_histogram: Counter[int] = Counter()
    forced_zero_rows: list[dict] = []
    univariate_rows: list[dict] = []
    constant_rows: list[dict] = []
    low_support_count = 0
    low_support_sample: list[dict] = []
    nonparameter_term_count = 0
    nonparameter_generator_counts: Counter[str] = Counter()

    for tag, polynomial in iter_json_array(text, constraint_pos):
        constraint_count += 1
        if not polynomial:
            zero_count += 1
            continue

        support: set[int] = set()
        collapsed: dict[tuple[tuple[int, int], ...], Fraction] = {}
        max_degree = 0
        affine = True
        has_nonparameter_generator = False

        for sparse, values in polynomial:
            term_count += 1
            coefficient = Fraction(values[0])
            monomial: list[tuple[int, int]] = []
            degree = 0
            for generator_index, power in sparse:
                if not power:
                    continue
                if generator_index in parameter_indices:
                    parameter = parameter_indices[generator_index]
                    support.add(parameter)
                    monomial.append((parameter, power))
                    degree += power
                else:
                    has_nonparameter_generator = True
                    nonparameter_term_count += 1
                    nonparameter_generator_counts[names[generator_index]] += 1
            max_degree = max(max_degree, degree)
            degree_histogram[degree] += 1
            key = tuple(monomial)
            collapsed[key] = collapsed.get(key, Fraction(0)) + coefficient
            if has_nonparameter_generator or degree > 1 or len(monomial) > 1:
                affine = False

        support_histogram[len(support)] += 1
        if len(support) <= 3 and not has_nonparameter_generator:
            low_support_count += 1
            if len(low_support_sample) < 40:
                low_support_sample.append({
                    "tag": tag,
                    "support": [f"p{i}" for i in sorted(support)],
                    "terms": len(polynomial),
                    "max_parameter_degree": max_degree,
                })

        if affine:
            row = [Fraction(0)] * (nparameters + 1)
            for monomial, coefficient in collapsed.items():
                if not monomial:
                    row[0] += coefficient
                elif len(monomial) == 1 and monomial[0][1] == 1:
                    row[monomial[0][0] + 1] += coefficient
                else:
                    affine = False
                    break
            if affine:
                affine_rows.append(row)

        if len(collapsed) == 1:
            monomial, coefficient = next(iter(collapsed.items()))
            if coefficient and not has_nonparameter_generator:
                if len(monomial) == 1:
                    parameter, power = monomial[0]
                    forced_zero_rows.append({
                        "parameter": f"p{parameter}",
                        "power": power,
                        "tag": tag,
                        "coefficient": str(coefficient),
                    })
                elif not monomial:
                    constant_rows.append({"tag": tag, "value": str(coefficient)})

        if len(support) == 1 and not has_nonparameter_generator:
            parameter = next(iter(support))
            univariate_rows.append({
                "tag": tag,
                "parameter": f"p{parameter}",
                "terms": len(polynomial),
                "terms_by_power": {
                    str(monomial[0][1]) if monomial else "0": str(coefficient)
                    for monomial, coefficient in sorted(collapsed.items())
                    if coefficient
                },
            })

    coefficient_rank = exact_rank([row[1:] for row in affine_rows])
    augmented_rank = exact_rank(affine_rows)
    return {
        "source_checkpoint": str(source),
        "source_sha256": sha256(source),
        "branch": metadata["branch"],
        "checkpoint_depth": depth,
        "carried_parameter_count": nparameters,
        "constraint_count": constraint_count,
        "term_count": term_count,
        "zero_constraint_count": zero_count,
        "nonparameter_term_count": nonparameter_term_count,
        "nonparameter_generator_counts": dict(sorted(nonparameter_generator_counts.items())),
        "support_histogram": {
            str(k): v for k, v in sorted(support_histogram.items())
        },
        "parameter_degree_histogram": {
            str(k): v for k, v in sorted(degree_histogram.items())
        },
        "affine_constraint_count": len(affine_rows),
        "affine_rank_coefficients_Q": coefficient_rank,
        "affine_rank_augmented_Q": augmented_rank,
        "affine_system_inconsistent": augmented_rank > coefficient_rank,
        "forced_coordinate_zero_rows": forced_zero_rows,
        "univariate_constraint_count": len(univariate_rows),
        "univariate_constraints": univariate_rows,
        "nonzero_constant_constraints": constant_rows,
        "constraints_with_support_at_most_3_count": low_support_count,
        "constraints_with_support_at_most_3_sample": low_support_sample,
        "scope": (
            "Exact profile of saved necessary equations only; not a Keller-pair "
            "or whole-branch exclusion certificate."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    result = profile(args.source)
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    summary_keys = (
        "source_sha256", "branch", "checkpoint_depth", "carried_parameter_count",
        "constraint_count", "term_count", "zero_constraint_count",
        "nonparameter_term_count", "nonparameter_generator_counts",
        "support_histogram", "parameter_degree_histogram", "affine_constraint_count",
        "affine_rank_coefficients_Q", "affine_rank_augmented_Q",
        "affine_system_inconsistent", "forced_coordinate_zero_rows",
        "univariate_constraint_count", "nonzero_constant_constraints",
        "constraints_with_support_at_most_3_count",
    )
    print(json.dumps({key: result[key] for key in summary_keys}, indent=2))
    print(f"report={args.report}")


if __name__ == "__main__":
    main()
