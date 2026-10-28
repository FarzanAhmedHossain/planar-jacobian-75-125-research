"""Profile exact sparse-ring constraints in a saved F2 lift checkpoint.

The script records polynomial-degree strata, detects constant obstructions,
and computes the exact rational rank of all constraints that are affine in
the carried parameters.  It is diagnostic only; it does not drop equations or
claim branch exclusion.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path


def exact_rank(rows: list[list[Fraction]]) -> tuple[int, list[int]]:
    basis: dict[int, list[Fraction]] = {}
    for source in rows:
        row = list(source)
        for pivot in sorted(basis):
            if row[pivot]:
                factor = row[pivot]
                row = [a - factor * b for a, b in zip(row, basis[pivot])]
        pivot = next((i for i, value in enumerate(row) if value), None)
        if pivot is not None:
            scale = row[pivot]
            basis[pivot] = [value / scale for value in row]
    return len(basis), sorted(basis)


def profile(path: Path) -> dict:
    raw = path.read_bytes()
    state = json.loads(raw)
    if state.get("format") != "f2-sparse-ring-v1":
        raise ValueError("unexpected checkpoint format")
    if state.get("quadratic"):
        raise ValueError("this profiler currently expects a rational coefficient field")

    names = state["generators"]
    nparameters = int(state["nparameters"])
    p_names = [f"p{i}" for i in range(nparameters)]
    p_indices = {names.index(name): i for i, name in enumerate(p_names)}
    term_count = 0
    degree_hist = Counter()
    mixed_degree = 0
    t_terms = 0
    other_generator_terms = 0
    constants = []
    identically_zero_constraints = 0
    affine_rows = []
    affine_tags = []
    constraint_profiles = []

    for tag, polynomial in state.get("constraints", []):
        if not polynomial:
            identically_zero_constraints += 1
        histogram = Counter()
        max_degree = 0
        affine = True
        vector = [Fraction(0)] * (nparameters + 1)
        for sparse, values in polynomial:
            term_count += 1
            coefficient = Fraction(values[0])
            degree = 0
            p_variable = None
            t_power = 0
            other_power = 0
            for generator_index, power in sparse:
                if generator_index == 0:
                    t_power += power
                elif generator_index in p_indices:
                    degree += power
                    if p_variable is None:
                        p_variable = p_indices[generator_index]
                    elif power:
                        p_variable = -1
                else:
                    other_power += power
            histogram[degree] += 1
            degree_hist[degree] += 1
            max_degree = max(max_degree, degree)
            t_terms += bool(t_power)
            other_generator_terms += bool(other_power)
            if t_power or other_power:
                affine = False
            if degree > 1:
                affine = False
            elif degree == 0:
                vector[0] += coefficient
            else:
                if p_variable is None or p_variable < 0:
                    affine = False
                else:
                    vector[p_variable + 1] += coefficient
        if len(histogram) > 1:
            mixed_degree += 1
        if len(polynomial) == 1 and not polynomial[0][0]:
            constants.append({"tag": tag, "value": polynomial[0][1][0]})
        if affine and polynomial:
            affine_rows.append(vector)
            affine_tags.append(tag)
        constraint_profiles.append({
            "tag": tag,
            "terms": len(polynomial),
            "max_parameter_degree": max_degree,
            "degree_histogram": {str(k): v for k, v in sorted(histogram.items())},
            "affine_parameter_equation": affine,
        })

    augmented_rank, _ = exact_rank(affine_rows)
    coefficient_rank, pivots = exact_rank([row[1:] for row in affine_rows])
    affine_basis = []
    for pivot in pivots:
        row_index = next(i for i, row in enumerate(affine_rows)
                         if row[pivot + 1])
        row = affine_rows[row_index]
        scale = row[pivot + 1]
        coefficients = {"constant": str(row[0] / scale)}
        for i, value in enumerate(row[1:]):
            if value:
                coefficients[p_names[i]] = str(value / scale)
        affine_basis.append({
            "tag": affine_tags[row_index],
            "normalized_pivot": p_names[pivot],
            "equation_coefficients": coefficients,
        })
    return {
        "source_checkpoint": str(path),
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "branch": state["branch"],
        "checkpoint_depth": len(state["p"]) - 1,
        "carried_parameter_count": nparameters,
        "constraint_count": len(constraint_profiles),
        "constraint_term_count": term_count,
        "parameter_degree_histogram": {str(k): v for k, v in sorted(degree_hist.items())},
        "mixed_parameter_degree_constraint_count": mixed_degree,
        "terms_with_t_powers": t_terms,
        "terms_with_nonparameter_generators": other_generator_terms,
        "nonzero_constant_obstructions": constants,
        "identically_zero_constraint_count": identically_zero_constraints,
        "affine_constraint_count": len(affine_rows),
        "affine_coefficient_rank_over_Q": coefficient_rank,
        "affine_augmented_rank_over_Q": augmented_rank,
        "affine_system_inconsistent": augmented_rank > coefficient_rank,
        "affine_pivot_parameters": [p_names[i] for i in pivots if 0 <= i < nparameters],
        "affine_basis_equations": affine_basis,
        "affine_tags": affine_tags,
        "constraint_profiles": constraint_profiles,
        "scope": "Exact diagnostic profile only; equations are not removed and no branch is excluded.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("receipt", type=Path)
    args = parser.parse_args()
    result = profile(args.checkpoint)
    result["receipt"] = str(args.receipt)
    args.receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "source_sha256", "branch", "checkpoint_depth", "carried_parameter_count",
        "constraint_count", "constraint_term_count", "parameter_degree_histogram",
        "mixed_parameter_degree_constraint_count", "terms_with_t_powers",
        "terms_with_nonparameter_generators", "nonzero_constant_obstructions",
        "identically_zero_constraint_count",
        "affine_constraint_count", "affine_coefficient_rank_over_Q",
        "affine_augmented_rank_over_Q", "affine_system_inconsistent",
        "affine_pivot_parameters", "affine_basis_equations", "receipt",
    )}, indent=2))


if __name__ == "__main__":
    main()
