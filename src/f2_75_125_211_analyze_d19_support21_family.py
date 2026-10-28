"""Exact Q-row compression and weight-space audit for a d19 support family."""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
import time

import sympy as sp
from flint import fmpq
from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import QQ


INPUT = Path("research/f2_75_125_211_d19_support21_family_input_20261025.json")
OUTPUT = Path("research/f2_75_125_211_d19_support21_analysis_20261025.json")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def add_scaled(target, source, scale):
    if not scale:
        return
    for key, value in source.items():
        new_value = target.get(key, fmpq(0)) + scale*value
        if new_value:
            target[key] = new_value
        else:
            target.pop(key, None)


def main() -> None:
    started = time.monotonic()
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    names = tuple(data["target_support"])
    if len(names) != 21 or data["equation_count"] != len(data["equations"]):
        raise ValueError("support-21 input metadata mismatch")
    name_index = {name: i for i, name in enumerate(names)}

    pivots = {}
    independent_indices = []
    dependencies = []
    input_terms = 0
    difference_rows = set()
    weighted_degrees = []

    for equation_index, equation in enumerate(data["equations"]):
        row = {}
        monomials = []
        for monomial, coefficient in equation["terms"]:
            exponents = [0] * len(names)
            for name, power in monomial:
                exponents[name_index[name]] += int(power)
            key = tuple(exponents)
            value = fmpq(str(coefficient))
            row[key] = row.get(key, fmpq(0)) + value
            monomials.append(key)
            input_terms += 1
        row = {key: value for key, value in row.items() if value}
        if not row:
            dependencies.append(equation_index)
            continue

        base = monomials[0]
        for monomial in monomials[1:]:
            delta = tuple(a-b for a, b in zip(monomial, base))
            if any(delta):
                difference_rows.add(delta)

        while row:
            pivot = max(row)
            known = pivots.get(pivot)
            if known is None:
                coefficient = row[pivot]
                inverse = fmpq(1)/coefficient
                row = {key: value*inverse for key, value in row.items() if value}
                pivots[pivot] = row
                independent_indices.append(equation_index)
                break
            coefficient = row[pivot]
            add_scaled(row, known, -coefficient)
        else:
            dependencies.append(equation_index)

        if (equation_index+1) % 20 == 0:
            print(f"equations={equation_index+1}; rank={len(independent_indices)}; "
                  f"weight_differences={len(difference_rows)}", flush=True)

    ordered_differences = sorted(difference_rows)
    matrix = DomainMatrix.from_Matrix(sp.Matrix(ordered_differences)).convert_to(QQ)
    reduced, pivot_columns = matrix.rref()
    rank = len(pivot_columns)
    reduced_matrix = reduced.to_Matrix()
    free_columns = [j for j in range(len(names)) if j not in pivot_columns]
    nullspace = []
    for free in free_columns:
        vector = [sp.S.Zero] * len(names)
        vector[free] = sp.S.One
        for row_index, pivot_column in enumerate(pivot_columns):
            vector[pivot_column] = -reduced_matrix[row_index, free]
        nullspace.append(vector)

    positive_weight = None
    positive_check = "not_assessed"
    if len(nullspace) == 1:
        vector = nullspace[0]
        if all(value > 0 for value in vector):
            positive_weight = vector
        elif all(value < 0 for value in vector):
            positive_weight = [-value for value in vector]
        positive_check = "complete_one_dimensional_nullspace_check"
    elif not nullspace:
        positive_check = "no_nonzero_weight_vector"

    if positive_weight is not None:
        lcm = sp.ilcm(*[sp.denom(value) for value in positive_weight])
        integer_weights = [int(value*lcm) for value in positive_weight]
        common = sp.igcd(*integer_weights)
        integer_weights = [value//common for value in integer_weights]
        all_homogeneous = True
        degree_histogram = defaultdict(int)
        for equation in data["equations"]:
            degrees = set()
            for monomial, _coefficient in equation["terms"]:
                exponents = dict(monomial)
                degrees.add(sum(integer_weights[name_index[name]]*int(power)
                                for name, power in exponents.items()))
            if len(degrees) != 1:
                all_homogeneous = False
                break
            degree_histogram[str(next(iter(degrees)))] += 1
    else:
        integer_weights = None
        all_homogeneous = None
        degree_histogram = {}

    result = {
        "scope": (
            "Exact Q-linear compression and exact positive-weight nullspace audit "
            "for the selected support-21 subsystem only. Neither is a Groebner "
            "exclusion certificate."
        ),
        "input": str(INPUT),
        "input_sha256": sha256(INPUT),
        "d19_checkpoint_sha256": data["source_sha256"],
        "branch": data["branch"],
        "selected_variables": list(names),
        "input_equations": len(data["equations"]),
        "input_sparse_terms": input_terms,
        "exact_Q_linear_rank": len(independent_indices),
        "dependent_equation_count": len(dependencies),
        "independent_equation_indices_zero_based": independent_indices,
        "weight_difference_row_count": len(ordered_differences),
        "weight_equation_rank_over_Q": rank,
        "weight_nullity": len(names)-rank,
        "rational_weight_nullspace_basis": [
            [str(value) for value in vector] for vector in nullspace
        ],
        "positive_weight_check": positive_check,
        "primitive_positive_integer_weights": integer_weights,
        "all_equations_weight_homogeneous": all_homogeneous,
        "weighted_degree_histogram": dict(degree_histogram),
        "elapsed_seconds": time.monotonic()-started,
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "input_equations", "input_sparse_terms", "exact_Q_linear_rank",
        "dependent_equation_count", "weight_difference_row_count",
        "weight_equation_rank_over_Q", "weight_nullity",
        "positive_weight_check", "primitive_positive_integer_weights",
        "all_equations_weight_homogeneous", "weighted_degree_histogram",
        "elapsed_seconds",
    )}, indent=2), flush=True)
    print(f"report={OUTPUT}", flush=True)


if __name__ == "__main__":
    main()
