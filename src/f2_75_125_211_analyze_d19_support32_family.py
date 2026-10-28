"""Exact Q-row compression and positive-weight audit for support 32."""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from pathlib import Path
import time

import sympy as sp
from flint import fmpq


INPUT = Path("research/f2_75_125_211_d19_support32_family_input_20261025.json")
OUTPUT = Path("research/f2_75_125_211_d19_support32_analysis_20261025.json")


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


def main(input_path: Path = INPUT, output_path: Path = OUTPUT) -> None:
    started = time.monotonic()
    data = json.loads(input_path.read_text(encoding="utf-8"))
    family_label = input_path.name.split("_support", 1)[-1].split("_family", 1)[0]
    names = tuple(sorted({
        name for equation in data["equations"]
        for name in equation["reduced_support"]
    }))
    if data["equation_count"] != len(data["equations"]) or not names:
        raise ValueError("support-32 input metadata mismatch or empty family")
    name_index = {name: i for i, name in enumerate(names)}

    pivots = {}
    independent_indices = []
    dependencies = []
    input_terms = 0
    difference_rows = set()
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
            add_scaled(row, known, -row[pivot])
        else:
            dependencies.append(equation_index)
        if (equation_index+1) % 20 == 0:
            print(f"equations={equation_index+1}; rank={len(independent_indices)}; "
                  f"weight_differences={len(difference_rows)}", flush=True)

    # Sparse exact elimination in at most 28 coordinates avoids materializing
    # a potentially large dense SymPy matrix of all monomial differences.
    weight_pivots = {}
    for delta in difference_rows:
        row = [fmpq(value) for value in delta]
        while True:
            pivot = next((i for i in range(len(row)-1, -1, -1) if row[i]), None)
            if pivot is None:
                break
            known = weight_pivots.get(pivot)
            if known is None:
                inverse = fmpq(1)/row[pivot]
                weight_pivots[pivot] = [value*inverse for value in row]
                break
            coefficient = row[pivot]
            row = [a-coefficient*b for a, b in zip(row, known)]

    free_columns = [j for j in range(len(names)) if j not in weight_pivots]
    nullspace = []
    for free in free_columns:
        vector = [fmpq(0)]*len(names)
        vector[free] = fmpq(1)
        for pivot in sorted(weight_pivots):
            row = weight_pivots[pivot]
            vector[pivot] = -sum((row[j]*vector[j] for j in range(pivot)), fmpq(0))
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

    integer_weights = None
    all_homogeneous = None
    degree_histogram = {}
    if positive_weight is not None:
        rational_weights = [sp.Rational(str(value)) for value in positive_weight]
        lcm = sp.ilcm(*[sp.denom(value) for value in rational_weights])
        integer_weights = [int(value*lcm) for value in rational_weights]
        common = sp.igcd(*integer_weights)
        integer_weights = [value//common for value in integer_weights]
        all_homogeneous = True
        histogram = defaultdict(int)
        for equation in data["equations"]:
            degrees = set()
            for monomial, _coefficient in equation["terms"]:
                exponents = dict(monomial)
                degrees.add(sum(integer_weights[name_index[name]]*int(power)
                                for name, power in exponents.items()))
            if len(degrees) != 1:
                all_homogeneous = False
                break
            histogram[str(next(iter(degrees)))] += 1
        degree_histogram = dict(histogram)

    result = {
        "scope": (
            f"Exact Q-linear compression and positive-weight nullspace audit "
            f"for the complete original support-{family_label} family after "
            "the certified d19 substitutions. This is not a Groebner "
            "exclusion certificate."
        ),
        "input": str(input_path),
        "input_sha256": sha256(input_path),
        "d19_checkpoint_sha256": data["source_sha256"],
        "branch": data["branch"],
        "selected_variables": list(names),
        "input_equations": len(data["equations"]),
        "input_sparse_terms": input_terms,
        "exact_Q_linear_rank": len(independent_indices),
        "dependent_equation_count": len(dependencies),
        "independent_equation_indices_zero_based": independent_indices,
        "weight_difference_row_count": len(difference_rows),
        "weight_equation_rank_over_Q": len(weight_pivots),
        "weight_nullity": len(names)-len(weight_pivots),
        "rational_weight_nullspace_basis": [
            [str(value) for value in vector] for vector in nullspace
        ],
        "positive_weight_check": positive_check,
        "primitive_positive_integer_weights": integer_weights,
        "all_equations_weight_homogeneous": all_homogeneous,
        "weighted_degree_histogram": degree_histogram,
        "elapsed_seconds": time.monotonic()-started,
    }
    output_path.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "input_equations", "input_sparse_terms", "exact_Q_linear_rank",
        "dependent_equation_count", "weight_difference_row_count",
        "weight_equation_rank_over_Q", "weight_nullity",
        "positive_weight_check", "primitive_positive_integer_weights",
        "all_equations_weight_homogeneous", "weighted_degree_histogram",
        "elapsed_seconds",
    )}, indent=2), flush=True)
    print(f"report={output_path}", flush=True)


if __name__ == "__main__":
    main()
