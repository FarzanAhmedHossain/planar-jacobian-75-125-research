#!/usr/bin/env python3
"""Find and verify the positive integer grading of the d18 support-12 ideal."""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from pathlib import Path

from sympy import Matrix


INPUT = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_input.json")
LINEAR = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_linear_reduce_verified.json")
EXPECTED_INPUT_SHA256 = "A9DF8405D43EAFFB98C724165C0C7A647C0952D7D8D1A8630CD26700898C6F56"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest().upper()


def main() -> None:
    input_hash = sha256(INPUT)
    if input_hash != EXPECTED_INPUT_SHA256:
        raise ValueError(f"support-12 input hash mismatch: {input_hash}")
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    linear = json.loads(LINEAR.read_text(encoding="utf-8"))
    if linear.get("input_sha256") != input_hash or not linear.get(
            "all_142_equations_span_exactly_4_Q_generators"):
        raise ValueError("exact linear reduction receipt is missing or stale")
    names = tuple(data["selected_variables"])
    selected = linear["independent_equation_indices"]
    if len(names) != 12 or len(selected) != 4:
        raise ValueError("unexpected support-12 generator metadata")

    differences = []
    for index in selected:
        exponents = []
        for monomial, _ in data["equations"][index]["terms"]:
            powers = dict(monomial)
            exponents.append(tuple(int(powers.get(name, 0)) for name in names))
        base = exponents[0]
        differences.extend(tuple(a - b for a, b in zip(exp, base))
                           for exp in exponents[1:])
    matrix = Matrix(differences)
    nullspace = matrix.nullspace()
    if len(nullspace) != 1:
        raise ValueError(f"expected a one-dimensional weight space, got {len(nullspace)}")
    vector = nullspace[0]
    denominator_lcm = math.lcm(*(int(value.q) for value in vector))
    integer_weights = [int(value * denominator_lcm) for value in vector]
    common_gcd = math.gcd(*integer_weights)
    weights = [value // common_gcd for value in integer_weights]
    if any(value <= 0 for value in weights):
        raise ValueError(f"weight vector is not positive: {weights}")

    def weighted_degrees(equation: dict) -> set[int]:
        result = set()
        for monomial, _ in equation["terms"]:
            powers = dict(monomial)
            result.add(sum(weights[i] * int(powers.get(name, 0))
                           for i, name in enumerate(names)))
        return result

    degrees_by_row = [sorted(weighted_degrees(row)) for row in data["equations"]]
    homogeneous = all(len(degrees) == 1 for degrees in degrees_by_row)
    degree_histogram = Counter(degrees[0] for degrees in degrees_by_row if len(degrees) == 1)
    result = {
        "input_sha256": input_hash,
        "linear_reduction_receipt_sha256": sha256(LINEAR),
        "branch": data["branch"],
        "selected_variables": list(names),
        "selected_generator_indices": selected,
        "weight_matrix_shape": list(matrix.shape),
        "weight_matrix_rank": matrix.rank(),
        "weight_space_dimension": len(nullspace),
        "primitive_positive_integer_weights": weights,
        "all_142_input_equations_weight_homogeneous": homogeneous,
        "weighted_degree_histogram": dict(sorted(degree_histogram.items())),
        "selected_generator_weight_degrees": [degrees_by_row[i] for i in selected],
        "scope": "Exact grading of the support-12 necessary subsystem only; not a branch exclusion.",
    }
    output = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_weighted_homogeneity.json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({**result, "output": str(output)}, indent=2), flush=True)


if __name__ == "__main__":
    main()
