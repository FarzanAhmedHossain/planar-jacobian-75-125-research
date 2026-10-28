#!/usr/bin/env python3
"""Find exact positive weight gradings on a reduced small-support core."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import lcm
from pathlib import Path

import sympy as sp

from f2_75_125_22_root10_block_elimination import N_PARAMETERS
from f2_75_125_22_root10_support_profile import constraint_array


def primitive_integer_vector(vector):
    denominators = [int(value.q) for value in vector]
    scale = 1
    for denominator in denominators:
        scale = lcm(scale, denominator)
    integers = [int(value * scale) for value in vector]
    common = 0
    for value in integers:
        common = int(sp.igcd(common, abs(value)))
    if common:
        integers = [value // common for value in integers]
    first_nonzero = next((value for value in integers if value), 1)
    if first_nonzero < 0:
        integers = [-value for value in integers]
    return integers


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("state", type=Path)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    certificate = json.loads(args.certificate.read_text(encoding="utf-8"))
    zeros = {int(name[1:]) for name in certificate["forced_zero_parameters_after_third_block"]}
    core = (3, 10, 11, 12, 13, 14, 15, 16, 17)
    differences = []
    included = 0
    for tag, encoded in constraint_array(args.state):
        reduced_coefficients = {}
        for sparse, values in encoded:
            exp_map = {int(i): int(e) for i, e in sparse if e}
            if 0 in exp_map:
                raise AssertionError("root-10 constraints must not retain t")
            exponent = tuple(exp_map.get(i + 1, 0) for i in range(N_PARAMETERS))
            if any(exponent[i] for i in zeros):
                continue
            if len(values) == 1:
                coefficient = (Fraction(values[0]), Fraction(0))
            elif len(values) == 2:
                coefficient = (Fraction(values[0]), Fraction(values[1]))
            else:
                raise AssertionError("unexpected Q(sqrt(5)) coefficient encoding")
            old = reduced_coefficients.get(exponent, (Fraction(0), Fraction(0)))
            reduced_coefficients[exponent] = (old[0] + coefficient[0], old[1] + coefficient[1])
        reduced = [exponent for exponent, coefficient in reduced_coefficients.items()
                   if coefficient != (0, 0)]
        if not reduced:
            continue
        support = {i for i in range(N_PARAMETERS)
                   if any(exponent[i] for exponent in reduced)}
        if not support or not support <= set(core):
            continue
        included += 1
        vectors = [sp.Matrix([exponent[i] for i in core]) for exponent in reduced]
        base = vectors[0]
        for vector in vectors[1:]:
            differences.append(list(vector - base))
    matrix = sp.Matrix(differences) if differences else sp.zeros(0, len(core))
    nullspace = matrix.nullspace()
    integer_basis = [primitive_integer_vector(vector) for vector in nullspace]
    positive_basis = [vector for vector in integer_basis if all(value > 0 for value in vector)]
    result = {
        "scope": "exact exponent-difference grading audit on the small-support core of one F3 root-10 necessary system",
        "variables": [f"p{i}" for i in core],
        "included_equation_count": included,
        "weight_difference_matrix_shape": list(matrix.shape),
        "grading_nullity": len(nullspace),
        "integer_nullspace_basis": integer_basis,
        "strictly_positive_basis_vectors": positive_basis,
        "interpretation": "Every included polynomial is homogeneous for each vector in the nullspace; a positive vector supplies a weighted-projective grading of this core.",
        "warning": "A grading describes the local necessary system and does not by itself exclude the full 2+2 branch.",
    }
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
