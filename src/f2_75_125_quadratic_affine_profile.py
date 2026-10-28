"""Audit affine necessary equations over Q(sqrt(5)) in a sparse F2 checkpoint.

The profiler streams no conclusions from nonlinear constraints.  If the affine
subsystem is inconsistent, it records an exact linear combination equal to 1,
which is already a unit-ideal certificate for the checkpoint's component.  If
it is consistent, the rank and pivot equations are diagnostic only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path


Elt = tuple[Fraction, Fraction]  # a + b*sqrt(5)
ZERO: Elt = (Fraction(0), Fraction(0))
ONE: Elt = (Fraction(1), Fraction(0))


def add(a: Elt, b: Elt) -> Elt:
    return a[0] + b[0], a[1] + b[1]


def sub(a: Elt, b: Elt) -> Elt:
    return a[0] - b[0], a[1] - b[1]


def mul(a: Elt, b: Elt) -> Elt:
    return a[0] * b[0] + 5 * a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def inv(a: Elt) -> Elt:
    norm = a[0] * a[0] - 5 * a[1] * a[1]
    if not norm:
        raise ZeroDivisionError("zero element in Q(sqrt(5))")
    return a[0] / norm, -a[1] / norm


def iszero(a: Elt) -> bool:
    return not a[0] and not a[1]


def from_json(values: list[str]) -> Elt:
    if not values:
        raise ValueError("missing coefficient")
    if len(values) > 2:
        raise ValueError("coefficient is not represented in Q(sqrt(5))")
    # AlgebraicField.to_list() stores coefficients in descending power order:
    # [b, a] represents a + b*sqrt(5), not [a, b].
    if len(values) == 1:
        return Fraction(values[0]), Fraction(0)
    return Fraction(values[1]), Fraction(values[0])


def show(a: Elt) -> str:
    rational, radical = a
    if not radical:
        return str(rational)
    radical_text = f"{radical}*sqrt(5)"
    if not rational:
        return radical_text
    sign = "+" if radical > 0 else "-"
    return f"{rational} {sign} {abs(radical)}*sqrt(5)"


def rank_affine(rows: list[list[Elt]]) -> tuple[int, list[dict[int, Elt]], dict | None]:
    """Row-reduce affine equations, retaining exact source-row combinations."""
    basis: dict[int, tuple[list[Elt], dict[int, Elt]]] = {}
    for source_index, source in enumerate(rows):
        row = list(source)
        witness: dict[int, Elt] = {source_index: ONE}

        for pivot in sorted(basis):
            factor = row[pivot + 1]
            if iszero(factor):
                continue
            base_row, base_witness = basis[pivot]
            row = [sub(value, mul(factor, base_value))
                   for value, base_value in zip(row, base_row)]
            for index, coefficient in base_witness.items():
                updated = sub(witness.get(index, ZERO), mul(factor, coefficient))
                if iszero(updated):
                    witness.pop(index, None)
                else:
                    witness[index] = updated

        pivot = next((j for j in range(len(row) - 1) if not iszero(row[j + 1])), None)
        if pivot is None:
            if not iszero(row[0]):
                scale = inv(row[0])
                certificate = {index: mul(scale, coefficient)
                               for index, coefficient in witness.items()}
                return len(basis), [], {"source_row_combination": certificate}
            continue

        scale = inv(row[pivot + 1])
        row = [mul(value, scale) for value in row]
        witness = {index: mul(scale, coefficient)
                   for index, coefficient in witness.items()}
        basis[pivot] = (row, witness)

    pivots = []
    for pivot in sorted(basis):
        row, _witness = basis[pivot]
        pivots.append({"pivot_column": pivot, "row": [show(value) for value in row]})
    return len(basis), pivots, None


def profile(source: Path) -> dict:
    raw = source.read_bytes()
    state = json.loads(raw)
    if state.get("format") != "f2-sparse-ring-v1":
        raise ValueError("unexpected checkpoint format")
    if not state.get("quadratic"):
        raise ValueError("expected a quadratic coefficient field")

    names = state["generators"]
    nparameters = int(state["nparameters"])
    parameter_names = [f"p{i}" for i in range(nparameters)]
    parameter_indices = {names.index(name): i for i, name in enumerate(parameter_names)}
    term_count = 0
    degree_histogram: Counter[int] = Counter()
    affine_rows: list[list[Elt]] = []
    affine_tags: list[str] = []
    affine_degrees: list[dict] = []
    constants = []
    zero_constraints = 0
    mixed_degree_constraints = 0
    max_degree_histogram: Counter[int] = Counter()
    profiles = []

    for tag, polynomial in state.get("constraints", []):
        if not polynomial:
            zero_constraints += 1
            continue
        row = [ZERO for _ in range(nparameters + 1)]
        histogram: Counter[int] = Counter()
        max_degree = 0
        affine = True
        for sparse, values in polynomial:
            term_count += 1
            coefficient = from_json(values)
            degree = 0
            p_variable = None
            t_power = 0
            other_power = 0
            for generator_index, power in sparse:
                if generator_index == 0:
                    t_power += power
                elif generator_index in parameter_indices:
                    degree += power
                    if p_variable is None:
                        p_variable = parameter_indices[generator_index]
                    elif power:
                        p_variable = -1
                else:
                    other_power += power
            histogram[degree] += 1
            degree_histogram[degree] += 1
            max_degree = max(max_degree, degree)
            if t_power or other_power or degree > 1:
                affine = False
            elif degree == 0:
                row[0] = add(row[0], coefficient)
            elif p_variable is None or p_variable < 0:
                affine = False
            else:
                row[p_variable + 1] = add(row[p_variable + 1], coefficient)

        if len(histogram) > 1:
            mixed_degree_constraints += 1
        max_degree_histogram[max_degree] += 1
        if len(polynomial) == 1 and not polynomial[0][0]:
            constants.append({"tag": tag, "value": show(from_json(polynomial[0][1]))})
        if affine:
            affine_rows.append(row)
            affine_tags.append(tag)
            affine_degrees.append({"tag": tag, "max_parameter_degree": max_degree})
        profiles.append({
            "tag": tag,
            "terms": len(polynomial),
            "max_parameter_degree": max_degree,
            "degree_histogram": {str(k): v for k, v in sorted(histogram.items())},
            "affine_parameter_equation": affine,
        })

    rank, pivots, inconsistency = rank_affine(affine_rows)
    certificate = None
    if inconsistency is not None:
        weights = inconsistency["source_row_combination"]
        combined = [ZERO for _ in range(nparameters + 1)]
        terms = []
        for row_index, coefficient in sorted(weights.items()):
            combined = [add(value, mul(coefficient, term))
                        for value, term in zip(combined, affine_rows[row_index])]
            terms.append({"tag": affine_tags[row_index], "coefficient": show(coefficient)})
        if combined[0] != ONE or any(not iszero(value) for value in combined[1:]):
            raise AssertionError("reconstructed affine combination is not exactly 1")
        certificate = {
            "identity": "1 = sum(coefficient_i * constraint_i)",
            "terms": terms,
            "exact_replay": True,
        }

    return {
        "source_checkpoint": str(source),
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "branch": state["branch"],
        "quadratic_field": "Q(sqrt(5))",
        "checkpoint_depth": len(state["p"]) - 1,
        "carried_parameter_count": nparameters,
        "constraint_count": len(profiles),
        "constraint_term_count": term_count,
        "parameter_degree_histogram": {str(k): v for k, v in sorted(degree_histogram.items())},
        "maximum_degree_histogram": {str(k): v for k, v in sorted(max_degree_histogram.items())},
        "mixed_degree_constraint_count": mixed_degree_constraints,
        "identically_zero_constraint_count": zero_constraints,
        "nonzero_constant_obstructions": constants,
        "affine_constraint_count": len(affine_rows),
        "affine_coefficient_rank_over_Qsqrt5": rank,
        "affine_system_inconsistent": certificate is not None,
        "affine_pivots": pivots,
        "affine_unit_certificate": certificate,
        "affine_constraint_profiles": affine_degrees,
        "constraint_profiles": profiles,
        "scope": (
            "If the affine subsystem is inconsistent, the saved identity proves the "
            "necessary checkpoint ideal is the unit ideal. Otherwise this is an exact "
            "profile only; nonlinear constraints remain to be decided."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("receipt", type=Path)
    args = parser.parse_args()
    if args.receipt.exists():
        raise FileExistsError(f"refusing to overwrite {args.receipt}")
    result = profile(args.checkpoint)
    args.receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "source_sha256", "branch", "quadratic_field", "checkpoint_depth",
        "carried_parameter_count", "constraint_count", "constraint_term_count",
        "parameter_degree_histogram", "affine_constraint_count",
        "affine_coefficient_rank_over_Qsqrt5", "affine_system_inconsistent",
        "affine_unit_certificate", "receipt",
    ) if key in result}, indent=2))
    print(f"receipt={args.receipt}")


if __name__ == "__main__":
    main()
