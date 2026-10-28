#!/usr/bin/env python3
"""Exact axis-gcd screen of the 673-equation depth-14 2+1+1 checkpoint."""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict, Counter
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_sparse_ring_lift_211_terminal_roots14_extension_state.json"
OUT = HERE / "f2_75_125_211_roots14_checkpoint_axis_screen.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def trim(poly: list[Fraction]) -> list[Fraction]:
    while poly and poly[-1] == 0:
        poly.pop()
    return poly


def monic(poly: list[Fraction]) -> list[Fraction]:
    poly = trim(poly[:])
    if not poly:
        return []
    lead = poly[-1]
    return [coefficient / lead for coefficient in poly]


def divrem(a: list[Fraction], b: list[Fraction]):
    remainder = trim(a[:])
    b = trim(b[:])
    if not b:
        raise ZeroDivisionError
    quotient = [Fraction(0)] * max(0, len(remainder) - len(b) + 1)
    while len(remainder) >= len(b) and remainder:
        shift = len(remainder) - len(b)
        factor = remainder[-1] / b[-1]
        quotient[shift] += factor
        for i, coefficient in enumerate(b):
            remainder[shift + i] -= factor * coefficient
        trim(remainder)
    return trim(quotient), remainder


def gcd_poly(a: list[Fraction], b: list[Fraction]) -> list[Fraction]:
    a, b = trim(a[:]), trim(b[:])
    while b:
        _, remainder = divrem(a, b)
        a, b = b, remainder
    return monic(a)


def encode(poly: list[Fraction]) -> list[str]:
    return [str(coefficient) for coefficient in poly]


def main() -> None:
    data = json.loads(SOURCE.read_text(encoding="utf-8"))
    if data.get("branch") != "2+1+1" or data.get("nparameters") != 108:
        raise AssertionError("unexpected checkpoint metadata")
    generators = data["generators"]
    parameter_names = [name for name in generators
                       if name.startswith("p") and name[1:].isdigit()]
    positions = {name: index for index, name in enumerate(generators)}
    axis_terms = defaultdict(lambda: defaultdict(lambda: defaultdict(Fraction)))
    selected_subspace = {"p46", "p47", "p48", "p49", "p50"}
    selected_restrictions = defaultdict(lambda: defaultdict(Fraction))
    row_count = 0
    term_count = 0
    max_degrees = Counter()
    constant_terms = []

    for row, (tag, terms) in enumerate(data["constraints"]):
        row_count += 1
        degree_set = set()
        for sparse, values in terms:
            term_count += 1
            if len(values) != 1:
                raise ValueError("expected exact rational coefficients")
            coefficient = Fraction(values[0])
            degree_set.add(sum(power for _, power in sparse))
            if not sparse:
                if coefficient:
                    constant_terms.append((row, tag, str(coefficient)))
                continue
            support_names = {generators[index] for index, _ in sparse}
            if support_names and support_names.issubset(selected_subspace):
                exponent_tuple = tuple(
                    next((power for index, power in sparse
                          if generators[index] == name), 0)
                    for name in sorted(selected_subspace)
                )
                selected_restrictions[row][exponent_tuple] += coefficient
            if len(sparse) != 1:
                continue
            generator_index, exponent = sparse[0]
            name = generators[generator_index]
            if name in parameter_names:
                axis_terms[name][row][exponent] += coefficient
        if degree_set:
            max_degrees[max(degree_set)] += 1

    if row_count != 673 or term_count != 327868:
        raise AssertionError((row_count, term_count))
    selected_rows_nonzero = [
        row for row, terms in selected_restrictions.items()
        if any(coefficient for coefficient in terms.values())
    ]

    axis_results = {}
    for name in parameter_names:
        common = None
        nonzero_restrictions = 0
        for row in range(row_count):
            by_degree = axis_terms[name].get(row, {})
            if not by_degree:
                continue
            polynomial = [by_degree.get(power, Fraction(0))
                          for power in range(max(by_degree) + 1)]
            trim(polynomial)
            if not polynomial:
                continue
            nonzero_restrictions += 1
            common = polynomial if common is None else gcd_poly(common, polynomial)

        if common is None:
            axis_results[name] = {
                "status": "entire_axis_survives_checkpoint",
                "nonzero_restricted_equations": 0,
            }
            continue
        common = monic(common)
        nonzero_part = common[:]
        while len(nonzero_part) > 1 and nonzero_part[0] == 0:
            nonzero_part = nonzero_part[1:]
        axis_results[name] = {
            "status": ("possible_nonzero_axis_points" if len(nonzero_part) > 1
                       else "axis_excluded"),
            "nonzero_restricted_equations": nonzero_restrictions,
            "common_gcd_ascending_coefficients_Q": encode(common),
            "nonzero_root_factor_ascending_coefficients_Q": encode(nonzero_part),
        }

    result = {
        "source": SOURCE.name,
        "source_sha256": sha256_file(SOURCE),
        "branch": "2+1+1",
        "checkpoint_parameters": data["nparameters"],
        "parameter_axes_screened": len(parameter_names),
        "ideal_equations": row_count,
        "serialized_terms": term_count,
        "ordinary_max_degree_profile": dict(sorted((str(k), v) for k, v in max_degrees.items())),
        "nonzero_constant_terms": constant_terms,
        "selected_subspace": sorted(selected_subspace),
        "selected_subspace_nonzero_restricted_equation_count": len(selected_rows_nonzero),
        "selected_subspace_restricted_equation_examples": selected_rows_nonzero[:20],
        "selected_subspace_contained_in_zero_locus": not selected_rows_nonzero,
        "scope": "exact axis intersections of the saved depth-14 necessary ideal only; not a full branch exclusion",
        "axis_results": axis_results,
        "axis_status_counts": {
            status: sum(value["status"] == status for value in axis_results.values())
            for status in sorted({value["status"] for value in axis_results.values()})
        },
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "source_sha256": result["source_sha256"],
        "checkpoint_parameters": result["checkpoint_parameters"],
        "parameter_axes_screened": result["parameter_axes_screened"],
        "ideal_equations": row_count,
        "serialized_terms": term_count,
        "ordinary_max_degree_profile": result["ordinary_max_degree_profile"],
        "nonzero_constant_terms": len(constant_terms),
        "selected_subspace_nonzero_restricted_equation_count": len(selected_rows_nonzero),
        "selected_subspace_contained_in_zero_locus": not selected_rows_nonzero,
        "axis_status_counts": result["axis_status_counts"],
        "surviving_axes": [name for name, value in axis_results.items()
                           if value["status"] != "axis_excluded"],
        "receipt": OUT.name,
    }, indent=2))


if __name__ == "__main__":
    main()
