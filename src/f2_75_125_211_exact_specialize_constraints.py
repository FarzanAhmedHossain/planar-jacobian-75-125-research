#!/usr/bin/env python3
"""Propagate exact rational assignments forced by univariate constraints.

This is an exact necessary-condition reducer for saved sparse F2 checkpoints.
At each round it substitutes forced values, detects nonzero constants, and
computes the gcd of all surviving univariate equations in each parameter.
Only a linear gcd over Q is converted into a new parameter assignment.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

from f2_75_125_211_d13_univariate_gcd import (
    add,
    encode_poly,
    mul,
    scale,
    xgcd,
)
from f2_75_125_quadratic_affine_profile_stream import (
    iter_json_array,
    top_metadata_and_constraints,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def specialize(source: Path, assignments: dict[str, Fraction]) -> dict:
    text = source.read_text(encoding="utf-8")
    metadata, constraint_pos = top_metadata_and_constraints(text)
    if metadata.get("format") != "f2-sparse-ring-v1":
        raise ValueError("unexpected checkpoint format")
    if metadata.get("branch") != "2+1+1" or metadata.get("quadratic"):
        raise ValueError("expected a rational 2+1+1 checkpoint")
    names = metadata["generators"]
    assigned_generators = {names.index(name): value for name, value in assignments.items()}

    univariate: dict[int, list[tuple[str, list[Fraction]]]] = defaultdict(list)
    zero_rows = 0
    constant_rows: list[dict] = []
    constraint_count = 0
    term_count_before = 0
    term_count_after = 0
    occurrences = [0] * int(metadata["nparameters"])
    parent = list(range(len(occurrences)))
    component_rank = [0] * len(parent)
    auxiliary_generators: set[int] = set()
    support_histogram: Counter[int] = Counter()
    supports_by_size: dict[int, Counter[tuple[int, ...]]] = defaultdict(Counter)
    term_count_histogram_by_support: dict[int, Counter[int]] = defaultdict(Counter)
    degree_histogram: Counter[int] = Counter()
    largest_polynomial_terms = 0

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(left: int, right: int) -> None:
        root_left, root_right = find(left), find(right)
        if root_left == root_right:
            return
        if component_rank[root_left] < component_rank[root_right]:
            root_left, root_right = root_right, root_left
        parent[root_right] = root_left
        if component_rank[root_left] == component_rank[root_right]:
            component_rank[root_left] += 1

    for tag, terms in iter_json_array(text, constraint_pos):
        constraint_count += 1
        polynomial: dict[tuple[tuple[int, int], ...], Fraction] = {}
        for sparse, values in terms:
            term_count_before += 1
            if len(values) != 1:
                raise ValueError("expected rational coefficients")
            coefficient = Fraction(values[0])
            remaining = []
            for generator, exponent in sparse:
                if generator in assigned_generators:
                    coefficient *= assigned_generators[generator] ** exponent
                else:
                    remaining.append((generator, exponent))
            if coefficient:
                monomial = tuple(remaining)
                polynomial[monomial] = polynomial.get(monomial, Fraction(0)) + coefficient
        polynomial = {m: c for m, c in polynomial.items() if c}
        term_count_after += len(polynomial)
        if not polynomial:
            zero_rows += 1
            support_histogram[0] += 1
            continue
        largest_polynomial_terms = max(largest_polynomial_terms, len(polynomial))
        total_degree = max(sum(power for _, power in monomial) for monomial in polynomial)
        degree_histogram[total_degree] += 1
        support = {generator for monomial in polynomial for generator, _ in monomial}
        parameter_support = []
        for generator in support:
            if 1 <= generator <= len(occurrences):
                parameter_support.append(generator - 1)
            else:
                auxiliary_generators.add(generator)
        for variable in parameter_support:
            occurrences[variable] += 1
        if parameter_support:
            for variable in parameter_support[1:]:
                union(parameter_support[0], variable)
        parameter_support_tuple = tuple(sorted(parameter_support))
        support_histogram[len(parameter_support_tuple)] += 1
        supports_by_size[len(parameter_support_tuple)][parameter_support_tuple] += 1
        term_count_histogram_by_support[len(parameter_support_tuple)][len(polynomial)] += 1
        if not support:
            constant_rows.append({"tag": tag, "constant": str(polynomial[()] )})
            continue
        if len(support) == 1:
            generator = next(iter(support))
            coefficients: dict[int, Fraction] = defaultdict(Fraction)
            for monomial, coefficient in polynomial.items():
                if any(index != generator for index, _ in monomial):
                    raise AssertionError("univariate support accounting failed")
                exponent = sum(power for _, power in monomial)
                coefficients[exponent] += coefficient
            equation = [Fraction(0)] * (max(coefficients) + 1)
            for exponent, coefficient in coefficients.items():
                equation[exponent] = coefficient
            while len(equation) > 1 and not equation[-1]:
                equation.pop()
            univariate[generator].append((tag, equation))

    gcd_results = []
    for generator, equations in sorted(univariate.items()):
        gcd_poly = equations[0][1]
        bezout = [[Fraction(1)]]
        for _tag, equation in equations[1:]:
            gcd_poly, left, right = xgcd(gcd_poly, equation)
            bezout = [mul(left, value) for value in bezout]
            bezout.append(right)
        if len(gcd_poly) == 1:
            normalization = 1 / gcd_poly[0]
            gcd_poly = [Fraction(1)]
            bezout = [scale(value, normalization) for value in bezout]
        identity = [Fraction(0)]
        for coefficient, (_tag, equation) in zip(bezout, equations):
            identity = add(identity, mul(coefficient, equation))
        if identity != gcd_poly:
            raise AssertionError(f"Bezout identity failed for {names[generator]}")
        gcd_results.append({
            "variable": names[generator],
            "generator_index": generator,
            "equation_count": len(equations),
            "gcd_coefficients_low_to_high": encode_poly(gcd_poly),
            "gcd_degree": len(gcd_poly) - 1,
            "equation_tags": [tag for tag, _ in equations],
            "bezout_coefficients_low_to_high": [encode_poly(value) for value in bezout]
                if len(gcd_poly) <= 1 else None,
        })

    components: dict[int, list[str]] = defaultdict(list)
    for index, count in enumerate(occurrences):
        if count:
            components[find(index)].append(f"p{index}")
    component_sets = sorted(
        (values for values in components.values()), key=lambda values: (-len(values), values)
    )

    return {
        "source": str(source),
        "source_sha256": sha256(source),
        "branch": metadata["branch"],
        "nparameters": int(metadata["nparameters"]),
        "assignments": {name: str(value) for name, value in sorted(assignments.items())},
        "constraint_count": constraint_count,
        "input_term_count": term_count_before,
        "specialized_term_count": term_count_after,
        "rows_becoming_zero": zero_rows,
        "nonzero_constant_rows": constant_rows,
        "univariate_gcds": gcd_results,
        "occurring_parameter_count": sum(bool(value) for value in occurrences),
        "parameter_incidence_component_sizes": [len(values) for values in component_sets],
        "parameter_incidence_components": component_sets,
        "auxiliary_generators_in_constraints": sorted(auxiliary_generators),
        "support_size_histogram": {
            str(size): count for size, count in sorted(support_histogram.items())
        },
        "distinct_support_sets_by_size": {
            str(size): len(supports) for size, supports in sorted(supports_by_size.items())
        },
        "term_count_histogram_by_support": {
            str(size): {str(terms): count for terms, count in sorted(histogram.items())}
            for size, histogram in sorted(term_count_histogram_by_support.items())
        },
        "most_repeated_supports": {
            str(size): [
                {"variables": [f"p{index}" for index in support], "constraints": count}
                for support, count in counts.most_common(5)
            ]
            for size, counts in sorted(supports_by_size.items())
        },
        "polynomial_total_degree_histogram": {
            str(degree): count for degree, count in sorted(degree_histogram.items())
        },
        "largest_specialized_constraint_term_count": largest_polynomial_terms,
        "status": "inconsistent" if constant_rows or any(r["gcd_degree"] == 0
                                                            for r in gcd_results)
                  else "necessary_constraints_reduced",
        "scope": (
            "Exact specialization of retained necessary equations only; a contradiction "
            "excludes the checkpoint chart, while absence of a contradiction is not survival."
        ),
    }


def parse_assignment(value: str) -> tuple[str, Fraction]:
    name, separator, scalar = value.partition("=")
    if not separator or not name.startswith("p") or not name[1:].isdigit():
        raise argparse.ArgumentTypeError("assignment must look like p79=0")
    try:
        return name, Fraction(scalar)
    except (ValueError, ZeroDivisionError) as error:
        raise argparse.ArgumentTypeError(f"invalid rational assignment {value!r}") from error


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("report", type=Path)
    parser.add_argument("--assign", action="append", type=parse_assignment, default=[])
    args = parser.parse_args()
    assignments = dict(args.assign)
    result = specialize(args.source, assignments)
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        key: result[key] for key in (
            "source_sha256", "assignments", "constraint_count", "input_term_count",
            "specialized_term_count", "rows_becoming_zero", "nonzero_constant_rows",
            "univariate_gcds", "occurring_parameter_count",
            "parameter_incidence_component_sizes", "parameter_incidence_components",
            "auxiliary_generators_in_constraints", "support_size_histogram",
            "distinct_support_sets_by_size", "term_count_histogram_by_support",
            "most_repeated_supports",
            "polynomial_total_degree_histogram", "largest_specialized_constraint_term_count",
            "status", "scope"
        )
    }, indent=2))
    print(f"report={args.report}")


if __name__ == "__main__":
    main()
