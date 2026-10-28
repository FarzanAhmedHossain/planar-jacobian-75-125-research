#!/usr/bin/env python3
"""Certify the weighted-homogeneous core split and exact row-space bases."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix

from f2_75_125_22_root10_block_elimination import N_PARAMETERS
from f2_75_125_22_root10_support_profile import constraint_array


CORE = (3, 10, 11, 12, 13, 14, 15, 16)
P17 = 17
WEIGHTS = (1, 2, 2, 2, 2, 2, 2, 2)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def coefficient_pair(values):
    if len(values) == 1:
        return Fraction(values[0]), Fraction(0)
    if len(values) == 2:
        return Fraction(values[0]), Fraction(values[1])
    raise ValueError(f"unexpected Q(sqrt(5)) coefficient tuple: {values!r}")


def matrix_rank_and_rows(rows, width, field):
    if not rows:
        return 0, []
    matrix = DomainMatrix(rows, (len(rows), width), field)
    _, pivots = matrix.rref()
    transpose = [[rows[i][j] for i in range(len(rows))] for j in range(width)]
    _, row_indices = DomainMatrix(transpose, (width, len(rows)), field).rref()
    return len(pivots), list(row_indices)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("state", type=Path)
    parser.add_argument("zero_certificate", type=Path)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    certificate = json.loads(args.zero_certificate.read_text(encoding="utf-8"))
    zero_indices = {
        int(name[1:])
        for name in certificate["forced_zero_parameters_after_third_block"]
    }
    sqrt5 = sp.sqrt(5)
    field = QQ.algebraic_field(sqrt5)
    weights_by_degree = {"p17_exponent_per_equation": defaultdict(int)}
    branch_polys = {"p17_zero": defaultdict(list), "p17_nonzero": defaultdict(list)}
    equation_tags = {"p17_zero": {}, "p17_nonzero": {}}
    included = 0
    for tag, encoded in constraint_array(args.state):
        reduced = {}
        for sparse, values in encoded:
            exp_map = {int(i): int(e) for i, e in sparse if e}
            if 0 in exp_map:
                raise AssertionError("root-10 source unexpectedly retains t")
            exponent = tuple(exp_map.get(i + 1, 0) for i in range(N_PARAMETERS))
            if any(exponent[i] for i in zero_indices):
                continue
            pair = coefficient_pair(values)
            old = reduced.get(exponent, (Fraction(0), Fraction(0)))
            reduced[exponent] = (old[0] + pair[0], old[1] + pair[1])
        reduced = {e: c for e, c in reduced.items() if c != (0, 0)}
        if not reduced:
            continue
        support = {
            i for i in range(N_PARAMETERS)
            if any(exponent[i] for exponent in reduced)
        }
        if not support <= set(CORE) | {P17}:
            continue
        if not support:
            continue
        included += 1
        p17_exponents = {exponent[P17] for exponent in reduced}
        if len(p17_exponents) != 1:
            raise AssertionError(f"p17 exponent varies within {tag}: {p17_exponents}")
        p17_exponent = next(iter(p17_exponents))
        weights_by_degree["p17_exponent_per_equation"][str(p17_exponent)] += 1
        core_poly = {}
        weighted_degrees = set()
        for exponent, coefficient in reduced.items():
            core_exp = tuple(exponent[i] for i in CORE)
            weighted_degree = sum(core_exp[j] * WEIGHTS[j] for j in range(len(CORE)))
            weighted_degrees.add(weighted_degree)
            core_poly[core_exp] = coefficient
        if len(weighted_degrees) != 1:
            raise AssertionError(f"core-weight degree is not constant in {tag}: {weighted_degrees}")
        weighted_degree = next(iter(weighted_degrees))
        weights_by_degree.setdefault("weighted_degree_histogram", defaultdict(int))
        weights_by_degree["weighted_degree_histogram"][str(weighted_degree)] += 1
        if p17_exponent == 0:
            branch = "p17_zero"
        elif p17_exponent == 1:
            branch = "p17_nonzero"
        else:
            raise AssertionError(f"unexpected common p17 exponent {p17_exponent} in {tag}")
        row = {}
        for exponent, (rational, radical) in core_poly.items():
            expression = sp.Rational(rational.numerator, rational.denominator)
            expression += sp.Rational(radical.numerator, radical.denominator) * sqrt5
            row[exponent] = field.convert(expression)
        branch_polys[branch][weighted_degree].append(row)
        equation_tags[branch].setdefault(weighted_degree, []).append(tag)

    span_report = {}
    for branch in ("p17_zero", "p17_nonzero"):
        selected_tags = []
        degree_report = {}
        for degree in sorted(branch_polys[branch]):
            polynomials = branch_polys[branch][degree]
            monomials = sorted({monomial for polynomial in polynomials for monomial in polynomial})
            monomial_index = {monomial: i for i, monomial in enumerate(monomials)}
            rows = []
            for polynomial in polynomials:
                row = [field.zero] * len(monomials)
                for monomial, coefficient in polynomial.items():
                    row[monomial_index[monomial]] = coefficient
                rows.append(row)
            rank, independent = matrix_rank_and_rows(rows, len(monomials), field)
            tags = equation_tags[branch][degree]
            selected = [tags[i] for i in independent]
            selected_tags.extend(selected)
            degree_report[str(degree)] = {
                "equation_count": len(polynomials),
                "monomial_count": len(monomials),
                "coefficient_matrix_rank": rank,
                "independent_source_tags": selected,
            }
        span_report[branch] = {
            "weighted_degree_slices": degree_report,
            "independent_equation_count": len(selected_tags),
            "independent_source_tags": selected_tags,
            "ideal_interpretation": (
                "p17=0 branch uses core rows only"
                if branch == "p17_zero"
                else "on p17!=0, divide p17 from each row and add the resulting core equations"
            ),
        }

    result = {
        "source": str(args.state),
        "source_sha256": sha256(args.state),
        "zero_certificate": str(args.zero_certificate),
        "zero_certificate_sha256": sha256(args.zero_certificate),
        "scope": "exact weighted-homogeneous split and coefficient-row-space reduction inside one F3 root-10 necessary system",
        "coefficient_field": "Q(sqrt(5))",
        "core_variables": [f"p{i}" for i in CORE],
        "core_weights": list(WEIGHTS),
        "p17_variable": "p17",
        "included_equation_count": included,
        "p17_exponent_histogram": dict(sorted(weights_by_degree["p17_exponent_per_equation"].items())),
        "weighted_degree_histogram": dict(sorted(weights_by_degree["weighted_degree_histogram"].items())),
        "branches": span_report,
        "warning": "The two p17 cases and row-space bases are local to one necessary root-10 sheet; they do not exclude the full 2+2 F2 branch.",
    }
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
