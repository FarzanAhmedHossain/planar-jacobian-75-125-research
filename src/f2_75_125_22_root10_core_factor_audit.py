#!/usr/bin/env python3
"""Inspect/factor the small-support core of the reduced F3 root-10 system."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

import sympy as sp

from f2_75_125_22_root10_block_elimination import (
    N_PARAMETERS,
    field_expr,
    parse_source_polynomial,
)
from f2_75_125_22_root10_support_profile import constraint_array


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("state", type=Path)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()

    certificate = json.loads(args.certificate.read_text(encoding="utf-8"))
    zero_names = certificate["forced_zero_parameters_after_third_block"]
    zero_indices = {int(name[1:]) for name in zero_names}
    core_indices = (3, 10, 11, 12, 13, 14, 15, 16, 17)
    symbols = sp.symbols(" ".join(f"p{i}" for i in core_indices))
    symbol_for = dict(zip(core_indices, symbols))
    sqrt5 = sp.sqrt(5)

    support_counts = Counter()
    factor_report = []
    included = 0
    for tag, encoded in constraint_array(args.state):
        poly = parse_source_polynomial(encoded)
        reduced = {}
        for exponent, coefficient in poly.items():
            if any(exponent[i] for i in zero_indices):
                continue
            reduced[exponent] = coefficient
        if not reduced:
            continue
        support = tuple(i for i in range(N_PARAMETERS)
                        if any(exponent[i] for exponent in reduced))
        if not set(support) <= set(core_indices):
            continue
        included += 1
        support_counts[len(support)] += 1
        expression = sp.S.Zero
        for exponent, coefficient in reduced.items():
            monomial = sp.prod(symbol_for[i] ** exponent[i]
                               for i in support if exponent[i])
            expression += field_expr(coefficient, sqrt5) * monomial
        expression = sp.expand(expression)
        factors = sp.factor_list(expression, *symbols, extension=sqrt5)
        original_degree = sp.Poly(expression, *symbols, extension=sqrt5).total_degree()
        factor_degrees = [
            sp.Poly(factor, *symbols, extension=sqrt5).total_degree()
            for factor, _multiplicity in factors[1]
        ]
        nontrivial = len(factor_degrees) > 1 or any(
            degree < original_degree for degree in factor_degrees
        )
        if nontrivial:
            factor_report.append({
                "tag": tag,
                "support": [f"p{i}" for i in support],
                "factorization": str(factors),
            })

    result = {
        "scope": "factor audit of constraints whose reduced support lies in p3,p10,...,p17 on one F3 root-10 necessary system",
        "zero_substitution_certificate": str(args.certificate),
        "variables": [f"p{i}" for i in core_indices],
        "included_equations": included,
        "support_size_histogram": dict(sorted(support_counts.items())),
        "nontrivially_factored_count": len(factor_report),
        "nontrivial_factorizations": factor_report,
        "warning": "Factorization and local core inspection do not alone exclude the full 2+2 branch.",
    }
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
