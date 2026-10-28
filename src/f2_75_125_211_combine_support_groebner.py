#!/usr/bin/env python3
"""Combine exact Groebner bases for overlapping checkpoint support slices."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.groebnertools import groebner as ring_groebner
from sympy.polys.rings import ring as polynomial_ring


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("report", type=Path)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--method", choices=("f5b", "buchberger"), default="f5b")
    parser.add_argument("--order", choices=("grevlex", "lex"), default="grevlex")
    args = parser.parse_args()

    sources = [json.loads(path.read_text(encoding="utf-8")) for path in args.inputs]
    hashes = {source["source_sha256"] for source in sources}
    if len(hashes) != 1:
        raise ValueError("input reports do not describe the same checkpoint")
    if any(source.get("unit_ideal") for source in sources):
        # A unit source basis already certifies the combined ideal is unit.
        variables = tuple(dict.fromkeys(
            name for source in sources for name in source["selected_variables"]
        ))
        result = {
            "source_sha256": next(iter(hashes)),
            "input_reports": [str(path) for path in args.inputs],
            "selected_variables": list(variables),
            "input_basis_lengths": [source["groebner_basis_length"]
                                    for source in sources],
            "groebner_basis": ["1"],
            "groebner_basis_length": 1,
            "unit_ideal": True,
            "zero_dimensional": True,
            "scope": "The sum of the selected exact support-slice ideals is unit.",
        }
        args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(result, indent=2), flush=True)
        return

    variables = tuple(dict.fromkeys(
        name for source in sources for name in source["selected_variables"]
    ))
    symbols = sp.symbols(" ".join(variables))
    ring, *_ = polynomial_ring(variables, QQ, order=args.order)
    input_polynomials = []
    for source in sources:
        for expression in source["groebner_basis"]:
            polynomial = sp.Poly(sp.sympify(expression), *symbols, domain=QQ)
            input_polynomials.append(
                ring.from_dict(dict(polynomial.terms()))
            )

    basis = ring_groebner(input_polynomials, ring, method=args.method)
    basis_strings = [str(poly) for poly in basis]
    leading_monomials = [poly.terms()[0][0] for poly in basis if poly]
    unit_ideal = any(poly == ring.one for poly in basis)
    zero_dimensional = (True if unit_ideal else all(
        any(monomial[i] > 0 and
            all(other == 0 for j, other in enumerate(monomial) if j != i)
            for monomial in leading_monomials)
        for i in range(len(variables))
    ))
    result = {
        "source_sha256": next(iter(hashes)),
        "input_reports": [str(path) for path in args.inputs],
        "selected_variables": list(variables),
        "input_basis_lengths": [source["groebner_basis_length"]
                                for source in sources],
        "input_generators": len(input_polynomials),
        "groebner_method": args.method,
        "monomial_order": args.order,
        "groebner_basis_length": len(basis_strings),
        "groebner_basis": basis_strings,
        "leading_monomials": [list(monomial) for monomial in leading_monomials],
        "unit_ideal": unit_ideal,
        "zero_dimensional": zero_dimensional,
        "scope": (
            "Exact sum of the reported complete-support necessary-system ideals. "
            "A unit basis excludes their joint selected support slice only."
        ),
    }
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "source_sha256", "selected_variables", "input_basis_lengths",
        "input_generators", "groebner_method", "monomial_order",
        "groebner_basis_length", "leading_monomials", "unit_ideal",
        "zero_dimensional", "scope",
    )}, indent=2), flush=True)
    print(f"report={args.report}", flush=True)


if __name__ == "__main__":
    main()
