#!/usr/bin/env python3
"""Exact F5B basis for the 12-variable subsystem after verified L=M=0."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.groebnertools import groebner as ring_groebner
from sympy.polys.rings import ring as polynomial_ring


EXPECTED_SOURCE_SHA256 = "270f0eb737463c75d3bdfbcdfb3a303a2e7fa0af2da13f63a9e2d01ec5e04b59"
EXPECTED_SUPPORT8_REPORT_SHA256 = "7a041a934a7df94ff5a3b12833a29dbe43f83adc62f6ec8823248202f901b89d"
EXPECTED_VARIABLES = (
    "p15", "p24", "p25", "p26", "p27", "p28",
    "p44", "p45", "p46", "p47", "p48", "p5",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("report", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--order", choices=("grevlex", "lex"), default="grevlex")
    args = parser.parse_args()

    data = json.loads(args.source.read_text(encoding="utf-8"))
    if data.get("source_sha256") != EXPECTED_SOURCE_SHA256:
        raise ValueError("support-12 equations are bound to another checkpoint")
    if data.get("support8_report_sha256") != EXPECTED_SUPPORT8_REPORT_SHA256:
        raise ValueError("support-12 equations are bound to another support8 report")
    if tuple(data.get("variables", ())) != EXPECTED_VARIABLES:
        raise ValueError("unexpected 12-variable support")
    if int(data.get("equation_count", -1)) != 142 or len(data.get("equations", ())) != 142:
        raise ValueError("unexpected support-12 equation count")
    report_hash = sha256(args.report)
    if report_hash != EXPECTED_SUPPORT8_REPORT_SHA256:
        raise ValueError("support8 report hash does not match the verified input")

    variables = sp.symbols(" ".join(EXPECTED_VARIABLES))
    positions = {name: index for index, name in enumerate(EXPECTED_VARIABLES)}
    input_polynomials: list[tuple[str, sp.Poly]] = []
    max_degree = 0
    max_terms = 0
    for row in data["equations"]:
        coefficients: dict[tuple[int, ...], Fraction] = defaultdict(Fraction)
        for monomial, coefficient in row["terms"]:
            exponents = [0] * len(EXPECTED_VARIABLES)
            for name, exponent in monomial:
                if name not in positions:
                    raise ValueError(f"equation {row['tag']} escaped the declared support")
                exponents[positions[name]] += int(exponent)
            coefficients[tuple(exponents)] += Fraction(coefficient)
        coefficients = {monomial: coefficient for monomial, coefficient in coefficients.items()
                        if coefficient}
        if not coefficients:
            continue
        polynomial = sp.Poly.from_dict(
            {monomial: sp.Rational(value.numerator, value.denominator)
             for monomial, value in coefficients.items()},
            variables,
            domain=QQ,
        )
        max_degree = max(max_degree, polynomial.total_degree())
        max_terms = max(max_terms, len(polynomial.terms()))
        input_polynomials.append((row["tag"], polynomial))

    unique: dict[tuple, tuple[str, sp.Poly]] = {}
    for tag, polynomial in input_polynomials:
        monic = polynomial.monic()
        unique.setdefault(tuple(sorted(monic.terms())), (tag, monic))
    polynomials = [entry[1] for entry in unique.values()]
    ring_basis, *_ = polynomial_ring(EXPECTED_VARIABLES, QQ, order=args.order)
    ring_inputs = [ring_basis.from_dict(dict(poly.terms())) for poly in polynomials]
    basis = ring_groebner(ring_inputs, ring_basis, method="f5b")
    basis_strings = [str(poly) for poly in basis]
    unit_ideal = any(poly == ring_basis.one for poly in basis)

    result = {
        "equations_source": str(args.source),
        "equations_source_sha256": sha256(args.source),
        "checkpoint_sha256": EXPECTED_SOURCE_SHA256,
        "support8_report": str(args.report),
        "support8_report_sha256": report_hash,
        "variables": list(EXPECTED_VARIABLES),
        "input_equations": len(input_polynomials),
        "unique_equations_up_to_nonzero_rational_scale": len(polynomials),
        "max_total_degree": max_degree,
        "max_terms_per_equation": max_terms,
        "method": "f5b",
        "order": args.order,
        "groebner_basis_length": len(basis_strings),
        "groebner_basis": basis_strings,
        "unit_ideal": unit_ideal,
        "scope": (
            "Exact basis of the 142 complete checkpoint equations that, after the "
            "verified support-eight consequences L=M=0, have precisely this support. "
            "A nonunit basis is not an exclusion; a unit basis excludes only this "
            "necessary chart subsystem."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "equations_source_sha256", "input_equations",
        "unique_equations_up_to_nonzero_rational_scale", "max_total_degree",
        "max_terms_per_equation", "method", "order",
        "groebner_basis_length", "groebner_basis", "unit_ideal", "scope",
    )}, indent=2), flush=True)
    print(f"report={args.output}", flush=True)


if __name__ == "__main__":
    main()
