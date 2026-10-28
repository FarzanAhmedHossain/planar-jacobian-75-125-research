#!/usr/bin/env python3
"""Exact Gröbner screen on equations supported only on the five-variable block."""

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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("report", type=Path)
    parser.add_argument("--method", choices=("f5b", "buchberger"), default="f5b")
    parser.add_argument("--order", choices=("grevlex", "lex"), default="grevlex")
    args = parser.parse_args()

    text = args.source.read_text(encoding="utf-8")
    metadata, constraint_pos = top_metadata_and_constraints(text)
    if (metadata.get("format") != "f2-sparse-ring-v1"
            or metadata.get("branch") != "2+1+1" or metadata.get("quadratic")):
        raise ValueError("expected a rational 2+1+1 sparse checkpoint")

    names = ("p1", "p7", "p15", "p16", "p17")
    symbols = sp.symbols(" ".join(names))
    symbol_by_generator = {
        metadata["generators"].index(name): symbol
        for name, symbol in zip(names, symbols)
    }
    p79_index = metadata["generators"].index("p79")
    input_polynomials: list[tuple[str, sp.Poly]] = []
    for tag, sparse_terms in iter_json_array(text, constraint_pos):
        specialized: dict[tuple[tuple[int, int], ...], Fraction] = defaultdict(Fraction)
        for sparse, values in sparse_terms:
            if len(values) != 1:
                raise ValueError("expected rational coefficients")
            if any(generator == p79_index and exponent for generator, exponent in sparse):
                continue
            monomial = tuple(sorted(
                (int(generator), int(exponent))
                for generator, exponent in sparse if exponent
            ))
            specialized[monomial] += Fraction(values[0])
        specialized = {
            monomial: coefficient for monomial, coefficient in specialized.items()
            if coefficient
        }
        # Retain an equation only if its entire p79-specialized polynomial lies
        # in QQ[p1,p7,p15,p16,p17]. Dropping mixed monomials would not preserve
        # the equation's zero locus for unrestricted values of the other
        # parameters.
        if any(
            generator not in symbol_by_generator
            for monomial in specialized
            for generator, _exponent in monomial
        ):
            continue
        terms: dict[tuple[int, ...], Fraction] = defaultdict(Fraction)
        for monomial, coefficient in specialized.items():
            exponents = [0] * len(names)
            for generator, exponent in monomial:
                name = metadata["generators"][generator]
                exponents[names.index(name)] += exponent
            terms[tuple(exponents)] += coefficient
        terms = {monomial: coefficient for monomial, coefficient in terms.items() if coefficient}
        if not terms:
            continue
        support = {
            i for i in range(len(names))
            if any(monomial[i] for monomial in terms)
        }
        if len(support) != len(names):
            continue
        polynomial = sp.Poly.from_dict(
            {monomial: sp.Rational(coefficient.numerator, coefficient.denominator)
             for monomial, coefficient in terms.items()},
            symbols,
            domain=QQ,
        )
        input_polynomials.append((tag, polynomial))

    # Remove scalar duplicates without changing the generated ideal.
    unique: dict[tuple, tuple[str, sp.Poly]] = {}
    for tag, polynomial in input_polynomials:
        monic = polynomial.monic()
        key = tuple(sorted(monic.terms()))
        unique.setdefault(key, (tag, monic))
    polynomials = [entry[1] for entry in unique.values()]

    common = polynomials[0] if polynomials else None
    for polynomial in polynomials[1:]:
        common = sp.gcd(common, polynomial)
        if common.total_degree() == 0:
            break

    groebner = None
    groebner_ring = None
    if polynomials:
        groebner_ring, *_ = polynomial_ring(names, QQ, order=args.order)
        ring_polynomials = [
            groebner_ring.from_dict(dict(polynomial.terms()))
            for polynomial in polynomials
        ]
        groebner = ring_groebner(ring_polynomials, groebner_ring, method=args.method)
    basis = [] if groebner is None else [str(poly) for poly in groebner]
    is_unit = bool(groebner is not None and any(poly == groebner_ring.one for poly in groebner))
    if groebner is None:
        zero_dimensional = None
    elif is_unit:
        zero_dimensional = True
    else:
        leading_monomials = [poly.terms()[0][0] for poly in groebner if poly]
        zero_dimensional = all(
            any(monomial[index] > 0 and all(other == 0 for j, other in enumerate(monomial) if j != index)
                for monomial in leading_monomials)
            for index in range(len(names))
        )
    result = {
        "source": str(args.source),
        "source_sha256": sha256(args.source),
        "branch": metadata["branch"],
        "checkpoint_parameters": int(metadata["nparameters"]),
        "forced_specialization": {"p79": "0"},
        "variables": list(names),
        "selection_semantics": (
            "Keep only equations whose complete polynomial after p79=0 has no "
            "monomial involving any parameter outside the five listed variables."
        ),
        "input_equations": len(input_polynomials),
        "unique_equations_up_to_nonzero_rational_scale": len(polynomials),
        "common_gcd": None if common is None else str(common.as_expr()),
        "groebner_method": args.method,
        "monomial_order": args.order,
        "groebner_basis_length": len(basis),
        "groebner_basis": basis,
        "unit_ideal": is_unit,
        "zero_dimensional": zero_dimensional,
        "scope": (
            "These are necessary equations for the checkpoint chart. An exact unit ideal "
            "excludes this chart; a nonunit basis does not resolve the remaining equations."
        ),
    }
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "source_sha256", "variables", "input_equations",
        "unique_equations_up_to_nonzero_rational_scale", "common_gcd",
        "groebner_method", "monomial_order", "groebner_basis_length",
        "unit_ideal", "zero_dimensional", "scope"
    )}, indent=2))
    print(f"report={args.report}")


if __name__ == "__main__":
    main()
