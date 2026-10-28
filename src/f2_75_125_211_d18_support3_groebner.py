#!/usr/bin/env python3
"""Exact Groebner screen on the d18 equations supported only on p13,p14,p15.

The input is the full saved sparse necessary-equation checkpoint. Only complete
equations whose entire polynomial support is exactly the named three-parameter
set are retained. A unit ideal therefore excludes the saved chart slice; a
nonunit ideal is only a partial screen.
"""

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

from f2_75_125_211_streaming_checkpoint_scan import JsonStream


NAMES = ("p13", "p14", "p15")


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

    source_hash = sha256(args.source)
    reader = JsonStream(args.source)
    try:
        metadata, _ = reader.top_level_constraints()
        if (metadata.get("format") != "f2-sparse-ring-v1"
                or metadata.get("branch") != "2+1+1"
                or metadata.get("quadratic")):
            raise ValueError("expected a rational 2+1+1 sparse checkpoint")
        generators = metadata["generators"]
        target_ids = {generators.index(name): i for i, name in enumerate(NAMES)}
        target_support = set(range(len(NAMES)))

        symbols = sp.symbols(" ".join(NAMES))
        input_rows: list[tuple[str, sp.Poly]] = []
        input_terms = rows_scanned = 0
        for tag, sparse_terms in reader.array_values():
            rows_scanned += 1
            coefficients: dict[tuple[int, ...], Fraction] = defaultdict(Fraction)
            support: set[int] = set()
            has_other_generator = False
            for sparse, values in sparse_terms:
                input_terms += 1
                if len(values) != 1:
                    raise ValueError("expected rational checkpoint coefficients")
                exponents = [0] * len(NAMES)
                for generator, exponent in sparse:
                    exponent = int(exponent)
                    if not exponent:
                        continue
                    index = target_ids.get(int(generator))
                    if index is None:
                        has_other_generator = True
                        continue
                    exponents[index] += exponent
                    support.add(index)
                coefficients[tuple(exponents)] += Fraction(values[0])

            if has_other_generator or support != target_support:
                continue
            coefficients = {monomial: value for monomial, value in coefficients.items()
                            if value}
            if not coefficients:
                continue
            actual_support = {
                index for index in range(len(NAMES))
                if any(monomial[index] for monomial in coefficients)
            }
            if actual_support != target_support:
                continue
            polynomial = sp.Poly.from_dict(
                {monomial: sp.Rational(value.numerator, value.denominator)
                 for monomial, value in coefficients.items()},
                symbols,
                domain=QQ,
            )
            input_rows.append((tag, polynomial))
    finally:
        reader.close()

    unique: dict[tuple, tuple[str, sp.Poly]] = {}
    for tag, polynomial in input_rows:
        monic = polynomial.monic()
        key = tuple(sorted(monic.terms()))
        unique.setdefault(key, (tag, monic))
    polynomials = [entry[1] for entry in unique.values()]

    common = None
    if polynomials:
        common = polynomials[0]
        for polynomial in polynomials[1:]:
            common = sp.gcd(common, polynomial)
            if common.total_degree() == 0:
                break

    groebner_ring = None
    basis = []
    if polynomials:
        groebner_ring, *_ = polynomial_ring(NAMES, QQ, order=args.order)
        ring_polynomials = [
            groebner_ring.from_dict(dict(polynomial.terms()))
            for polynomial in polynomials
        ]
        basis = ring_groebner(ring_polynomials, groebner_ring,
                              method=args.method)
    basis_strings = [str(poly) for poly in basis]
    unit_ideal = bool(groebner_ring is not None
                      and any(poly == groebner_ring.one for poly in basis))
    leading_monomials = [poly.terms()[0][0] for poly in basis if poly]
    zero_dimensional = (True if unit_ideal else
                        all(any(monomial[i] > 0 and
                                all(other == 0 for j, other in enumerate(monomial)
                                    if j != i)
                                for monomial in leading_monomials)
                            for i in range(len(NAMES))) if basis else None)

    result = {
        "source": str(args.source),
        "source_sha256": source_hash,
        "branch": metadata["branch"],
        "checkpoint_parameters": int(metadata["nparameters"]),
        "selected_variables": list(NAMES),
        "selection_semantics": (
            "Keep complete checkpoint equations whose full polynomial uses only "
            "p13,p14,p15 and uses all three. No mixed monomials are projected away."
        ),
        "checkpoint_rows_scanned": rows_scanned,
        "input_equations": len(input_rows),
        "unique_equations_up_to_nonzero_rational_scale": len(polynomials),
        "input_sparse_terms_scanned": input_terms,
        "common_gcd": None if common is None else str(common.as_expr()),
        "groebner_method": args.method,
        "monomial_order": args.order,
        "groebner_basis_length": len(basis_strings),
        "groebner_basis": basis_strings,
        "unit_ideal": unit_ideal,
        "zero_dimensional": zero_dimensional,
        "scope": (
            "An exact unit ideal excludes the saved necessary-system chart slice. "
            "A nonunit basis does not exclude the U1=0 subcase or the whole F2 branch."
        ),
    }
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "source_sha256", "selected_variables", "input_equations",
        "unique_equations_up_to_nonzero_rational_scale", "common_gcd",
        "groebner_method", "monomial_order", "groebner_basis_length",
        "groebner_basis", "unit_ideal", "zero_dimensional", "scope",
    )}, indent=2), flush=True)
    print(f"report={args.report}", flush=True)


if __name__ == "__main__":
    main()
