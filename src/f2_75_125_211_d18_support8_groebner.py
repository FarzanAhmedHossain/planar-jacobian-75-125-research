#!/usr/bin/env python3
"""Exact d18 screen on the repeated eight-variable support after line reduction."""

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

from f2_75_125_211_d18_support3_reduction_scan import (
    TARGET_RELATIONS,
    TARGET_NAMES,
    verify_support3_relations,
)
from f2_75_125_211_streaming_checkpoint_scan import JsonStream


SOURCE = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_state.json")
SUPPORT3_REPORT = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support3_gb.json")
SUPPORT8 = ("p5", "p34", "p35", "p36", "p37", "p38", "p39", "r")
EXPECTED_SOURCE_SHA256 = "e2e35061e883437806d6486d8f61fdc4cc6b963c49655366a6795db210d874e7"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("report", type=Path)
    parser.add_argument("--method", choices=("f5b", "buchberger"), default="f5b")
    parser.add_argument("--order", choices=("grevlex", "lex"), default="grevlex")
    args = parser.parse_args()

    source_hash = sha256(SOURCE)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise ValueError(f"degree-18 checkpoint SHA-256 mismatch: {source_hash}")
    basis_hash, basis_claim = verify_support3_relations(SUPPORT3_REPORT, source_hash)
    reader = JsonStream(SOURCE)
    try:
        metadata, _ = reader.top_level_constraints()
        generators = metadata["generators"]
        target_ids = {generators.index(name): name for name in TARGET_RELATIONS}
        other_ids = {
            generators.index(f"p{i}"): f"p{i}"
            for i in range(int(metadata["nparameters"]))
            if f"p{i}" not in TARGET_RELATIONS
        }
        names = tuple(name for name in SUPPORT8 if name != "r")
        symbols = sp.symbols(" ".join(names))
        variables = (*names, "r")
        selected = []
        rows_scanned = term_count = 0
        origin_failures = []
        for tag, sparse_terms in reader.array_values():
            rows_scanned += 1
            full_coefficients: dict[tuple[tuple[str, int], ...], Fraction] = defaultdict(Fraction)
            has_auxiliary = False
            for sparse, values in sparse_terms:
                term_count += 1
                if len(values) != 1:
                    raise ValueError("expected rational checkpoint coefficients")
                coefficient = Fraction(values[0])
                powers: dict[str, int] = defaultdict(int)
                for generator, exponent in sparse:
                    exponent = int(exponent)
                    if not exponent:
                        continue
                    if generator in target_ids:
                        name = target_ids[generator]
                        powers["r"] += exponent
                        coefficient *= TARGET_RELATIONS[name] ** exponent
                    elif generator in other_ids:
                        powers[other_ids[generator]] += exponent
                    else:
                        has_auxiliary = True
                full_monomial = tuple(sorted(
                    (name, exponent) for name, exponent in powers.items()
                    if exponent
                ))
                full_coefficients[full_monomial] += coefficient

            if has_auxiliary:
                raise ValueError(f"unexpected auxiliary generator in {tag}")
            full_coefficients = {
                monomial: coefficient
                for monomial, coefficient in full_coefficients.items()
                if coefficient
            }
            if not full_coefficients:
                continue
            # Do not project mixed equations into the selected ring. Keep only
            # complete equations whose entire support is exactly the target.
            actual_support = {
                name for monomial in full_coefficients for name, _ in monomial
            }
            if actual_support != set(variables):
                continue
            coefficients = {}
            for monomial, coefficient in full_coefficients.items():
                exponents = dict(monomial)
                coefficients[tuple(exponents.get(name, 0) for name in variables)] = coefficient
            constant = coefficients.get((0,)*len(variables), Fraction(0))
            if constant:
                origin_failures.append({"tag": tag, "constant": str(constant)})
            poly_dict = {
                monomial: sp.Rational(value.numerator, value.denominator)
                for monomial, value in coefficients.items()
            }
            polynomial = sp.Poly.from_dict(poly_dict, (*symbols, sp.Symbol("r")), domain=QQ)
            selected.append((tag, polynomial))
    finally:
        reader.close()

    unique: dict[tuple, tuple[str, sp.Poly]] = {}
    for tag, polynomial in selected:
        monic = polynomial.monic()
        unique.setdefault(tuple(monic.terms()), (tag, monic))
    polynomials = [row[1] for row in unique.values()]

    common = None
    if polynomials and origin_failures:
        common = polynomials[0]
        for polynomial in polynomials[1:]:
            common = sp.gcd(common, polynomial)
            if common.total_degree() == 0:
                break

    basis = []
    ring = None
    if polynomials and not origin_failures:
        ring_variables = (*symbols, sp.Symbol("r"))
        ring, *_ = polynomial_ring(ring_variables, QQ, order=args.order)
        inputs = [ring.from_dict(dict(poly.terms())) for poly in polynomials]
        basis = ring_groebner(inputs, ring, method=args.method)
    basis_strings = [str(poly) for poly in basis]
    unit_ideal = bool(ring is not None and any(poly == ring.one for poly in basis))
    result = {
        "source": str(SOURCE),
        "source_sha256": source_hash,
        "support3_basis_report_sha256": basis_hash,
        "support3_relations": basis_claim,
        "branch": metadata["branch"],
        "forced_specialization": {"p14": "-14*r", "p13": "(452/5)*r", "p15": "r"},
        "selected_variables": [*names, "r"],
        "selection_semantics": "Complete equations after the exact support-three line substitution with full eight-variable support; no mixed terms are discarded.",
        "checkpoint_rows_scanned": rows_scanned,
        "checkpoint_sparse_terms_scanned": term_count,
        "input_equations": len(selected),
        "unique_equations_up_to_nonzero_rational_scale": len(polynomials),
        "common_gcd": None if common is None else str(common.as_expr()),
        "origin_witness": None if origin_failures else {name: 0 for name in (*names, "r")},
        "origin_nonzero_equations": origin_failures,
        "groebner_method": None if origin_failures else args.method,
        "monomial_order": None if origin_failures else args.order,
        "groebner_basis_length": len(basis_strings),
        "groebner_basis": basis_strings,
        "unit_ideal": unit_ideal,
        "scope": "This exact necessary subsystem cannot exclude the chart if the origin witness is verified; it is only one component of the full reduced system.",
    }
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "source_sha256", "selected_variables", "input_equations",
        "unique_equations_up_to_nonzero_rational_scale", "common_gcd",
        "origin_witness", "origin_nonzero_equations", "groebner_method",
        "groebner_basis_length", "unit_ideal", "scope",
    )}, indent=2), flush=True)
    print(f"report={args.report}", flush=True)


if __name__ == "__main__":
    main()
