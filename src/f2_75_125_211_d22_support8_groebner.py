#!/usr/bin/env python3
"""Exact Groebner screen on a repeated eight-parameter checkpoint support."""

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


NAMES = ("p5", "p15", "p34", "p35", "p36", "p37", "p38", "p39")


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
        rows_scanned = input_terms = 0
        skipped_auxiliary = 0
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

            if has_other_generator:
                skipped_auxiliary += 1
                continue
            if support != target_support:
                continue
            coefficients = {monomial: value for monomial, value in coefficients.items()
                            if value}
            if not coefficients:
                continue
            polynomial = sp.Poly.from_dict(
                {monomial: sp.Rational(value.numerator, value.denominator)
                 for monomial, value in coefficients.items()},
                symbols,
                domain=QQ,
            )
            input_rows.append((tag, polynomial))
            if rows_scanned % 200 == 0:
                print(f"scanned={rows_scanned} selected={len(input_rows)} "
                      f"sparse_terms={input_terms}", flush=True)
    finally:
        reader.close()

    unique: dict[tuple, tuple[str, sp.Poly]] = {}
    for tag, polynomial in input_rows:
        monic = polynomial.monic()
        unique.setdefault(tuple(sorted(monic.terms())), (tag, monic))
    polynomials = [entry[1] for entry in unique.values()]

    ring_basis = None
    basis = []
    if polynomials:
        ring_basis, *_ = polynomial_ring(NAMES, QQ, order=args.order)
        ring_inputs = [ring_basis.from_dict(dict(poly.terms()))
                       for poly in polynomials]
        basis = ring_groebner(ring_inputs, ring_basis, method=args.method)

    basis_strings = [str(poly) for poly in basis]
    unit_ideal = bool(ring_basis is not None
                      and any(poly == ring_basis.one for poly in basis))
    result = {
        "source": str(args.source),
        "source_sha256": source_hash,
        "branch": metadata["branch"],
        "checkpoint_parameters": int(metadata["nparameters"]),
        "selected_variables": list(NAMES),
        "selection_semantics": (
            "Keep complete checkpoint equations whose full polynomial uses only "
            "p5,p15,p34,p35,p36,p37,p38,p39 and uses every one of them. "
            "Terms involving t or any other generator are excluded, not projected away."
        ),
        "checkpoint_rows_scanned": rows_scanned,
        "checkpoint_sparse_terms_scanned": input_terms,
        "skipped_rows_with_other_generators": skipped_auxiliary,
        "input_equations": len(input_rows),
        "unique_equations_up_to_nonzero_rational_scale": len(polynomials),
        "selected_tags": [tag for tag, _ in input_rows],
        "origin_nonzero_equations": [
            tag for tag, poly in input_rows
            if poly.as_expr().subs({symbol: 0 for symbol in symbols}) != 0
        ],
        "groebner_method": args.method,
        "monomial_order": args.order,
        "groebner_basis_length": len(basis_strings),
        "groebner_basis": basis_strings,
        "unit_ideal": unit_ideal,
        "scope": (
            "This is an exact basis of a complete low-variable necessary subsystem. "
            "A nonunit basis is not a branch exclusion; a unit basis would exclude "
            "the represented U1=0 chart slice only."
        ),
    }
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "source_sha256", "selected_variables", "checkpoint_rows_scanned",
        "checkpoint_sparse_terms_scanned", "skipped_rows_with_other_generators",
        "input_equations", "unique_equations_up_to_nonzero_rational_scale",
        "origin_nonzero_equations", "groebner_method", "monomial_order",
        "groebner_basis_length", "groebner_basis", "unit_ideal", "scope",
    )}, indent=2), flush=True)
    print(f"report={args.report}", flush=True)


if __name__ == "__main__":
    main()
