#!/usr/bin/env python3
"""Exact Groebner screen on complete checkpoint equations of one support set.

Only equations whose full parameter support is exactly the requested set are
included; no mixed monomials are projected away. A unit ideal excludes that
specific necessary-system slice, not automatically the full F2 branch.
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
    parser.add_argument("--variables", nargs="+", required=True,
                        help="full parameter support of the equations to keep")
    parser.add_argument("--method", choices=("f5b", "buchberger"), default="f5b")
    parser.add_argument("--order", choices=("grevlex", "lex"), default="grevlex")
    args = parser.parse_args()

    names = tuple(args.variables)
    if len(names) < 2 or len(set(names)) != len(names):
        raise ValueError("provide at least two distinct parameter names")

    source_hash = sha256(args.source)
    reader = JsonStream(args.source)
    try:
        metadata, _ = reader.top_level_constraints()
        if (metadata.get("format") != "f2-sparse-ring-v1"
                or metadata.get("branch") != "2+1+1"
                or metadata.get("quadratic")):
            raise ValueError("expected a rational 2+1+1 sparse checkpoint")
        generators = metadata["generators"]
        missing_names = sorted(set(names) - set(generators))
        if missing_names:
            raise ValueError(f"checkpoint lacks selected parameters: {missing_names}")
        target_ids = {generators.index(name): i for i, name in enumerate(names)}
        target_support = set(range(len(names)))
        symbols = sp.symbols(" ".join(names))

        input_rows: list[tuple[str, sp.Poly]] = []
        input_terms = rows_scanned = 0
        for tag, sparse_terms in reader.array_values():
            rows_scanned += 1
            if rows_scanned % 200 == 0:
                print(f"checkpoint scan {rows_scanned} rows; "
                      f"selected={len(input_rows)}", flush=True)
            coefficients: dict[tuple[int, ...], Fraction] = defaultdict(Fraction)
            support: set[int] = set()
            has_other_generator = False
            for sparse, values in sparse_terms:
                input_terms += 1
                if len(values) != 1:
                    raise ValueError("expected rational checkpoint coefficients")
                exponents = [0] * len(names)
                for generator, exponent in sparse:
                    exponent = int(exponent)
                    if not exponent:
                        continue
                    index = target_ids.get(int(generator))
                    if index is None:
                        has_other_generator = True
                    else:
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
                index for index in range(len(names))
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

    # Exact Q-linear row reduction preserves the generated ideal and removes
    # redundant checkpoint equations before the more expensive Groebner step.
    echelon: dict[tuple[int, ...], dict[tuple[int, ...], Fraction]] = {}
    independent_rows: list[tuple[str, dict[tuple[int, ...], Fraction]]] = []
    # Scale duplicates do not change the Q-row span. Reduce only unique rows;
    # this matters for support groups where dozens of checkpoint constraints
    # are scalar copies of a small number of large sparse polynomials.
    for row_index, (tag, polynomial) in enumerate(unique.values(), start=1):
        row = {monomial: Fraction(coefficient)
               for monomial, coefficient in polynomial.terms() if coefficient}
        while row:
            pivot = max(row)
            coefficient = row[pivot]
            known = echelon.get(pivot)
            if known is not None:
                for monomial, basis_coefficient in known.items():
                    value = row.get(monomial, Fraction(0)) - coefficient*basis_coefficient
                    if value:
                        row[monomial] = value
                    else:
                        row.pop(monomial, None)
                continue
            inverse = 1 / coefficient
            row = {monomial: value*inverse for monomial, value in row.items() if value}
            echelon[pivot] = row
            independent_rows.append((tag, row))
            break
        if row_index % 10 == 0 or row_index == len(unique):
            print(f"exact-Q row reduction {row_index}/{len(unique)}; "
                  f"rank={len(independent_rows)}", flush=True)

    polynomials = [
        sp.Poly.from_dict(
            {monomial: sp.Rational(value.numerator, value.denominator)
             for monomial, value in row.items()},
            symbols,
            domain=QQ,
        )
        for _tag, row in independent_rows
    ]

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
        groebner_ring, *_ = polynomial_ring(names, QQ, order=args.order)
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
                            for i in range(len(names))) if basis else None)

    result = {
        "source": str(args.source),
        "source_sha256": source_hash,
        "branch": metadata["branch"],
        "checkpoint_parameters": int(metadata["nparameters"]),
        "selected_variables": list(names),
        "selection_semantics": (
            "Keep complete checkpoint equations whose full polynomial uses only "
            f"{', '.join(names)} and uses all selected variables; no mixed terms "
            "are discarded."
        ),
        "checkpoint_rows_scanned": rows_scanned,
        "input_equations": len(input_rows),
        "unique_equations_up_to_nonzero_rational_scale": len(unique),
        "exact_Q_linear_rank": len(polynomials),
        "input_sparse_terms_scanned": input_terms,
        "common_gcd": None if common is None else str(common.as_expr()),
        "groebner_method": args.method,
        "monomial_order": args.order,
        "groebner_basis_length": len(basis_strings),
        "groebner_basis": basis_strings,
        "unit_ideal": unit_ideal,
        "zero_dimensional": zero_dimensional,
        "scope": (
            "An exact unit ideal excludes this complete-support slice of the saved "
            "necessary system. A nonunit basis does not exclude the entire F2 branch."
        ),
    }
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "source_sha256", "selected_variables", "input_equations",
        "unique_equations_up_to_nonzero_rational_scale", "exact_Q_linear_rank",
        "input_sparse_terms_scanned",
        "common_gcd", "groebner_method", "monomial_order",
        "groebner_basis_length", "groebner_basis", "unit_ideal",
        "zero_dimensional", "scope",
    )}, indent=2), flush=True)
    print(f"report={args.report}", flush=True)


if __name__ == "__main__":
    main()
