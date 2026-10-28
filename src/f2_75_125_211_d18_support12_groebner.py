#!/usr/bin/env python3
"""Exact QQ Groebner computation for the extracted 12-variable subsystem."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.groebnertools import groebner as ring_groebner
from sympy.polys.rings import ring as polynomial_ring


INPUT = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_input.json")
EXPECTED_INPUT_SHA256 = "A9DF8405D43EAFFB98C724165C0C7A647C0952D7D8D1A8630CD26700898C6F56"
WEIGHTS = (4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 2, 3)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest().upper()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("report", type=Path)
    ap.add_argument("--order", choices=("grevlex", "lex", "wgrevlex"), default="grevlex")
    ap.add_argument("--method", choices=("f5b", "buchberger"), default="f5b")
    ap.add_argument("--use-linear-reduction", action="store_true",
                    help="use the independently verified exact Q-linear generator basis")
    args = ap.parse_args()

    in_hash = sha256(INPUT)
    if in_hash != EXPECTED_INPUT_SHA256:
        raise ValueError(f"support-12 input hash mismatch: {in_hash}")
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    names = tuple(data["selected_variables"])
    if len(names) != 12 or data["equation_count"] != len(data["equations"]):
        raise ValueError("support-12 extraction metadata mismatch")
    name_index = {name: i for i, name in enumerate(names)}

    start = time.perf_counter()
    polys = []
    term_count = max_degree = 0
    origin_violations = []
    for row in data["equations"]:
        coeffs = {}
        used = set()
        for monomial, coefficient in row["terms"]:
            exponents = [0] * len(names)
            for name, exponent in monomial:
                if name not in name_index:
                    raise ValueError(f"unexpected variable {name} in {row['tag']}")
                exponent = int(exponent)
                exponents[name_index[name]] += exponent
                used.add(name)
            key = tuple(exponents)
            coeffs[key] = coeffs.get(key, sp.Rational(0)) + sp.Rational(coefficient)
            term_count += 1
            max_degree = max(max_degree, sum(key))
        coeffs = {m: c for m, c in coeffs.items() if c}
        if used != set(names):
            raise ValueError(f"row {row['tag']} does not have complete support")
        if coeffs.get((0,) * len(names), 0):
            origin_violations.append(row["tag"])
        polys.append(sp.Poly.from_dict(coeffs, sp.symbols(" ".join(names)), domain=QQ))
    if origin_violations:
        raise ValueError(f"claimed origin witness fails: {origin_violations[:5]}")

    selected_indices = list(range(len(polys)))
    reduction_verified = False
    if args.use_linear_reduction:
        reduction_path = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_linear_reduce.json")
        verified_path = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_linear_reduce_verified.json")
        reduction = json.loads(reduction_path.read_text(encoding="utf-8"))
        verified = json.loads(verified_path.read_text(encoding="utf-8"))
        if reduction.get("input_sha256") != in_hash or verified.get("input_sha256") != in_hash:
            raise ValueError("linear reduction is not bound to the support-12 input")
        if verified.get("reduction_report_sha256") != sha256(reduction_path):
            raise ValueError("linear reduction verification receipt is stale")
        if not verified.get("all_142_equations_span_exactly_4_Q_generators"):
            raise ValueError("exact Q-linear reduction was not verified")
        selected_indices = list(verified["independent_equation_indices"])
        if len(selected_indices) != verified.get("independent_equations"):
            raise ValueError("verified independent-row metadata mismatch")
        if len(selected_indices) + verified.get("verified_exact_dependencies", -1) != len(polys):
            raise ValueError("verified dependency count does not cover all input equations")
        polys = [polys[i] for i in selected_indices]
        reduction_verified = True
    print(json.dumps({
        "stage": "parsed",
        "input_sha256": in_hash,
        "branch": data["branch"],
        "variables": names,
        "equations": len(polys),
        "generator_selection": selected_indices,
        "exact_linear_reduction_verified": reduction_verified,
        "input_terms": term_count,
        "max_total_degree": max_degree,
        "origin_witness": {name: 0 for name in names},
        "elapsed_seconds": round(time.perf_counter() - start, 3),
    }), flush=True)

    rows_for_polys = [data["equations"][i] for i in selected_indices]
    unique = {}
    for row, poly in zip(rows_for_polys, polys):
        monic = poly.monic()
        unique.setdefault(tuple(monic.terms()), (row["tag"], monic))
    inputs = [poly for _, poly in unique.values()]
    print(json.dumps({"stage": "groebner_start", "unique_generators": len(inputs),
                      "method": args.method, "order": args.order,
                      "elapsed_seconds": round(time.perf_counter() - start, 3)}), flush=True)

    symbols = sp.symbols(" ".join(names))
    if args.order == "wgrevlex":
        weight_path = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_weighted_homogeneity.json")
        weight_data = json.loads(weight_path.read_text(encoding="utf-8"))
        if weight_data.get("input_sha256") != in_hash or not weight_data.get(
                "all_142_input_equations_weight_homogeneous"):
            raise ValueError("weighted-homogeneity receipt is missing or stale")
        if tuple(weight_data.get("primitive_positive_integer_weights", ())) != WEIGHTS:
            raise ValueError("unexpected weighted-order vector")
        def weighted_grevlex(monomial):
            return (sum(weight * exponent for weight, exponent in zip(WEIGHTS, monomial)),
                    tuple(reversed([-exponent for exponent in monomial])))
        order = weighted_grevlex
    else:
        order = args.order
    ring, *_ = polynomial_ring(symbols, QQ, order=order)
    ring_inputs = [ring.from_dict(dict(poly.terms())) for poly in inputs]
    basis = ring_groebner(ring_inputs, ring, method=args.method)
    basis_records = []
    for poly in basis:
        terms = poly.terms()
        basis_records.append({
            "polynomial": str(poly),
            "total_degree": max((sum(m) for m, _ in terms), default=0),
            "term_count": len(terms),
        })
    unit = any(poly == ring.one for poly in basis)
    result = {
        "input": str(INPUT),
        "input_sha256": in_hash,
        "upstream_checkpoint_sha256": data["source_sha256"],
        "upstream_support8_report_sha256": data["support8_report_sha256"],
        "upstream_support8_line_section_report_sha256": data["support8_line_section_report_sha256"],
        "branch": data["branch"],
        "selected_variables": list(names),
        "selection_semantics": data["selection_semantics"],
        "equation_count": len(polys),
        "source_equation_count": data["equation_count"],
        "selected_equation_indices": selected_indices,
        "exact_linear_reduction_verified": reduction_verified,
        "unique_generators_up_to_QQ_scaling": len(inputs),
        "input_terms": term_count,
        "max_total_degree": max_degree,
        "origin_witness": {name: 0 for name in names},
        "groebner_method": args.method,
        "monomial_order": args.order,
        "monomial_weights": list(WEIGHTS) if args.order == "wgrevlex" else None,
        "basis_length": len(basis_records),
        "unit_ideal": unit,
        "basis": basis_records,
        "elapsed_seconds": round(time.perf_counter() - start, 3),
        "scope": "Exact ideal of the selected support-12 subsystem only. If the verified linear reduction was used, its four selected generators span all 142 extracted equations over QQ. A nonunit basis does not exclude the full U1=0 slice; a unit basis excludes this support-12 ideal, subject to the recorded upstream exact substitutions and extraction semantics.",
    }
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "input_sha256", "selected_variables", "equation_count",
        "unique_generators_up_to_QQ_scaling", "input_terms", "max_total_degree",
        "origin_witness", "groebner_method", "monomial_order", "basis_length",
        "unit_ideal", "elapsed_seconds", "scope",
    )}, indent=2), flush=True)
    print(f"report={args.report}", flush=True)


if __name__ == "__main__":
    main()
