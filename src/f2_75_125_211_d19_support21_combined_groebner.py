"""Exact weighted Groebner screen for the support-12 + support-21 slice."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import time

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.groebnertools import groebner as ring_groebner
from sympy.polys.rings import ring as polynomial_ring


ROOT = Path("research")
FAMILY = ROOT / "f2_75_125_211_d19_support21_family_input_20261025.json"
ANALYSIS = ROOT / "f2_75_125_211_d19_support21_analysis_20261025.json"
SUPPORT12_GB = ROOT / "f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_gb_4gen_weighted_flint.json"
SUPPORT12_INPUT = ROOT / "f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_input.json"
SUPPORT12_REUSE = ROOT / "f2_75_125_211_verify_d19_support12_reuse_20261025.json"
OUTPUT = ROOT / "f2_75_125_211_d19_support21_combined_gb_20261025.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    started = time.monotonic()
    family = json.loads(FAMILY.read_text(encoding="utf-8"))
    analysis = json.loads(ANALYSIS.read_text(encoding="utf-8"))
    old_report = json.loads(SUPPORT12_GB.read_text(encoding="utf-8"))
    old_input = json.loads(SUPPORT12_INPUT.read_text(encoding="utf-8"))
    reuse = json.loads(SUPPORT12_REUSE.read_text(encoding="utf-8"))

    names = tuple(family["target_support"])
    weights = tuple(int(value)
                    for value in analysis["primitive_positive_integer_weights"])
    if len(names) != 21 or len(weights) != len(names):
        raise ValueError("support-21 variable/weight metadata mismatch")
    if not analysis["all_equations_weight_homogeneous"]:
        raise ValueError("support-21 equations lack a verified positive grading")
    if not reuse.get("same_142_equations_and_full_support12_group"):
        raise ValueError("d19 support-12 group is not verified as the d18 system")
    if old_input.get("equation_count") != 142:
        raise ValueError("unexpected support-12 input size")
    if old_report.get("unit_ideal"):
        raise ValueError("unexpected old support-12 unit basis")

    symbols = sp.symbols(" ".join(names))
    symbol_map = dict(zip(names, symbols))
    name_index = {name: i for i, name in enumerate(names)}

    def weighted_grevlex(monomial):
        return (
            sum(weight*exponent for weight, exponent in zip(weights, monomial)),
            tuple(reversed([-exponent for exponent in monomial])),
        )

    polynomial_ring_21, *_ = polynomial_ring(
        symbols, QQ, order=weighted_grevlex
    )

    selected_indices = analysis["independent_equation_indices_zero_based"]
    if len(selected_indices) != analysis["exact_Q_linear_rank"]:
        raise ValueError("support-21 exact row-rank receipt mismatch")
    inputs = []
    input_labels = []
    input_weight_degrees = []

    for index in selected_indices:
        equation = family["equations"][index]
        coefficients = {}
        for monomial, coefficient in equation["terms"]:
            exponents = [0] * len(names)
            for name, power in monomial:
                exponents[name_index[name]] += int(power)
            key = tuple(exponents)
            coefficients[key] = coefficients.get(key, sp.Rational(0)) + sp.Rational(coefficient)
        coefficients = {key: value for key, value in coefficients.items() if value}
        poly = polynomial_ring_21.from_dict(coefficients)
        inputs.append(poly)
        input_labels.append(equation["tag"])
        input_weight_degrees.append(sorted({
            sum(w*e for w, e in zip(weights, monomial))
            for monomial in coefficients
        }))

    old_names = old_report["selected_variables"]
    old_symbol_map = {name: symbol_map[name] for name in old_names}
    for index, record in enumerate(old_report["basis"]):
        expression = sp.sympify(record["polynomial"], locals=old_symbol_map)
        poly_sympy = sp.Poly(expression, *symbols, domain=QQ)
        coefficients = {
            tuple(monomial): coefficient
            for monomial, coefficient in poly_sympy.terms()
            if coefficient
        }
        poly = polynomial_ring_21.from_dict(coefficients)
        inputs.append(poly)
        input_labels.append(f"support12-G{index}")
        input_weight_degrees.append(sorted({
            sum(w*e for w, e in zip(weights, monomial))
            for monomial in coefficients
        }))

    if len(inputs) != 8:
        raise ValueError("expected four support-21 rows and four support-12 basis rows")
    origin_monomial = tuple([0]*len(names))
    if any(dict(poly.terms()).get(origin_monomial, 0) for poly in inputs):
        raise ValueError("origin witness does not satisfy every generator")

    print(json.dumps({
        "stage": "groebner_start",
        "variables": len(names),
        "generators": len(inputs),
        "labels": input_labels,
        "weight_degrees": input_weight_degrees,
        "weights": weights,
    }), flush=True)
    basis = ring_groebner(inputs, polynomial_ring_21, method="f5b")
    basis_records = []
    pure_powers = {}
    for poly in basis:
        terms = poly.terms()
        leading_monomial = terms[0][0] if terms else None
        basis_records.append({
            "polynomial": str(poly),
            "term_count": len(terms),
            "total_degree": max((sum(monomial) for monomial, _ in terms), default=0),
            "weighted_degree": max((sum(w*e for w, e in zip(weights, monomial))
                                    for monomial, _ in terms), default=0),
            "leading_monomial": list(leading_monomial) if leading_monomial else None,
        })
        if len(terms) == 1:
            monomial = terms[0][0]
            nonzero = [i for i, exponent in enumerate(monomial) if exponent]
            if len(nonzero) == 1:
                name = names[nonzero[0]]
                pure_powers[name] = int(monomial[nonzero[0]])

    unit_ideal = any(poly == polynomial_ring_21.one for poly in basis)
    result = {
        "scope": (
            "Exact weighted Groebner basis for the ideal generated by the four "
            "support-12 basis elements and four Q-linearly independent generators "
            "spanning the complete 139-row support-21 family. This is one necessary "
            "slice only; a nonunit result does not exclude the whole F2 branch."
        ),
        "support21_input_sha256": sha256(FAMILY),
        "support21_analysis_sha256": sha256(ANALYSIS),
        "support12_groebner_report_sha256": sha256(SUPPORT12_GB),
        "support12_input_sha256": sha256(SUPPORT12_INPUT),
        "support12_reuse_receipt_sha256": sha256(SUPPORT12_REUSE),
        "variables": list(names),
        "weights": list(weights),
        "weighted_monomial_order": "positive-weight grevlex, lex tie-break inherited on the support-12 subring",
        "support21_equations": 139,
        "support21_exact_Q_linear_rank": analysis["exact_Q_linear_rank"],
        "support21_independent_tags": input_labels[:4],
        "support12_basis_generators": 4,
        "input_generators": len(inputs),
        "input_labels": input_labels,
        "input_weight_degrees": input_weight_degrees,
        "origin_common_zero_verified": True,
        "groebner_method": "f5b",
        "basis_length": len(basis_records),
        "unit_ideal": unit_ideal,
        "basis": basis_records,
        "pure_power_generators": pure_powers,
        "krull_dimension_lower_bound": len(names)-len(inputs),
        "elapsed_seconds": time.monotonic()-started,
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "input_generators", "input_labels", "input_weight_degrees",
        "basis_length", "unit_ideal", "pure_power_generators",
        "krull_dimension_lower_bound", "elapsed_seconds",
    )}, indent=2), flush=True)
    print(f"report={OUTPUT}", flush=True)


if __name__ == "__main__":
    main()
