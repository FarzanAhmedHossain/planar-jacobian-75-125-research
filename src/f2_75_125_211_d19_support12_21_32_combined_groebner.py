"""Exact weighted F5B basis of the support-12/21/32 necessary subsystem."""

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
FAMILY21 = ROOT / "f2_75_125_211_d19_support21_family_input_20261025.json"
ANALYSIS21 = ROOT / "f2_75_125_211_d19_support21_analysis_20261025.json"
FAMILY32 = ROOT / "f2_75_125_211_d19_support32_family_input_20261025.json"
ANALYSIS32 = ROOT / "f2_75_125_211_d19_support32_analysis_20261025.json"
SUPPORT12_GB = ROOT / "f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_gb_4gen_weighted_flint.json"
SUPPORT12_INPUT = ROOT / "f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_input.json"
SUPPORT12_REUSE = ROOT / "f2_75_125_211_verify_d19_support12_reuse_20261025.json"
LOW_SUPPORT_REDUNDANCY = ROOT / "f2_75_125_211_new_small_support_redundancy_20261025.json"
OUTPUT = ROOT / "f2_75_125_211_d19_support12_21_32_combined_gb_20261025.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def family_polynomial(record, names, ring):
    index = {name: i for i, name in enumerate(names)}
    coefficients = {}
    for monomial, coefficient in record["terms"]:
        exponents = [0]*len(names)
        for name, power in monomial:
            exponents[index[name]] += int(power)
        key = tuple(exponents)
        coefficients[key] = coefficients.get(key, sp.Rational(0)) + sp.Rational(coefficient)
    return ring.from_dict({key: value for key, value in coefficients.items() if value})


def main() -> None:
    started = time.monotonic()
    family21 = json.loads(FAMILY21.read_text(encoding="utf-8"))
    analysis21 = json.loads(ANALYSIS21.read_text(encoding="utf-8"))
    family32 = json.loads(FAMILY32.read_text(encoding="utf-8"))
    analysis32 = json.loads(ANALYSIS32.read_text(encoding="utf-8"))
    old_report = json.loads(SUPPORT12_GB.read_text(encoding="utf-8"))
    old_input = json.loads(SUPPORT12_INPUT.read_text(encoding="utf-8"))
    reuse = json.loads(SUPPORT12_REUSE.read_text(encoding="utf-8"))
    redundancy = json.loads(LOW_SUPPORT_REDUNDANCY.read_text(encoding="utf-8"))

    names = tuple(analysis32["selected_variables"])
    weights = tuple(int(value) for value in analysis32["primitive_positive_integer_weights"])
    if len(names) != 22 or len(weights) != 22:
        raise ValueError("support-32 variable/weight metadata mismatch")
    if analysis32["exact_Q_linear_rank"] != 4 or analysis21["exact_Q_linear_rank"] != 4:
        raise ValueError("family compression metadata mismatch")
    if not analysis32["all_equations_weight_homogeneous"]:
        raise ValueError("support-32 grading is not verified")
    if not reuse.get("same_142_equations_and_full_support12_group"):
        raise ValueError("d19 support-12 group is not verified as d18")
    if old_input.get("equation_count") != 142 or old_report.get("unit_ideal"):
        raise ValueError("unexpected support-12 source system")
    if not redundancy.get("all_new_low_support_equations_redundant"):
        raise ValueError("support-17 redundancy certificate is not verified")
    if not set(old_report["selected_variables"]).issubset(names):
        raise ValueError("support-12 variables escape the support-32 ring")

    symbols = sp.symbols(" ".join(names))
    symbol_map = dict(zip(names, symbols))
    name_index = {name: i for i, name in enumerate(names)}

    def weighted_grevlex(monomial):
        return (
            sum(weight*exponent for weight, exponent in zip(weights, monomial)),
            tuple(reversed([-exponent for exponent in monomial])),
        )

    ring, *_ = polynomial_ring(symbols, QQ, order=weighted_grevlex)
    inputs = []
    labels = []
    input_degrees = []

    for family, analysis, tag in (
        (family21, analysis21, "support21"),
        (family32, analysis32, "support32"),
    ):
        selected = analysis["independent_equation_indices_zero_based"]
        if len(selected) != analysis["exact_Q_linear_rank"]:
            raise ValueError(f"{tag} row-space receipt mismatch")
        if tag == "support32":
            redundant_tags = {
                row["tag"] for row in redundancy["equations"]
                if row.get("all_coefficient_components_divisible")
            }
            redundant_selected = [
                family["equations"][index]["tag"] for index in selected
                if family["equations"][index]["tag"] in redundant_tags
            ]
            expected_redundant = {
                "corrected-power-distance17:t^9", "17:0:t^4"
            }
            if set(redundant_selected) != expected_redundant:
                raise ValueError("support-17 divisibility rows do not match the exact receipt")
            selected = [
                index for index in selected
                if family["equations"][index]["tag"] not in redundant_tags
            ]
        for index in selected:
            record = family["equations"][index]
            poly = family_polynomial(record, names, ring)
            inputs.append(poly)
            labels.append(f"{tag}:{record['tag']}")
            input_degrees.append(sorted({
                sum(w*e for w, e in zip(weights, monomial))
                for monomial, _ in poly.terms()
            }))

    for index, record in enumerate(old_report["basis"]):
        expression = sp.sympify(record["polynomial"], locals=symbol_map)
        poly_sympy = sp.Poly(expression, *symbols, domain=QQ)
        poly = ring.from_dict({tuple(monomial): coefficient
                               for monomial, coefficient in poly_sympy.terms()
                               if coefficient})
        inputs.append(poly)
        labels.append(f"support12-G{index}")
        input_degrees.append(sorted({
            sum(w*e for w, e in zip(weights, monomial))
            for monomial, _ in poly.terms()
        }))

    if len(inputs) != 10:
        raise ValueError("expected four support-21, two residual support-32, and four support-12 generators")
    origin = tuple([0]*len(names))
    if any(dict(poly.terms()).get(origin, 0) for poly in inputs):
        raise ValueError("origin is not a common zero of the compressed system")
    if any(len(degrees) != 1 for degrees in input_degrees):
        raise ValueError("combined generators are not weighted homogeneous")

    print(json.dumps({
        "stage": "groebner_start",
        "variables": len(names),
        "generators": len(inputs),
        "labels": labels,
        "weight_degrees": input_degrees,
        "weights": weights,
    }), flush=True)
    basis = ring_groebner(inputs, ring, method="f5b")

    basis_records = []
    pure_powers = {}
    for poly in basis:
        terms = poly.terms()
        leading = terms[0][0] if terms else None
        basis_records.append({
            "polynomial": str(poly),
            "term_count": len(terms),
            "total_degree": max((sum(m) for m, _ in terms), default=0),
            "weighted_degree": max((sum(w*e for w, e in zip(weights, m))
                                    for m, _ in terms), default=0),
            "leading_monomial": list(leading) if leading else None,
        })
        if len(terms) == 1:
            monomial = terms[0][0]
            nonzero = [i for i, exponent in enumerate(monomial) if exponent]
            if len(nonzero) == 1:
                pure_powers[names[nonzero[0]]] = int(monomial[nonzero[0]])

    result = {
        "scope": (
            "Exact weighted F5B basis of the combined support-12, support-21, "
            "and support-32 necessary subsystem after the d19 substitutions. "
            "The inputs vanish at the origin; a nonunit result alone is not a "
            "whole-chart or whole-branch exclusion."
        ),
        "source_hashes": {
            "support21_family": sha256(FAMILY21),
            "support21_analysis": sha256(ANALYSIS21),
            "support32_family": sha256(FAMILY32),
            "support32_analysis": sha256(ANALYSIS32),
            "support12_groebner": sha256(SUPPORT12_GB),
            "support12_input": sha256(SUPPORT12_INPUT),
            "support12_reuse": sha256(SUPPORT12_REUSE),
            "low_support_redundancy": sha256(LOW_SUPPORT_REDUNDANCY),
        },
        "variables": list(names),
        "weights": list(weights),
        "weighted_order": "positive-weight grevlex",
        "support21_equations": family21["equation_count"],
        "support21_rank": analysis21["exact_Q_linear_rank"],
        "support32_equations": family32["equation_count"],
        "support32_rank": analysis32["exact_Q_linear_rank"],
        "support32_rank_mod_support12_ideal": 2,
        "support32_redundant_selected_tags": sorted(expected_redundant),
        "support12_generators": len(old_report["basis"]),
        "input_generators": len(inputs),
        "input_labels": labels,
        "input_weight_degrees": input_degrees,
        "origin_common_zero_verified": True,
        "groebner_method": "f5b",
        "basis_length": len(basis_records),
        "unit_ideal": any(poly == ring.one for poly in basis),
        "pure_power_generators": pure_powers,
        "krull_dimension_lower_bound": len(names)-len(inputs),
        "basis": basis_records,
        "elapsed_seconds": time.monotonic()-started,
    }
    OUTPUT.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "input_generators", "input_labels", "input_weight_degrees",
        "basis_length", "unit_ideal", "pure_power_generators",
        "krull_dimension_lower_bound", "elapsed_seconds",
    )}, indent=2), flush=True)
    print(f"report={OUTPUT}", flush=True)


if __name__ == "__main__":
    main()
