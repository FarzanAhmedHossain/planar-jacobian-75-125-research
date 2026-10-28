"""Exact weighted F5B basis of the compressed support-12/21/32/38 system."""

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
FAMILY38 = ROOT / "f2_75_125_211_d19_support38_family_input_20261025.json"
ANALYSIS38 = ROOT / "f2_75_125_211_d19_support38_analysis_20261025.json"
SUPPORT32_REDUNDANCY = ROOT / "f2_75_125_211_support32_support17_redundancy_verify_20261025.json"
SUPPORT38_LOW_REDUCTION = ROOT / "f2_75_125_211_support38_low_rows_mod_support12_20261025.json"
SUPPORT12_GB = ROOT / "f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_gb_4gen_weighted_flint.json"
SUPPORT12_INPUT = ROOT / "f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_input.json"
SUPPORT12_REUSE = ROOT / "f2_75_125_211_verify_d19_support12_reuse_20261025.json"
OUTPUT = ROOT / "f2_75_125_211_d19_support12_21_32_38_combined_gb_20261025.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def coefficient_for_ring(value, ring):
    """Coerce a rational coefficient, including by exact reduction to GF(p)."""
    modulus = getattr(ring.domain, "mod", None)
    if modulus is None:
        return value
    rational = sp.Rational(value)
    prime = int(modulus)
    residue = (int(rational.p) * pow(int(rational.q), -1, prime)) % prime
    return ring.domain.convert(residue)


def family_polynomial(record, names, ring):
    index = {name: i for i, name in enumerate(names)}
    coefficients = {}
    for monomial, coefficient in record["terms"]:
        exponents = [0]*len(names)
        for name, power in monomial:
            exponents[index[name]] += int(power)
        key = tuple(exponents)
        coefficients[key] = coefficients.get(key, sp.Rational(0)) + sp.Rational(coefficient)
    return ring.from_dict({
        key: coefficient_for_ring(value, ring)
        for key, value in coefficients.items() if value
    })


def main() -> None:
    started = time.monotonic()
    family21 = json.loads(FAMILY21.read_text(encoding="utf-8"))
    analysis21 = json.loads(ANALYSIS21.read_text(encoding="utf-8"))
    family32 = json.loads(FAMILY32.read_text(encoding="utf-8"))
    analysis32 = json.loads(ANALYSIS32.read_text(encoding="utf-8"))
    family38 = json.loads(FAMILY38.read_text(encoding="utf-8"))
    analysis38 = json.loads(ANALYSIS38.read_text(encoding="utf-8"))
    redundancy32 = json.loads(SUPPORT32_REDUNDANCY.read_text(encoding="utf-8"))
    reduction38 = json.loads(SUPPORT38_LOW_REDUCTION.read_text(encoding="utf-8"))
    support12 = json.loads(SUPPORT12_GB.read_text(encoding="utf-8"))
    support12_input = json.loads(SUPPORT12_INPUT.read_text(encoding="utf-8"))
    support12_reuse = json.loads(SUPPORT12_REUSE.read_text(encoding="utf-8"))

    names = tuple(analysis38["selected_variables"])
    weights = tuple(int(v) for v in analysis38["primitive_positive_integer_weights"])
    index = {name: i for i, name in enumerate(names)}
    if len(names) != 28 or len(weights) != 28:
        raise ValueError("support-38 grading metadata mismatch")
    if not analysis38["all_equations_weight_homogeneous"]:
        raise ValueError("support-38 weights are not verified")
    if analysis21["exact_Q_linear_rank"] != 4 or analysis32["exact_Q_linear_rank"] != 4:
        raise ValueError("earlier family rank receipts mismatch")
    if not redundancy32.get("both_divide_exactly"):
        raise ValueError("support-32 support-17 divisibility is not verified")
    if not reduction38.get("support38_support12_rows_exactly_match_old_input"):
        raise ValueError("support-38 support-12 equalities are not verified")
    expected_low_tags = {
        "corrected-power-distance18:t^9", "18:0:t^4"
    }
    low_rows = reduction38.get("rows", [])
    if (reduction38.get("support38_support22_rows_reduced") != 2
            or {row.get("tag") for row in low_rows} != expected_low_tags
            or any(row.get("normal_form_zero") is not False
                   or int(row.get("normal_form_term_count", 0)) <= 0
                   for row in low_rows)):
        raise ValueError(
            "support-38 support-22 reductions are missing, stale, or not "
            "verified nonzero; retain both rows as generators"
        )
    for key, path in (
        ("support38_family", FAMILY38),
        ("support38_analysis", ANALYSIS38),
        ("support12_groebner", SUPPORT12_GB),
        ("support12_input", SUPPORT12_INPUT),
        ("support12_reuse", SUPPORT12_REUSE),
    ):
        if reduction38.get("input_hashes", {}).get(key) != sha256(path):
            raise ValueError(f"support-38 low-row receipt has stale {key} input")
    if not support12_reuse.get("same_142_equations_and_full_support12_group"):
        raise ValueError("full support-12 reuse receipt is not verified")
    if not set(support12["selected_variables"]).issubset(names):
        raise ValueError("support-12 variables escape support-38 ring")

    for analysis in (analysis21, analysis32):
        for name, weight in zip(analysis["selected_variables"],
                                analysis["primitive_positive_integer_weights"]):
            if weights[index[name]] != int(weight):
                raise ValueError("nested support gradings disagree")

    symbols = sp.symbols(" ".join(names))
    symbol_map = dict(zip(names, symbols))

    def weighted_grevlex(monomial):
        return (
            sum(w*e for w, e in zip(weights, monomial)),
            tuple(reversed([-e for e in monomial])),
        )

    ring, *_ = polynomial_ring(symbols, QQ, order=weighted_grevlex)
    inputs = []
    labels = []
    degrees = []

    def append_family(family, analysis, label, omit_tags=frozenset()):
        selected = analysis["independent_equation_indices_zero_based"]
        if len(selected) != analysis["exact_Q_linear_rank"]:
            raise ValueError(f"{label} row-rank receipt mismatch")
        retained = []
        for row_index in selected:
            record = family["equations"][row_index]
            if record["tag"] in omit_tags:
                continue
            poly = family_polynomial(record, names, ring)
            inputs.append(poly)
            labels.append(f"{label}:{record['tag']}")
            degrees.append(sorted({
                sum(w*e for w, e in zip(weights, monomial))
                for monomial, _ in poly.terms()
            }))
            retained.append(record["tag"])
        return retained

    append_family(family21, analysis21, "support21")
    support32_omit = {
        "corrected-power-distance17:t^9", "17:0:t^4"
    }
    retained32 = append_family(family32, analysis32, "support32", support32_omit)
    if len(retained32) != 2:
        raise ValueError("expected two support-32 rows after ideal reduction")

    support38_support12_tags = {
        row["tag"] for row in family38["equations"]
        if len(row["reduced_support"]) == 12
    }
    old_support12_tags = {row["tag"] for row in support12_input["equations"]}
    if len(support38_support12_tags) != 34 or not support38_support12_tags.issubset(old_support12_tags):
        raise ValueError("support-38 support-12 tags are not all in the old input")
    selected38 = analysis38["independent_equation_indices_zero_based"]
    if len(selected38) != analysis38["exact_Q_linear_rank"]:
        raise ValueError("support-38 row-rank receipt mismatch")
    retained38 = []
    for row_index in selected38:
        record = family38["equations"][row_index]
        if record["tag"] in support38_support12_tags:
            continue
        poly = family_polynomial(record, names, ring)
        inputs.append(poly)
        labels.append(f"support38:{record['tag']}")
        degrees.append(sorted({
            sum(w*e for w, e in zip(weights, monomial))
            for monomial, _ in poly.terms()
        }))
        retained38.append(record["tag"])
    if len(retained38) != 8:
        raise ValueError(f"expected eight support-38 rows modulo support-12, got {len(retained38)}")

    for index, record in enumerate(support12["basis"]):
        expression = sp.sympify(record["polynomial"], locals=symbol_map)
        poly_sympy = sp.Poly(expression, *symbols, domain=sp.QQ)
        poly = ring.from_dict({
            tuple(m): coefficient_for_ring(c, ring)
            for m, c in poly_sympy.terms() if c
        })
        inputs.append(poly)
        labels.append(f"support12-G{index}")
        degrees.append(sorted({
            sum(w*e for w, e in zip(weights, monomial))
            for monomial, _ in poly.terms()
        }))

    if len(inputs) != 18:
        raise ValueError(f"expected 18 compressed generators, got {len(inputs)}")
    origin = tuple([0]*len(names))
    if any(dict(poly.terms()).get(origin, 0) for poly in inputs):
        raise ValueError("origin is not a common zero of the compressed system")
    if any(len(ds) != 1 for ds in degrees):
        raise ValueError("a combined input is not weight homogeneous")

    print(json.dumps({
        "stage": "groebner_start",
        "variables": len(names),
        "generators": len(inputs),
        "labels": labels,
        "weight_degrees": degrees,
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
            "Exact weighted F5B basis for the compressed support-12/21/32/38 "
            "necessary subsystem after d19 substitutions. It still contains "
            "the origin and is not by itself a branch exclusion."
        ),
        "source_hashes": {
            "support21_family": sha256(FAMILY21),
            "support21_analysis": sha256(ANALYSIS21),
            "support32_family": sha256(FAMILY32),
            "support32_analysis": sha256(ANALYSIS32),
            "support38_family": sha256(FAMILY38),
            "support38_analysis": sha256(ANALYSIS38),
            "support32_redundancy": sha256(SUPPORT32_REDUNDANCY),
            "support38_low_reduction": sha256(SUPPORT38_LOW_REDUCTION),
            "support12_groebner": sha256(SUPPORT12_GB),
            "support12_input": sha256(SUPPORT12_INPUT),
            "support12_reuse": sha256(SUPPORT12_REUSE),
        },
        "variables": list(names),
        "weights": list(weights),
        "weighted_order": "positive-weight grevlex",
        "support21_equations": family21["equation_count"],
        "support21_rank": analysis21["exact_Q_linear_rank"],
        "support32_equations": family32["equation_count"],
        "support32_rank": analysis32["exact_Q_linear_rank"],
        "support32_generators_after_support12_reduction": len(retained32),
        "support38_equations": family38["equation_count"],
        "support38_rank": analysis38["exact_Q_linear_rank"],
        "support38_generators_after_support12_reduction": len(retained38),
        "support12_generators": len(support12["basis"]),
        "input_generators": len(inputs),
        "input_labels": labels,
        "input_weight_degrees": degrees,
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
