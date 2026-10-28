"""Exact weighted F5B basis of the four-generator support-21 subsystem."""

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
OUTPUT = ROOT / "f2_75_125_211_d19_support21_only_gb_20261025.json"


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
    names = tuple(family["target_support"])
    weights = tuple(int(value)
                    for value in analysis["primitive_positive_integer_weights"])
    if len(names) != 21 or len(weights) != 21:
        raise ValueError("support-21 metadata mismatch")
    if not analysis.get("all_equations_weight_homogeneous"):
        raise ValueError("the positive grading is not verified")

    symbols = sp.symbols(" ".join(names))
    name_index = {name: i for i, name in enumerate(names)}

    def weighted_grevlex(monomial):
        return (
            sum(weight*exponent for weight, exponent in zip(weights, monomial)),
            tuple(reversed([-exponent for exponent in monomial])),
        )

    ring, *_ = polynomial_ring(symbols, QQ, order=weighted_grevlex)
    independent = analysis["independent_equation_indices_zero_based"]
    if len(independent) != 4:
        raise ValueError("expected four exact Q-linear generators")
    inputs = []
    tags = []
    for index in independent:
        equation = family["equations"][index]
        coefficients = {}
        for monomial, coefficient in equation["terms"]:
            powers = [0]*len(names)
            for name, exponent in monomial:
                powers[name_index[name]] += int(exponent)
            key = tuple(powers)
            coefficients[key] = coefficients.get(key, sp.Rational(0)) + sp.Rational(coefficient)
        coefficients = {key: value for key, value in coefficients.items() if value}
        inputs.append(ring.from_dict(coefficients))
        tags.append(equation["tag"])

    origin = tuple([0]*len(names))
    if any(dict(poly.terms()).get(origin, 0) for poly in inputs):
        raise ValueError("support-21 generator does not vanish at origin")
    print(json.dumps({
        "stage": "groebner_start", "variables": len(names),
        "generators": len(inputs), "tags": tags, "weights": weights,
    }), flush=True)

    basis = ring_groebner(inputs, ring, method="f5b")
    basis_records = []
    pure_powers = {}
    for poly in basis:
        terms = poly.terms()
        basis_records.append({
            "polynomial": str(poly),
            "term_count": len(terms),
            "total_degree": max((sum(m) for m, _ in terms), default=0),
            "weighted_degree": max((sum(w*e for w, e in zip(weights, m))
                                    for m, _ in terms), default=0),
        })
        if len(terms) == 1:
            monomial = terms[0][0]
            nonzero = [i for i, exponent in enumerate(monomial) if exponent]
            if len(nonzero) == 1:
                pure_powers[names[nonzero[0]]] = int(monomial[nonzero[0]])

    result = {
        "scope": (
            "Exact weighted F5B Groebner basis of the four Q-linearly independent "
            "rows spanning the complete support-21 family. This subsystem alone "
            "is homogeneous and has the origin as a zero, so it cannot exclude "
            "the entire F2 branch."
        ),
        "family_sha256": sha256(FAMILY),
        "analysis_sha256": sha256(ANALYSIS),
        "variables": list(names),
        "weights": list(weights),
        "weighted_monomial_order": "positive-weight grevlex",
        "equations": len(family["equations"]),
        "exact_Q_linear_rank": analysis["exact_Q_linear_rank"],
        "independent_tags": tags,
        "generator_count": len(inputs),
        "origin_common_zero_verified": True,
        "groebner_method": "f5b",
        "basis_length": len(basis_records),
        "unit_ideal": any(poly == ring.one for poly in basis),
        "pure_power_generators": pure_powers,
        "basis": basis_records,
        "krull_dimension_lower_bound": len(names)-len(inputs),
        "elapsed_seconds": time.monotonic()-started,
    }
    OUTPUT.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "generator_count", "basis_length", "unit_ideal",
        "pure_power_generators", "krull_dimension_lower_bound",
        "elapsed_seconds",
    )}, indent=2), flush=True)
    print(f"report={OUTPUT}", flush=True)


if __name__ == "__main__":
    main()
