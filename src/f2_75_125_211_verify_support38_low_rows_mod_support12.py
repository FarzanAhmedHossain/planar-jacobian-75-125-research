"""Reduce all support-12/22 support-38 rows modulo the support-12 GB."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import time

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.rings import ring as polynomial_ring


ROOT = Path("research")
FAMILY = ROOT / "f2_75_125_211_d19_support38_family_input_20261025.json"
ANALYSIS = ROOT / "f2_75_125_211_d19_support38_analysis_20261025.json"
SUPPORT12 = ROOT / "f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_gb_4gen_weighted_flint.json"
SUPPORT12_INPUT = ROOT / "f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_input.json"
SUPPORT12_REUSE = ROOT / "f2_75_125_211_verify_d19_support12_reuse_20261025.json"
OUTPUT = ROOT / "f2_75_125_211_support38_low_rows_mod_support12_20261025.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def family_poly(record, names, ring):
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
    family = json.loads(FAMILY.read_text(encoding="utf-8"))
    analysis = json.loads(ANALYSIS.read_text(encoding="utf-8"))
    support12 = json.loads(SUPPORT12.read_text(encoding="utf-8"))
    support12_input = json.loads(SUPPORT12_INPUT.read_text(encoding="utf-8"))
    reuse = json.loads(SUPPORT12_REUSE.read_text(encoding="utf-8"))
    if not reuse.get("same_142_equations_and_full_support12_group"):
        raise ValueError("full d19/d18 support-12 reuse receipt is missing")
    old_rows = {row["tag"]: row for row in support12_input["equations"]}
    support12_rows = [row for row in family["equations"]
                      if len(row["reduced_support"]) == 12]

    def canonical_terms(record):
        return {
            tuple(sorted((name, int(power)) for name, power in monomial if power)):
            sp.Rational(coefficient)
            for monomial, coefficient in record["terms"]
        }

    support12_mismatches = []
    for row in support12_rows:
        old = old_rows.get(row["tag"])
        if old is None or canonical_terms(row) != canonical_terms(old):
            support12_mismatches.append(row["tag"])
    if len(support12_rows) != 34 or support12_mismatches:
        raise ValueError("support-38 support-12 rows do not match the verified old group")
    names = tuple(analysis["selected_variables"])
    weights = tuple(int(value) for value in analysis["primitive_positive_integer_weights"])
    symbols = sp.symbols(" ".join(names))
    symbol_map = dict(zip(names, symbols))

    def weighted_grevlex(monomial):
        return (
            sum(w*e for w, e in zip(weights, monomial)),
            tuple(reversed([-e for e in monomial])),
        )

    ring, *_ = polynomial_ring(symbols, QQ, order=weighted_grevlex)
    basis = []
    for row in support12["basis"]:
        expression = sp.sympify(row["polynomial"], locals=symbol_map)
        poly_sympy = sp.Poly(expression, *symbols, domain=QQ)
        basis.append(ring.from_dict({tuple(m): c for m, c in poly_sympy.terms() if c}))

    targets = [row for row in family["equations"]
               if len(row["reduced_support"]) == 22]
    if len(targets) != 2:
        raise ValueError(f"expected two support-22 rows, found {len(targets)}")
    results = []
    for record in targets:
        polynomial = family_poly(record, names, ring)
        quotients, remainder = polynomial.div(tuple(basis))
        results.append({
            "tag": record["tag"],
            "support_size": len(record["reduced_support"]),
            "term_count": record["term_count"],
            "normal_form_zero": remainder == ring.zero,
            "normal_form_term_count": len(remainder.terms()),
            "normal_form": str(remainder),
        })
    result = {
        "scope": (
            "Exact weighted-grevlex reduction of support-38's two support-22 "
            "rows modulo the support-12 Groebner basis, extended to the "
            "28-variable ring. The 34 support-12 rows are separately verified "
            "identical to rows in the old support-12 input."
        ),
        "input_hashes": {
            "support38_family": sha256(FAMILY),
            "support38_analysis": sha256(ANALYSIS),
            "support12_groebner": sha256(SUPPORT12),
            "support12_input": sha256(SUPPORT12_INPUT),
            "support12_reuse": sha256(SUPPORT12_REUSE),
        },
        "support38_support12_rows": len(support12_rows),
        "support38_support12_rows_exactly_match_old_input": not support12_mismatches,
        "support12_reuse_receipt_verified": True,
        "support38_support22_rows_reduced": len(results),
        "support12_groebner_basis_length": len(basis),
        "rows": results,
        "all_support22_normal_forms_zero": all(row["normal_form_zero"] for row in results),
        "elapsed_seconds": time.monotonic()-started,
    }
    OUTPUT.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({
        "support38_support12_rows": len(support12_rows),
        "support38_support22_rows_reduced": len(results),
        "zero_normal_forms": sum(row["normal_form_zero"] for row in results),
        "nonzero_normal_forms": sum(not row["normal_form_zero"] for row in results),
        "nonzero_rows": [
            {key: row[key] for key in ("tag", "support_size", "normal_form_term_count")}
            for row in results if not row["normal_form_zero"]
        ],
        "elapsed_seconds": result["elapsed_seconds"],
        "output": str(OUTPUT),
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
