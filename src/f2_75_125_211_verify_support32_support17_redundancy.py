"""Independently verify support-32 rows 81 and 82 lie in the support-12 ideal."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import time

import sympy as sp
from sympy.polys.domains import QQ


ROOT = Path("research")
FAMILY = ROOT / "f2_75_125_211_d19_support32_family_input_20261025.json"
ANALYSIS = ROOT / "f2_75_125_211_d19_support32_analysis_20261025.json"
SUPPORT12 = ROOT / "f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_gb_4gen_weighted_flint.json"
PRIOR_RECEIPT = ROOT / "f2_75_125_211_new_small_support_redundancy_20261025.json"
OUTPUT = ROOT / "f2_75_125_211_support32_support17_redundancy_verify_20261025.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def record_poly(record, symbols, index):
    expression = sp.S.Zero
    for monomial, coefficient in record["terms"]:
        term = sp.Rational(coefficient)
        for name, power in monomial:
            term *= symbols[index[name]]**int(power)
        expression += term
    return sp.Poly(expression, *symbols, domain=QQ)


def main() -> None:
    started = time.monotonic()
    family = json.loads(FAMILY.read_text(encoding="utf-8"))
    analysis = json.loads(ANALYSIS.read_text(encoding="utf-8"))
    basis_data = json.loads(SUPPORT12.read_text(encoding="utf-8"))
    prior = json.loads(PRIOR_RECEIPT.read_text(encoding="utf-8"))
    names = tuple(analysis["selected_variables"])
    symbols = sp.symbols(" ".join(names))
    symbol_index = dict(zip(names, range(len(names))))
    divisor_index = int(prior["divisor_basis_index_zero_based"])
    divisor_expr = sp.sympify(
        basis_data["basis"][divisor_index]["polynomial"],
        locals=dict(zip(names, symbols)),
    )
    divisor = sp.Poly(divisor_expr, *symbols, domain=QQ)
    if len(divisor.terms()) != int(prior["divisor_term_count"]):
        raise ValueError("support-12 divisor term count mismatch")

    target_tags = {
        "corrected-power-distance17:t^9", "17:0:t^4"
    }
    rows = [record for record in family["equations"] if record["tag"] in target_tags]
    if {record["tag"] for record in rows} != target_tags:
        raise ValueError("the two support-17 rows are missing")
    results = []
    for record in rows:
        polynomial = record_poly(record, symbols, symbol_index)
        quotient, remainder = polynomial.div(divisor)
        if remainder != 0 or quotient*divisor != polynomial:
            raise AssertionError(f"exact divisibility failed for {record['tag']}")
        results.append({
            "tag": record["tag"],
            "term_count": len(polynomial.terms()),
            "divisor_basis_index_zero_based": divisor_index,
            "divisor_term_count": len(divisor.terms()),
            "quotient": str(quotient.as_expr()),
            "quotient_term_count": len(quotient.terms()),
            "remainder_zero": True,
        })

    result = {
        "scope": (
            "Direct exact-QQ multivariate division of both nonzero support-17 "
            "support-32 equations by the designated support-12 Groebner basis "
            "element. This independently verifies the ideal redundancy used "
            "to remove these two rows from the joint Groebner input."
        ),
        "input_hashes": {
            "support32_family": sha256(FAMILY),
            "support32_analysis": sha256(ANALYSIS),
            "support12_groebner": sha256(SUPPORT12),
            "prior_redundancy_receipt": sha256(PRIOR_RECEIPT),
        },
        "support32_rank": analysis["exact_Q_linear_rank"],
        "support32_selected_indices": analysis["independent_equation_indices_zero_based"],
        "support17_rows": results,
        "both_divide_exactly": all(row["remainder_zero"] for row in results),
        "elapsed_seconds": time.monotonic()-started,
    }
    OUTPUT.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
