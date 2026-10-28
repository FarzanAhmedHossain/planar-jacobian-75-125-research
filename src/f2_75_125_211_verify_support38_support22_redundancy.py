"""Exact-QQ divisibility check for support-38's two support-22 rows."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import time

import sympy as sp
from sympy.polys.domains import QQ


ROOT = Path("research")
FAMILY = ROOT / "f2_75_125_211_d19_support38_family_input_20261025.json"
ANALYSIS = ROOT / "f2_75_125_211_d19_support38_analysis_20261025.json"
SUPPORT12 = ROOT / "f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_gb_4gen_weighted_flint.json"
OUTPUT = ROOT / "f2_75_125_211_support38_support22_redundancy_verify_20261025.json"


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
    names = tuple(analysis["selected_variables"])
    symbols = sp.symbols(" ".join(names))
    symbol_index = dict(zip(names, range(len(names))))
    divisor_index = 3
    divisor = sp.Poly(
        sp.sympify(basis_data["basis"][divisor_index]["polynomial"],
                   locals=dict(zip(names, symbols))),
        *symbols, domain=QQ,
    )
    selected_indices = analysis["independent_equation_indices_zero_based"]
    rows = [family["equations"][i] for i in selected_indices
            if len(family["equations"][i]["reduced_support"]) == 22]
    if len(rows) != 2:
        raise ValueError("expected exactly two independent support-22 rows")
    results = []
    for record in rows:
        polynomial = record_poly(record, symbols, symbol_index)
        quotient, remainder = polynomial.div(divisor)
        results.append({
            "tag": record["tag"],
            "term_count": len(polynomial.terms()),
            "divisor_index_zero_based": divisor_index,
            "divisor_term_count": len(divisor.terms()),
            "quotient": str(quotient.as_expr()),
            "quotient_term_count": len(quotient.terms()),
            "remainder_zero": remainder == 0,
        })
        if remainder != 0 or quotient*divisor != polynomial:
            raise AssertionError(f"support-12 divisibility failed for {record['tag']}")
    result = {
        "scope": (
            "Direct exact-QQ multivariate division of support-38's two "
            "support-22 independent rows by support-12 basis element G."
        ),
        "input_hashes": {
            "support38_family": sha256(FAMILY),
            "support38_analysis": sha256(ANALYSIS),
            "support12_groebner": sha256(SUPPORT12),
        },
        "support38_rank": analysis["exact_Q_linear_rank"],
        "support38_independent_indices": selected_indices,
        "rows": results,
        "both_divide_exactly": all(row["remainder_zero"] for row in results),
        "elapsed_seconds": time.monotonic()-started,
    }
    OUTPUT.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
