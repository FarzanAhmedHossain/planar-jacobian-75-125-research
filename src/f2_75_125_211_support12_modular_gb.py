#!/usr/bin/env python3
"""Hash-bound finite-field Groebner reconnaissance for the d18 support-12 slice.

The result is a solver aid only. It is not a characteristic-zero exclusion.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from fractions import Fraction
from pathlib import Path

from sympy.polys.domains import GF
from sympy.polys.groebnertools import groebner as ring_groebner
from sympy.polys.rings import ring as polynomial_ring


INPUT = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_input.json")
EXPECTED_INPUT_SHA256 = "A9DF8405D43EAFFB98C724165C0C7A647C0952D7D8D1A8630CD26700898C6F56"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest().upper()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("report", type=Path)
    ap.add_argument("--prime", type=int, default=65521)
    ap.add_argument("--order", choices=("grevlex", "lex"), default="grevlex")
    ap.add_argument("--method", choices=("f5b", "buchberger"), default="f5b")
    args = ap.parse_args()

    input_hash = sha256(INPUT)
    if input_hash != EXPECTED_INPUT_SHA256:
        raise ValueError(f"support-12 input hash mismatch: {input_hash}")
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    names = tuple(data["selected_variables"])
    if len(names) != 12 or data["equation_count"] != len(data["equations"]):
        raise ValueError("support-12 extraction metadata mismatch")
    if args.prime < 2:
        raise ValueError("prime must be at least 2")

    p = args.prime
    domain = GF(p)
    index = {name: i for i, name in enumerate(names)}
    symbols = __import__("sympy").symbols(" ".join(names))
    raw_polys = []
    terms_read = 0
    for row in data["equations"]:
        coeffs: dict[tuple[int, ...], int] = {}
        for monomial, coefficient in row["terms"]:
            exp = [0] * len(names)
            for name, power in monomial:
                if name not in index:
                    raise ValueError(f"unexpected variable {name} in {row['tag']}")
                exp[index[name]] += int(power)
            q = Fraction(coefficient)
            if q.denominator % p == 0:
                raise ValueError(f"prime {p} divides a coefficient denominator")
            value = q.numerator * pow(q.denominator, -1, p) % p
            key = tuple(exp)
            coeffs[key] = (coeffs.get(key, 0) + value) % p
            terms_read += 1
        coeffs = {m: c for m, c in coeffs.items() if c % p}
        if (0,) * len(names) in coeffs:
            raise ValueError(f"claimed origin witness fails modulo {p}: {row['tag']}")
        raw_polys.append(coeffs)

    ring, *_ = polynomial_ring(symbols, domain, order=args.order)
    unique = {}
    for coeffs in raw_polys:
        poly = ring.from_dict(coeffs)
        if not poly:
            continue
        monic = poly.monic()
        unique.setdefault(tuple(monic.terms()), monic)
    inputs = list(unique.values())

    start = time.perf_counter()
    basis = ring_groebner(inputs, ring, method=args.method)
    basis_records = []
    for poly in basis:
        terms = poly.terms()
        basis_records.append({
            "leading_monomial": list(terms[0][0]) if terms else None,
            "total_degree": max((sum(m) for m, _ in terms), default=0),
            "term_count": len(terms),
        })
    result = {
        "input": str(INPUT),
        "input_sha256": input_hash,
        "branch": data["branch"],
        "selected_variables": list(names),
        "equation_count": len(raw_polys),
        "unique_nonzero_generators_up_to_scaling_mod_p": len(inputs),
        "input_terms": terms_read,
        "field": f"GF({p})",
        "groebner_method": args.method,
        "monomial_order": args.order,
        "basis_length": len(basis_records),
        "unit_ideal": any(poly == ring.one for poly in basis),
        "origin_witness": {name: 0 for name in names},
        "basis_leading_data": basis_records,
        "elapsed_seconds": round(time.perf_counter() - start, 3),
        "scope": "Finite-field reconnaissance on the exact hash-bound support-12 equations only; not a characteristic-zero ideal certificate or a branch exclusion.",
    }
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
