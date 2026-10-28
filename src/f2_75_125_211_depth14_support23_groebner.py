#!/usr/bin/env python3
"""Groebner reconnaissance on the 35 depth-14 compatibility equations
supported on p0,...,p24 (excluding p5,p13).

The modular mode is a structural probe only.  It never claims a Q-exclusion.
If it finds only the origin, an exact QQ replay is still required before that
conclusion can be used.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from fractions import Fraction
from pathlib import Path

import ijson
import sympy as sp
from sympy.polys.domains import GF, QQ


SOURCE = Path("research/f2_75_125_211_depth14_global_compatibility_reduction_20261021.json")
OUT = Path("research/f2_75_125_211_depth14_support23_groebner.json")
ROWS = set(range(86, 121))
SUPPORT = [f"p{i}" for i in range(25) if i not in (5, 13)]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--field", choices=("modular", "QQ"), default="modular")
    parser.add_argument("--prime", type=int, default=32749)
    parser.add_argument("--nilpotence-cap", type=int, default=32)
    args = parser.parse_args()

    source_hash = sha256_file(SOURCE)
    symbols = sp.symbols(" ".join(SUPPORT))
    support_pos = {int(name[1:]): j for j, name in enumerate(SUPPORT)}
    polynomials = []
    tags = []
    term_counts = []
    with SOURCE.open("rb") as stream:
        for row in ijson.items(stream, "global_compatibility_equations.item"):
            if row["constraint_row_zero_based"] not in ROWS:
                continue
            if row["support_variables"] != SUPPORT:
                raise AssertionError(f"support mismatch in {row['constraint_tag']}")
            terms = {}
            for sparse, coefficient in row["terms"]:
                exponents = [0] * len(SUPPORT)
                for variable, power in sparse:
                    name = f"p{variable}"
                    if variable not in support_pos:
                        raise AssertionError(f"unexpected variable {name}")
                    exponents[support_pos[variable]] += power
                q = Fraction(coefficient)
                if args.field == "modular":
                    denominator = q.denominator % args.prime
                    if denominator == 0:
                        raise ValueError(f"prime divides denominator in {row['constraint_tag']}")
                    value = (q.numerator % args.prime) * pow(denominator, -1, args.prime) % args.prime
                else:
                    value = sp.Rational(q.numerator, q.denominator)
                terms[tuple(exponents)] = value
            domain = GF(args.prime) if args.field == "modular" else QQ
            polynomial = sp.Poly.from_dict(terms, *symbols, domain=domain)
            if not polynomial or polynomial.total_degree() != 14:
                raise AssertionError(f"unexpected polynomial in {row['constraint_tag']}")
            polynomials.append(polynomial.as_expr())
            tags.append(row["constraint_tag"])
            term_counts.append(len(terms))

    if len(polynomials) != 35:
        raise AssertionError(f"expected 35 support-23 equations, found {len(polynomials)}")

    start = time.monotonic()
    print(f"field={args.field} prime={args.prime if args.field == 'modular' else 'QQ'} source_sha256={source_hash}", flush=True)
    print(f"selected={len(polynomials)} variables={len(symbols)} input_terms={sum(term_counts)} starting Groebner", flush=True)
    domain = GF(args.prime) if args.field == "modular" else QQ
    basis = sp.groebner(polynomials, *symbols, domain=domain, order="grevlex", method="f5b")
    elapsed = time.monotonic() - start
    zero_dimensional = basis.is_zero_dimensional
    print(f"basis_size={len(basis.polys)} zero_dimensional={zero_dimensional} elapsed={elapsed:.2f}s", flush=True)

    powers = {}
    if zero_dimensional:
        for symbol in symbols:
            found = None
            for power in range(1, args.nilpotence_cap + 1):
                if basis.reduce(symbol**power)[1] == 0:
                    found = power
                    break
            powers[str(symbol)] = found
            print(f"nilpotence={symbol}^{found if found else 'not-found-through-'+str(args.nilpotence_cap)}", flush=True)

    result = {
        "source": str(SOURCE),
        "source_sha256": source_hash,
        "scope": "35 support-23 equations from the global depth-14 compatibility reduction; modular mode is reconnaissance, not exclusion",
        "field": "GF(prime)" if args.field == "modular" else "Q",
        "prime": args.prime if args.field == "modular" else None,
        "support_variables": SUPPORT,
        "selected_constraints": len(polynomials),
        "selected_tags": tags,
        "input_term_count": sum(term_counts),
        "monomial_order": "grevlex",
        "algorithm": "SymPy F5B",
        "basis_size": len(basis.polys),
        "zero_dimensional": zero_dimensional,
        "variable_nilpotence_powers_through_cap": powers,
        "groebner_seconds": elapsed,
        "warning": "The equations form only a necessary subideal of the 2+1+1 branch. A modular result is not a characteristic-zero theorem; even an exact origin-only subideal would only set these 23 parameters to zero and would leave the remaining variables to analyze.",
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"receipt={OUT}", flush=True)


if __name__ == "__main__":
    main()
