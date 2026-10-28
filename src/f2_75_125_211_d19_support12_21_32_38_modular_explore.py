#!/usr/bin/env python3
"""Modular reconnaissance for the exact d19 support-12/21/32/38 ideal.

This is only a prime-field exploration.  It cannot certify a characteristic-zero
branch exclusion; any useful leading-monomial pattern must be replayed over QQ.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import sympy as sp
from sympy.polys.domains import GF


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import f2_75_125_211_d19_support12_21_32_38_combined_groebner as builder


PRIME = 32003
OUTPUT = HERE / "f2_75_125_211_d19_support12_21_32_38_mod32003_20261028.json"


def main() -> None:
    if not sp.isprime(PRIME):
        raise ValueError("selected modulus is not prime")
    builder.QQ = GF(PRIME)
    builder.OUTPUT = OUTPUT
    builder.main()

    result = json.loads(OUTPUT.read_text(encoding="utf-8"))
    result["scope"] = (
        "prime-field F5B exploration of the support-12/21/32/38 necessary subsystem; "
        "not a characteristic-zero certificate and not a branch exclusion"
    )
    result["coefficient_field"] = f"GF({PRIME})"
    result["modular_caveat"] = (
        "The ideal contains the origin. A modular basis or projective dimension result "
        "is reconnaissance only; exact characteristic-zero replay is required."
    )
    result["annotated_report_payload_sha256"] = hashlib.sha256(
        json.dumps(result, indent=2).encode("utf-8")
    ).hexdigest()
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result.get(key) for key in (
        "scope", "coefficient_field", "variables", "input_generators",
        "basis_length", "pure_power_generators", "krull_dimension_lower_bound", "unit_ideal",
        "modular_caveat"
    )}, indent=2), flush=True)
    print(f"report={OUTPUT}", flush=True)


if __name__ == "__main__":
    main()
