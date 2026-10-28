#!/usr/bin/env python3
"""Exact F5B Groebner computation for the smallest-support U1=0 slice.

The input is the immutable set of all 143 complete equations supported on
{p36,p37,p38,p39,p40,p41,p7,r}, extracted from the exact d17 checkpoint
after the certified specializations p79=p0=p1=0 and the exact strict5
set-theoretic reductions p15=452*r/5, p16=-14*r, p17=r.

A unit ideal excludes this slice of the U1=0 subcase. A nonunit result is
only a diagnostic on this necessary subsystem; it does not exclude the full
slice or the whole 2+1+1 branch.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

from sympy.polys.domains import QQ
from sympy.polys.groebnertools import groebner as ring_groebner
from sympy.polys.rings import ring as polynomial_ring


INPUT = Path(
    "research/f2_75_125_sparse_ring_lift_211_roots17_dd_20261024_"
    "d17_u1zero_support8_inputs.json"
)
OUTPUT = Path(
    "research/f2_75_125_sparse_ring_lift_211_roots17_dd_20261024_"
    "d17_u1zero_support8_gb.json"
)
EXPECTED_INPUT_SHA256 = (
    "1676730e3ab0714c51202ff67ee16a1a6368bb5095b79c1156549310e2833678"
)
VARIABLES = ("p36", "p37", "p38", "p39", "p40", "p41", "p7", "r")
GROUP_KEY = ",".join(VARIABLES)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    input_hash = sha256(INPUT)
    if input_hash != EXPECTED_INPUT_SHA256:
        raise ValueError(f"support8 input SHA-256 mismatch: {input_hash}")
    document = json.loads(INPUT.read_text(encoding="utf-8"))
    group = document["smallest_support_groups"][GROUP_KEY]
    if tuple(group["variables"]) != VARIABLES:
        raise ValueError("support8 variable ordering changed")
    equations = group["equations"]
    if len(equations) != 143:
        raise ValueError(f"expected 143 complete equations, got {len(equations)}")

    ring, *gens = polynomial_ring(VARIABLES, QQ, order="grevlex")
    unique: dict[tuple, tuple[str, object]] = {}
    total_terms = 0
    for equation in equations:
        encoded = {}
        for monomial_data, coefficient_text in equation["terms"]:
            exponents = [0] * len(VARIABLES)
            for variable, exponent in monomial_data:
                exponents[VARIABLES.index(variable)] = int(exponent)
            coefficient = Fraction(coefficient_text)
            if coefficient:
                encoded[tuple(exponents)] = (
                    QQ.convert(coefficient.numerator)
                    / QQ.convert(coefficient.denominator)
                )
        polynomial = ring.from_dict(encoded)
        if not polynomial:
            raise ValueError(f"unexpected zero input equation: {equation['tag']}")
        total_terms += len(polynomial)
        monic = polynomial.monic()
        unique.setdefault(tuple(sorted(monic.terms())), (equation["tag"], monic))

    inputs = [item[1] for item in unique.values()]
    zero_monomial = (0,) * len(VARIABLES)
    origin_is_common_zero = all(
        dict(poly.terms()).get(zero_monomial, 0) == 0 for poly in inputs
    )
    # These are upper-band residual equations. If the origin is a common zero,
    # an ordinary GB run on this subsystem cannot exclude the admissible slice;
    # it can only rediscover a degenerate solution of the truncated equations.
    if origin_is_common_zero:
        basis = []
        unit = False
        basis_strings = []
        scope = (
            "The exact origin witness (all eight variables zero) satisfies "
            "every equation in this smallest-support subsystem. Thus this "
            "subsystem is nonunit and cannot exclude the U1=0 slice. A valid "
            "exclusion must attach the omitted nonzero Jacobian/normalization "
            "condition or add further equations."
        )
    else:
        basis = ring_groebner(inputs, ring, method="f5b")
        unit = any(poly == ring.one for poly in basis)
        basis_strings = [str(poly) for poly in basis]
        scope = (
            "A unit ideal proves that this U1=0 slice is empty. A nonunit basis "
            "does not exclude this slice or the full 2+1+1 branch."
        )
    result = {
        "input": str(INPUT),
        "input_sha256": input_hash,
        "branch": document["branch"],
        "specialization": document["specialization"],
        "additional_set_theoretic_reduction": document[
            "additional_set_theoretic_reduction"
        ],
        "variables": list(VARIABLES),
        "equations_before_scalar_deduplication": len(equations),
        "equations_after_scalar_deduplication": len(inputs),
        "terms_after_scalar_deduplication": sum(len(poly) for poly in inputs),
        "terms_before_scalar_deduplication": total_terms,
        "groebner_method": "F5B skipped" if origin_is_common_zero else "F5B",
        "monomial_order": "grevlex",
        "origin_is_common_zero": origin_is_common_zero,
        "origin_witness": {name: "0" for name in VARIABLES}
        if origin_is_common_zero else None,
        "basis_length": len(basis) if not origin_is_common_zero else None,
        "unit_ideal": unit,
        "basis": basis_strings,
        "scope": scope,
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "input_sha256",
                    "equations_before_scalar_deduplication",
                    "equations_after_scalar_deduplication",
                    "terms_before_scalar_deduplication",
                    "terms_after_scalar_deduplication",
                    "origin_is_common_zero",
                    "origin_witness",
                    "basis_length",
                    "unit_ideal",
                    "scope",
                )
            },
            indent=2,
        ),
        flush=True,
    )
    print(f"report={OUTPUT}", flush=True)


if __name__ == "__main__":
    main()
