"""Exact symbolic replay of the first-two-band divisor in F2 2+1+1.

For h=(t-1)^2 R(t), with R monic, squarefree of degree two, and R(1) != 0,
the band equations give h | U1^2.  The rootwise terminal condition gives
ord_{t=1}(U1) >= 5.  Together these imply (t-1)^5 R | U1.  Since deg U1 <= 11,
write U1=(t-1)^5 R A with deg A <= 4.  This leaves five U1 coefficients for
each fixed quartic h, including the zero case A=0.

The script checks the parameterized identities exactly in Q[t] with symbolic
coefficients.  It does not exclude the 2+1+1 branch.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "F2-75-125-band-model.md"
OUTPUT = HERE / "f2_75_125_211_u1_divisor_reduction_20261021.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite {OUTPUT}")

    t, c1, c0 = sp.symbols("t c1 c0")
    a = sp.symbols("a0:5")
    R = t**2 + c1 * t + c0
    A = sum(a[i] * t**i for i in range(5))
    h = sp.expand((t - 1) ** 2 * R)
    divisor = sp.expand((t - 1) ** 5 * R)
    U1 = sp.expand(divisor * A)

    quotient, remainder = sp.div(sp.Poly(U1**2, t), sp.Poly(h, t))
    if remainder.as_expr() != 0:
        raise AssertionError("h does not divide U1^2 under the parametrization")
    if sp.Poly(U1, t).degree() > 11:
        raise AssertionError("the parameterized U1 exceeds its support degree")
    if sp.Poly(A, t).degree() > 4:
        raise AssertionError("the quotient exceeds degree four")
    if any(sp.diff(U1, t, k).subs(t, 1) != 0 for k in range(5)):
        raise AssertionError("the terminal order at t=1 is below five")

    exact_quotient = sp.cancel(U1**2 / h)
    expected_quotient = sp.expand((t - 1) ** 8 * R * A**2)
    if sp.expand(exact_quotient - expected_quotient) != 0:
        raise AssertionError("the exact h-divisibility quotient changed")

    receipt = {
        "source": str(SOURCE),
        "source_sha256": sha256(SOURCE),
        "branch": "2+1+1",
        "field": "Q(c0,c1,a0,a1,a2,a3,a4)",
        "quartic_factor": "h=(t-1)^2 R(t)",
        "R_assumptions": ["monic quadratic", "squarefree", "R(1) != 0"],
        "band_constraint": "h divides U1^2",
        "terminal_constraint": "ord_(t=1)(U1) >= 5",
        "deduced_divisor": "(t-1)^5 R(t) divides U1",
        "support_bound": "degree(U1) <= 11",
        "parameterization": "U1=(t-1)^5 R(t) A(t), degree(A) <= 4",
        "free_U1_coefficients_for_fixed_R": 5,
        "zero_U1_case_included": True,
        "checked_identities": {
            "h_divides_U1_squared": True,
            "U1_squared_over_h": "(t-1)^8 R(t) A(t)^2",
            "terminal_order_at_t1_at_least_5": True,
            "degree_U1_at_most_11": True,
        },
        "scope_warning": (
            "This is a necessary first-band reduction only. It does not exclude "
            "the 2+1+1 branch or the degree pair."
        ),
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "EXACT_DIVISOR_REPLAY_PASSED",
        "branch": receipt["branch"],
        "parameterization": receipt["parameterization"],
        "free_U1_coefficients_for_fixed_R": receipt["free_U1_coefficients_for_fixed_R"],
        "output": str(OUTPUT),
        "output_sha256": sha256(OUTPUT),
    }, indent=2))


if __name__ == "__main__":
    main()
