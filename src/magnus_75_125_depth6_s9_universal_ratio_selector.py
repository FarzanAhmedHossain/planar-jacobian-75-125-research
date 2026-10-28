"""Exact universal s=9 Magnus selector and algebraic branch cross-check.

This recomputes the next [Z,R] slice using only the already-forced rho_0=0
condition; it does not set any coefficients of S_60 or S_59 to zero.  It then
compares the resulting ratio polynomial with the two cubics arising from the
previous nonzero-vector branches.  The resulting eliminations are local
pruning statements only, not an exclusion of the full (75,125) F2 case.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "f2_75_125_flint_deps"))
os.environ.setdefault("SYMPY_GROUND_TYPES", "flint")

import sympy as sp


def bracket(f: sp.Expr, g: sp.Expr, x: sp.Symbol, y: sp.Symbol) -> sp.Expr:
    return sp.diff(f, x) * sp.diff(g, y) - sp.diff(f, y) * sp.diff(g, x)


def ycoeff(poly: sp.Expr, y: sp.Symbol, degree: int) -> sp.Expr:
    return sp.Poly(sp.expand(poly), y).coeff_monomial(y**degree)


def laurent_coeff(poly: sp.Expr, x: sp.Symbol, degree: int, clear: int = 40) -> sp.Expr:
    shifted = sp.Poly(sp.expand(poly * x**clear), x)
    return shifted.nth(degree + clear)


def main() -> None:
    x, y = sp.symbols("x y")
    mu = sp.symbols("mu", nonzero=True)
    eta = sp.symbols("eta9")
    a3, a4 = sp.symbols("a3 a4")
    rho = sp.symbols("rho0:5")
    sigma = sp.symbols("sigma0:5")
    tau = sp.symbols("tau0:5")
    r15c = sp.symbols("r15_0:5")

    s, e = 9, 36
    gamma = sp.Rational(9, 5)

    # rho_0=0 is the only branch substitution used in the [Z,R] calculation.
    r19 = x**3 * (a3 + a4*x)
    r18 = sum(rho[j]*x**j for j in range(5)).subs(rho[0], 0)
    r17 = sum(sigma[j]*x**j for j in range(5))
    r16 = sum(tau[j]*x**j for j in range(5))
    r15 = sum(r15c[j]*x**j for j in range(5))

    z0 = mu*x**s
    z1 = mu*gamma*x**(s-5)*r19
    z2 = mu*(gamma*x**(s-5)*r18
             + sp.binomial(gamma, 2)*x**(s-10)*r19**2)
    z3 = mu*(gamma*x**(s-5)*r17
             + 2*sp.binomial(gamma, 2)*x**(s-10)*r19*r18
             + sp.binomial(gamma, 3)*x**(s-15)*r19**3)
    z4 = sp.expand(mu*(
        gamma*x**(s-5)*r16
        + sp.binomial(gamma, 2)*x**(s-10)*(2*r19*r17+r18**2)
        + 3*sp.binomial(gamma, 3)*x**(s-15)*r19**2*r18
        + sp.binomial(gamma, 4)*x**(s-20)*r19**4
    ) + eta*x**(s-1))

    z5_coeffs = sp.symbols("z5_0:10")
    z5 = sum(z5_coeffs[i]*x**i for i in range(10))
    Zpre = z0*y**36+z1*y**35+z2*y**34+z3*y**33+z4*y**32
    Rtrunc = (x**5*y**20+r19*y**19+r18*y**18+r17*y**17
              +r16*y**16+r15*y**15)
    slice_rhs = ycoeff(bracket(Zpre, Rtrunc, x, y), y, e+14)
    slice_lhs = ycoeff(
        bracket(z5*y**(e-5), x**5*y**20, x, y), y, e+14
    )
    commutator_slice = sp.expand(slice_rhs+slice_lhs)
    comm_poly = sp.Poly(sp.expand(commutator_slice*x**40), x)

    z5_solve: dict[sp.Symbol, sp.Expr] = {}
    for degree in range(4, 18):
        coeff_expr = sp.factor_terms(comm_poly.nth(degree+40).subs(z5_solve))
        matching = [v for v in z5_coeffs if sp.diff(coeff_expr, v) != 0]
        if matching:
            v = matching[0]
            pivot = sp.factor_terms(sp.diff(coeff_expr, v))
            if pivot.is_number and pivot != 0:
                z5_solve[v] = sp.factor_terms(-coeff_expr.subs(v, 0)/pivot)

    residual = sp.expand(commutator_slice.subs(z5_solve))
    residual_coeffs = {
        degree: sp.factor(laurent_coeff(residual, x, degree))
        for degree in range(-40, 10)
    }
    nonzero_residuals = {k: v for k, v in residual_coeffs.items() if v != 0}
    expected = -sp.Rational(63, 625)*mu*a3*(
        11*a3**4-100*a3**2*rho[1]+250*rho[1]**2
    )
    assert nonzero_residuals == {3: expected}, nonzero_residuals

    t = sp.symbols("t")
    ratio = 250*t**2-100*t+11
    cubic_u2c2 = 70750*t**3-26100*t**2+2925*t-172
    cubic_u3c3 = 69250*t**3-24300*t**2+2115*t-52
    results = {
        "ratio_discriminant": int(sp.discriminant(ratio, t)),
        "u2c2_gcd": str(sp.gcd(ratio, cubic_u2c2)),
        "u2c2_resultant": int(sp.resultant(ratio, cubic_u2c2, t)),
        "u3c3_gcd": str(sp.gcd(ratio, cubic_u3c3)),
        "u3c3_resultant": int(sp.resultant(ratio, cubic_u3c3, t)),
    }
    assert results == {
        "ratio_discriminant": -1000,
        "u2c2_gcd": "1",
        "u2c2_resultant": 295619000000,
        "u3c3_gcd": "1",
        "u3c3_resultant": 221699000000,
    }, results

    prior = HERE / "magnus_75_125_depth5_s9_selector_audit.py"
    receipt = {
        "claim": "The universal s=9 [Z,R] slice is incompatible with either previous nonzero-vector ratio cubic when a3 != 0.",
        "assumptions": ["characteristic zero", "s=9", "rho0=0", "the two stated nonzero-vector cubic branches"],
        "universal_slice_residual": str(expected),
        "ratio_polynomial": str(ratio),
        "previous_u2c2_cubic": str(cubic_u2c2),
        "previous_u3c3_cubic": str(cubic_u3c3),
        **results,
        "branch_consequence": "On either nonzero-vector cubic branch, a3!=0 is impossible; a3=0 then forces rho1=0. The zero-vector leaf is not excluded.",
        "scope": "local necessary-condition pruning only; not a whole-branch exclusion",
        "source_depth5_script_sha256": hashlib.sha256(prior.read_bytes()).hexdigest(),
    }
    out = HERE / "magnus_75_125_depth6_s9_universal_ratio_selector_20261019.json"
    out.write_text(json.dumps(receipt, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
