"""Exact z_6 edge-pin audit for the s=9 Magnus vertex (9,36).

This replays only the [Z,R] slice recurrence through z_6, then imposes the
F2 boundary condition [x^9]Z_30=0. It tests the two already-known quartic
edge values, without claiming a full F2 branch exclusion.
"""

from __future__ import annotations

import json
import os
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


def solve_slice(
    equation: sp.Expr,
    x: sp.Symbol,
    coefficients: tuple[sp.Symbol, ...],
    degrees: range,
) -> tuple[sp.Expr, dict[sp.Symbol, sp.Expr]]:
    poly = sp.Poly(sp.expand(equation * x**40), x)
    solved: dict[sp.Symbol, sp.Expr] = {}
    for degree in degrees:
        row = sp.factor_terms(poly.nth(degree + 40).subs(solved))
        matching = [v for v in coefficients if sp.diff(row, v) != 0]
        if matching:
            variable = matching[0]
            pivot = sp.factor_terms(sp.diff(row, variable))
            if not (pivot.is_number and pivot != 0):
                raise AssertionError(f"unexpected slice pivot {pivot}")
            solved[variable] = sp.factor_terms(-row.subs(variable, 0) / pivot)
    return sp.expand(equation.subs(solved)), solved


def main() -> None:
    x, y = sp.symbols("x y")
    mu = sp.Symbol("mu", nonzero=True)
    alpha, xi = sp.symbols("alpha xi")
    a4 = sp.Symbol("a4")
    rho2, rho3, rho4 = sp.symbols("rho2 rho3 rho4")
    sigma1, sigma2, sigma3, sigma4 = sp.symbols("sigma1 sigma2 sigma3 sigma4")
    tau = sp.symbols("tau0:5")
    r15lo = sp.symbols("r15_0:4")
    r14c = sp.symbols("r14_0:4")

    s, e = 9, 36
    gamma = sp.Rational(9, 5)
    r19 = a4 * x**4
    r18 = rho2*x**2 + rho3*x**3 + rho4*x**4
    r17 = sigma1*x + sigma2*x**2 + sigma3*x**3 + sigma4*x**4
    r16 = sum(tau[j]*x**j for j in range(5))
    r15 = sum(r15lo[j]*x**j for j in range(4)) + xi*x**4
    r14 = sum(r14c[j]*x**j for j in range(4))

    z0 = mu*x**s
    z1 = mu*gamma*x**(s-5)*r19
    z2 = mu*(gamma*x**(s-5)*r18
             + sp.binomial(gamma, 2)*x**(s-10)*r19**2)
    z3 = mu*(gamma*x**(s-5)*r17
             + 2*sp.binomial(gamma, 2)*x**(s-10)*r19*r18
             + sp.binomial(gamma, 3)*x**(s-15)*r19**3)
    z3 = sp.expand(z3)
    eta = sp.Symbol("eta9")
    z4 = sp.expand(mu*(
        gamma*x**(s-5)*r16
        + sp.binomial(gamma, 2)*x**(s-10)*(2*r19*r17+r18**2)
        + 3*sp.binomial(gamma, 3)*x**(s-15)*r19**2*r18
        + sp.binomial(gamma, 4)*x**(s-20)*r19**4
    ) + eta*x**(s-1))

    z5_coeffs = sp.symbols("z5_0:10")
    z5 = sum(z5_coeffs[i]*x**i for i in range(10))
    Zpre = z0*y**36 + z1*y**35 + z2*y**34 + z3*y**33 + z4*y**32
    Rtrunc = (x**5*y**20 + r19*y**19 + r18*y**18 + r17*y**17
              + r16*y**16 + r15*y**15)
    z5_rhs = ycoeff(bracket(Zpre, Rtrunc, x, y), y, e+14)
    z5_equation = sp.expand(z5_rhs + ycoeff(
        bracket(z5*y**(e-5), x**5*y**20, x, y), y, e+14
    ))
    z5_residual, solved_z5 = solve_slice(
        z5_equation, x, z5_coeffs, range(4, 18)
    )
    z5_solved = sp.expand(z5.subs(solved_z5))
    z5_residual_poly = sp.Poly(sp.expand(z5_residual*x**40), x)
    z5_residual_rows = [
        (k, sp.factor_terms(z5_residual_poly.nth(k+40)))
        for k in range(-40, 10)
        if z5_residual_poly.nth(k+40) != 0
    ]

    z6_coeffs = sp.symbols("z6_0:10")
    z6 = sum(z6_coeffs[i]*x**i for i in range(10))
    Zpre6 = Zpre + z5_solved*y**31
    Rtrunc6 = Rtrunc + r14*y**14
    z6_rhs = ycoeff(bracket(Zpre6, Rtrunc6, x, y), y, e+13)
    z6_equation = sp.expand(z6_rhs + ycoeff(
        bracket(z6*y**(e-6), x**5*y**20, x, y), y, e+13
    ))
    z6_residual, solved_z6 = solve_slice(
        z6_equation, x, z6_coeffs, range(4, 14)
    )
    z6_solved = sp.expand(z6.subs(solved_z6))
    edge_row = sp.factor_terms(sp.Poly(z6_solved, x).coeff_monomial(x**9))
    z6_residual_poly = sp.Poly(sp.expand(z6_residual*x**40), x)
    z6_residual_rows = [
        (k, sp.factor_terms(z6_residual_poly.nth(k+40)))
        for k in range(-40, 14)
        if z6_residual_poly.nth(k+40) != 0
    ]

    results = []
    for name, edge_value in (
        ("2+1+1", sp.Integer(-5)),
        ("2+2", -2*(alpha+1)),
    ):
        expr = sp.factor_terms(edge_row.subs(xi, edge_value))
        if name == "2+2":
            expr = sp.rem(sp.together(expr).as_numer_denom()[0],
                          alpha**2-3*alpha+1, alpha)
        results.append({
            "branch": name,
            "r15_4": str(edge_value),
            "edge_condition_x9_z6": str(expr),
            "is_identically_zero": bool(sp.simplify(expr) == 0),
            "operation_count": int(sp.count_ops(expr)),
        })

    output = {
        "vertex": [9, 36],
        "zero_leaf_specialization": ["a3=0", "rho0=rho1=0", "sigma0=0"],
        "solved_z5_coefficients": len(solved_z5),
        "z5_commutator_residual_rows": [
            {"x_degree": k, "equation": str(v)} for k, v in z5_residual_rows
        ],
        "solved_z6_coefficients": len(solved_z6),
        "z6_commutator_residual_rows": [
            {"x_degree": k, "equation": str(v)} for k, v in z6_residual_rows
        ],
        "z6_edge_row": results,
        "scope": "necessary edge-row conditions only; no F2 branch exclusion",
    }
    target = HERE / "magnus_75_125_s9_z6_edge_pin_audit_20261021.json"
    target.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
