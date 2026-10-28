"""Audit the next exact Magnus remainder row on the s=9 a3=rho1=0 leaf.

This continues the exact quotient recurrence from the audited y^(e+54) row
to y^(e+53), after the universal selector and previous ratio cubics have
forced a3=rho1=0 on either nonzero-vector leaf. It retains u2,c2,u3,c3 as
free coordinates. New Laurent rows and commutator obstructions are kept as
necessary equations; no branch parameter is divided out.
"""

from __future__ import annotations

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


def laurent_coeff(poly: sp.Expr, x: sp.Symbol, degree: int, clear: int = 40) -> sp.Expr:
    shifted = sp.Poly(sp.expand(poly * x**clear), x)
    return sp.factor_terms(shifted.nth(degree + clear))


def laurent_coeffs(
    poly: sp.Expr, x: sp.Symbol, degrees: range, clear: int = 40
) -> dict[int, sp.Expr]:
    """Expand once, then read a whole consecutive Laurent coefficient block."""
    shifted = sp.Poly(sp.expand(poly * x**clear), x)
    return {degree: shifted.nth(degree + clear) for degree in degrees}


def main() -> None:
    x, y = sp.symbols("x y")
    mu = sp.symbols("mu", nonzero=True)
    eta = sp.symbols("eta9")
    a3, a4 = sp.symbols("a3 a4")
    rho = sp.symbols("rho0:5")
    sigma = sp.symbols("sigma0:5")
    tau = sp.symbols("tau0:5")
    r15c = sp.symbols("r15_0:5")
    u = sp.symbols("u0:16")
    c = sp.symbols("c0:17")
    d = sp.symbols("d0:18")
    f = sp.symbols("f0:18")
    g = sp.symbols("g0:18")
    h = sp.symbols("h0:18")

    s, e = 9, 36
    gamma = sp.Rational(9, 5)
    branch = {
        a3: 0, rho[0]: 0, rho[1]: 0,
        u[0]: 0, u[1]: 0,
        c[0]: 0, c[1]: 0,
    }

    r19 = x**3 * (a3 + a4*x)
    r18 = sum(rho[j]*x**j for j in range(5))
    r17 = sum(sigma[j]*x**j for j in range(5))
    r16 = sum(tau[j]*x**j for j in range(5))
    r15 = sum(r15c[j]*x**j for j in range(5))
    r14c = sp.symbols("r14_0:4")
    r14 = sum(r14c[j]*x**j for j in range(4))

    S60 = x**2 * sum(u[j]*x**j for j in range(16))
    S59 = sum(c[j]*x**j for j in range(17))
    S58 = sum(d[j]*x**j for j in range(18))
    S57 = sum(f[j]*x**j for j in range(18))
    S56 = sum(g[j]*x**j for j in range(18))
    S55 = sum(h[j]*x**j for j in range(18))
    s54c = sp.symbols("s54_0:8")
    S54 = sum(s54c[j]*x**j for j in range(8))

    # Apply the complete zero-vector specialization before continuing.
    S60 = sp.expand(S60.subs(branch))
    S59 = sp.expand(S59.subs(branch))
    r19, r18, r17, r16 = (
        sp.expand(q.subs(branch)) for q in (r19, r18, r17, r16)
    )

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

    # Previous rows solve the x^8,x^9 coefficients of S58 and S57.
    q1 = 2*x**5*r19
    q2 = r19**2 + 2*x**5*r18
    q3 = 2*x**5*r17 + 2*r19*r18
    q4 = 2*x**5*r16 + 2*r19*r17 + r18**2
    q5 = 2*x**5*r15 + 2*r19*r16 + 2*r18*r17

    C0 = 60*sp.diff(z0, x)*S60 - e*z0*sp.diff(S60, x)
    C1 = (59*sp.diff(z0, x)*S59 - e*z0*sp.diff(S59, x)
          + 60*sp.diff(z1, x)*S60 - (e-1)*z1*sp.diff(S60, x))
    C2 = (58*sp.diff(z0, x)*S58 - e*z0*sp.diff(S58, x)
          + 59*sp.diff(z1, x)*S59 - (e-1)*z1*sp.diff(S59, x)
          + 60*sp.diff(z2, x)*S60 - (e-2)*z2*sp.diff(S60, x))

    T0 = sp.cancel(-C0/x**10)
    T1 = sp.cancel(-(C1+q1*T0)/x**10)
    N2 = sp.expand(C2+q1*T1+q2*T0)
    solved_d = {}
    for degree in (8, 9):
        j = degree-8
        equation = sp.factor_terms(laurent_coeff(N2.subs(solved_d), x, degree))
        pivot = sp.factor_terms(sp.diff(equation, d[j]))
        assert sp.simplify(pivot-mu*s*(58-4*j)) == 0
        solved_d[d[j]] = sp.factor_terms(-equation.subs(d[j], 0)/pivot)
    S58 = sp.expand(S58.subs(branch).subs(solved_d))
    N2 = sp.expand(N2.subs(solved_d))
    N2_low = laurent_coeffs(N2, x, range(8))
    for k, value in N2_low.items():
        assert value == 0, ("prior K4 row", k)
    T2 = sp.cancel(-N2/x**10)

    C3fixed = (58*sp.diff(z1, x)*S58-(e-1)*z1*sp.diff(S58, x)
               + 59*sp.diff(z2, x)*S59-(e-2)*z2*sp.diff(S59, x)
               + 60*sp.diff(z3, x)*S60-(e-3)*z3*sp.diff(S60, x))
    N3fixed = sp.expand(C3fixed+q1*T2+q2*T1+q3*T0)
    solved_f = {}
    for degree in (8, 9):
        j = degree-8
        equation = sp.factor_terms(laurent_coeff(
            (N3fixed+57*sp.diff(z0, x)*S57-e*z0*sp.diff(S57, x)).subs(solved_f),
            x, degree,
        ))
        pivot = sp.factor_terms(sp.diff(equation, f[j]))
        assert sp.simplify(pivot-mu*s*(57-4*j)) == 0
        solved_f[f[j]] = sp.factor_terms(-equation.subs(f[j], 0)/pivot)
    S57 = sp.expand(S57.subs(solved_f))
    N3 = sp.expand(N3fixed+57*sp.diff(z0, x)*S57-e*z0*sp.diff(S57, x))
    N3_low = laurent_coeffs(N3, x, range(8))
    prior_k5 = [(k, value) for k, value in N3_low.items() if value != 0]
    prior_k5_residual = sp.expand(sum(
        value*x**degree for degree, value in prior_k5
    ))
    # These low coefficients are equations on the current nonzero-vector
    # branch. Retain them explicitly, then quotient only the part known to
    # vanish on that branch's solution set.
    N3_for_quotient = sp.expand(N3-prior_k5_residual)
    T3 = sp.cancel(-N3_for_quotient/x**10)

    C4fixed = (57*sp.diff(z1, x)*S57-(e-1)*z1*sp.diff(S57, x)
               + 58*sp.diff(z2, x)*S58-(e-2)*z2*sp.diff(S58, x)
               + 59*sp.diff(z3, x)*S59-(e-3)*z3*sp.diff(S59, x)
               + 60*sp.diff(z4, x)*S60-(e-4)*z4*sp.diff(S60, x))
    N4fixed = sp.expand(C4fixed+q1*T3+q2*T2+q3*T1+q4*T0)
    previous_fixed = laurent_coeffs(N4fixed, x, range(8))
    # The e+55 row also contains [z0*y^e,S56*y^56]. Its x^8,x^9
    # equations solve g0,g1 with the constant pivots mu*s*(56-4*j).
    solved_g = {}
    for degree in (8, 9):
        j = degree-8
        equation = sp.factor_terms(laurent_coeff(
            (N4fixed+56*sp.diff(z0, x)*S56-e*z0*sp.diff(S56, x)).subs(solved_g),
            x, degree,
        ))
        pivot = sp.factor_terms(sp.diff(equation, g[j]))
        assert sp.simplify(pivot-mu*s*(56-4*j)) == 0
        solved_g[g[j]] = sp.factor_terms(-equation.subs(g[j], 0)/pivot)
    S56 = sp.expand(S56.subs(solved_g))
    N4 = sp.expand(N4fixed+56*sp.diff(z0, x)*S56-e*z0*sp.diff(S56, x))
    # The x^6,x^7 coefficients are residual constraints from the preceding
    # (e+55) row. A genuine branch point satisfies them, so remove those
    # constrained low terms before taking the quotient that enters the next
    # remainder row. Do not silently assert that they vanish identically.
    previous_residual = sp.expand(sum(
        value*x**degree for degree, value in previous_fixed.items()
    ))
    N4_for_quotient = sp.expand(N4-previous_residual)
    T4 = sp.cancel(-N4_for_quotient/x**10)

    # Solve [Z,R]=0 at the next slice for z5, retaining eta. The coefficient
    # of z5 has only rational integer pivots (20*i-31), never a3 or another
    # branch parameter, so this step does not discard zero cases.
    z5_coeffs = sp.symbols("z5_0:10")
    z5 = sum(z5_coeffs[i]*x**i for i in range(10))
    Zpre = z0*y**36+z1*y**35+z2*y**34+z3*y**33+z4*y**32
    Rtrunc = (x**5*y**20+r19*y**19+r18*y**18+r17*y**17
              +r16*y**16+r15*y**15)
    z5_rhs = ycoeff(bracket(Zpre, Rtrunc, x, y), y, e+14)
    commutator_slice = sp.expand(z5_rhs+ycoeff(
        bracket(z5*y**(e-5), x**5*y**20, x, y), y, e+14
    ))
    comm_poly = sp.Poly(sp.expand(commutator_slice*x**40), x)
    z5_solve = {}
    for degree in range(4, 18):
        coeff_expr = sp.factor_terms(comm_poly.nth(degree+40).subs(z5_solve))
        matching = [v for v in z5_coeffs if sp.diff(coeff_expr, v) != 0]
        if matching:
            v = matching[0]
            pivot = sp.factor_terms(sp.diff(coeff_expr, v))
            if pivot.is_number and pivot != 0:
                z5_solve[v] = sp.factor_terms(-coeff_expr.subs(v, 0)/pivot)
    z5_solved = sp.expand(z5.subs(z5_solve))
    residual_z5 = sp.expand(commutator_slice.subs(z5_solve))
    z5_degrees = range(-40, 10)
    z5_residual_coeffs = laurent_coeffs(residual_z5, x, z5_degrees)
    z5_solved_coeffs = laurent_coeffs(z5_solved, x, z5_degrees)
    z5_obstructions = [(k, value) for k, value in z5_residual_coeffs.items()
                       if value != 0]
    z5_support = [k for k, value in z5_solved_coeffs.items() if value != 0]
    z5_nonpolynomial = [(k, z5_solved_coeffs[k])
                        for k in z5_support if k < 0 or k > 9]

    # The sixth remainder numerator is the y^(e+54) coefficient. S55 and
    # S56 are kept general; no pivot/nonzero assumptions are made here.
    Ztrunc = Zpre+z5_solved*y**31
    Strunc = (S60*y**60+S59*y**59+S58*y**58+S57*y**57
              +S56*y**56+S55*y**55)
    C5 = ycoeff(bracket(Ztrunc, Strunc, x, y), y, e+54)
    directC5 = (55*sp.diff(z0, x)*S55-e*z0*sp.diff(S55, x)
                +56*sp.diff(z1, x)*S56-(e-1)*z1*sp.diff(S56, x)
                +57*sp.diff(z2, x)*S57-(e-2)*z2*sp.diff(S57, x)
                +58*sp.diff(z3, x)*S58-(e-3)*z3*sp.diff(S58, x)
                +59*sp.diff(z4, x)*S59-(e-4)*z4*sp.diff(S59, x)
                +60*sp.diff(z5_solved, x)*S60-(e-5)*z5_solved*sp.diff(S60, x))
    assert sp.expand(C5-directC5) == 0
    direct_q5 = ycoeff(Rtrunc**2, y, 35)
    assert sp.expand(direct_q5-q5) == 0

    N5 = sp.expand(C5+q1*T4+q2*T3+q3*T2+q4*T1+q5*T0)
    # The x^8,x^9 coefficients determine h0,h1 using the fixed nonzero
    # pivots mu*s*(55-4j); retaining this exact graph makes T5 available.
    solved_h = {}
    for degree in (8, 9):
        j = degree-8
        equation = sp.factor_terms(laurent_coeff(N5.subs(solved_h), x, degree))
        pivot = sp.factor_terms(sp.diff(equation, h[j]))
        assert sp.simplify(pivot-mu*s*(55-4*j)) == 0
        solved_h[h[j]] = sp.factor_terms(-equation.subs(h[j], 0)/pivot)
    S55 = sp.expand(S55.subs(solved_h))
    N5 = sp.expand(N5.subs(solved_h))
    N5_low = laurent_coeffs(N5, x, range(8))
    fixed = list(N5_low.items())
    fixed = [(k, value) for k, value in fixed if value != 0]

    # Retain the depth-six fixed Laurent rows as equations, then take the
    # quotient on their zero set to define T5.  This mirrors the earlier
    # N3/N4 quotient handling and does not assert those rows vanish identically.
    depth6_residual = sp.expand(sum(value*x**degree for degree, value in fixed))
    T5 = sp.cancel(-(N5-depth6_residual)/x**10)

    # Next solve [Z,R]=0 at y^(e+13)=y^49 for z6.  Include the full required
    # R14 slice; its four coefficients remain free and no pivot depends on
    # branch parameters.
    z6_coeffs = sp.symbols("z6_0:10")
    z6 = sum(z6_coeffs[i]*x**i for i in range(10))
    Zpre6 = Zpre+z5_solved*y**31
    Rtrunc6 = (Rtrunc+r14*y**14)
    z6_rhs = ycoeff(bracket(Zpre6, Rtrunc6, x, y), y, e+13)
    z6_equation = sp.expand(z6_rhs+ycoeff(
        bracket(z6*y**(e-6), x**5*y**20, x, y), y, e+13
    ))
    z6_poly = sp.Poly(sp.expand(z6_equation*x**40), x)
    solved_z6 = {}
    for degree in range(4, 14):
        equation = sp.factor_terms(z6_poly.nth(degree+40).subs(solved_z6))
        matching = [v for v in z6_coeffs if sp.diff(equation, v) != 0]
        if matching:
            v = matching[0]
            pivot = sp.factor_terms(sp.diff(equation, v))
            assert sp.simplify(pivot-(20*(degree-4)-150)) == 0
            if pivot.is_number and pivot != 0:
                solved_z6[v] = sp.factor_terms(-equation.subs(v, 0)/pivot)
    z6_solved = sp.expand(z6.subs(solved_z6))
    z6_residual = sp.expand(z6_equation.subs(solved_z6))
    z6_residual_poly = sp.Poly(sp.expand(z6_residual*x**40), x)
    z6_obstructions = [(degree-40, value)
                       for (degree,), value in z6_residual_poly.terms()
                       if value != 0]

    # Degree-seven remainder numerator: C6 plus q1*T5+...+q6*T0.  The first
    # eight x-coefficients are independent of S54 because the z0 commutator
    # starts at x^8.  We keep the S54 symbols in the source expression and
    # assert this independence rather than silently truncating a term.
    q6 = ycoeff(Rtrunc6**2, y, 34)
    zlist = [z0, z1, z2, z3, z4, z5_solved, z6_solved]
    slist = [S60, S59, S58, S57, S56, S55, S54]
    C6 = sp.expand(sum(
        (60-j)*sp.diff(zlist[i], x)*slist[j]
        -(e-i)*zlist[i]*sp.diff(slist[j], x)
        for i in range(7) for j in (6-i,)
    ))
    q_values = [x**10, q1, q2, q3, q4, q5, q6]
    T_values = [T0, T1, T2, T3, T4, T5]
    N6 = sp.expand(C6+sum(q_values[i]*T_values[6-i]
                           for i in range(1, 7)))
    N6_low = laurent_coeffs(N6, x, range(8))
    if any(value.has(*s54c) for value in N6_low.values()):
        raise AssertionError("S54 unexpectedly affects x^0..x^7")
    depth7_fixed = [(k, value) for k, value in N6_low.items() if value != 0]

    print("branch=s9,a3=rho1=u0=u1=c0=c1=0;u2,c2,u3,c3 free")
    prior = [(k, v) for k, v in previous_fixed.items() if v != 0]
    print(f"depth5_fixed_row_count={len(prior)}")
    for degree, equation in prior:
        print(f"depth5_fixed_x{degree}={sp.factor_terms(equation)}")
    print(f"prior_k5_row_count={len(prior_k5)}")
    for degree, equation in prior_k5:
        print(f"prior_k5_fixed_x{degree}={sp.factor_terms(equation)}")
    print(f"z5_support={z5_support}")
    print(f"z5_nonpolynomial_conditions={z5_nonpolynomial}")
    print(f"z5_remaining_commutator_rows={z5_obstructions}")
    print(f"depth6_fixed_row_count={len(fixed)}")
    for degree, equation in fixed:
        print(f"fixed_x{degree}={sp.factor_terms(equation)}")
    print(f"solved_S55_low_coefficients={{{', '.join(f'{key}:{sp.factor_terms(value)}' for key, value in solved_h.items())}}}")
    print(f"z6_solved_coefficient_count={len(solved_z6)}")
    print(f"z6_remaining_commutator_rows={[(k, sp.factor_terms(v)) for k, v in z6_obstructions]}")
    print(f"depth7_fixed_row_count={len(depth7_fixed)}")
    for degree, equation in depth7_fixed:
        print(f"depth7_fixed_x{degree}={sp.factor_terms(equation)}")
    print("scope=exact next necessary rows only; not a branch exclusion")


if __name__ == "__main__":
    main()
