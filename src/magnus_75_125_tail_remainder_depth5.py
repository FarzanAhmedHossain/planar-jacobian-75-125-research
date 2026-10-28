"""Continue the exact R^2-remainder expansion through y^(e+55).

The preceding rows determine S_58 and the low coefficients of S_57.  This
checker then forms the next quotient numerator, retaining the resonant
fourth Z-slice.  For e=24 that slice is free at the [Z,R] cutoff; for
e=28,32,36 it is the forced fractional-power slice plus its x^(s-1)
kernel.  The reported fixed rows are necessary conditions only.
"""

from __future__ import annotations

import argparse
import sympy as sp


def bracket(f, g, x, y):
    return sp.diff(f, x) * sp.diff(g, y) - sp.diff(f, y) * sp.diff(g, x)


def coeff(poly, x, degree):
    return sp.Poly(sp.expand(poly), x).nth(degree)


def ycoeff(poly, y, degree):
    return sp.Poly(sp.expand(poly), y).coeff_monomial(y**degree)


def laurent_coeff(poly, x, degree):
    """Read an x^degree coefficient after clearing the known x^-10 tail."""
    shifted = sp.Poly(sp.expand(poly*x**10), x)
    return shifted.nth(degree+10)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s", type=int, choices=(6, 7, 8, 9))
    parser.add_argument("--leaf")
    parser.add_argument("--compact", action="store_true")
    args = parser.parse_args()

    x, y = sp.symbols("x y")
    mu = sp.symbols("mu", nonzero=True)
    a3, a4 = sp.symbols("a3 a4")
    rho = sp.symbols("rho0:5")
    sigma = sp.symbols("sigma0:5")
    tau = sp.symbols("tau0:5")
    u = sp.symbols("u0:16")
    c = sp.symbols("c0:17")
    d = sp.symbols("d0:18")
    f = sp.symbols("f0:18")

    r19 = x**3 * (a3 + a4*x)
    r18 = sum(rho[j]*x**j for j in range(5))
    r17 = sum(sigma[j]*x**j for j in range(5))
    r16 = sum(tau[j]*x**j for j in range(5))
    S59_full = sum(c[j]*x**j for j in range(17))
    S58_full = sum(d[j]*x**j for j in range(18))
    S57_full = sum(f[j]*x**j for j in range(18))

    c_pins = {
        6: {c[0]: 0, c[1]: 0, c[2]: 0,
            c[3]: sp.Rational(55, 47)*a3*u[0],
            c[4]: sp.Rational(6, 5)*a3*u[1]+a4*u[0]},
        7: {c[0]: 0, c[1]: 0,
            c[2]: sp.Rational(248, 255)*a3*u[0],
            c[3]: a3*u[1]+sp.Rational(4, 5)*a4*u[0]},
        8: {c[0]: 0,
            c[1]: sp.Rational(213, 275)*a3*u[0],
            c[2]: sp.Rational(4, 5)*a3*u[1]+sp.Rational(3, 5)*a4*u[0]},
        9: {c[0]: sp.Rational(34, 59)*a3*u[0],
            c[1]: sp.Rational(3, 5)*a3*u[1]+sp.Rational(2, 5)*a4*u[0]},
    }
    leaves = {
        6: [("a3=0", {a3: 0}), ("rho0=0", {rho[0]: 0})],
        7: [("rho1=a3^2/5", {rho[0]: 0, rho[1]: a3**2/5}),
            ("rho1=2a3^2/5", {rho[0]: 0, rho[1]: 2*a3**2/5})],
        8: [("u0=0", {rho[0]: 0, u[0]: 0}),
            ("u0!=0", {rho[0]: 0, rho[1]: sp.Rational(2529, 15125)*a3**2})],
        9: [("u0!=0", {rho[0]: 0,
                        rho[1]: sp.Rational(59, 240)*a3**2,
                        rho[2]: sp.Rational(57, 95)*a3*a4
                        - sp.Rational(133, 9120)*a3**2*u[1]/u[0]}),
            ("u0=0,u1=0", {rho[0]: 0, u[0]: 0, u[1]: 0}),
            ("u0=0,u1!=0", {rho[0]: 0, u[0]: 0,
                              rho[1]: sp.Rational(11, 47)*a3**2})],
    }

    for s in (6, 7, 8, 9):
        if args.s is not None and s != args.s:
            continue
        e = 4*s
        gamma = sp.Rational(s, 5)
        for leaf_name, branch in leaves[s]:
            if args.leaf is not None and leaf_name != args.leaf:
                continue
            S60 = x**(11-s)*sum(u[j]*x**j for j in range(s+7))
            S59 = sp.expand(S59_full.subs(c_pins[s]).subs(branch))
            S58 = sp.expand(S58_full.subs(branch))
            S57 = sp.expand(S57_full.subs(branch))
            r19b = sp.expand(r19.subs(branch))
            r18b = sp.expand(r18.subs(branch))
            r17b = sp.expand(r17.subs(branch))
            r16b = sp.expand(r16.subs(branch))
            S60b = sp.expand(S60.subs(branch))
            S59b = sp.expand(S59)

            z0 = mu*x**s
            z1 = mu*gamma*x**(s-5)*r19b
            z2 = mu*(gamma*x**(s-5)*r18b
                     + sp.binomial(gamma, 2)*x**(s-10)*r19b**2)
            z3 = mu*(gamma*x**(s-5)*r17b
                     + 2*sp.binomial(gamma, 2)*x**(s-10)*r19b*r18b
                     + sp.binomial(gamma, 3)*x**(s-15)*r19b**3)

            C0 = 60*sp.diff(z0, x)*S60b - e*z0*sp.diff(S60b, x)
            C1 = (59*sp.diff(z0, x)*S59b - e*z0*sp.diff(S59b, x)
                  + 60*sp.diff(z1, x)*S60b - (e-1)*z1*sp.diff(S60b, x))
            C2 = (58*sp.diff(z0, x)*S58 - e*z0*sp.diff(S58, x)
                  + 59*sp.diff(z1, x)*S59b - (e-1)*z1*sp.diff(S59b, x)
                  + 60*sp.diff(z2, x)*S60b - (e-2)*z2*sp.diff(S60b, x))

            q1 = 2*x**5*r19b
            q2 = r19b**2 + 2*x**5*r18b
            q3 = 2*x**5*r17b + 2*r19b*r18b
            q4 = 2*x**5*r16b + 2*r19b*r17b + r18b**2
            T0 = sp.cancel(-C0/x**10)
            T1 = sp.cancel(-(C1 + q1*T0)/x**10)
            N2 = sp.expand(C2 + q1*T1 + q2*T0)

            # Solve the S_58 coefficients used by the preceding quotient row.
            solved_d = {}
            for degree in range(s-1, 10):
                j = degree-(s-1)
                eq = sp.factor(coeff(N2.subs(solved_d), x, degree))
                pivot = sp.factor(sp.diff(eq, d[j]))
                assert sp.simplify(pivot - mu*s*(58-4*j)) == 0
                solved_d[d[j]] = sp.factor(-eq.subs(d[j], 0)/pivot)
            S58b = sp.expand(S58.subs(solved_d))
            N2 = sp.expand(N2.subs(solved_d))
            T2 = sp.cancel(-N2/x**10)
            assert sp.Poly(T2, x).is_univariate
            assert all(coeff(N2, x, k) == 0 for k in range(10))

            # The next row determines the low S_57 coefficients, with the
            # already-audited fixed rows below x^(s-1) retained as checks.
            C3_fixed = (58*sp.diff(z1, x)*S58b-(e-1)*z1*sp.diff(S58b, x)
                        + 59*sp.diff(z2, x)*S59b-(e-2)*z2*sp.diff(S59b, x)
                        + 60*sp.diff(z3, x)*S60b-(e-3)*z3*sp.diff(S60b, x))
            N3_fixed = sp.expand(C3_fixed + q1*T2 + q2*T1 + q3*T0)
            # Coefficients below x^(s-1) are the preceding row's fixed
            # compatibility equations.  Keep them as branch guards; they
            # need not vanish identically on the parent leaf.
            solved_f = {}
            for degree in range(s-1, 10):
                j = degree-(s-1)
                eq = sp.factor(coeff((N3_fixed + 57*sp.diff(z0, x)*S57
                                      - e*z0*sp.diff(S57, x)).subs(solved_f),
                                     x, degree))
                pivot = sp.factor(sp.diff(eq, f[j]))
                assert sp.simplify(pivot - mu*s*(57-4*j)) == 0
                solved_f[f[j]] = sp.factor(-eq.subs(f[j], 0)/pivot)
            S57b = sp.expand(S57.subs(solved_f))
            N3 = sp.expand(N3_fixed + 57*sp.diff(z0, x)*S57b
                           - e*z0*sp.diff(S57b, x))
            T3 = sp.cancel(-N3/x**10)
            # Lower coefficients of N3 are the prior fixed equations.  Keep
            # them unevaluated; on any candidate branch they vanish, making
            # this Laurent quotient the required polynomial T_3.

            # At e=24 the fourth Z-slice lies on the allowed [Z,R] cutoff,
            # so retain it as an arbitrary polynomial.  At the other three
            # vertices [Z,R] forces the fractional-power slice, with its
            # one-dimensional resonant x^(s-1) kernel kept free.
            if s == 6:
                eta = sp.symbols("eta6_0:7")
                z4 = sum(eta[j]*x**j for j in range(7))
            else:
                z4_base = mu*(
                    gamma*x**(s-5)*r16b
                    + sp.binomial(gamma, 2)*x**(s-10)
                    * (2*r19b*r17b + r18b**2)
                    + 3*sp.binomial(gamma, 3)*x**(s-15)*r19b**2*r18b
                    + sp.binomial(gamma, 4)*x**(s-20)*r19b**4
                )
                eta = sp.symbols(f"eta{s}")
                z4 = sp.expand(z4_base + eta*x**(s-1))
                assert sp.Poly(z4, x).is_univariate, (s, leaf_name, z4)

            C4 = (57*sp.diff(z1, x)*S57b-(e-1)*z1*sp.diff(S57b, x)
                  + 58*sp.diff(z2, x)*S58b-(e-2)*z2*sp.diff(S58b, x)
                  + 59*sp.diff(z3, x)*S59b-(e-3)*z3*sp.diff(S59b, x)
                  + 60*sp.diff(z4, x)*S60b-(e-4)*z4*sp.diff(S60b, x))
            N4 = sp.expand(C4 + q1*T3 + q2*T2 + q3*T1 + q4*T0)

            Ztrunc = (z0*y**e + z1*y**(e-1) + z2*y**(e-2)
                      + z3*y**(e-3) + z4*y**(e-4))
            Strunc = (S60b*y**60 + S59b*y**59 + S58b*y**58 + S57b*y**57)
            assert sp.expand(ycoeff(bracket(Ztrunc, Strunc, x, y), y, e+55)-C4) == 0
            Rtrunc = (x**5*y**20 + r19b*y**19 + r18b*y**18
                      + r17b*y**17 + r16b*y**16)
            if s >= 7:
                # The y^(e+15) slice is above the established [Z,R]
                # degree cap.  This verifies the forced fractional slice;
                # eta*x^(s-1) is in the resonant kernel.
                assert sp.expand(ycoeff(bracket(Ztrunc, Rtrunc, x, y),
                                        y, e+15)) == 0
            assert sp.expand(ycoeff(Rtrunc**2, y, 36)-q4) == 0

            fixed_map = {
                k: sp.factor(laurent_coeff(N4, x, k))
                for k in range(s-1)
            }
            fixed = [(k, value) for k, value in fixed_map.items() if value != 0]

            if s == 8 and leaf_name == "u0=0":
                expected_x3 = sp.Rational(96, 625)*mu*u[1]*(
                    203*a3**4-990*a3**2*rho[1]+25*rho[1]**2
                )
                assert sp.factor(fixed_map[3]-expected_x3) == 0
                reduced_x5 = sp.factor(fixed_map[5].subs({
                    a3: 0, rho[1]: 0, sigma[0]: 0
                }))
                assert sp.factor(
                    reduced_x5 - sp.Rational(16, 16875)*mu
                    * 55300*rho[2]**2*u[1]
                ) == 0
                reduced_x6 = sp.factor(fixed_map[6].subs({
                    a3: 0, rho[1]: 0, rho[2]: 0, sigma[0]: 0
                }))
                expected_x6 = sp.Rational(24, 25)*mu*(
                    552*a4*sigma[1]*u[1]
                    -535*c[3]*sigma[1]
                    -620*tau[0]*u[1]
                )
                assert sp.factor(reduced_x6-expected_x6) == 0

            if s == 8 and leaf_name == "u0!=0":
                assert sp.factor(
                    fixed_map[2]
                    - sp.Rational(42996609432, 5719140625)
                    * a3**4*mu*u[0]
                ) == 0
                reduced_x4 = sp.factor(fixed_map[4].subs({
                    a3: 0, rho[1]: 0, sigma[0]: 0
                }))
                expected_x4 = sp.Rational(4, 30883359375)*mu*(
                    504016425000*rho[2]**2*u[0]
                )
                assert sp.factor(reduced_x4-expected_x4) == 0
                reduced_x5 = sp.factor(fixed_map[5].subs({
                    a3: 0, rho[1]: 0, rho[2]: 0,
                    sigma[0]: 0, sigma[1]: 0
                }))
                expected_x5 = -sp.Rational(16, 30883359375)*mu*(
                    1185921000000*tau[0]*u[0]
                )
                assert sp.factor(reduced_x5-expected_x5) == 0

            if s == 7:
                if leaf_name == "rho1=a3^2/5":
                    expected_x2 = -sp.Rational(114436, 31875)*mu*a3**4*u[0]
                    expected_x3 = sp.Rational(56, 95625)*mu*a3*(
                        -5185*a3**3*u[1] + 13788*a3**2*a4*u[0]
                        - 101875*a3*rho[2]*u[0]
                        + 1022550*sigma[0]*u[0]
                    )
                else:
                    expected_x2 = -sp.Rational(388136, 31875)*mu*a3**4*u[0]
                    expected_x3 = -sp.Rational(112, 57375)*mu*a3*(
                        3553*a3**3*u[1] + 11454*a3**2*a4*u[0]
                        + 33512*a3*rho[2]*u[0]
                        - 306765*sigma[0]*u[0]
                    )
                assert sp.factor(fixed_map[2]-expected_x2) == 0
                assert sp.factor(fixed_map[3]-expected_x3) == 0
                reduced_x4 = sp.factor(fixed_map[4].subs(a3, 0))
                assert sp.factor(
                    reduced_x4-sp.Rational(28, 135)*mu*u[0]
                    *(2835*a4*sigma[0]+32*rho[2]**2)
                ) == 0
                assert sp.factor(
                    reduced_x4.subs(sigma[0], 0)
                    - sp.Rational(896, 135)*mu*rho[2]**2*u[0]
                ) == 0

            if args.compact:
                print(f"s={s} leaf={leaf_name} fixed_degrees={[k for k, _ in fixed]}", flush=True)
                for k, equation in fixed:
                    print(f"fixed_factor_x{k}={sp.factor(equation)}", flush=True)
            else:
                print(f"s={s} leaf={leaf_name} fixed_rows={fixed}", flush=True)

    print("scope=next exact quotient row with resonant Z slice; necessary only")


if __name__ == "__main__":
    main()
