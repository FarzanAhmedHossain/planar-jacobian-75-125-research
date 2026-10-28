"""Next low-x obstructions from the y^(e+56) remainder row.

Starting from the exact y^(e+57) zero/nonzero leaves, solve the coefficients
of S_58 needed to make the preceding quotient slice polynomial.  Then compute
the fixed x-degrees below the image of the new S_57 slice.  This is necessary
condition arithmetic only, not a branch exclusion.
"""

from __future__ import annotations

import sympy as sp


def bracket(f: sp.Expr, g: sp.Expr, x: sp.Symbol, y: sp.Symbol) -> sp.Expr:
    return sp.diff(f, x)*sp.diff(g, y)-sp.diff(f, y)*sp.diff(g, x)


def ycoeff(expr: sp.Expr, y: sp.Symbol, degree: int) -> sp.Expr:
    return sp.Poly(sp.expand(expr), y).coeff_monomial(y**degree)


def xcoeff(expr: sp.Expr, x: sp.Symbol, degree: int) -> sp.Expr:
    return sp.Poly(sp.expand(expr), x).nth(degree)


def main() -> None:
    x, y = sp.symbols("x y")
    mu = sp.symbols("mu", nonzero=True)
    a3, a4 = sp.symbols("a3 a4")
    rho0, rho1, rho2, rho3, rho4 = sp.symbols(
        "rho0 rho1 rho2 rho3 rho4"
    )
    sigma = sp.symbols("sigma0:5")
    u = sp.symbols("u0:16")
    c = sp.symbols("c0:17")
    d = sp.symbols("d0:18")

    r19 = x**3*(a3+a4*x)
    r18 = rho0+rho1*x+rho2*x**2+rho3*x**3+rho4*x**4
    r17 = sum(sigma[j]*x**j for j in range(5))
    S59_full = sum(c[j]*x**j for j in range(17))
    S58_full = sum(d[j]*x**j for j in range(18))

    # Coefficients c_i previously pinned by x^8 and x^9 at y^(e+58).
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

    # Full zero/nonzero leaves of the already verified first two quotient rows.
    leaves = {
        6: [
            ("a3=0", {a3: 0}),
            ("rho0=0", {rho0: 0}),
        ],
        7: [
            ("rho1=a3^2/5", {rho0: 0, rho1: a3**2/sp.Integer(5)}),
            ("rho1=2a3^2/5", {rho0: 0, rho1: 2*a3**2/sp.Integer(5)}),
        ],
        8: [
            ("u0=0", {rho0: 0, u[0]: 0}),
            ("u0!=0", {rho0: 0, rho1: sp.Rational(2529, 15125)*a3**2}),
        ],
        9: [
            ("u0!=0", {
                rho0: 0,
                rho1: sp.Rational(59, 240)*a3**2,
                rho2: sp.Rational(57, 95)*a3*a4
                     -sp.Rational(133, 9120)*a3**2*u[1]/u[0],
            }),
            ("u0=0,u1=0", {rho0: 0, u[0]: 0, u[1]: 0}),
            ("u0=0,u1!=0", {
                rho0: 0, u[0]: 0,
                rho1: sp.Rational(11, 47)*a3**2,
            }),
        ],
    }

    for s in (6, 7, 8, 9):
        e = 4*s
        gamma = sp.Rational(s, 5)
        for leaf_name, branch_subs in leaves[s]:
            S60 = x**(11-s)*sum(u[j]*x**j for j in range(s+7))
            S59 = sp.expand(S59_full.subs(c_pins[s]))
            S58 = S58_full
            r18_b = sp.expand(r18.subs(branch_subs))
            r19_b = sp.expand(r19.subs(branch_subs))
            r17_b = sp.expand(r17.subs(branch_subs))
            S60_b = sp.expand(S60.subs(branch_subs))
            S59_b = sp.expand(S59.subs(branch_subs))
            S58_b = sp.expand(S58.subs(branch_subs))

            z0 = sp.expand((mu*x**s).subs(branch_subs))
            z1 = sp.expand((mu*gamma*x**(s-5)*r19).subs(branch_subs))
            z2 = sp.expand((mu*(gamma*x**(s-5)*r18
                         +sp.binomial(gamma, 2)*x**(s-10)*r19**2))
                         .subs(branch_subs))
            z3 = sp.expand((mu*(gamma*x**(s-5)*r17
                         +2*sp.binomial(gamma, 2)*x**(s-10)*r19*r18
                         +sp.binomial(gamma, 3)*x**(s-15)*r19**3))
                         .subs(branch_subs))
            r19_b = sp.expand(r19_b)

            # Rows m=0,1,2 of R^2*T=1-[Z,S], with m=2 using S_58.
            C0 = 60*sp.diff(z0, x)*S60_b-e*z0*sp.diff(S60_b, x)
            C1 = (59*sp.diff(z0, x)*S59_b-e*z0*sp.diff(S59_b, x)
                  +60*sp.diff(z1, x)*S60_b-(e-1)*z1*sp.diff(S60_b, x))
            C2 = (58*sp.diff(z0, x)*S58_b-e*z0*sp.diff(S58_b, x)
                  +59*sp.diff(z1, x)*S59_b-(e-1)*z1*sp.diff(S59_b, x)
                  +60*sp.diff(z2, x)*S60_b-(e-2)*z2*sp.diff(S60_b, x))

            q1 = sp.expand(2*x**5*r19_b)
            q2 = sp.expand(r19_b**2+2*x**5*r18_b)
            q3 = sp.expand(2*x**5*r17_b+2*r19_b*r18_b)
            T0 = sp.cancel(-C0/x**10)
            T1 = sp.cancel(-(C1+q1*T0)/x**10)
            assert sp.Poly(T0, x).is_univariate
            assert sp.Poly(T1, x).is_univariate

            N2 = sp.expand(C2+q1*T1+q2*T0)
            # Earlier low rows must vanish below the S_58 image.  The
            # remaining coefficients x^(s-1),...,x^9 uniquely solve
            # d_0,...,d_(10-s), since their multipliers are nonzero.
            for degree in range(s-1):
                assert sp.factor(xcoeff(N2, x, degree)) == 0, (
                    s, leaf_name, degree, sp.factor(xcoeff(N2, x, degree))
                )
            solved_d = {}
            for degree in range(s-1, 10):
                j = degree-(s-1)
                assert j <= 10-s
                eq = sp.factor(xcoeff(N2, x, degree).subs(solved_d))
                pivot = sp.factor(sp.diff(eq, d[j]))
                assert sp.simplify(pivot-mu*s*(58-4*j)) == 0
                value = sp.factor(-eq.subs(d[j], 0)/pivot)
                solved_d[d[j]] = value
            N2 = sp.expand(N2.subs(solved_d))
            S58_solved = sp.expand(S58_b.subs(solved_d))
            T2 = sp.cancel(-N2/x**10)
            assert sp.Poly(T2, x).is_univariate
            assert all(xcoeff(N2, x, degree) == 0 for degree in range(10))

            # The next commutator slice below S_57 is [z1,S58]+[z2,S59]+[z3,S60].
            C3 = (58*sp.diff(z1, x)*S58_solved-(e-1)*z1*sp.diff(S58_solved, x)
                  +59*sp.diff(z2, x)*S59_b-(e-2)*z2*sp.diff(S59_b, x)
                  +60*sp.diff(z3, x)*S60_b-(e-3)*z3*sp.diff(S60_b, x))

            # Independently verify the collected commutator and R^2 slices.
            Zfull = z0*y**e+z1*y**(e-1)+z2*y**(e-2)+z3*y**(e-3)
            Sfull = S60_b*y**60+S59_b*y**59+S58_solved*y**58
            direct_C3 = ycoeff(bracket(Zfull, Sfull, x, y), y, e+56)
            Rtrunc = x**5*y**20+r19_b*y**19+r18_b*y**18+r17_b*y**17
            direct_q3 = ycoeff(Rtrunc**2, y, 37)
            assert sp.expand(direct_C3-C3) == 0
            assert sp.expand(direct_q3-q3) == 0

            residual = sp.expand(C3+q1*T2+q2*T1+q3*T0)
            fixed = []
            for degree in range(s-1):
                value = sp.factor(xcoeff(residual, x, degree))
                if value != 0:
                    fixed.append((degree, value))

            by_degree = dict(fixed)
            if s == 6 and leaf_name == "a3=0":
                assert fixed == []
            elif s == 6:
                expected = sp.Rational(49288, 3525)*a3*mu*u[0]*(a3**2-5*rho1)
                assert sp.factor(by_degree[4]-expected) == 0
            elif s == 7:
                expected4 = (
                    sp.Rational(4004, 1275)*a3**3*mu*u[0]
                    if leaf_name == "rho1=a3^2/5"
                    else sp.Rational(2912, 255)*a3**3*mu*u[0]
                )
                assert sp.factor(by_degree[4]-expected4) == 0
                expected5_a3zero = (
                    -sp.Rational(28*240975, 11475)*mu*sigma[0]*u[0]
                )
                assert sp.factor(by_degree[5].subs(a3, 0)-expected5_a3zero) == 0
            elif s == 8 and leaf_name == "u0!=0":
                expected4 = -sp.Rational(123872664, 20796875)*a3**3*mu*u[0]
                assert sp.factor(by_degree[4]-expected4) == 0
                expected5_a3zero = (
                    -sp.Rational(8*29266875, 378125)*mu*sigma[0]*u[0]
                )
                assert sp.factor(by_degree[5].subs(a3, 0)-expected5_a3zero) == 0
                expected6_a3sigma0zero = (
                    -sp.Rational(4*1335386250, 10209375)*mu*u[0]
                    *(sigma[1]-sp.Rational(2, 5)*a4*rho2)
                )
                assert sp.factor(
                    by_degree[6].subs({a3: 0, sigma[0]: 0})
                    -expected6_a3sigma0zero
                ) == 0
            elif s == 8:
                expected5 = -sp.Rational(3488, 125)*a3*mu*u[1]*(a3**2-5*rho1)
                assert sp.factor(by_degree[5]-expected5) == 0
            elif s == 9 and leaf_name == "u0!=0":
                expected4 = sp.Rational(2409, 1000)*a3**3*mu*u[0]
                expected5_a3zero = -sp.Rational(3*415200, 2000)*mu*sigma[0]*u[0]
                expected6_a3sigma0zero = (
                    -sp.Rational(3*343200, 2000)*mu*sigma[1]*u[0]
                )
                expected7_a3sigma01zero = (
                    sp.Rational(3, 2000)*mu*u[0]
                    *(-43392*a4**3+162720*a4*rho3-271200*sigma[2])
                )
                assert sp.factor(by_degree[4]-expected4) == 0
                assert sp.factor(by_degree[5].subs(a3, 0)-expected5_a3zero) == 0
                assert sp.factor(
                    by_degree[6].subs({a3: 0, sigma[0]: 0})
                    -expected6_a3sigma0zero
                ) == 0
                assert sp.factor(
                    by_degree[7].subs({a3: 0, sigma[0]: 0, sigma[1]: 0})
                    -expected7_a3sigma01zero
                ) == 0
            elif s == 9 and leaf_name == "u0=0,u1!=0":
                expected5 = -sp.Rational(459, 5875)*a3**3*mu*u[1]
                expected6_a3zero = (
                    -sp.Rational(9*398325, 5875)*mu*sigma[0]*u[1]
                )
                expected7_a3sigma0zero = (
                    -sp.Rational(9, 5875)*mu
                    *(-289050*a4*rho2*u[1]+263200*c[2]*rho2
                      +327825*sigma[1]*u[1])
                )
                assert sp.factor(by_degree[5]-expected5) == 0
                assert sp.factor(by_degree[6].subs(a3, 0)-expected6_a3zero) == 0
                assert sp.factor(
                    by_degree[7].subs({a3: 0, sigma[0]: 0})
                    -expected7_a3sigma0zero
                ) == 0
            print(f"s={s} leaf={leaf_name} fixed_rows={fixed}")

    print("scope=next necessary quotient row; zero cases are separate leaves; no branch exclusion")


if __name__ == "__main__":
    main()
