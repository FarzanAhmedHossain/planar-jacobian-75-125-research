"""First coefficient obstructions in the next R^2-remainder slice.

This continues the y^(e+58) quotient row to y^(e+57) for the four high
Magnus vertices.  It checks all x-degrees below the image of the new S_58
slice, retaining the exact K4/K5 branch pins available at each vertex.
"""

from __future__ import annotations

import sympy as sp


def bracket(f: sp.Expr, g: sp.Expr, x: sp.Symbol, y: sp.Symbol) -> sp.Expr:
    return sp.diff(f, x)*sp.diff(g, y)-sp.diff(f, y)*sp.diff(g, x)


def main() -> None:
    x = sp.symbols("x")
    y = sp.symbols("y")
    mu = sp.symbols("mu", nonzero=True)
    a3, a4 = sp.symbols("a3 a4")
    rho0, rho1, rho2, rho3, rho4 = sp.symbols(
        "rho0 rho1 rho2 rho3 rho4"
    )
    u = sp.symbols("u0:16")
    c = sp.symbols("c0:17")
    d = sp.symbols("d0:17")

    r19 = x**3*(a3+a4*x)
    r18 = rho0+rho1*x+rho2*x**2+rho3*x**3+rho4*x**4
    S59_full = sum(c[i]*x**i for i in range(17))
    S58 = sum(d[i]*x**i for i in range(17))

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

    for s in (6, 7, 8, 9):
        e = 4*s
        S60 = x**(11-s)*sum(u[j]*x**j for j in range(s+7))
        S59 = sp.expand(S59_full.subs(c_pins[s]))
        Z0 = mu*x**s
        Z1 = mu*sp.Rational(s, 5)*x**(s-5)*r19
        Z2 = mu*(
            sp.Rational(s, 5)*x**(s-5)*r18
            + sp.binomial(sp.Rational(s, 5), 2)*x**(s-10)*r19**2
        )

        C_top = sp.expand(60*sp.diff(Z0, x)*S60-e*Z0*sp.diff(S60, x))
        R = x**5*y**20+r19*y**19+r18*y**18
        S = S60*y**60+S59*y**59+S58*y**58
        Z = Z0*y**e+Z1*y**(e-1)+Z2*y**(e-2)

        # Check the hand-collected slices against direct bivariate brackets.
        direct_top = sp.Poly(
            sp.expand(bracket(Z0*y**e, S60*y**60, x, y)), y
        ).coeff_monomial(y**(e+59))
        direct_next = sp.Poly(
            sp.expand(bracket(Z0*y**e+Z1*y**(e-1), S, x, y)), y
        ).coeff_monomial(y**(e+58))
        direct_low = sp.Poly(sp.expand(bracket(Z, S, x, y)), y).coeff_monomial(
            y**(e+57)
        )
        direct_R2_38 = sp.Poly(sp.expand(R**2), y).coeff_monomial(y**38)
        assert sp.expand(direct_top-C_top) == 0

        Ttop = sp.cancel(-C_top/x**10)
        C_next = sp.expand(
            59*sp.diff(Z0, x)*S59-e*Z0*sp.diff(S59, x)
            +60*sp.diff(Z1, x)*S60-(e-1)*Z1*sp.diff(S60, x)
        )
        assert sp.expand(direct_next-C_next) == 0
        Tnext = sp.cancel(-(C_next+2*x**5*r19*Ttop)/x**10)
        assert sp.Poly(Tnext, x).is_univariate

        C_low = sp.expand(
            58*sp.diff(Z0, x)*S58-e*Z0*sp.diff(S58, x)
            +59*sp.diff(Z1, x)*S59-(e-1)*Z1*sp.diff(S59, x)
            +60*sp.diff(Z2, x)*S60-(e-2)*Z2*sp.diff(S60, x)
        )
        R2_38 = sp.expand(r19**2+2*x**5*r18)
        assert sp.expand(direct_R2_38-R2_38) == 0
        residual = sp.expand(
            C_low+2*x**5*r19*Tnext+R2_38*Ttop
        )
        assert sp.expand(direct_low-C_low) == 0

        substitutions = {}
        if s in (7, 8, 9):
            substitutions[rho0] = 0  # x | r_18 from the literal [Z,R] rows.
        residual = sp.expand(residual.subs(substitutions))

        nonzero = []
        for degree in range(s-1):
            value = sp.factor(sp.expand(sp.Poly(residual, x).nth(degree)))
            if value != 0:
                nonzero.append((degree, value))

        # The only rows below the S_58 image start occur at the last two
        # branches.  These are the expected exact split selectors.
        if s in (6, 7):
            assert nonzero == [], (s, nonzero)
        elif s == 8:
            expected = sp.Rational(48, 1375)*mu*u[0]*(2529*a3**2-15125*rho1)
            assert nonzero == [(6, sp.factor(expected))], nonzero
            eq8 = 2529*a3**2-15125*rho1
            assert sp.factor(nonzero[0][1]-sp.Rational(48, 1375)*mu*u[0]*eq8) == 0
        else:
            assert [degree for degree, _ in nonzero] == [6, 7], nonzero
            eq9_6 = 59*a3**2-240*rho1
            eq9_7 = (
                55*a3**2*u[1]+114*a3*a4*u[0]
                -235*rho1*u[1]-190*rho2*u[0]
            )
            assert sp.factor(nonzero[0][1]-sp.Rational(54, 25)*mu*u[0]*eq9_6) == 0
            assert sp.factor(nonzero[1][1]-sp.Rational(54, 25)*mu*eq9_7) == 0

            # Complete the zero/nonzero split without treating u0 or u1 as
            # globally invertible.  The divisions occur only on stated leaves.
            rho1_u0 = sp.Rational(59, 240)*a3**2
            rho2_u0 = (
                sp.Rational(57, 95)*a3*a4
                - sp.Rational(133, 9120)*a3**2*u[1]/u[0]
            )
            assert sp.factor(eq9_6.subs(rho1, rho1_u0)) == 0
            assert sp.factor(eq9_7.subs({rho1: rho1_u0, rho2: rho2_u0})) == 0
            assert sp.factor(eq9_7.subs(u[0], 0)-u[1]*(55*a3**2-235*rho1)) == 0
            assert sp.Rational(55, 235) == sp.Rational(11, 47)

        print(f"s={s} fixed_rows_below_S58_image={nonzero}")

    print("scope=necessary remainder constraints and exact zero/nonzero split; no branch exclusion")


if __name__ == "__main__":
    main()
