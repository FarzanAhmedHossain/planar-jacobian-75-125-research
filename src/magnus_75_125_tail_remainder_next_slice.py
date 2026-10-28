"""Exact next-slice constraints from R^2 | 1-[Z,S] in high Magnus branches.

For e=4s in {24,28,32,36}, this replays the y^(e+58) coefficient
after the already-derived top quotient slice.  It solves the x^8 and x^9
divisibility conditions for two coefficients of S_59.  These are necessary
conditions, not a branch exclusion.
"""

from __future__ import annotations

import sympy as sp


def coeff(poly: sp.Expr, x: sp.Symbol, degree: int) -> sp.Expr:
    return sp.Poly(sp.expand(poly), x).nth(degree)


def main() -> None:
    x = sp.symbols("x")
    mu, a3, a4, b0, b1 = sp.symbols("mu a3 a4 b0 b1", nonzero=True)
    b_rest = sp.symbols("b2:16")
    c = sp.symbols("c0:17")

    r19 = x**3 * (a3 + a4*x)
    for s in (6, 7, 8, 9):
        e = 4*s
        S60 = x**(11-s) * (b0 + b1*x + sum(b_rest[j-2]*x**j for j in range(2, s+7)))
        # The support caps give deg_x(r_19)<=4 and deg_x(S_60)<=17.
        assert sp.Poly(r19, x).degree() <= 4
        assert sp.Poly(S60, x).degree() <= 17
        S59 = sum(c[i]*x**i for i in range(17))

        Ztop = mu*x**s
        Zprev = mu*sp.Rational(s, 5)*x**(s-5)*r19

        top_bracket = sp.expand(60*sp.diff(Ztop, x)*S60 - e*Ztop*sp.diff(S60, x))
        assert sp.cancel(top_bracket / x**10).is_polynomial(x)
        Ttop = sp.cancel(-top_bracket / x**10)

        next_bracket = sp.expand(
            59*sp.diff(Ztop, x)*S59 - e*Ztop*sp.diff(S59, x)
            + 60*sp.diff(Zprev, x)*S60 - (e-1)*Zprev*sp.diff(S60, x)
        )
        # From R^2 T=1-[Z,S], the next quotient row is
        # x^10*T_(e+18) + 2*x^5*r19*T_(e+19) = -[Z,S]_(e+58).
        # Therefore x^8 and x^9 in this expression must vanish.
        residual = sp.expand(next_bracket + 2*x**5*r19*Ttop)

        first_expected = (
            mu*s*(23+4*s)*c[9-s]
            + sp.Rational(1, 5)*mu*s*(4*s**2-25*s-269)*a3*b0
        )
        second_expected = (
            mu*s*(19+4*s)*c[10-s]
            + sp.Rational(1, 5)*mu*s*(4*s**2-29*s-228)*a3*b1
            + sp.Rational(1, 5)*mu*s*(4*s**2-25*s-209)*a4*b0
        )
        assert sp.expand(coeff(residual, x, 8)-first_expected) == 0
        assert sp.expand(coeff(residual, x, 9)-second_expected) == 0

        solved_c8 = sp.solve(first_expected, c[9-s])[0]
        solved_c9 = sp.solve(second_expected, c[10-s])[0]
        assert solved_c8 == sp.cancel(
            -sp.Rational(4*s**2-25*s-269, 5*(23+4*s))*a3*b0
        )
        assert solved_c9 == sp.cancel(
            -sp.Rational(4*s**2-29*s-228, 5*(19+4*s))*a3*b1
            -sp.Rational(4*s**2-25*s-209, 5*(19+4*s))*a4*b0
        )
        print(
            f"e={e} S59_x{9-s}={sp.factor(solved_c8)} "
            f"S59_x{10-s}={sp.factor(solved_c9)}"
        )

    print("scope=exact quotient-row constraints; no high-branch exclusion")


if __name__ == "__main__":
    main()
