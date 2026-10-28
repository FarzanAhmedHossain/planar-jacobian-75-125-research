"""Exact exclusion of ord_x(p_59)=9 in the pure fifth-power y-carrier.

This extends the same top-y recurrence used by the order-8 certificate.
It keeps generic higher coefficients in every slice that could contribute
to the tested low x-degrees and verifies the final guarded Bezout identity.
"""

from __future__ import annotations

import sympy as sp


x = sp.symbols("x")


def coeff(poly: sp.Expr, degree: int) -> sp.Expr:
    return sp.Poly(sp.expand(poly), x).nth(degree)


def main() -> None:
    c = sp.symbols("c")
    a1_tail = sp.symbols("a1_10:15")
    A1 = c * x**9 + sum(a1_tail[k - 10] * x**k for k in range(10, 15))
    B1 = sp.Rational(5, 3) * x**10 * A1

    a2 = sp.symbols("a2_0:15")
    A2 = sum(a2[k] * x**k for k in range(15))
    C2 = sp.expand(
        99 * sp.diff(A1, x) * B1 - 59 * A1 * sp.diff(B1, x)
    )

    # Depth 2 fixes B2_10..B2_13 from rows x^24..x^27.
    b2 = sp.symbols("b2_10:14")
    B2 = sum(b2[k - 10] * x**k for k in range(10, 14))
    for k in range(14, 25):
        B2 += sp.symbols(f"b2_{k}") * x**k
    E2 = sp.expand(
        C2
        + 98 * 15 * x**14 * B2 - 60 * x**15 * sp.diff(B2, x)
        + 100 * sp.diff(A2, x) * x**25 - 58 * A2 * 25 * x**24
    )
    forced_b2 = {
        b2[0]: sp.Rational(5, 3) * a2[0],
        b2[1]: sp.Rational(5, 3) * a2[1],
        b2[2]: sp.Rational(5, 3) * a2[2],
        b2[3]: sp.Rational(5, 3) * a2[3] + sp.Rational(5, 9) * c**2,
    }
    assert all(coeff(E2, degree) == 0 for degree in range(14, 24))
    for degree in range(24, 28):
        assert sp.expand(coeff(E2, degree).subs(forced_b2)) == 0

    B2 = sp.expand(B2.subs(forced_b2))
    C3 = sp.expand(
        98 * sp.diff(A1, x) * B2 - 59 * A1 * sp.diff(B2, x)
        + 99 * sp.diff(A2, x) * B1 - 58 * A2 * sp.diff(B1, x)
    )
    b3low = []
    for ell in range(10):
        b3low.append(
            sp.cancel(-coeff(C3, ell + 14) / (15 * (97 - 4 * ell)))
        )
    B3 = sum(b3low[ell] * x**ell for ell in range(10))
    for k in range(10, 25):
        B3 += sp.symbols(f"b3_{k}") * x**k

    a3 = sp.symbols("a3_0:15")
    A3 = sum(a3[k] * x**k for k in range(15))
    C4 = sp.expand(
        97 * sp.diff(A1, x) * B3 - 59 * A1 * sp.diff(B3, x)
        + 98 * sp.diff(A2, x) * B2 - 58 * A2 * sp.diff(B2, x)
        + 99 * sp.diff(A3, x) * B1 - 57 * A3 * sp.diff(B1, x)
    )

    # The first three sub-14 rows successively force a2_0=a2_1=a2_2=0.
    assert sp.expand(coeff(C4, 9)) == -sp.Rational(2900, 3) * a2[0] ** 2
    C4_a20 = sp.expand(C4.subs(a2[0], 0))
    assert sp.expand(coeff(C4_a20, 11) + 900 * a2[1] ** 2) == 0
    C4_a21 = sp.expand(C4_a20.subs(a2[1], 0))
    assert sp.expand(coeff(C4_a21, 13) + sp.Rational(2500, 3) * a2[2] ** 2) == 0
    C4_small = sp.expand(C4_a21.subs(a2[2], 0))

    # With those coefficients zero, the depth-4 x^15 selector is F=0.
    u = a2[3]
    F = 27 * u**2 - 9 * c**2 * u + c**4
    C4_x15 = sp.factor(coeff(C4_small, 15))
    assert sp.expand(C4_x15 + sp.Rational(2300, 81) * F) == 0

    # B4_0=0 and B4_1=-C4_x15/1380 from the depth-4 rows x^14,x^15.
    assert coeff(C4_small, 14) == 0
    B4_1 = sp.cancel(-C4_x15 / 1380)
    B4 = B4_1 * x
    for k in range(2, 25):
        B4 += sp.symbols(f"b4_{k}") * x**k

    a4 = sp.symbols("a4_0:15")
    A4 = sum(a4[k] * x**k for k in range(15))
    C5 = sp.expand(
        96 * sp.diff(A1, x) * B4 - 59 * A1 * sp.diff(B4, x)
        + 97 * sp.diff(A2, x) * B3 - 58 * A2 * sp.diff(B3, x)
        + 98 * sp.diff(A3, x) * B2 - 57 * A3 * sp.diff(B2, x)
        + 99 * sp.diff(A4, x) * B1 - 56 * A4 * sp.diff(B1, x)
    )
    C5_x9 = sp.factor(coeff(C5, 9).subs({a2[0]: 0, a2[1]: 0, a2[2]: 0}))
    expected = sp.Rational(575, 81) * c * u * (c**2 - 18 * u)
    assert sp.expand(C5_x9 - expected - sp.Rational(4025, 243) * c * F) == 0

    G = u * (c**2 - 18 * u)
    bezout = (90 * u + 7 * c**2) * F + (135 * u - 27 * c**2) * G
    assert sp.expand(bezout - 7 * c**6) == 0

    print("pure_fifth_y_carrier_ord_x(p59)=9=excluded")
    print("depth4_selector=27*u^2-9*c^2*u+c^4")
    print("depth5_residual_on_selector=(575/81)*c*u*(c^2-18*u)")
    print("guard=c!=0; Bezout_identity=7*c^6, hence contradiction")
    print("scope=exact valuation-branch exclusion; orders>=10 and p59=0 remain")


if __name__ == "__main__":
    main()
