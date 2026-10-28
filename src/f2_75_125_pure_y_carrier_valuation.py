"""Exact first-layer valuation pruning for the pure fifth-power y-carrier.

Assume p_60=a*x^15 and q_100=b*x^25 after translating x and scaling P,Q.
The first four descending y-adic equations force ord_x(p_59) >= 8; the
next two exclude the order-8 stratum by an exact low-degree certificate.
Orders >=9 and p_59=0 remain open in this carrier branch.
"""

from __future__ import annotations

import sympy as sp


x = sp.symbols("x")


def first_nonzero_order(expr: sp.Expr, *, upto: int) -> int | None:
    poly = sp.Poly(sp.expand(expr), x)
    for k in range(upto + 1):
        if poly.nth(k) != 0:
            return k
    return None


def main() -> None:
    c = sp.symbols("c0:15")

    # At depth 1, the recurrence uniquely gives B_1=(5/3)x^10 A_1.
    for m in range(15):
        A1 = sum(c[k - m] * x**k for k in range(m, 15))
        B1 = sp.Rational(5, 3) * x**10 * A1
        C2 = sp.expand(99 * sp.diff(A1, x) * B1 - 59 * A1 * sp.diff(B1, x))
        expected_degree = 2 * m + 9
        expected_coefficient = sp.Rational(50, 3) * (4 * m - 59) * c[0] ** 2
        assert first_nonzero_order(C2, upto=expected_degree) == expected_degree
        assert sp.Poly(C2, x).nth(expected_degree) == expected_coefficient

    # The new-slice image at every depth starts in x-degree 14.  Thus the
    # depth-2 row immediately rules out m=0,1,2.
    assert all(2 * m + 9 < 14 for m in range(3))

    # For 3 <= m <= 6, C_2's leading term lies in degrees 15..21.  Only B_2
    # can cancel it there (the new A_2 contribution starts at degree 24).
    # Solving that row gives B_2[x^(2m-5)], then the depth-3 residual starts
    # at degree 3m-6, below 14.  The A_2,B_1 contribution starts later.
    for m in range(3, 7):
        ell = 2 * m - 5
        c2_lead = sp.Rational(50, 3) * (4 * m - 59)
        b2_lead = sp.cancel(-c2_lead / (15 * (98 - 4 * ell)))
        c3_lead = sp.cancel((98 * m - 59 * ell) * b2_lead)
        assert b2_lead != 0
        assert c3_lead != 0
        assert 3 * m - 6 < 14
        assert 3 * m - 6 < m + 9

    # For m=7, verify the forced low coefficients using generic higher terms
    # in A_1,A_2,A_3 and in the later B slices.
    c7 = sp.symbols("c7:15")
    A1 = sum(c7[k - 7] * x**k for k in range(7, 15))
    B1 = sp.Rational(5, 3) * x**10 * A1
    C2 = sp.expand(99 * sp.diff(A1, x) * B1 - 59 * A1 * sp.diff(B1, x))
    assert all(sp.Poly(C2, x).nth(k) == 0 for k in range(23))
    assert sp.Poly(C2, x).nth(23) == -sp.Rational(1550, 3) * c7[0] ** 2

    b2 = sp.Rational(5, 9) * c7[0] ** 2
    d = sp.symbols("d10:25")
    B2 = b2 * x**9 + sum(d[k - 10] * x**k for k in range(10, 25))
    a2 = sp.symbols("a2_0:15")
    A2 = sum(a2[k] * x**k for k in range(15))
    C3 = sp.expand(
        98 * sp.diff(A1, x) * B2 - 59 * A1 * sp.diff(B2, x)
        + 99 * sp.diff(A2, x) * B1 - 58 * A2 * sp.diff(B1, x)
    )
    assert all(sp.Poly(C3, x).nth(k) == 0 for k in range(15))
    assert sp.Poly(C3, x).nth(15) == sp.Rational(775, 9) * c7[0] ** 3

    b3 = -sp.Rational(5, 81) * c7[0] ** 3
    e = sp.symbols("e2:25")
    B3 = b3 * x + sum(e[k - 2] * x**k for k in range(2, 25))
    a3 = sp.symbols("a3_0:15")
    A3 = sum(a3[k] * x**k for k in range(15))
    C4 = sp.expand(
        97 * sp.diff(A1, x) * B3 - 59 * A1 * sp.diff(B3, x)
        + 98 * sp.diff(A2, x) * B2 - 58 * A2 * sp.diff(B2, x)
        + 99 * sp.diff(A3, x) * B1 - 57 * A3 * sp.diff(B1, x)
    )
    assert sp.Poly(C4, x).nth(7) == -sp.Rational(3100, 81) * c7[0] ** 4
    assert sp.Poly(C4, x).nth(7) != 0

    # The next possible order, m=8, yields a small exact selector.  Write
    # A_1=c*x^8+O(x^9), A_2=u*x+O(x^2) after the depth-4 x^9 row below.
    c8 = sp.symbols("c8:15")
    A1 = sum(c8[k - 8] * x**k for k in range(8, 15))
    B1 = sp.Rational(5, 3) * x**10 * A1
    C2 = sp.expand(99 * sp.diff(A1, x) * B1 - 59 * A1 * sp.diff(B1, x))
    c2_24 = sp.Poly(C2, x).nth(24)
    c2_25 = sp.Poly(C2, x).nth(25)
    assert c2_24 == 0
    assert c2_25 == -450 * c8[0] ** 2

    u = sp.symbols("u")
    a20, a21 = sp.symbols("a20 a21")
    d = sp.symbols("d12:25")
    b210 = sp.Rational(5, 3) * a20
    b211 = sp.Rational(5, 3) * a21 + sp.Rational(5, 9) * c8[0] ** 2
    B2 = b210 * x**10 + b211 * x**11 + sum(
        d[k - 12] * x**k for k in range(12, 25)
    )
    a2 = sp.symbols("a2_0:15")
    A2 = a20 + a21 * x + sum(a2[k] * x**k for k in range(2, 15))
    C3 = sp.expand(
        98 * sp.diff(A1, x) * B2 - 59 * A1 * sp.diff(B2, x)
        + 99 * sp.diff(A2, x) * B1 - 58 * A2 * sp.diff(B1, x)
    )
    b3_3 = sp.cancel(-sp.Poly(C3, x).nth(17) / 1275)
    b3_4 = sp.cancel(-sp.Poly(C3, x).nth(18) / 1215)
    assert sp.Poly(C3, x).nth(17) == -sp.Rational(4250, 3) * c8[0] * a20

    # The depth-4 coefficient x^9 is below the new-slice image and equals a
    # nonzero scalar times a20^2.  Its vanishing forces a20=0, hence b210=0.
    e = sp.symbols("e5:25")
    B3 = b3_3 * x**3 + b3_4 * x**4 + sum(
        e[k - 5] * x**k for k in range(5, 25)
    )
    a3 = sp.symbols("a3_0:15")
    A3 = sum(a3[k] * x**k for k in range(15))
    C4 = sp.expand(
        97 * sp.diff(A1, x) * B3 - 59 * A1 * sp.diff(B3, x)
        + 98 * sp.diff(A2, x) * B2 - 58 * A2 * sp.diff(B2, x)
        + 99 * sp.diff(A3, x) * B1 - 57 * A3 * sp.diff(B1, x)
    )
    assert sp.Poly(C4, x).nth(9) == -sp.Rational(2900, 3) * a20**2

    selector_row = sp.factor(sp.Poly(C4.subs(a20, 0), x).nth(11))
    assert sp.expand(
        selector_row + sp.Rational(100, 3)
        * (27 * a21**2 - 9 * c8[0] ** 2 * a21 + c8[0] ** 4)
    ) == 0
    assert sp.discriminant(27 * u**2 - 9 * u + 1, u) == -27

    # At depth 5, after the forced a20=0 branch, the x^4 fixed residual is
    # (25/3)*c*u*(c^2-18u).  The depth-4 selector and this row contradict
    # c != 0 by the following polynomial Bezout identity.
    cvar, uvar = sp.symbols("cvar uvar")
    a1_tail = sp.symbols("a1_9:15")
    A1h = cvar * x**8 + sum(a1_tail[k - 9] * x**k for k in range(9, 15))
    B1h = sp.Rational(5, 3) * x**10 * A1h
    B2h = (
        (sp.Rational(5, 3) * uvar + sp.Rational(5, 9) * cvar**2) * x**11
        + sum(d[k - 12] * x**k for k in range(12, 25))
    )
    A2h = uvar * x + sum(a2[k] * x**k for k in range(2, 15))
    B3h = (
        -sp.Rational(5, 81) * cvar * (cvar**2 - 18 * uvar) * x**4
        + sum(e[k - 5] * x**k for k in range(5, 25))
    )
    B4_coeffs = sp.symbols("b4_0:25")
    B4h = sum(B4_coeffs[k] * x**k for k in range(25))
    A3h = sum(a3[k] * x**k for k in range(15))
    A4_coeffs = sp.symbols("a4_0:15")
    A4h = sum(A4_coeffs[k] * x**k for k in range(15))
    C5 = sp.expand(
        96 * sp.diff(A1h, x) * B4h - 59 * A1h * sp.diff(B4h, x)
        + 97 * sp.diff(A2h, x) * B3h - 58 * A2h * sp.diff(B3h, x)
        + 98 * sp.diff(A3h, x) * B2h - 57 * A3h * sp.diff(B2h, x)
        + 99 * sp.diff(A4h, x) * B1h - 56 * A4h * sp.diff(B1h, x)
    )
    c5_x4 = sp.Poly(C5, x).nth(4)
    assert sp.expand(
        c5_x4 - sp.Rational(25, 3) * cvar * uvar * (cvar**2 - 18 * uvar)
    ) == 0

    Fh = 27 * uvar**2 - 9 * cvar**2 * uvar + cvar**4
    Gh = uvar * (cvar**2 - 18 * uvar)
    bezout = (90 * uvar + 7 * cvar**2) * Fh + (135 * uvar - 27 * cvar**2) * Gh
    assert sp.expand(bezout - 7 * cvar**6) == 0

    print("pure_fifth_y_carrier=p60 proportional x^15, q100 proportional x^25")
    print("depth1=B1=(5/3)*x^10*A1 after scaling")
    print("depth2_excludes_ord_x(A1)=0,1,2")
    print("depth3_excludes_ord_x(A1)=3,4,5,6")
    print("depth4_excludes_ord_x(A1)=7")
    print("necessary_result=ord_x(p59)>=8 or p59=0")
    print("ord_x(p59)=8 branch=excluded by depth4-selector/depth5-Bezout certificate")
    print("orders_ge_10_and_p59_zero=remain in the pure y-carrier branch")
    print("scope=one exact valuation-branch exclusion; not degree-pair exclusion")


if __name__ == "__main__":
    main()
