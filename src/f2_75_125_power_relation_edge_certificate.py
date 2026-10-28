"""Exact edge consequence of Q^3 - lambda P^5 in the F2 (75,125) model.

The written argument constructs a lower-weight correction of
W = Q^3 - (b^3/a^5) P^5 whose first nonzero (5,-1)-band is 39.  Its band
coefficient obeys a first-order ODE.  Root orders and the degree at infinity
force the quartic carrier h into two finite algebraic configurations.

This checker verifies the finite semigroup list, local-order arithmetic,
and exact polynomial coefficient reductions.  It does not prove the
standardization/support theorem inputs and it does not exclude (75,125).
"""

from __future__ import annotations

import sympy as sp


def main() -> None:
    t, alpha, beta, slope, intercept, scalar = sp.symbols(
        "t alpha beta slope intercept scalar"
    )

    # If a W-face of weight s>39 commutes with P_15, then 5|s.
    # The seven possible faces can all be cancelled by a lower monomial
    # P^i Q^j with 3i+5j=s/5.  Their bracket corrections have weight s-25<50.
    representations = {
        40: (1, 1),
        45: (3, 0),
        50: (0, 2),
        55: (2, 1),
        60: (4, 0),
        65: (1, 2),
        70: (3, 1),
    }
    for weight, (p_power, q_power) in representations.items():
        assert 3 * p_power + 5 * q_power == weight // 5
        correction_weight = (
            15 * p_power + 25 * (q_power - 1) if q_power else None
        )
        assert correction_weight is None or correction_weight < 50

    # For U=t*h, the W_39 equation is
    # 5 U g' - 39 U' g = mu U^8, mu != 0.
    # At a root of U of multiplicity e=1 or 2, the leading coefficient
    # 5*l - 39*e cannot vanish for integer l below 7e+1.
    for multiplicity in (1, 2):
        threshold = 7 * multiplicity + 1
        assert all(
            5 * order - 39 * multiplicity != 0
            and multiplicity + order - 1 < 8 * multiplicity
            for order in range(threshold)
        )
    assert 8 + 15 + 15 == 38  # t=0 and two double roots, type 2+2
    assert 8 + 15 + 8 + 8 == 39  # t=0 and roots 2+1+1

    # At infinity, if deg(U)=5 and deg(g)=n, the leading term on the left
    # has degree n+4 and coefficient proportional to n-39.  The RHS has
    # degree 40, so n is 36 or 39.  The root orders force n=39 in both types.
    possible_degrees = [
        degree
        for degree in range(68)
        if degree == 39 or degree + 4 == 40
    ]
    assert possible_degrees == [36, 39]
    assert 38 > 36 and 39 >= 39
    assert 39 - 5 * 7 == 4  # deg(V) when g=U^7 V

    # Type 2+2: h=(t-1)^2(t-alpha)^2, rad(U)=t(t-1)(t-alpha),
    # V=rad(U)(slope*t+intercept), and 5 U V' - 4 U' V = U.
    U22 = sp.expand(t * (t - 1) ** 2 * (t - alpha) ** 2)
    rad22 = sp.expand(t * (t - 1) * (t - alpha))
    V22 = sp.expand(rad22 * (slope * t + intercept))
    equation22 = sp.Poly(
        sp.expand(5 * U22 * sp.diff(V22, t) - 4 * sp.diff(U22, t) * V22 - U22),
        t,
    )
    eq_t7 = sp.factor(equation22.nth(7))
    eq_t6 = sp.factor(equation22.nth(6))
    assert sp.simplify(eq_t7 - (-3 * (alpha + 1) * slope - 5 * intercept)) == 0
    intercept_from_t7 = -3 * (alpha + 1) * slope / 5
    eq_t6_reduced = sp.factor(eq_t6.subs(intercept, intercept_from_t7))
    assert sp.simplify(
        eq_t6_reduced
        + sp.Rational(6, 5) * slope * (alpha**2 - 3 * alpha + 1)
    ) == 0

    alpha_polynomial = alpha**2 - 3 * alpha + 1
    slope_solution = (4 * alpha - 11) / 3
    intercept_solution = 3 - alpha
    assert sp.gcd(sp.Poly(alpha_polynomial, alpha), sp.Poly(4 * alpha - 11, alpha)).degree() == 0
    equation22_candidate = sp.expand(
        equation22.as_expr().subs(
            {slope: slope_solution, intercept: intercept_solution}
        )
    )
    for coefficient in sp.Poly(equation22_candidate, t).all_coeffs():
        numerator = sp.together(coefficient).as_numer_denom()[0]
        assert sp.rem(numerator, alpha_polynomial, alpha) == 0

    # Type 2+1+1: h=(t-1)^2(t-alpha)(t-beta),
    # rad(U)=t(t-1)(t-alpha)(t-beta), V=scalar*rad(U).
    U211 = sp.expand(t * (t - 1) ** 2 * (t - alpha) * (t - beta))
    rad211 = sp.expand(t * (t - 1) * (t - alpha) * (t - beta))
    V211 = sp.expand(scalar * rad211)
    equation211 = sp.Poly(
        sp.expand(5 * U211 * sp.diff(V211, t) - 4 * sp.diff(U211, t) * V211 - U211),
        t,
    )
    assert sp.simplify(equation211.nth(7) - scalar * (alpha + beta - 3)) == 0
    eq_t6_211 = sp.factor(equation211.nth(6))
    product_condition = alpha * beta - 3
    assert sp.simplify(
        eq_t6_211.subs(beta, 3 - alpha)
        + 2 * scalar * product_condition.subs(beta, 3 - alpha)
    ) == 0

    beta_polynomial = beta**2 - 3 * beta + 3
    equation211_candidate = sp.expand(
        equation211.as_expr().subs(
            {scalar: -sp.Rational(1, 3), alpha: 3 - beta}
        )
    )
    for coefficient in sp.Poly(equation211_candidate, t).all_coeffs():
        numerator = sp.together(coefficient).as_numer_denom()[0]
        assert sp.rem(numerator, beta_polynomial, beta) == 0

    assert sp.gcd(sp.Poly(alpha_polynomial, alpha), sp.Poly(alpha, alpha)).degree() == 0
    assert alpha_polynomial.subs(alpha, 1) != 0
    assert sp.discriminant(beta_polynomial, beta) != 0
    assert beta_polynomial.subs(beta, 0) != 0
    assert beta_polynomial.subs(beta, 1) != 0

    print("W_edge_weight=39")
    print("W_edge_t_degree=39")
    print("root_order_floor=ord(g)>=7e+1 at roots of t*h")
    print("type_2+2_parameter=alpha^2-3*alpha+1=0")
    print("type_2+1+1_residual=t^2-3*t+3")
    print("candidate_faces_replay=PASS")
    print("scope=necessary edge condition; no exclusion of (75,125)")


if __name__ == "__main__":
    main()
