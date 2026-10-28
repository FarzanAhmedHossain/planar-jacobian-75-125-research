"""Exact positive control for independent F2 edge and upper-y carriers.

For a quartic h and degree-five R with matching leading coefficient, construct
a polynomial A with those two faces, then take P=A^3 and Q=A^5. The script
checks the finite support arithmetic and the chain-rule Jacobian identity
without expanding the large fifth power. It is a hostile control for
homogeneous/face-only eliminations, not a Keller pair.
"""

from __future__ import annotations

import json

import sympy as sp


def support(poly: sp.Expr, x: sp.Symbol, y: sp.Symbol) -> set[tuple[int, int]]:
    return {
        monomial
        for monomial, _coeff in sp.Poly(sp.expand(poly), x, y).terms()
    }


def face(poly: sp.Expr, x: sp.Symbol, y: sp.Symbol, weight: int) -> sp.Expr:
    result = sp.Integer(0)
    for (i, j), coeff in sp.Poly(sp.expand(poly), x, y).terms():
        if 5 * i - j == weight:
            result += coeff * x**i * y**j
    return sp.expand(result)


def coefficient_slice(poly: sp.Expr, x: sp.Symbol, y: sp.Symbol, degree: int) -> sp.Expr:
    result = sp.Integer(0)
    for (i, j), coeff in sp.Poly(sp.expand(poly), x, y).terms():
        if j == degree:
            result += coeff * x**i
    return sp.expand(result)


def total_degree_face(
    poly: sp.Expr, x: sp.Symbol, y: sp.Symbol, degree: int
) -> sp.Expr:
    result = sp.Integer(0)
    for (i, j), coeff in sp.Poly(sp.expand(poly), x, y).terms():
        if i + j == degree:
            result += coeff * x**i * y**j
    return sp.expand(result)


def minkowski_power(
    exponents: set[tuple[int, int]], power: int
) -> set[tuple[int, int]]:
    result = {(0, 0)}
    for _ in range(power):
        result = {
            (i + u, j + v)
            for i, j in result
            for u, v in exponents
        }
    return result


def in_envelope(
    monomials: set[tuple[int, int]], i_max: int, j_max: int, weight_max: int
) -> bool:
    return all(
        0 <= i <= i_max and 0 <= j <= j_max and 5 * i - j <= weight_max
        for i, j in monomials
    )


def main() -> None:
    x, y, t = sp.symbols("x y t")
    s0, s1, s2 = sp.symbols("s0 s1 s2", nonzero=True)
    r0, r1, r2, r3, r4 = sp.symbols("r0 r1 r2 r3 r4")

    # h has the double root at 1 forced by the F2 complete-chain successor.
    # The lower coefficients of R are independent of the residual quadratic
    # in h; only the common-corner coefficient s2 is shared.
    h = sp.expand((t - 1) ** 2 * (s0 + s1 * t + s2 * t**2))
    R = s2 * x**5 + r4 * x**4 + r3 * x**3 + r2 * x**2 + r1 * x + r0
    edge_carrier = sp.expand(x * h.subs(t, x * y**5))
    A = sp.expand(edge_carrier + y**20 * (R - s2 * x**5))

    assert sp.degree(h, t) == 4
    assert sp.Poly(h, t).nth(0) == s0
    assert sp.Poly(h, t).nth(4) == s2
    assert sp.simplify(h.subs(t, 1)) == 0
    assert sp.simplify(sp.diff(h, t).subs(t, 1)) == 0
    assert sp.degree(R, x) == 5

    support_a = support(A, x, y)
    assert in_envelope(support_a, 5, 20, 5)
    assert face(A, x, y, 5) == edge_carrier
    assert coefficient_slice(A, x, y, 20) == R
    assert sp.Poly(A, x, y).total_degree() == 25
    assert total_degree_face(A, x, y, 25) == s2 * x**5 * y**20

    support_p = minkowski_power(support_a, 3)
    support_q = minkowski_power(support_a, 5)
    assert in_envelope(support_p, 15, 60, 15)
    assert in_envelope(support_q, 25, 100, 25)
    assert max(i + j for i, j in support_p) == 75
    assert max(i + j for i, j in support_q) == 125

    # The highest weighted face and highest y-slice of a product are the
    # corresponding products of the highest faces/slices of its factors.
    # This proves the carrier identities for P=A^3 and Q=A^5 without a
    # costly full symbolic expansion.
    p_face = sp.expand(edge_carrier**3)
    q_face = sp.expand(edge_carrier**5)
    p_y_top = sp.expand(R**3)
    q_y_top = sp.expand(R**5)

    # Exact chain-rule certificate: d(A^3)=3A^2 dA and d(A^5)=5A^4 dA.
    p_x, p_y = 3 * A**2 * sp.diff(A, x), 3 * A**2 * sp.diff(A, y)
    q_x, q_y = 5 * A**4 * sp.diff(A, x), 5 * A**4 * sp.diff(A, y)
    jacobian = sp.expand(p_x * q_y - p_y * q_x)
    assert jacobian == 0

    print(
        json.dumps(
            {
                "A_support_count": len(support_a),
                "P_possible_support_count": len(support_p),
                "Q_possible_support_count": len(support_q),
                "edge_carrier": "x*h(x*y^5), h=(t-1)^2*(s0+s1*t+s2*t^2)",
                "upper_y_carrier": "R=s2*x^5+r4*x^4+r3*x^3+r2*x^2+r1*x+r0",
                "carrier_parameters_independent_except_shared_corner": True,
                "P_Q_support_in_exhaustive_outer_envelopes": True,
                "P_edge": "edge_carrier^3",
                "Q_edge": "edge_carrier^5",
                "P_y60_slice": "R^3",
                "Q_y100_slice": "R^5",
                "degrees": [75, 125],
                "jacobian": 0,
                "scope": (
                    "exact positive control for support and homogeneous face rows; "
                    "not a Keller pair, since the required Jacobian is nonzero"
                ),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
