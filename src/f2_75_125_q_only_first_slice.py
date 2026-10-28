"""Exact audit of the first Q-only equation in the descending y-edge block.

At top-edge depth 61, the P slice has ended but q_39 is still present.  This
script checks the support degree of the residual and the full-rank structure
of its 13-dimensional new-layer operator for a generic degree-five carrier.
The symbolic leading-coefficient argument in the research note proves the
rank uniformly for every carrier with nonzero leading coefficient.
"""

from __future__ import annotations

import sympy as sp


x, a = sp.symbols("x a")
r0, r1, r2, r3, r4, r5 = sp.symbols("r0 r1 r2 r3 r4 r5")
R = r5 * x**5 + r4 * x**4 + r3 * x**3 + r2 * x**2 + r1 * x + r0


def support(i_max: int, j_max: int, edge_max: int) -> list[tuple[int, int]]:
    return [
        (i, j)
        for i in range(i_max + 1)
        for j in range(j_max + 1)
        if 5 * i - j <= edge_max
    ]


def image_basis(carrier: sp.Expr, scale: sp.Expr) -> list[sp.Poly]:
    return [
        sp.Poly(
            sp.expand(
                carrier**2
                * (
                    39 * sp.diff(carrier, x) * (3 * scale * x**d)
                    - 20 * carrier * sp.diff(3 * scale * x**d, x)
                )
            ),
            x,
        )
        for d in range(13)
    ]


def reduce_mod_image(
    polynomial: sp.Expr, carrier: sp.Expr, scale: sp.Expr
) -> sp.Poly:
    """Return the canonical degree-<=13 remainder modulo the image span."""
    remainder = sp.Poly(sp.expand(polynomial), x)
    basis = image_basis(carrier, scale)
    for d in reversed(range(13)):
        target_degree = d + 14
        pivot = sp.factor(basis[d].coeff_monomial(x**target_degree))
        assert pivot != 0
        quotient = sp.cancel(remainder.coeff_monomial(x**target_degree) / pivot)
        remainder = sp.Poly(
            sp.cancel(remainder.as_expr() - quotient * basis[d].as_expr()), x
        )
    assert remainder.is_zero or remainder.degree() <= 13
    return remainder


def main() -> None:
    p_support = support(15, 60, 15)
    q_support = support(25, 100, 25)
    assert len(p_support) == 586
    assert len(q_support) == 1576

    # The first Q-only equation in the descending top-y recurrence is the
    # coefficient of y^98.  Its fixed terms have y-exponents adding to 99.
    fixed_degrees = []
    for i, j in p_support:
        for u, v in q_support:
            if j + v != 99:
                continue
            if i * v - j * u == 0:
                continue
            fixed_degrees.append(i + u - 1)
    assert fixed_degrees
    assert max(fixed_degrees) == 26

    # At depth 61 the new P slice is zero and the new Q slice is q_39.
    # The support envelope gives deg(q_39) <= floor((39+25)/5) = 12.
    assert min(u for u, v in q_support if v == 39) == 0
    assert max(u for u, v in q_support if v == 39) == 12

    # T(x^d) = R^2 D_61(3a x^d), where
    # D_61(w) = 39 R' w - 20 R w'.  Its leading term is
    # 3 a r5^3 (195 - 20d) x^(d+14), which never vanishes for integral d.
    for d in range(13):
        monomial = x**d
        image = sp.Poly(
            sp.expand(
                R**2
                * (
                    39 * sp.diff(R, x) * (3 * a * monomial)
                    - 20 * R * sp.diff(3 * a * monomial, x)
                )
            ),
            x,
        )
        assert image.degree() == d + 14
        assert sp.factor(image.coeff_monomial(x ** (d + 14))) == (
            3 * a * r5**3 * (195 - 20 * d)
        )

    sample_R = x**5 + 2 * x + 1
    sample_a = sp.Integer(3)
    sample_basis = image_basis(sample_R, sample_a)
    sample_image = 2 * sample_basis[12] - 3 * sample_basis[4] + sample_basis[0]
    assert reduce_mod_image(sample_image.as_expr(), sample_R, sample_a).is_zero
    assert reduce_mod_image(1, sample_R, sample_a).as_expr() == 1

    print("P_support_count=586")
    print("Q_support_count=1576")
    print("fixed_y98_max_x_degree=26")
    print("q39_degree_bound=12")
    print("operator_domain_dimension=13")
    print("operator_rank=13")
    print("target_dimension=27")
    print("compatibility_codimension=14")
    print("image_reduction_controls=PASS")
    print("PASS: symbolic leading terms certify injectivity uniformly when r5 != 0")


if __name__ == "__main__":
    main()
