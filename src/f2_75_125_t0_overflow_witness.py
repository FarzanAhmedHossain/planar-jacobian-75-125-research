"""Exact t=0 support-face witness for homogeneous bands 30--35.

This demonstrates that the first forbidden Q pole terms at band distances
30--35 can cancel in the homogeneous first-integral block.  It is deliberately
only a local support-face certificate: it does not check polynomiality at the
roots of h, the inhomogeneous band 36, or the full Keller equation.
"""

from __future__ import annotations

from fractions import Fraction


def binomial_rational(a: Fraction, n: int) -> Fraction:
    """Generalized binomial coefficient (a choose n), exactly."""
    value = Fraction(1)
    for j in range(n):
        value *= (a - j) / (j + 1)
    return value


def poly_mul(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    out = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return out


def poly_pow(poly: list[Fraction], exponent: int) -> list[Fraction]:
    out = [Fraction(1)]
    for _ in range(exponent):
        out = poly_mul(out, poly)
    return out


def scalar_mul(poly: list[Fraction], scalar: Fraction) -> list[Fraction]:
    return [scalar * value for value in poly]


def homogeneous_face_coefficient(power: int, kappa: dict[int, Fraction]) -> Fraction:
    """[z^power] sum_r kappa_r z^r (1+z)^((5-r)/3)."""
    value = Fraction(0)
    for r, scalar in kappa.items():
        if r <= power:
            exponent = Fraction(5 - r, 3)
            value += scalar * binomial_rational(exponent, power - r)
    return value


def main() -> None:
    # Quartic edge factors in the two GGV-allowed root patterns.  Both are
    # normalized by h(0)=1 and have h'(0)=0, so h(tau^5)^3=1+O(tau^10).
    h_22 = [Fraction(1), Fraction(0), Fraction(-2), Fraction(0), Fraction(1)]
    h_211 = [
        Fraction(1),
        Fraction(0),
        Fraction(-17, 4),
        Fraction(9, 2),
        Fraction(-5, 4),
    ]
    factor_1 = [Fraction(-1), Fraction(1)]  # t - 1
    factor_minus_1 = [Fraction(1), Fraction(1)]  # t + 1
    factor_2 = [Fraction(-2), Fraction(1)]  # t - 2
    factor_minus_2_over_5 = [Fraction(2, 5), Fraction(1)]  # t + 2/5
    assert poly_pow(factor_1, 2) == [Fraction(1), Fraction(-2), Fraction(1)]
    assert poly_mul(poly_pow(factor_1, 2), poly_pow(factor_minus_1, 2)) == h_22
    h_211_from_roots = scalar_mul(
        poly_mul(
            poly_mul(poly_pow(factor_1, 2), factor_2),
            factor_minus_2_over_5,
        ),
        Fraction(-5, 4),
    )
    assert h_211_from_roots == h_211
    assert h_22[0] == h_211[0] == 1
    assert h_22[1] == h_211[1] == 0
    assert len(h_22) == len(h_211) == 5

    # Put tau^5=t and w=tau*v, and write xi=v^5.  For
    # F=h(t)^3+t^-1*w^5, its tau-face through order 5 is 1+xi for either h.
    for h in (h_22, h_211):
        h_cubed = poly_pow(h, 3)
        assert h_cubed[0] == 1
        assert h_cubed[1] == 0

    # Set the removable target-shear modes kappa_2=kappa_5 to zero, and set
    # kappa_1=kappa_3=kappa_4=0.  The two over-cap coefficients determine
    # kappa_6 and kappa_7 exactly.
    b6 = binomial_rational(Fraction(5, 3), 6)
    b7 = binomial_rational(Fraction(5, 3), 7)
    kappa6 = -b6
    kappa7 = -b7 - kappa6 * binomial_rational(Fraction(-1, 3), 1)
    assert b6 == Fraction(35, 6561)
    assert b7 == Fraction(-65, 19683)
    assert kappa6 == Fraction(-35, 6561)
    assert kappa7 == Fraction(10, 6561)

    kappa = {
        0: Fraction(1),
        1: Fraction(0),
        3: Fraction(0),
        4: Fraction(0),
        6: kappa6,
        7: kappa7,
    }
    assert homogeneous_face_coefficient(6, kappa) == 0
    assert homogeneous_face_coefficient(7, kappa) == 0

    # The seven forbidden Q pole terms in this block are:
    #   s=30..34: depth 6 (tau orders 0..4);
    #   s=35: depths 7 and 6 (tau orders 0 and 5).
    # The tau-face has no tau^1,...,tau^5 terms, and its xi^6,xi^7
    # coefficients have just been checked to vanish.
    forbidden = [
        {"distance": s, "pole_depth": 6, "tau_order": s - 30}
        for s in range(30, 35)
    ] + [
        {"distance": 35, "pole_depth": 7, "tau_order": 0},
        {"distance": 35, "pole_depth": 6, "tau_order": 5},
    ]
    assert len(forbidden) == 7

    print("Exact t=0 support-face cancellation at distances 30--35: PASS")
    print(f"kappa_6 = {kappa6}; kappa_7 = {kappa7}")
    print("Both the 2+2 and 2+1+1 quartic root patterns admit this face witness.")
    print("Scope: no conclusion about root polynomiality, band 36, or a Keller pair.")


if __name__ == "__main__":
    main()
