"""Replay the local double-root jet elimination for the F2 (75,125) model.

Scope: verify the exact successor-corner arithmetic, the polynomial identity
that excludes ord_alpha(h)=2, ord_alpha(U_1)=3, and the leading-coefficient
relation in the q=1, r=5 branch. This checker does not prove the band
equations or the GGV reduction; those are stated and derived in the
accompanying F2-75-125-band-model.md note.

Polynomials are tuples of integer coefficients in ascending degree order.
No CAS or floating-point arithmetic is used.
"""

from __future__ import annotations

from fractions import Fraction


def add(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    size = max(len(left), len(right))
    result = [0] * size
    for i, coefficient in enumerate(left):
        result[i] += coefficient
    for i, coefficient in enumerate(right):
        result[i] += coefficient
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result)


def multiply(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result)


def main() -> None:
    # GGV F2 successor data: m=3, l=1, k=3, (rho,sigma)=(5,-1),
    # A1=(7/5,2). The type-II.b formula gives A1=(k/(m*l),0)
    # + (m_lambda/m)*(-sigma/rho,1), so m_lambda=6.
    m = 3
    l = 1
    k = 3
    rho = 5
    sigma = -1
    a1 = (Fraction(7, 5), Fraction(2, 1))
    multiplicity_ratio = a1[1]
    assert a1[0] == Fraction(k, m * l) + multiplicity_ratio * Fraction(-sigma, rho)
    edge_root_multiplicity = m * multiplicity_ratio
    assert edge_root_multiplicity == 6
    # The edge polynomial is h(z^5)^3. At z=lambda != 0, z^5 is
    # unramified, so the selected root of h has multiplicity 6/3 = 2.
    assert edge_root_multiplicity / m == 2

    # In the q=1, r=5 double-root branch, order-19 cancellation in N5 gives
    # E = A*D/(6*c^3); substituting this cancels 810*c^9*D*E against
    # -135*c^6*A*D^2 exactly.
    assert 810 * Fraction(1, 6) - 135 == 0

    # f(x) = 27*x^2 - 9*x + 1 and g(x) = -135*x^2 + 60*x - 7.
    f = (1, -9, 27)
    g = (-7, 60, -135)
    bezout_f = (70, -225)
    bezout_g = (9, -45)

    certificate = add(multiply(bezout_f, f), multiply(bezout_g, g))
    assert certificate == (7,), certificate
    print("GGV successor arithmetic verified: selected h-root multiplicity = 2")
    print("exact certificate verified: Bezout_f*f + Bezout_g*g = 7")


if __name__ == "__main__":
    main()
