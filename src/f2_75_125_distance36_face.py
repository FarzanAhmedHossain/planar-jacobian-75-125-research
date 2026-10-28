"""Exact first tau-face of the inhomogeneous distance-36 equation.

The checker replays a compact polynomial identity obtained after setting
t=tau^5, w=tau*v in the ratio PDE.  It validates a compatible face jet; it
does not solve or exclude the full F2 coefficient system.
"""

from __future__ import annotations

from fractions import Fraction


Poly = list[Fraction]  # coefficients in ascending powers of z


def add(left: Poly, right: Poly) -> Poly:
    out = [Fraction(0)] * max(len(left), len(right))
    for i, value in enumerate(left):
        out[i] += value
    for i, value in enumerate(right):
        out[i] += value
    return trim(out)


def scale(poly: Poly, scalar: Fraction) -> Poly:
    return trim([scalar * value for value in poly])


def mul(left: Poly, right: Poly) -> Poly:
    out = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return trim(out)


def derivative(poly: Poly) -> Poly:
    if len(poly) <= 1:
        return [Fraction(0)]
    return trim([i * poly[i] for i in range(1, len(poly))])


def z_times(poly: Poly) -> Poly:
    return [Fraction(0), *poly]


def trim(poly: Poly) -> Poly:
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def main() -> None:
    # In the ramified face, z=v^5 and write
    # Fhat=A(z)+tau*v*B(z)+..., Ghat=G0(z)+tau*v*C(z)+...
    # The support bounds give deg(A)<=3, deg(B)<=3, deg(G0)<=5,
    # deg(C)<=5.  Define D=15A-5zA' and E=zG0'-5G0.
    A = [Fraction(1)]
    B = [Fraction(0), Fraction(0), Fraction(0), Fraction(1)]  # z^3
    C = [Fraction(0), Fraction(0), Fraction(0), Fraction(5, 3)]  # 5z^3/3
    G0 = [Fraction(1), Fraction(0), Fraction(0), Fraction(0), Fraction(-1)]

    D = add(scale(A, Fraction(15)), scale(z_times(derivative(A)), Fraction(-5)))
    E = add(z_times(derivative(G0)), scale(G0, Fraction(-5)))
    lhs = add(mul(D, C), scale(mul(B, E), Fraction(5)))
    rhs = [Fraction(0)] * 8
    rhs[7] = Fraction(5)
    assert lhs == rhs
    assert D == [Fraction(15)]
    assert E == [Fraction(-5), Fraction(0), Fraction(0), Fraction(0), Fraction(1)]

    # This tau^0 face comes from the homogeneous first integral with
    # kappa_0=1, kappa_4=-1 and all other low-face kappa_r=0:
    # G0=1-z^4.  Its tau^1 homogeneous prediction is
    # B*(5/3-z^4/3)=5z^3/3-z^7/3.  The z^7 term is at distance 36, where
    # the inhomogeneous source replaces the homogeneous recurrence.  The
    # actual support-compatible C has degree 3 and agrees through distance 35.
    homogeneous_tau1 = [Fraction(0)] * 8
    homogeneous_tau1[3] = Fraction(5, 3)
    homogeneous_tau1[7] = Fraction(-1, 3)
    actual_tau1 = [Fraction(0)] * 8
    actual_tau1[3] = Fraction(5, 3)
    assert homogeneous_tau1[:7] == actual_tau1[:7]
    assert actual_tau1[7] == 0

    # Translate the same face data back to the restrictions of P and Q on
    # y=0: p0=x^3, p1=1, q0=x^5-x, q1=5*x^2/3.  The source-face identity is
    # exactly p0' q1-p1 q0'=1.
    p0 = [Fraction(0), Fraction(0), Fraction(0), Fraction(1)]
    p1 = [Fraction(1)]
    q0 = [Fraction(0), Fraction(-1), Fraction(0), Fraction(0), Fraction(0), Fraction(1)]
    q1 = [Fraction(0), Fraction(0), Fraction(5, 3)]
    axis_jacobian = add(mul(derivative(p0), q1), scale(mul(p1, derivative(q0)), Fraction(-1)))
    assert axis_jacobian == [Fraction(1)]

    # The coefficient of z^7 in the general identity is exactly the
    # determinant of the linear part of (P,Q), normalized to 1:
    # a_2*c_5 - b_3*g_4 = 1.
    # Here a_2=0, c_5=0, b_3=1, g_4=-1.
    assert Fraction(0) * Fraction(0) - Fraction(1) * Fraction(-1) == 1

    print("Exact distance-36 first-face identity: PASS")
    print("D*C + 5*B*E = 5*z^7 for A=1, B=z^3, G0=1-z^4, C=5*z^3/3.")
    print("Equivalent axis restriction: p0'*q1-p1*q0'=1.")
    print("The top coefficient is the normalized linear Jacobian determinant.")
    print("Scope: this is a compatible face jet, not a full Keller pair or exclusion.")


if __name__ == "__main__":
    main()
