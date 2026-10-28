"""Exact third tau-face of the U1=0 distance-36 source equation.

This replay expands the ratio PDE in the ramified chart used in
F2-75-125-band-model.md. It first derives the tau^3 coefficient caps from
the exact F2 support envelopes, then asks whether the compatible first two
faces extend through the coefficient of tau^2. A solution is local evidence
only, not a lift at the other roots of h or an exclusion/constructive Keller
pair.
"""

from __future__ import annotations

import sympy as sp


def tau3_z_degrees(x_max: int, y_max: int, edge_cap: int, x_factor: int) -> list[int]:
    """Return z-powers in the tau^3 face induced by one exact support envelope.

    For a supported monomial x^i*y^j in x^x_factor*F(t,w), where
    t=x*y^5 and w=y, its F-coordinate is t^(i-x_factor)*w^(j-5*(i-x_factor)).
    Under t=tau^5,w=tau*v its tau-order is exactly j. If the Laurent
    t-exponent is -r, the tau^3 coefficient is tau^3*v^3*z^r.
    """
    powers: list[int] = []
    for i in range(x_max + 1):
        for j in range(y_max + 1):
            if 5 * i - j > edge_cap or j != 3:
                continue
            t_exp = i - x_factor
            w_exp = j - 5 * t_exp
            assert w_exp >= 0
            assert 5 * t_exp + w_exp == 3
            assert t_exp <= 0, (i, j, t_exp)
            z_degree = -t_exp
            assert w_exp == 3 + 5 * z_degree
            powers.append(z_degree)
    return sorted(powers)


def main() -> None:
    # The exact (75,125) F2 outer envelopes have edge caps 15 and 25.
    # Pulling their j=3 slices into the ramified chart proves, rather than
    # assumes, the coefficient degrees used below.
    f3_support_degrees = tau3_z_degrees(15, 60, 15, 3)
    g3_support_degrees = tau3_z_degrees(25, 100, 25, 5)
    assert f3_support_degrees == list(range(4))
    assert g3_support_degrees == list(range(6))

    tau, v = sp.symbols("tau v")
    z = v**5
    f3c = sp.symbols("f30:34")
    g3c = sp.symbols("g30:36")

    # The two already verified source faces.
    A = sp.Integer(1)
    B = z**3
    F2 = sp.Rational(5, 3) * z**2
    G0 = 1 - z**4
    C = sp.Rational(5, 3) * z**3
    G2 = sp.Rational(25, 9) * z**2

    F3 = sum(f3c[i] * z**i for i in range(4))
    G3 = sum(g3c[i] * z**i for i in range(6))
    F = A + tau * v * B + tau**2 * v**2 * F2 + tau**3 * v**3 * F3
    G = G0 + tau * v * C + tau**2 * v**2 * G2 + tau**3 * v**3 * G3

    F_minus_5_3 = sp.series(F ** sp.Rational(-5, 3), tau, 0, 4).removeO()
    F_minus_8_3 = sp.series(F ** sp.Rational(-8, 3), tau, 0, 3).removeO()
    R = sp.series(G * F_minus_5_3, tau, 0, 4).removeO()
    pde = (
        (15 - v * sp.diff(F, v) / F) * sp.diff(R, tau)
        + v * sp.diff(F, tau) / F * sp.diff(R, v)
        - 5 * v**36 * F_minus_8_3
    )
    pde = sp.series(pde, tau, 0, 3).removeO().expand()

    # Coefficients tau^0 and tau^1 are the two existing exact face checks.
    assert sp.expand(pde.coeff(tau, 0)) == 0
    assert sp.expand(pde.coeff(tau, 1)) == 0

    face = sp.Poly(sp.expand(pde.coeff(tau, 2)), v)
    powers = [monomial[0] for monomial, _ in face.terms()]
    assert all(power >= 3 and (power - 3) % 5 == 0 for power in powers)
    max_z_degree = max(((power - 3) // 5 for power in powers), default=0)
    equations = [
        sp.expand(face.coeff_monomial(v ** (3 + 5 * degree)))
        for degree in range(max_z_degree + 1)
    ]
    equations = [equation for equation in equations if equation != 0]
    variables = (*f3c, *g3c)
    matrix, rhs = sp.linear_eq_to_matrix(equations, variables)
    rank = matrix.rank()
    augmented_rank = matrix.row_join(rhs).rank()
    assert (len(equations), len(variables), rank, augmented_rank) == (8, 10, 8, 8)

    # An exact two-parameter family of solutions.  Taking a=b=0 gives the
    # especially small witness F3=(50/9)z, G3=(250/27)z.
    a, b = sp.symbols("a b")
    solution = {
        f3c[0]: -3 * a,
        f3c[1]: sp.Rational(50, 9) - 3 * b,
        f3c[2]: 0,
        f3c[3]: 0,
        g3c[0]: -5 * a,
        g3c[1]: sp.Rational(250, 27) - 5 * b,
        g3c[2]: 0,
        g3c[3]: 0,
        g3c[4]: a,
        g3c[5]: b,
    }
    assert all(sp.expand(equation.subs(solution)) == 0 for equation in equations)

    print("Exact distance-36 tau^2 face: compatible")
    print("support-derived tau^3 coefficient caps: deg(F3)<=3, deg(G3)<=5")
    print(f"coefficient equations={len(equations)}, unknown coefficients={len(variables)}")
    print(f"exact rational rank={rank}, augmented rank={augmented_rank}, free parameters=2")
    print("F3(z) = -3*a + (50/9 - 3*b)*z")
    print("G3(z) = -5*a + (250/27 - 5*b)*z + a*z^4 + b*z^5")
    print("Scope: local jet at t=0 only; other-root polynomiality and the full branch remain open.")


if __name__ == "__main__":
    main()
