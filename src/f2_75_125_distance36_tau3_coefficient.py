"""Exact next ramified source-face coefficient for the F2 (75,125) slice.

This continues the U1=0-compatible t=0 jet with the small witness a=b=0
from f2_75_125_distance36_tau3_face.py. The tau^4 coefficient bounds are
again read from the exact j=4 slices of the P and Q support envelopes.
Compatibility here would remain local t=0 evidence only.
"""

from __future__ import annotations

import sympy as sp


def tau4_z_degrees(x_max: int, y_max: int, edge_cap: int, x_factor: int) -> list[int]:
    """Project supported monomials at y-degree 4 into the ramified chart."""
    powers: list[int] = []
    for i in range(x_max + 1):
        for j in range(y_max + 1):
            if 5 * i - j > edge_cap or j != 4:
                continue
            t_exp = i - x_factor
            w_exp = j - 5 * t_exp
            assert w_exp >= 0
            assert 5 * t_exp + w_exp == 4
            assert t_exp <= 0, (i, j, t_exp)
            z_degree = -t_exp
            assert w_exp == 4 + 5 * z_degree
            powers.append(z_degree)
    return sorted(powers)


def main() -> None:
    assert tau4_z_degrees(15, 60, 15, 3) == list(range(4))
    assert tau4_z_degrees(25, 100, 25, 5) == list(range(6))

    tau, v = sp.symbols("tau v")
    z = v**5
    f4c = sp.symbols("f40:44")
    g4c = sp.symbols("g40:46")
    A = sp.Integer(1)
    B = z**3
    F2 = sp.Rational(5, 3) * z**2
    F3 = sp.Rational(50, 9) * z
    F4 = sum(f4c[i] * z**i for i in range(4))
    G0 = 1 - z**4
    C = sp.Rational(5, 3) * z**3
    G2 = sp.Rational(25, 9) * z**2
    G3 = sp.Rational(250, 27) * z
    G4 = sum(g4c[i] * z**i for i in range(6))

    F = A + tau * v * B + tau**2 * v**2 * F2 + tau**3 * v**3 * F3 + tau**4 * v**4 * F4
    G = G0 + tau * v * C + tau**2 * v**2 * G2 + tau**3 * v**3 * G3 + tau**4 * v**4 * G4
    F_minus_5_3 = sp.series(F ** sp.Rational(-5, 3), tau, 0, 5).removeO()
    F_minus_8_3 = sp.series(F ** sp.Rational(-8, 3), tau, 0, 4).removeO()
    R = sp.series(G * F_minus_5_3, tau, 0, 5).removeO()
    pde = (
        (15 - v * sp.diff(F, v) / F) * sp.diff(R, tau)
        + v * sp.diff(F, tau) / F * sp.diff(R, v)
        - 5 * v**36 * F_minus_8_3
    )
    pde = sp.series(pde, tau, 0, 4).removeO().expand()
    for degree in range(3):
        assert sp.expand(pde.coeff(tau, degree)) == 0

    face = sp.Poly(sp.expand(pde.coeff(tau, 3)), v)
    powers = [monomial[0] for monomial, _ in face.terms()]
    assert all(power >= 4 and (power - 4) % 5 == 0 for power in powers)
    max_z_degree = max(((power - 4) // 5 for power in powers), default=0)
    equations = [
        sp.expand(face.coeff_monomial(v ** (4 + 5 * degree)))
        for degree in range(max_z_degree + 1)
    ]
    equations = [equation for equation in equations if equation != 0]
    variables = (*f4c, *g4c)
    matrix, rhs = sp.linear_eq_to_matrix(equations, variables)
    rank = matrix.rank()
    augmented_rank = matrix.row_join(rhs).rank()
    assert rank == augmented_rank, (rank, augmented_rank)
    solution_set = sp.linsolve((matrix, rhs), variables)
    solution = next(iter(solution_set))
    assert all(sp.expand(e.subs(dict(zip(variables, solution)))) == 0 for e in equations)

    print("Exact distance-36 tau^3 PDE coefficient: compatible for the a=b=0 jet")
    print(f"support-derived tau^4 coefficient caps: deg(F4)<=3, deg(G4)<=5")
    print(f"coefficient equations={len(equations)}, unknown coefficients={len(variables)}")
    print(f"exact rational rank={rank}, augmented rank={augmented_rank}")
    print("one exact solution:")
    for name, value in zip(variables, solution):
        print(f"  {name} = {sp.factor(value)}")
    print("Scope: local jet at t=0 only; this is not a root lift or branch exclusion.")


if __name__ == "__main__":
    main()
