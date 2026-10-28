"""Stack the first non-top bands of [P,W*] and [Q,W*].

The weight-49 row of [P,W*] and weight-59 row of [Q,W*] share the same
support-bounded W*_38 coefficient.  This exact rank audit tests whether their
combined necessary linear system forces U1=0 on either surviving quartic
carrier configuration.  It is conditional on the F2 support/edge model and
does not by itself exclude (75,125).
"""

from __future__ import annotations

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import QQ

from f2_75_125_power_relation_second_band import build_matrix, rank_over_field


def q_side_matrix(h: sp.Expr, v: sp.Expr) -> sp.Matrix:
    t = sp.symbols("t")
    s = sp.symbols("s0:4")
    f_coeffs = sp.symbols("f8:68")
    unknowns = (*s, *f_coeffs)
    S = sum(s[k] * t**k for k in range(4))
    U = sp.expand(t * h)
    U1 = sp.expand(h**2 * S)
    A = sp.expand(t**3 * U1)
    q24 = sp.expand(sp.Rational(5, 3) * t**5 * h**2 * U1)
    g = sp.expand(U**7 * v)
    f = sum(f_coeffs[k - 8] * t**k for k in range(8, 68))

    # The weight-59 row of [Q,W*]=5 lambda J P^4, after dividing
    # by the common nonzero scalar J*b^3/a.
    equation = sp.Poly(
        sp.expand(
            5 * U**4 * (5 * U * sp.diff(f, t) - 38 * sp.diff(U, t) * f)
            + 24 * q24 * sp.diff(g, t)
            - 39 * sp.diff(q24, t) * g
            - 20 * U**9 * A
        ),
        t,
    )
    rows = []
    for degree in range(equation.degree() + 1):
        coefficient = sp.expand(equation.nth(degree))
        if not coefficient:
            continue
        if sp.expand(coefficient.subs({x: 0 for x in unknowns})) != 0:
            raise AssertionError("the Q-side band must be homogeneous")
        rows.append([sp.diff(coefficient, x) for x in unknowns])
    return sp.Matrix(rows)


def audit(tag: str, h: sp.Expr, v: sp.Expr, field) -> tuple[int, int, int, int]:
    M49 = build_matrix(h, v)
    M59 = q_side_matrix(h, v)
    if M49.cols != 64 or M59.cols != 64:
        raise AssertionError("both edge rows must use the same 64 coefficients")
    r49 = rank_over_field(M49, field)
    r59 = rank_over_field(M59, field)
    stacked = sp.Matrix.vstack(M49, M59)
    rstack = rank_over_field(stacked, field)
    data = (M49.rows, M59.rows, r59, rstack)
    print(f"branch={tag} P_row_rank={r49} Q_row_rank={r59} stacked_rank={rstack}")
    print(f"P_rows={M49.rows} Q_rows={M59.rows} stacked_kernel_dimension={64-rstack}")
    return data


def main() -> None:
    t = sp.symbols("t")
    sqrt5 = sp.sqrt(5)
    alpha = (3 + sqrt5) / 2
    h22 = sp.expand((t - 1) ** 2 * (t - alpha) ** 2)
    m, n = (4 * alpha - 11) / 3, 3 - alpha
    v22 = sp.expand(t * (t - 1) * (t - alpha) * (m * t + n))
    assert sp.simplify(5 * (t*h22) * sp.diff(v22,t)
                       - 4 * sp.diff(t*h22,t) * v22 - t*h22) == 0
    data22 = audit("2+2", h22, v22, QQ.algebraic_field(sqrt5))

    r = t**2 - 3*t + 3
    h211 = sp.expand((t - 1) ** 2 * r)
    v211 = sp.expand(-sp.Rational(1, 3) * t * (t - 1) * r)
    assert sp.expand(5 * (t*h211) * sp.diff(v211,t)
                     - 4 * sp.diff(t*h211,t) * v211 - t*h211) == 0
    data211 = audit("2+1+1", h211, v211, QQ)

    print(f"2+2_data={data22}")
    print(f"2+1+1_data={data211}")
    print("scope=two exact necessary rows; not a pair exclusion")


if __name__ == "__main__":
    main()
