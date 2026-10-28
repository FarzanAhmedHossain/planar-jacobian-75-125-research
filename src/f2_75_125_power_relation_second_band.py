"""Exact rank audit of the next [P,W*] edge band for the (75,125) F2 model.

After the weight-50 equation for W* = Q^3 - P^5 corrected by lower
P^i Q^j terms, this checks weight 49.  It relaxes all still-unknown
coefficient data to the support-allowed W*_38 polynomial and the first
adjacent P/Q bands, then computes the exact linear rank over the relevant
number fields.  It is conditional on the F2 support and edge hypotheses and
is not a pair exclusion.
"""

from __future__ import annotations

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import QQ


def rank_over_field(matrix: sp.Matrix, field) -> int:
    dm = DomainMatrix.from_Matrix(matrix).convert_to(field)
    return dm.rank()


def build_matrix(h: sp.Expr, v: sp.Expr) -> sp.Matrix:
    t = sp.symbols("t")
    s = sp.symbols("s0:4")
    f_coeffs = sp.symbols("f8:68")
    unknowns = (*s, *f_coeffs)

    S = sum(s[k] * t**k for k in range(4))
    if sp.degree(h, t) != 4:
        raise AssertionError("h must be quartic")
    if sp.degree(v, t) != 4:
        raise AssertionError("the normalized V must be quartic")

    # K5/K15 root-contact conditions give the degree-eight factor of U1;
    # identify it as h^2 divided by its scalar and the forced double root.
    # The caller passes h monic with h=(t-1)^2*r(t).
    r = sp.cancel(h / (t - 1) ** 2)
    if sp.denom(r) != 1:
        raise AssertionError("h must contain the forced double root at t=1")
    u1_factor = sp.expand((t - 1) ** 4 * r**2)
    u1 = sp.expand(u1_factor * S)

    U = sp.expand(t * h)
    A = sp.expand(t**3 * u1)  # P_14 coefficient
    q24 = sp.expand(sp.Rational(5, 3) * t**5 * h**2 * u1)
    g = sp.expand(U**7 * v)  # W*_39 after normalizing the nonzero Jacobian scalar

    f = sum(f_coeffs[k - 8] * t**k for k in range(8, 68))
    # [P_15,W*_38] + [P_14,W*_39] = 3(Q^2)_49 = 6 Q_25 Q_24.
    equation = sp.Poly(
        sp.expand(
            15 * U**3 * sp.diff(f, t)
            - 114 * U**2 * sp.diff(U, t) * f
            + 14 * A * sp.diff(g, t)
            - 39 * sp.diff(A, t) * g
            - 6 * U**5 * q24
        ),
        t,
    )

    rows = []
    for degree in range(equation.degree() + 1):
        coefficient = sp.expand(equation.nth(degree))
        if not coefficient:
            continue
        if sp.expand(coefficient.subs({x: 0 for x in unknowns})) != 0:
            raise AssertionError("the relaxed band equation must be homogeneous")
        rows.append([sp.diff(coefficient, x) for x in unknowns])
    return sp.Matrix(rows)


def main() -> None:
    t = sp.symbols("t")
    s5 = sp.sqrt(5)
    golden_roots = ((3 + s5) / 2, (3 - s5) / 2)
    outputs = []

    for alpha in golden_roots:
        h22 = sp.expand((t - 1) ** 2 * (t - alpha) ** 2)
        slope = (4 * alpha - 11) / 3
        intercept = 3 - alpha
        v22 = sp.expand(t * (t - 1) * (t - alpha) * (slope * t + intercept))
        assert sp.simplify(5 * (t * h22) * sp.diff(v22, t)
                           - 4 * sp.diff(t * h22, t) * v22 - t * h22) == 0
        matrix = build_matrix(h22, v22)
        field = QQ.algebraic_field(s5)
        rank = rank_over_field(matrix, field)
        outputs.append(("2+2", str(alpha), matrix.cols, rank))

    r211 = t**2 - 3 * t + 3
    h211 = sp.expand((t - 1) ** 2 * r211)
    v211 = sp.expand(-sp.Rational(1, 3) * t * (t - 1) * r211)
    assert sp.expand(5 * (t * h211) * sp.diff(v211, t)
                     - 4 * sp.diff(t * h211, t) * v211 - t * h211) == 0
    matrix211 = build_matrix(h211, v211)
    rank211 = rank_over_field(matrix211, QQ)
    outputs.append(("2+1+1", "rational h", matrix211.cols, rank211))

    expected = [
        ("2+2", str(golden_roots[0]), 64, 61),
        ("2+2", str(golden_roots[1]), 64, 61),
        ("2+1+1", "rational h", 64, 62),
    ]
    if outputs != expected:
        raise AssertionError(f"unexpected exact rank data: {outputs}")

    for branch, parameter, variables, rank in outputs:
        print(f"branch={branch} parameter={parameter} variables={variables} rank={rank}")
        print(f"kernel_dimension={variables - rank}")
    print("interpretation=rank is of a relaxed homogeneous necessary system")
    print("scope=conditional on F2 edge support; not a pair exclusion")
    print("rank_audit=PASS")


if __name__ == "__main__":
    main()
