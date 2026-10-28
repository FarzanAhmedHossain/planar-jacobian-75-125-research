"""Compute the parameter-independent matrix in the 46/56 band pair.

This extracts only the coefficients of the new W*_35 and U_4 slices from
the weight-46 [P,W*] and weight-56 [Q,W*] equations. The forcing terms from
the earlier affine blocks are intentionally omitted; this file reports rank
and fibre dimension, not compatibility or exclusion.
"""

from __future__ import annotations

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import QQ


def branch_matrix(tag: str, h: sp.Expr, V: sp.Expr, field) -> None:
    t = sp.symbols("t")
    U = sp.expand(t * h)
    g = sp.expand(U**7 * V)
    W_coeffs = sp.symbols("W35_7:68")
    U4_coeffs = sp.symbols("U4_0:12")
    W35 = sum(W_coeffs[i - 7] * t**i for i in range(7, 68))
    U4 = sum(U4_coeffs[i] * t**i for i in range(12))
    A4 = sp.expand(t**3 * U4)
    q21_linear = sp.expand(sp.Rational(5, 3) * t**5 * h**2 * U4)

    p46_linear = (
        15 * U**3 * sp.diff(W35, t)
        - 105 * U**2 * sp.diff(U, t) * W35
        + 11 * A4 * sp.diff(g, t)
        - 39 * sp.diff(A4, t) * g
        - 6 * U**5 * q21_linear
    )
    q56_linear = (
        5 * U**4 * (5 * U * sp.diff(W35, t) - 35 * sp.diff(U, t) * W35)
        + 21 * q21_linear * sp.diff(g, t)
        - 39 * sp.diff(q21_linear, t) * g
        - 20 * U**9 * A4
    )

    rows = []
    row_counts = {}
    for expression, label in ((p46_linear, "P46"), (q56_linear, "Q56")):
        polynomial = sp.Poly(sp.expand(expression), t)
        row_start = len(rows)
        for degree in range(polynomial.degree() + 1):
            coefficient = sp.expand(polynomial.nth(degree))
            if coefficient == 0:
                continue
            row = [sp.expand(sp.diff(coefficient, var))
                   for var in (*W_coeffs, *U4_coeffs)]
            if any(sp.sympify(entry).free_symbols.intersection((*W_coeffs, *U4_coeffs))
                   for entry in row):
                raise AssertionError(f"{label} matrix row is not linear")
            if sp.expand(coefficient - sum(row[i] * (*W_coeffs, *U4_coeffs)[i]
                                           for i in range(len(row)))) != 0:
                raise AssertionError(f"{label} contains a nonhomogeneous matrix term")
            rows.append(row)
        row_counts[label] = len(rows) - row_start

    matrix = sp.Matrix(rows)
    rank = DomainMatrix.from_Matrix(matrix).convert_to(field).rank()
    print(
        f"branch={tag} rows={matrix.rows} variables={matrix.cols} rank={rank} "
        f"left_cokernel_dimension={matrix.rows-rank} fibre_dimension={matrix.cols-rank} "
        f"P46_rows={row_counts['P46']} Q56_rows={row_counts['Q56']}",
        flush=True,
    )


def main() -> None:
    t = sp.symbols("t")
    sqrt5 = sp.sqrt(5)
    alpha = (3 + sqrt5) / 2
    h22 = sp.expand((t - 1) ** 2 * (t - alpha) ** 2)
    m, n = (4 * alpha - 11) / 3, 3 - alpha
    v22 = sp.expand(t * (t - 1) * (t - alpha) * (m * t + n))
    branch_matrix("2+2", h22, v22, QQ.algebraic_field(sqrt5))

    r = t**2 - 3 * t + 3
    h211 = sp.expand((t - 1) ** 2 * r)
    v211 = sp.expand(-sp.Rational(1, 3) * t * (t - 1) * r)
    branch_matrix("2+1+1", h211, v211, QQ)


if __name__ == "__main__":
    main()
