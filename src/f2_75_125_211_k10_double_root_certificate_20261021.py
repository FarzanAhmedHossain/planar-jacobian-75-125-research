#!/usr/bin/env python3
"""Local t=1 replay of the exact K10 obstruction for one F2 witness.

This independently checks the K10 polynomiality equations in the completed
local ring at the distinguished double root.  After multiplication by
`t*h^25`, regularity requires the coefficients of u^0,...,u^49 to vanish,
where u=t-1.  The supported U6 slice and kappa_1 are left completely free.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

import sympy as sp

from f2_75_125_211_exact_k9_tied_residue_check_20261021 import (
    coefficient_numerator,
    sha256,
)


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_exact_k9_tied_residue_check_20261021.json"
OUT = HERE / "f2_75_125_211_k10_double_root_certificate_20261021.json"


def coefficients_at_t1(expr: sp.Expr, t: sp.Symbol, u: sp.Symbol,
                       degree: int, field_generator: sp.Expr) -> list[sp.Expr]:
    shifted = sp.Poly(sp.expand(expr.subs(t, u + 1)), u,
                      extension=field_generator)
    return [sp.simplify(shifted.nth(i)) for i in range(degree)]


def main() -> None:
    t, u = sp.symbols("t u")
    s = sp.sqrt(-3)
    h = (t - 1) ** 2 * (t**2 - 3 * t + 3)
    R = t**2 - 3 * t + 3
    U1 = (t - 1) ** 5 * R**3
    U2 = (t - 1) ** 3 * R**2
    U3 = (t - 1) ** 2 * R**2
    d0 = sp.Rational(1, 6) - s / 18
    d1 = sp.Rational(1, 9) + 8 * s / 27
    U4 = R * (d0 + d1 * (t - 1))
    F5 = R  # U5/t, the actual w^5 coefficient of F.
    u6_names = sp.symbols("u6_0:12")
    U6 = sum(value * t**i for i, value in enumerate(u6_names))
    slices = [U1, U2, U3, U4, F5, U6 / t]

    # After clearing by t*h^25, K10 has the polynomial numerator
    # t*h^25*[w^10]F^(5/3) + kappa_1*h^25*[w^5]F^(4/3).
    main = coefficient_numerator(10, Fraction(5, 3), slices, h, t,
                                 clear_power=25, clear_t=True)
    correction = coefficient_numerator(
        5, Fraction(4, 3), [U1, U2, U3, U4, F5], h, t,
        clear_power=25, clear_t=False
    )
    base = sp.expand(main.subs({value: 0 for value in u6_names}))
    columns = [sp.expand(sp.diff(main, value)) for value in u6_names]
    columns.append(correction)

    # h has order 2 at t=1, so h^25 has order 50. Regularity at this root
    # forces exactly these first 50 Taylor coefficients of the cleared
    # numerator to vanish.  The omitted K10 terms are regular and hence
    # contribute order at least 50 after the same clearing.
    row_count = 50
    base_coeffs = coefficients_at_t1(base, t, u, row_count, s)
    column_coeffs = [coefficients_at_t1(column, t, u, row_count, s)
                     for column in columns]
    matrix = sp.Matrix(row_count, len(columns),
                       lambda row, col: column_coeffs[col][row])
    rhs = sp.Matrix([-value for value in base_coeffs])
    augmented = matrix.row_join(rhs)

    direct_rows = [i for i in range(row_count)
                   if all(matrix[i, j] == 0 for j in range(matrix.cols))
                   and rhs[i] != 0]
    rank = matrix.rank()
    augmented_rank = augmented.rank()
    if augmented_rank <= rank:
        raise AssertionError("the t=1 local system unexpectedly has a solution")

    witness = None
    for vector in matrix.T.nullspace():
        pairing = sp.simplify((vector.T * rhs)[0])
        if pairing:
            column_pairings = [sp.simplify(value) for value in vector.T * matrix]
            if any(column_pairings):
                raise AssertionError("local left-functional fails to annihilate a column")
            witness = {
                "nonzero_taylor_row_weights": [
                    [i, str(value)] for i, value in enumerate(vector) if value != 0
                ],
                "pairings_with_all_U6_and_kappa1_columns": [str(x) for x in column_pairings],
                "pairing_with_forcing": str(sp.simplify((vector.T * sp.Matrix(base_coeffs))[0])),
                "pairing_with_rhs": str(pairing),
            }
            break
    if witness is None:
        raise AssertionError("no exact local left-kernel certificate found")

    result = {
        "source_receipt": SOURCE.name,
        "source_receipt_sha256": sha256(SOURCE),
        "field": "Q(sqrt(-3))",
        "root_chart": "u=t-1; h=u^2*(t^2-3t+3)",
        "cleared_K10_numerator": "t*h^25*[w^10]F^(5/3) + kappa_1*h^25*[w^5]F^(4/3)",
        "regularity_equations": row_count,
        "unknown_columns": [*[str(x) for x in u6_names], "kappa_1"],
        "support_cap": "U6 is an arbitrary polynomial of degree <=11; F6=U6/t",
        "coefficient_matrix_rank": rank,
        "augmented_matrix_rank": augmented_rank,
        "rows_with_zero_unknown_coefficients_and_nonzero_forcing": direct_rows,
        "left_functional_inconsistency_certificate": witness,
        "conclusion": "the explicit K3-K8 witness has no supported K10 lift, already at t=1",
        "scope": "one explicit local coefficient assignment and its conjugate only; not a branch exclusion",
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "EXACT_T1_K10_LOCAL_NO_LIFT",
        "coefficient_matrix_rank": rank,
        "augmented_matrix_rank": augmented_rank,
        "direct_inconsistent_rows": direct_rows,
        "certificate": witness,
        "receipt": str(OUT),
        "script_sha256": sha256(Path(__file__)),
        "receipt_sha256": sha256(OUT),
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
