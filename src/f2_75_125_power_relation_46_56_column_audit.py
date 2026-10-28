"""Cross-check the direct 46/56 column assembler against full expressions."""

from __future__ import annotations

import sympy as sp

from f2_75_125_power_relation_fifthbands import (
    add,
    coefficient_rows,
    deriv,
    direct_46_56_system,
    mul,
    pow_poly,
    scale,
    shift,
    trim,
    verify_new_band_formulas,
)


def assert_same_rows(actual, expected, tag: str) -> None:
    actual_matrix, actual_rhs = actual
    expected_matrix, expected_rhs = expected
    if actual_matrix != expected_matrix:
        for i in range(max(actual_matrix.rows, expected_matrix.rows)):
            left = actual_matrix.row(i) if i < actual_matrix.rows else None
            right = expected_matrix.row(i) if i < expected_matrix.rows else None
            if left != right:
                raise AssertionError(f"{tag} matrix row {i}: {left} != {right}")
        raise AssertionError(f"{tag} matrix shapes differ")
    if len(actual_rhs) != len(expected_rhs):
        raise AssertionError(f"{tag} rhs lengths differ")
    for i, (left, right) in enumerate(zip(actual_rhs, expected_rhs)):
        if sp.expand(left - right) != 0:
            raise AssertionError(f"{tag} rhs row {i}: {left} != {right}")


def main() -> None:
    verify_new_band_formulas()
    U = [sp.Integer(1), sp.Integer(2), sp.Integer(1)]
    h = [sp.Integer(1), sp.Integer(2)]
    S = [sp.Integer(3), sp.Integer(1)]
    B = [sp.Integer(2)]
    U3 = [sp.Integer(1), sp.Integer(0), sp.Integer(1)]
    A = [sp.Integer(1), sp.Integer(2)]
    A2 = [sp.Integer(2), sp.Integer(1)]
    A3 = [sp.Integer(2), sp.Integer(0), sp.Integer(1)]
    q24 = [sp.Integer(1), sp.Integer(3)]
    q23 = [sp.Integer(2), sp.Integer(1)]
    q22 = [sp.Integer(-1), sp.Integer(2)]
    g = [sp.Integer(2), sp.Integer(1), sp.Integer(0), sp.Integer(1)]
    f38 = [sp.Integer(1), sp.Integer(1)]
    f37 = [sp.Integer(2), sp.Integer(3)]
    ell = [sp.Integer(1), sp.Integer(0), sp.Integer(2)]
    W35_coeffs = sp.symbols("w0:61")
    U4_coeffs = sp.symbols("u0:12")
    variables = (*W35_coeffs, *U4_coeffs)

    direct_matrix, direct_rhs, direct_variables, direct_p_count = direct_46_56_system(
        U, h, S, B, U3, A, A2, A3, q24, q23, q22, g, f38, f37, ell,
        W35_coeffs, U4_coeffs,
    )
    if direct_variables != variables:
        raise AssertionError("direct assembler changed the new-variable order")

    W35 = [sp.Integer(0)] * 7 + list(W35_coeffs)
    U4 = list(U4_coeffs)
    q21_0 = scale(
        sp.Rational(5, 243),
        shift(
            mul(
                h,
                add(
                    scale(54, mul(S, U3)),
                    scale(27, pow_poly(B, 2)),
                    scale(-9, mul(pow_poly(S, 2), B)),
                    pow_poly(S, 4),
                ),
            ),
            5,
        ),
    )
    q21 = add(
        scale(sp.Rational(5, 3), shift(mul(pow_poly(h, 2), U4), 5)),
        q21_0,
    )
    A4 = shift(U4, 3)

    p46 = add(
        scale(15, mul(pow_poly(U, 3), deriv(W35))),
        scale(-105, mul(pow_poly(U, 2), mul(deriv(U), W35))),
        scale(14, mul(A, deriv(ell))),
        scale(-36, mul(deriv(A), ell)),
        scale(13, mul(A2, deriv(f37))),
        scale(-37, mul(deriv(A2), f37)),
        scale(12, mul(A3, deriv(f38))),
        scale(-38, mul(deriv(A3), f38)),
        scale(11, mul(A4, deriv(g))),
        scale(-39, mul(deriv(A4), g)),
        scale(-6, mul(pow_poly(U, 5), q21)),
        scale(-6, mul(q24, q22)),
        scale(-3, pow_poly(q23, 2)),
    )
    q56 = add(
        scale(25, mul(pow_poly(U, 5), deriv(W35))),
        scale(-175, mul(pow_poly(U, 4), mul(deriv(U), W35))),
        scale(24, mul(q24, deriv(ell))),
        scale(-36, mul(deriv(q24), ell)),
        scale(23, mul(q23, deriv(f37))),
        scale(-37, mul(deriv(q23), f37)),
        scale(22, mul(q22, deriv(f38))),
        scale(-38, mul(deriv(q22), f38)),
        scale(21, mul(q21, deriv(g))),
        scale(-39, mul(deriv(q21), g)),
        scale(-20, mul(pow_poly(U, 9), A4)),
        scale(-60, mul(pow_poly(U, 6), mul(A, A3))),
        scale(-30, mul(pow_poly(U, 6), pow_poly(A2, 2))),
        scale(-60, mul(pow_poly(U, 3), mul(pow_poly(A, 2), A2))),
        scale(-5, pow_poly(A, 4)),
    )

    p_rows, p_rhs = coefficient_rows(p46, variables)
    q_rows, q_rhs = coefficient_rows(q56, variables)
    expected_matrix = sp.Matrix([*p_rows, *q_rows])
    expected_rhs = [*p_rhs, *q_rhs]
    print(
        f"row_split direct_p={direct_p_count} expected_p={len(p_rows)} "
        f"direct_q={direct_matrix.rows-direct_p_count} expected_q={len(q_rows)}",
        flush=True,
    )
    assert_same_rows(
        (direct_matrix[:direct_p_count, :], direct_rhs[:direct_p_count]),
        (sp.Matrix(p_rows), p_rhs),
        "46",
    )
    assert_same_rows(
        (direct_matrix[direct_p_count:, :], direct_rhs[direct_p_count:]),
        (sp.Matrix(q_rows), q_rhs),
        "56",
    )
    print(
        f"direct_46_56_reconstruction=PASS rows={direct_matrix.rows} "
        f"variables={direct_matrix.cols}",
        flush=True,
    )


if __name__ == "__main__":
    main()
