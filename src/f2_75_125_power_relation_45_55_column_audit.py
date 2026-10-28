"""Exact symbolic audit of the optimized 45/55 column assembler."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from f2_75_125_power_relation_fifthbands import (  # noqa: E402
    add,
    coefficient_rows,
    deriv,
    direct_45_55_system,
    mul,
    pow_poly,
    scale,
    shift,
)


def coefficients(expression, variable):
    polynomial = sp.Poly(sp.expand(expression), variable)
    return [polynomial.nth(i) for i in range(polynomial.degree() + 1)]


def main():
    t = sp.symbols("t")
    a, b, kappa = sp.symbols("a b kappa")
    data = {
        "U": coefficients(t + t**2, t),
        "h": coefficients(1 + t, t),
        "S": coefficients(1 + a*t, t),
        "B": coefficients(2 + b*t, t),
        "U3": coefficients(1 + a + t, t),
        "U4": coefficients(2 + b*t, t),
        "A": coefficients(1 + a*t, t),
        "A2": coefficients(1 + b*t, t),
        "A3": coefficients(2 + a*t, t),
        "A4": coefficients(1 + b*t + t**2, t),
        "q24": coefficients(1 + a*t, t),
        "q23": coefficients(2 + b*t, t),
        "q22": coefficients(1 + a + b*t, t),
        "q21": coefficients(2 + a*t, t),
        "g": coefficients(1 + t + t**2, t),
        "W38": coefficients(1 + a*t, t),
        "W37": coefficients(2 + b*t, t),
        "W36": coefficients(1 + a + b*t, t),
        "W35": coefficients(2 + a*t, t),
    }
    U, h = data["U"], data["h"]
    S, B, U3, U4 = (data[name] for name in ("S", "B", "U3", "U4"))
    A, A2, A3, A4 = (data[name] for name in ("A", "A2", "A3", "A4"))
    q24, q23, q22, q21 = (data[name] for name in ("q24", "q23", "q22", "q21"))
    g = data["g"]
    W38, W37, W36, W35 = (data[name] for name in ("W38", "W37", "W36", "W35"))

    U5_coeffs = sp.symbols("U5_0:12")
    W34_coeffs = sp.symbols("W34_7:67")
    new_variables = (*U5_coeffs, *W34_coeffs)
    U5 = list(U5_coeffs)
    W34 = [sp.Integer(0)] * 7 + list(W34_coeffs)

    q20_binomial = add(
        scale(sp.Rational(5, 3), mul(pow_poly(h, 2), U5)),
        scale(sp.Rational(10, 9), add(mul(mul(h, S), U4), mul(B, U3))),
        scale(sp.Rational(-5, 27), add(mul(pow_poly(S, 2), U3),
                                        mul(S, pow_poly(B, 2)))),
        scale(sp.Rational(20, 243), mul(pow_poly(S, 3), B)),
        scale(sp.Rational(-7, 729), pow_poly(S, 5)),
    )
    q20 = add(shift(q20_binomial, 5), scale(kappa, pow_poly(U, 4)))
    A5 = shift(U5, 3)
    p45 = add(
        scale(15, mul(pow_poly(U, 3), deriv(W34))),
        scale(-102, mul(pow_poly(U, 2), mul(deriv(U), W34))),
        scale(14, mul(A, deriv(W35))), scale(-35, mul(deriv(A), W35)),
        scale(13, mul(A2, deriv(W36))), scale(-36, mul(deriv(A2), W36)),
        scale(12, mul(A3, deriv(W37))), scale(-37, mul(deriv(A3), W37)),
        scale(11, mul(A4, deriv(W38))), scale(-38, mul(deriv(A4), W38)),
        scale(10, mul(A5, deriv(g))), scale(-39, mul(deriv(A5), g)),
        scale(-6, mul(pow_poly(U, 5), q20)),
        scale(-6, mul(q24, q21)), scale(-6, mul(q23, q22)),
        scale(3*kappa, pow_poly(U, 9)),
    )
    p4_drop5 = add(
        scale(4, mul(pow_poly(U, 9), A5)),
        scale(12, mul(pow_poly(U, 6), add(mul(A, A4), mul(A2, A3)))),
        scale(12, mul(pow_poly(U, 3), add(mul(pow_poly(A, 2), A3),
                                          mul(A, pow_poly(A2, 2))))),
        scale(4, mul(pow_poly(A, 3), A2)),
    )
    q55 = add(
        scale(25, mul(pow_poly(U, 5), deriv(W34))),
        scale(-170, mul(pow_poly(U, 4), mul(deriv(U), W34))),
        scale(24, mul(q24, deriv(W35))), scale(-35, mul(deriv(q24), W35)),
        scale(23, mul(q23, deriv(W36))), scale(-36, mul(deriv(q23), W36)),
        scale(22, mul(q22, deriv(W37))), scale(-37, mul(deriv(q22), W37)),
        scale(21, mul(q21, deriv(W38))), scale(-38, mul(deriv(q21), W38)),
        scale(20, mul(q20, deriv(g))), scale(-39, mul(deriv(q20), g)),
        scale(-5, p4_drop5), scale(-9*kappa, pow_poly(U, 11)),
    )

    p_rows, p_forcing = coefficient_rows(p45, new_variables)
    q_rows, q_forcing = coefficient_rows(q55, new_variables)
    expanded_matrix = sp.Matrix([*p_rows, *q_rows])
    expanded_forcing = [*p_forcing, *q_forcing]
    direct_matrix, direct_forcing, _ = direct_45_55_system(data, kappa)

    if expanded_matrix != direct_matrix or len(expanded_forcing) != len(direct_forcing):
        raise AssertionError("direct and expanded 45/55 system shapes differ")
    if any(sp.expand(expanded_matrix[i, j] - direct_matrix[i, j]) != 0
           for i in range(expanded_matrix.rows)
           for j in range(expanded_matrix.cols)):
        raise AssertionError("direct and expanded 45/55 coefficient matrices differ")
    if any(sp.expand(expanded_forcing[i] - direct_forcing[i]) != 0
           for i in range(len(expanded_forcing))):
        raise AssertionError("direct and expanded 45/55 forcing vectors differ")
    if any(value.free_symbols.intersection({a, b, kappa})
           for row in direct_matrix.tolist() for value in row):
        raise AssertionError("direct coefficient matrix depends on carried parameters")
    print(f"45/55_direct_column_audit=PASS rows={direct_matrix.rows} "
          f"variables={direct_matrix.cols}")


if __name__ == "__main__":
    main()
