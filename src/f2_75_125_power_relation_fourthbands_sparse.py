"""Coefficient-wise exact replay of the 46/56 forcing compatibility.

The original checker expands the complete expressions in t before extracting
their coefficients.  This version keeps each expression as a univariate
coefficient list and performs exact convolutions, so the same compatibility
test can be replayed without expanding the full 20-parameter expression at
once.  It uses the same edge hypotheses and first three band lifts as the
original checker; it is not an exclusion of (75,125) unless its necessary
model is attached to every hypothetical pair.
"""

from __future__ import annotations

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import QQ

from f2_75_125_power_relation_fourthbands import lift_first_two_blocks, verify_band_formulae
from f2_75_125_power_relation_pq_thirdbands import (
    coefficient_zero_checker,
    nullspace_from_rref,
)


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def from_expr(expr, t):
    p = sp.Poly(expr, t)
    return trim([p.nth(i) for i in range(p.degree() + 1)])


def add(*polys):
    n = max((len(p) for p in polys), default=1)
    out = [sp.Integer(0)] * n
    for p in polys:
        for i, c in enumerate(p):
            out[i] += c
    return trim([sp.expand(c) for c in out])


def scale(c, p):
    return trim([sp.expand(c * x) for x in p])


def neg(p):
    return scale(-1, p)


def sub(a, b):
    return add(a, neg(b))


def mul(a, b):
    out = [sp.Integer(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            if y != 0:
                out[i + j] += x * y
    return trim([sp.expand(c) for c in out])


def deriv(p):
    if len(p) <= 1:
        return [sp.Integer(0)]
    return trim([sp.expand(i * p[i]) for i in range(1, len(p))])


def shift(p, n):
    return [sp.Integer(0)] * n + list(p)


def pow_poly(p, e):
    out = [sp.Integer(1)]
    for _ in range(e):
        out = mul(out, p)
    return out


def coefficient_rows(poly, variables):
    rows = []
    constants = []
    for c in poly:
        row = [sp.diff(c, v) for v in variables]
        constant = c.subs({v: 0 for v in variables})
        reconstruction = constant + sum(row[i] * variables[i]
                                        for i in range(len(variables)))
        if sp.expand(c - reconstruction) != 0:
            raise AssertionError("coefficient is not affine-linear in new slices")
        if any(sp.sympify(x).free_symbols.intersection(variables) for x in row):
            raise AssertionError("new-slice matrix depends on carried parameters")
        if constant != 0 or any(x != 0 for x in row):
            rows.append(row)
            constants.append(sp.expand(constant))
    return rows, constants


def branch(tag: str, h, V, field):
    print(f"sparse_branch_start={tag}", flush=True)
    t = sp.symbols("t")
    data = lift_first_two_blocks(tag, h, V, field)
    parameters = data["parameters"]

    # Freeze all earlier slices as coefficient lists in t.
    U = from_expr(data["U"], t)
    hpoly = from_expr(h, t)
    S = from_expr(data["S"], t)
    B = from_expr(data["B"], t)
    U3 = from_expr(data["U3"], t)
    A = from_expr(data["A"], t)
    A2 = from_expr(data["A2"], t)
    A3 = from_expr(data["A3"], t)
    q24 = from_expr(data["q24"], t)
    q23 = from_expr(data["q23"], t)
    q22 = from_expr(data["q22"], t)
    g = from_expr(data["g"], t)
    f38 = from_expr(data["f38"], t)
    f37 = from_expr(data["f37"], t)
    ell = from_expr(data["ell"], t)

    # The known part of q_21; the new U_4 coefficient is kept separate.
    q21_force_inner = add(
        scale(54, mul(S, U3)),
        scale(27, pow_poly(B, 2)),
        scale(-9, mul(pow_poly(S, 2), B)),
        pow_poly(S, 4),
    )
    q21_force = scale(sp.Rational(5, 243),
                      shift(mul(hpoly, q21_force_inner), 5))

    # Forcing polynomials in the two bands, built by exact t-convolution.
    p_force = add(
        scale(14, mul(A, deriv(ell))),
        scale(-36, mul(deriv(A), ell)),
        scale(13, mul(A2, deriv(f37))),
        scale(-37, mul(deriv(A2), f37)),
        scale(12, mul(A3, deriv(f38))),
        scale(-38, mul(deriv(A3), f38)),
        scale(-6, mul(pow_poly(U, 5), q21_force)),
        scale(-6, mul(q24, q22)),
        scale(-3, pow_poly(q23, 2)),
    )
    q_force = add(
        scale(24, mul(q24, deriv(ell))),
        scale(-36, mul(deriv(q24), ell)),
        scale(23, mul(q23, deriv(f37))),
        scale(-37, mul(deriv(q23), f37)),
        scale(22, mul(q22, deriv(f38))),
        scale(-38, mul(deriv(q22), f38)),
        scale(21, mul(q21_force, deriv(g))),
        scale(-39, mul(deriv(q21_force), g)),
        scale(-60, mul(pow_poly(U, 6), mul(A, A3))),
        scale(-30, mul(pow_poly(U, 6), pow_poly(A2, 2))),
        scale(-60, mul(pow_poly(U, 3), mul(pow_poly(A, 2), A2))),
        scale(-5, pow_poly(A, 4)),
    )
    print(f"{tag}_forcing_lists_ready pdeg={len(p_force)-1} qdeg={len(q_force)-1}", flush=True)

    W_coeffs = sp.symbols("W35_7:68")
    U4_coeffs = sp.symbols("U4_0:12")
    W35 = [sp.Integer(0)] * 7 + list(W_coeffs)
    U4 = list(U4_coeffs)
    A4 = shift(U4, 3)
    q21_linear = scale(sp.Rational(5, 3),
                       shift(mul(pow_poly(hpoly, 2), U4), 5))
    p_linear = add(
        scale(15, mul(pow_poly(U, 3), deriv(W35))),
        scale(-105, mul(pow_poly(U, 2), mul(deriv(U), W35))),
        scale(11, mul(A4, deriv(g))),
        scale(-39, mul(deriv(A4), g)),
        scale(-6, mul(pow_poly(U, 5), q21_linear)),
    )
    q_linear = add(
        scale(25, mul(pow_poly(U, 5), deriv(W35))),
        scale(-175, mul(pow_poly(U, 4), mul(deriv(U), W35))),
        scale(21, mul(q21_linear, deriv(g))),
        scale(-39, mul(deriv(q21_linear), g)),
        scale(-20, mul(pow_poly(U, 9), A4)),
    )

    p_rows, p_constants = coefficient_rows(p_linear, (*W_coeffs, *U4_coeffs))
    q_rows, q_constants = coefficient_rows(q_linear, (*W_coeffs, *U4_coeffs))
    M = sp.Matrix([*p_rows, *q_rows])
    forcing = [
        *(p_force[i] if i < len(p_force) else 0 for i in range(len(p_linear))),
        *(q_force[i] if i < len(q_force) else 0 for i in range(len(q_linear))),
    ]
    # Keep exactly the nonzero coefficient rows selected above.
    forcing = [*(forcing[i] for i, c in enumerate(p_linear) if c != 0),
               *(forcing[len(p_linear) + i] for i, c in enumerate(q_linear) if c != 0)]
    if len(forcing) != M.rows:
        raise AssertionError("forcing and coefficient row indexing mismatch")

    transposed_domain = DomainMatrix.from_Matrix(M.T).convert_to(field)
    left_rref, pivots = transposed_domain.rref()
    rank = len(pivots)
    left = nullspace_from_rref(left_rref, pivots, M.rows)
    print(f"{tag}_matrix rows={M.rows} vars={M.cols} rank={rank} cokernel={left.rows}", flush=True)
    residuals = []
    coefficient_is_zero = coefficient_zero_checker(field)
    for i in range(left.rows):
        residual = sp.expand(sum(left[i, j] * forcing[j] for j in range(M.rows)))
        polynomial = sp.Poly(residual, *parameters, domain=sp.EX)
        if not all(coefficient_is_zero(c) for c in polynomial.coeffs()):
            residuals.append(sp.factor_terms(polynomial.as_expr()))
    print(f"{tag}_nonzero_residuals={len(residuals)}", flush=True)
    if residuals:
        for i, r in enumerate(residuals):
            print(f"{tag}_residual_{i}={r}", flush=True)
    return residuals


def main():
    verify_band_formulae()
    t = sp.symbols("t")
    sqrt5 = sp.sqrt(5)
    alpha = (3 + sqrt5) / 2
    h22 = sp.expand((t - 1) ** 2 * (t - alpha) ** 2)
    m, n = (4 * alpha - 11) / 3, 3 - alpha
    v22 = sp.expand(t * (t - 1) * (t - alpha) * (m * t + n))
    branch("2+2", h22, v22, QQ.algebraic_field(sqrt5))

    r = t**2 - 3 * t + 3
    h211 = sp.expand((t - 1) ** 2 * r)
    v211 = sp.expand(-sp.Rational(1, 3) * t * (t - 1) * r)
    branch("2+1+1", h211, v211, QQ)
    print("scope=exact 46/56 forcing compatibility via coefficient convolutions; edge-model scope only")


if __name__ == "__main__":
    main()
