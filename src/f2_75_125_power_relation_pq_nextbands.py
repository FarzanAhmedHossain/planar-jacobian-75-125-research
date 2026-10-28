"""Couple the next [P,W*] and [Q,W*] bands exactly.

This is a support-relaxed necessary system conditional on the F2 edge model.
It carries the weight-49 kernel into weight 48 of [P,W*] and weight 58 of
[Q,W*], then eliminates the shared new coefficients W*_37 and U2=hB by
exact left-kernel compatibility over the relevant number fields.
"""

from __future__ import annotations

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import QQ

from f2_75_125_power_relation_second_band import build_matrix


def weight49_parameters(h: sp.Expr, V: sp.Expr, field, tag: str):
    t = sp.symbols("t")
    M49 = build_matrix(h, V)
    D49 = DomainMatrix.from_Matrix(M49).convert_to(field)
    kernel = D49.nullspace()
    dimension = kernel.shape[0]
    if kernel.shape[1] != 64:
        raise AssertionError(f"unexpected weight-49 kernel shape in {tag}: {kernel.shape}")
    s_projection = kernel.extract(list(range(dimension)), list(range(4)))
    _, pivots = s_projection.rref()
    if len(pivots) != dimension:
        raise AssertionError(f"weight-49 kernel projection to S is not injective in {tag}")
    pivot_block = s_projection.extract(list(range(dimension)), list(pivots))
    K = pivot_block.inv().matmul(kernel).to_Matrix()

    c = sp.symbols(f"c0:{dimension}")
    S = sp.expand(sum(c[i] * sum(K[i, j] * t**j for j in range(4))
                      for i in range(dimension)))
    f38 = sp.expand(sum(c[i] * sum(K[i, 4 + (j - 8)] * t**j
                                   for j in range(8, 68))
                          for i in range(dimension)))
    return c, S, f38, K


def coefficient_system(expression: sp.Expr, t, c, linear_variables, tag: str):
    polynomial = sp.Poly(sp.expand(expression), t)
    zero_linear = {x: 0 for x in linear_variables}
    rows, forcing = [], []
    for degree in range(polynomial.degree() + 1):
        coefficient = sp.expand(polynomial.nth(degree))
        if not coefficient:
            continue
        constant = sp.expand(coefficient.subs(zero_linear))
        row = [sp.expand(sp.diff(coefficient, x).subs(zero_linear))
               for x in linear_variables]
        reconstruction = constant + sum(row[i] * linear_variables[i]
                                        for i in range(len(linear_variables)))
        if sp.expand(coefficient - reconstruction) != 0:
            raise AssertionError(f"{tag} is not linear in the new slices")
        if any(sp.sympify(entry).free_symbols.intersection(c) for entry in row):
            raise AssertionError(f"{tag} matrix depends on old parameters")
        rows.append(row)
        forcing.append(constant)
    return sp.Matrix(rows), forcing


def audit_branch(tag: str, h: sp.Expr, V: sp.Expr, field) -> None:
    t = sp.symbols("t")
    U = sp.expand(t * h)
    c, S, f38, _ = weight49_parameters(h, V, field, tag)

    B_coeffs = sp.symbols("B0:8")
    f37_coeffs = sp.symbols("F8:68")
    B = sum(B_coeffs[i] * t**i for i in range(8))
    f37 = sum(f37_coeffs[i - 8] * t**i for i in range(8, 68))

    U1 = sp.expand(h**2 * S)
    U2 = sp.expand(h * B)
    A = sp.expand(t**3 * U1)
    A2 = sp.expand(t**3 * U2)
    q24 = sp.expand(sp.Rational(5, 3) * t**5 * h**2 * U1)
    q23 = sp.expand(sp.Rational(5, 3) * t**5 * h**2 * U2
                     + sp.Rational(5, 9) * t**5 * h**3 * S**2)
    g = sp.expand(U**7 * V)

    # Weight 48 in [P,W*] = 3J Q^2.
    p48 = (
        15 * U**3 * sp.diff(f37, t)
        - 111 * U**2 * sp.diff(U, t) * f37
        + 14 * A * sp.diff(f38, t)
        - 38 * sp.diff(A, t) * f38
        + 13 * A2 * sp.diff(g, t)
        - 39 * sp.diff(A2, t) * g
        - 6 * U**5 * q23
        - 3 * q24**2
    )

    # Weight 58 in [Q,W*] = 5 lambda J P^4. Corrections to W have
    # commutator weight at most 55 and so do not enter this band.
    q58 = (
        5 * U**4 * (5 * U * sp.diff(f37, t) - 37 * sp.diff(U, t) * f37)
        + 24 * q24 * sp.diff(f38, t)
        - 38 * sp.diff(q24, t) * f38
        + 23 * q23 * sp.diff(g, t)
        - 39 * sp.diff(q23, t) * g
        - 20 * U**9 * A2
        - 30 * U**6 * A**2
    )

    new_variables = (*f37_coeffs, *B_coeffs)
    Mp, bp = coefficient_system(p48, t, c, new_variables, f"P48 {tag}")
    Mq, bq = coefficient_system(q58, t, c, new_variables, f"Q58 {tag}")
    M = sp.Matrix.vstack(Mp, Mq)
    b = [*bp, *bq]

    M_domain = DomainMatrix.from_Matrix(M).convert_to(field)
    rank = M_domain.rank()
    left_kernel = DomainMatrix.from_Matrix(M.T).convert_to(field).nullspace().to_Matrix()
    residuals = []
    for i in range(left_kernel.rows):
        residual = sp.expand(sum(left_kernel[i, j] * b[j] for j in range(M.rows)))
        residual_poly = sp.Poly(residual, *c, domain=field)
        if not residual_poly.is_zero:
            residuals.append(sp.factor(residual_poly.as_expr()))

    print(
        f"branch={tag} P_rows={Mp.rows} Q_rows={Mq.rows} variables={M.cols} "
        f"rank={rank} left_kernel={left_kernel.rows} residuals={len(residuals)}",
        flush=True,
    )
    for i, residual in enumerate(residuals):
        print(f"residual_{i}={residual}", flush=True)
    if residuals:
        basis = sp.groebner(residuals, *c, domain=field)
        print(f"residual_groebner_basis={list(basis)}", flush=True)
        print(f"residual_ideal_unit={basis.contains(1)}", flush=True)
    else:
        print("residual_ideal=zero", flush=True)


def main() -> None:
    t = sp.symbols("t")
    sqrt5 = sp.sqrt(5)
    alpha = (3 + sqrt5) / 2
    h22 = sp.expand((t - 1) ** 2 * (t - alpha) ** 2)
    m, n = (4 * alpha - 11) / 3, 3 - alpha
    v22 = sp.expand(t * (t - 1) * (t - alpha) * (m * t + n))
    audit_branch("2+2", h22, v22, QQ.algebraic_field(sqrt5))

    r = t**2 - 3 * t + 3
    h211 = sp.expand((t - 1) ** 2 * r)
    v211 = sp.expand(-sp.Rational(1, 3) * t * (t - 1) * r)
    audit_branch("2+1+1", h211, v211, QQ)
    print("scope=two coupled necessary bands in relaxed F2 edge model; not a pair exclusion")


if __name__ == "__main__":
    main()
