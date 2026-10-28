"""Exact weight-48 compatibility audit for both power-relation branches.

This carries the weight-49 kernel into the next [P,W*] row, allowing the
support-bounded W*_37 band and U2=hB band to vary freely.  It eliminates
those new linear variables by the left kernel and checks the residual ideal
in the weight-49 parameters.  The calculation is a necessary-condition
relaxation under the F2 support and edge hypotheses, not a pair exclusion.
"""

from __future__ import annotations

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import QQ

from f2_75_125_power_relation_second_band import build_matrix


def audit_branch(tag, h, V, field, expected_rank49, expected_rank48, expected_residuals):
    t = sp.symbols("t")
    U = sp.expand(t * h)
    if sp.simplify(5 * U * sp.diff(V, t) - 4 * sp.diff(U, t) * V - U) != 0:
        raise AssertionError(f"normalized W_39 profile fails its ODE in {tag}")
    M49 = build_matrix(h, V)
    D49 = DomainMatrix.from_Matrix(M49).convert_to(field)
    K49_domain = D49.nullspace()
    kernel_dimension = K49_domain.shape[0]
    if kernel_dimension != 64 - expected_rank49 or K49_domain.shape[1] != 64:
        raise AssertionError(f"unexpected weight-49 kernel shape in {tag}: {K49_domain.shape}")

    S_domain = K49_domain.extract(list(range(kernel_dimension)), list(range(4)))
    _, pivots = S_domain.rref()
    if len(pivots) != kernel_dimension:
        raise AssertionError("weight-49 kernel projection to S is not injective")
    pivot_block = S_domain.extract(list(range(kernel_dimension)), list(pivots))
    transform = pivot_block.inv()
    K49 = transform.matmul(K49_domain).to_Matrix()

    c = sp.symbols(f"c0:{K49.rows}")
    B_coeffs = sp.symbols("B0:8")
    f37_coeffs = sp.symbols("F8:68")
    S = sp.expand(sum(c[row] * sum(K49[row, k] * t**k for k in range(4))
                      for row in range(K49.rows)))
    f38 = sp.expand(sum(c[row] * sum(K49[row, 4 + k - 8] * t**k
                                    for k in range(8, 68))
                      for row in range(K49.rows)))

    B = sum(B_coeffs[k] * t**k for k in range(8))
    f37 = sum(f37_coeffs[k - 8] * t**k for k in range(8, 68))
    U1 = sp.expand(h**2 * S)
    U2 = sp.expand(h * B)
    A = sp.expand(t**3 * U1)
    A2 = sp.expand(t**3 * U2)
    q24 = sp.expand(sp.Rational(5, 3) * t**5 * h**2 * U1)
    q23 = sp.expand(sp.Rational(5, 3) * t**5 * h**2 * U2
                     + sp.Rational(5, 9) * t**5 * h**3 * S**2)
    g = sp.expand(U**7 * V)

    equation = sp.Poly(
        sp.expand(
            15 * U**3 * sp.diff(f37, t)
            - 111 * U**2 * sp.diff(U, t) * f37
            + 14 * A * sp.diff(f38, t)
            - 38 * sp.diff(A, t) * f38
            + 13 * A2 * sp.diff(g, t)
            - 39 * sp.diff(A2, t) * g
            - 6 * U**5 * q23
            - 3 * q24**2
        ),
        t,
    )

    linear_variables = (*f37_coeffs, *B_coeffs)
    zero_linear = {x: 0 for x in linear_variables}
    rows, forcing = [], []
    for degree in range(equation.degree() + 1):
        coefficient = sp.expand(equation.nth(degree))
        if not coefficient:
            continue
        constant = sp.expand(coefficient.subs(zero_linear))
        row = [sp.expand(sp.diff(coefficient, x).subs(zero_linear))
               for x in linear_variables]
        reconstruction = constant + sum(row[k] * linear_variables[k]
                                        for k in range(len(linear_variables)))
        if sp.expand(coefficient - reconstruction) != 0:
            raise AssertionError(f"weight-48 row is not linear in new bands ({tag})")
        rows.append(row)
        forcing.append(constant)

    M48 = sp.Matrix(rows)
    D48 = DomainMatrix.from_Matrix(M48.T).convert_to(field)
    left_kernel = D48.nullspace().to_Matrix()
    residuals = []
    for i in range(left_kernel.rows):
        residual = sp.expand(sum(left_kernel[i, j] * forcing[j]
                                 for j in range(len(forcing))))
        residual = sp.Poly(residual, *c, domain=field)
        if not residual.is_zero:
            residuals.append(sp.factor(residual.as_expr()))

    rank48 = D48.transpose().rank()
    actual = (M48.rows, M48.cols, rank48, left_kernel.rows, len(residuals))
    print(f"branch={tag} raw_weight48_data={actual}", flush=True)
    for index, residual in enumerate(residuals):
        print(f"raw_residual_{index}={residual}", flush=True)
    if (M48.cols, rank48, left_kernel.rows, len(residuals)) != (
        68, expected_rank48, M48.rows - expected_rank48, expected_residuals
    ):
        raise AssertionError(f"unexpected weight-48 rank data in {tag}: {actual}")

    print(f"branch={tag} weight49_kernel={K49.rows}")
    print(f"weight48_rows={M48.rows} linear_unknowns={M48.cols} rank={rank48}")
    print(f"left_kernel={left_kernel.rows} residual_polynomials={len(residuals)}")
    if residuals:
        basis = sp.groebner(residuals, *c, domain=field)
        print(f"residual_groebner_basis={list(basis)}")
        print(f"residual_ideal_unit={basis.contains(1)}")
    else:
        print("residual_ideal=zero")


def main() -> None:
    t = sp.symbols("t")
    sqrt5 = sp.sqrt(5)
    # The second golden root is the Galois conjugate of this calculation.
    for alpha in ((3 + sqrt5) / 2,):
        h = sp.expand((t - 1) ** 2 * (t - alpha) ** 2)
        slope, intercept = (4 * alpha - 11) / 3, 3 - alpha
        V = sp.expand(t * (t - 1) * (t - alpha) * (slope * t + intercept))
        audit_branch(
            f"2+2 alpha={alpha}", h, V, QQ.algebraic_field(sqrt5),
            expected_rank49=61, expected_rank48=61, expected_residuals=0,
        )

    r = t**2 - 3 * t + 3
    h = sp.expand((t - 1) ** 2 * r)
    V = sp.expand(-sp.Rational(1, 3) * t * (t - 1) * r)
    audit_branch(
        "2+1+1", h, V, QQ,
        expected_rank49=62, expected_rank48=62, expected_residuals=0,
    )
    print("weight48_audit=PASS")
    print("2+2_other_conjugate=covered by Q(sqrt(5)) automorphism")
    print("scope=relaxed necessary bands; no exclusion of (75,125)")


if __name__ == "__main__":
    main()
