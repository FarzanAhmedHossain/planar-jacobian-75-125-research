"""Test the next coupled commutator bands after the weight-49 selector.

The script first solves the compatible weight-48/58 linear block over the full
weight-49 kernel, parametrizes its exact affine fibres, and then imposes
weight 47 of [P,W*] with weight 57 of [Q,W*].  It eliminates the next shared
coefficients W*_36 and U3 by an exact left-kernel calculation.

All systems are necessary-condition relaxations conditional on the F2 edge
and support model; even a branch exclusion here would not by itself prove the
degree-pair exclusion.
"""

from __future__ import annotations

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import QQ

from f2_75_125_power_relation_pq_nextbands import (
    coefficient_system,
    weight49_parameters,
)


def exact_affine_solution(M: sp.Matrix, forcing: list[sp.Expr], field,
                          variables: tuple[sp.Symbol, ...], parameters: tuple[sp.Symbol, ...],
                          tag: str):
    """Solve M*x + forcing(parameters)=0 without dividing by a variable."""
    domain = DomainMatrix.from_Matrix(M).convert_to(field)
    rref, pivot_columns = domain.rref()
    transposed = DomainMatrix.from_Matrix(M.T).convert_to(field)
    transposed_rref, pivot_rows = transposed.rref()
    rank = len(pivot_columns)
    if len(pivot_rows) != rank:
        raise AssertionError(f"inconsistent exact rank profile for {tag}")
    left_kernel = nullspace_from_rref(transposed_rref, pivot_rows, M.rows)
    print(f"{tag}_pivot_profile_ready rank={rank}", flush=True)
    nonpivot_columns = [j for j in range(M.cols) if j not in pivot_columns]
    free = sp.symbols(f"d0:{len(nonpivot_columns)}")
    pivot_block = M.extract(pivot_rows, pivot_columns)
    free_block = M.extract(pivot_rows, nonpivot_columns)
    forcing_block = sp.Matrix([forcing[i] for i in pivot_rows])
    pivot_domain = DomainMatrix.from_Matrix(pivot_block).convert_to(field)
    inverse_numerator, inverse_denominator = pivot_domain.inv_den()
    denominator = field.to_sympy(inverse_denominator)
    pivot_values = -(inverse_numerator.to_Matrix() * (
        forcing_block + free_block * sp.Matrix(free)
    )) / denominator
    print(f"{tag}_pivot_solve_ready", flush=True)
    solution = [sp.Integer(0)] * M.cols
    for column, symbol in zip(nonpivot_columns, free):
        solution[column] = symbol
    for column, value in zip(pivot_columns, pivot_values):
        solution[column] = sp.expand(value)

    # Exact left-kernel compatibility guarantees that the selected pivot-row
    # solution satisfies every omitted equation as a polynomial identity.
    compatibility = residual_polynomials(
        M, forcing, field, parameters, f"compatibility {tag}",
        left_kernel=left_kernel, rank=rank,
    )
    if compatibility:
        raise AssertionError(f"incompatible forcing in {tag}")
    print(f"{tag}_affine_solution_verified", flush=True)
    return solution, free, rank


def coefficient_zero_checker(field):
    """Build an exact zero checker for coefficients in QQ or QQ(sqrt(5))."""
    if field == sp.QQ:
        return lambda value: sp.cancel(value) == 0
    alpha = sp.Dummy("alpha")
    generator = field.ext.as_expr()
    modulus = sp.Poly(sp.minimal_polynomial(generator, alpha), alpha, domain=sp.QQ)

    def check(value):
        formal = sp.cancel(value).xreplace({generator: alpha})
        numerator, denominator = sp.fraction(sp.cancel(formal))
        numerator_poly = sp.Poly(numerator, alpha, domain=sp.QQ)
        denominator_poly = sp.Poly(denominator, alpha, domain=sp.QQ)
        if sp.rem(denominator_poly, modulus).is_zero:
            raise ZeroDivisionError("coefficient has zero denominator in the number field")
        return sp.rem(numerator_poly, modulus).is_zero

    return check


def nullspace_from_rref(rref: DomainMatrix, pivot_columns, ncolumns: int):
    """Build a row-basis of a right kernel from an already computed RREF."""
    reduced = rref.to_Matrix()
    free_columns = [j for j in range(ncolumns) if j not in pivot_columns]
    basis = sp.zeros(len(free_columns), ncolumns)
    for i, free_column in enumerate(free_columns):
        basis[i, free_column] = 1
        for row, pivot_column in enumerate(pivot_columns):
            basis[i, pivot_column] = -reduced[row, free_column]
    return basis


def coefficient_is_zero_over_field(value: sp.Expr, field) -> bool:
    """Decide exact zero for one coefficient in QQ or QQ(sqrt(5))."""
    return coefficient_zero_checker(field)(value)


def parameter_polynomial_is_zero(value: sp.Expr, parameters, field) -> bool:
    """Test a parameter polynomial without SymPy's EX/algebraic autodetection.

    In the quadratic branch, regard sqrt(5) as a formal symbol, clear the
    parameter-free field denominator, and reduce each coefficient modulo
    alpha**2-5. This avoids expression-domain zero assumptions that can try
    to factor enormous exact integers through a floating logarithm.
    """
    if field == sp.QQ:
        polynomial = sp.Poly(value, *parameters, domain=sp.QQ)
        return polynomial.is_zero

    alpha = sp.Dummy("alpha")
    generator = field.ext.as_expr()
    formal = sp.together(sp.expand(value).xreplace({generator: alpha}))
    numerator, denominator = sp.fraction(formal)
    if denominator.free_symbols.intersection(parameters):
        raise AssertionError("parameter polynomial unexpectedly has a parameter denominator")
    modulus = sp.Poly(alpha**2 - 5, alpha, domain=sp.QQ)
    denominator_poly = sp.Poly(denominator, alpha, domain=sp.QQ)
    if sp.rem(denominator_poly, modulus).is_zero:
        raise ZeroDivisionError("parameter coefficient denominator vanishes in QQ(sqrt(5))")

    numerator_poly = sp.Poly(numerator, *parameters, alpha, domain=sp.QQ)
    reduced: dict[tuple[int, ...], list[sp.Rational]] = {}
    for powers, coefficient in numerator_poly.terms():
        parameter_powers, alpha_power = powers[:-1], powers[-1]
        pair = reduced.setdefault(parameter_powers, [sp.Rational(0), sp.Rational(0)])
        pair[alpha_power % 2] += coefficient * 5 ** (alpha_power // 2)
    return all(rational_part == 0 and sqrt5_part == 0
               for rational_part, sqrt5_part in reduced.values())


def residual_polynomials(M: sp.Matrix, forcing: list[sp.Expr], field,
                         parameters: tuple[sp.Symbol, ...], tag: str,
                         left_kernel: sp.Matrix | None = None,
                         rank: int | None = None):
    if left_kernel is None:
        left_kernel = DomainMatrix.from_Matrix(M.T).convert_to(field).nullspace().to_Matrix()
    Dleft = left_kernel
    if rank is None:
        rank = DomainMatrix.from_Matrix(M).convert_to(field).rank()
    residuals = []
    for i in range(Dleft.rows):
        residual = sp.expand(sum(Dleft[i, j] * forcing[j] for j in range(M.rows)))
        if not parameter_polynomial_is_zero(residual, parameters, field):
            residuals.append(sp.factor_terms(residual))
    print(
        f"{tag}_rows={M.rows} variables={M.cols} rank={rank} "
        f"left_kernel={Dleft.rows} nonzero_residuals={len(residuals)}",
        flush=True,
    )
    return residuals


def audit_branch(tag: str, h: sp.Expr, V: sp.Expr, field) -> None:
    print(f"branch_start={tag}", flush=True)
    t = sp.symbols("t")
    U = sp.expand(t * h)
    c, S, f38, _ = weight49_parameters(h, V, field, tag)
    old_parameters = c

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
    q58 = (
        5 * U**4 * (5 * U * sp.diff(f37, t) - 37 * sp.diff(U, t) * f37)
        + 24 * q24 * sp.diff(f38, t)
        - 38 * sp.diff(q24, t) * f38
        + 23 * q23 * sp.diff(g, t)
        - 39 * sp.diff(q23, t) * g
        - 20 * U**9 * A2
        - 30 * U**6 * A**2
    )

    old_variables = (*f37_coeffs, *B_coeffs)
    Mp, bp = coefficient_system(p48, t, c, old_variables, f"P48 {tag}")
    Mq, bq = coefficient_system(q58, t, c, old_variables, f"Q58 {tag}")
    M_old = sp.Matrix.vstack(Mp, Mq)
    b_old = [*bp, *bq]
    old_solution, fibre, rank_old = exact_affine_solution(
        M_old, b_old, field, old_variables, old_parameters, f"48/58 {tag}"
    )
    parameter_tuple = (*old_parameters, *fibre)
    substitution = dict(zip(old_variables, old_solution))
    B_solved = sp.expand(B.subs(substitution))
    f37_solved = sp.expand(f37.subs(substitution))

    U2_solved = sp.expand(h * B_solved)
    A2_solved = sp.expand(t**3 * U2_solved)
    U3_coeffs = sp.symbols("C0:12")
    U3 = sum(U3_coeffs[i] * t**i for i in range(12))
    A3 = sp.expand(t**3 * U3)
    q23_solved = sp.expand(
        sp.Rational(5, 3) * t**5 * h**2 * U2_solved
        + sp.Rational(5, 9) * t**5 * h**3 * S**2
    )
    q22 = sp.expand(
        t**5 * h**2 * (
            sp.Rational(5, 3) * U3
            + sp.Rational(10, 9) * S * B_solved
            - sp.Rational(5, 81) * S**3
        )
    )
    L_coeffs = sp.symbols("L8:68")
    ell = sum(L_coeffs[i - 8] * t**i for i in range(8, 68))

    # Weight 47 in [P,W*] = 3 J Q^2.
    p47 = (
        15 * U**3 * sp.diff(ell, t)
        - 108 * U**2 * sp.diff(U, t) * ell
        + 14 * A * sp.diff(f37_solved, t)
        - 37 * sp.diff(A, t) * f37_solved
        + 13 * A2_solved * sp.diff(f38, t)
        - 38 * sp.diff(A2_solved, t) * f38
        + 12 * A3 * sp.diff(g, t)
        - 39 * sp.diff(A3, t) * g
        - 6 * U**5 * q22
        - 6 * q24 * q23_solved
    )

    # Weight 57 in [Q,W*] = 5 lambda J P^4.
    q57 = (
        5 * U**4 * (5 * U * sp.diff(ell, t) - 36 * sp.diff(U, t) * ell)
        + 24 * q24 * sp.diff(f37_solved, t)
        - 37 * sp.diff(q24, t) * f37_solved
        + 23 * q23_solved * sp.diff(f38, t)
        - 38 * sp.diff(q23_solved, t) * f38
        + 22 * q22 * sp.diff(g, t)
        - 39 * sp.diff(q22, t) * g
        - 20 * U**9 * A3
        - 60 * U**6 * A * A2_solved
        - 20 * U**3 * A**3
    )

    new_variables = (*L_coeffs, *U3_coeffs)
    M1, b1 = coefficient_system(p47, t, parameter_tuple, new_variables,
                                f"P47 {tag}")
    M2, b2 = coefficient_system(q57, t, parameter_tuple, new_variables,
                                f"Q57 {tag}")
    M_next = sp.Matrix.vstack(M1, M2)
    b_next = [*b1, *b2]
    residuals = residual_polynomials(
        M_next, b_next, field, parameter_tuple, f"{tag}_nextbands"
    )
    rank_next = DomainMatrix.from_Matrix(M_next).convert_to(field).rank()
    print(
        f"branch={tag} old_rank={rank_old} old_kernel_parameters={len(parameter_tuple)} "
        f"new_fibre_dimension={len(new_variables) - rank_next}",
        flush=True,
    )
    for i, residual in enumerate(residuals):
        print(f"{tag}_residual_{i}={residual}", flush=True)
    if residuals:
        basis = sp.groebner(residuals, *parameter_tuple, domain=field)
        print(f"{tag}_groebner_basis={list(basis)}", flush=True)
        print(f"{tag}_residual_ideal_unit={basis.contains(1)}", flush=True)
    else:
        print(f"{tag}_residual_ideal=zero", flush=True)


def main() -> None:
    t = sp.symbols("t")
    sqrt5 = sp.sqrt(5)
    r = t**2 - 3 * t + 3
    h211 = sp.expand((t - 1) ** 2 * r)
    v211 = sp.expand(-sp.Rational(1, 3) * t * (t - 1) * r)
    audit_branch("2+1+1", h211, v211, QQ)

    alpha = (3 + sqrt5) / 2
    h22 = sp.expand((t - 1) ** 2 * (t - alpha) ** 2)
    m, n = (4 * alpha - 11) / 3, 3 - alpha
    v22 = sp.expand(t * (t - 1) * (t - alpha) * (m * t + n))
    audit_branch("2+2", h22, v22, QQ.algebraic_field(sqrt5))
    print("scope=three coupled necessary band levels in relaxed F2 model; not a pair exclusion")


if __name__ == "__main__":
    main()
