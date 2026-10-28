"""Exact quadratic-field affine solver using doubled rational linear algebra.

For the F2 2+2 lift, all coefficient-matrix entries lie in Q(sqrt(5)).
This module represents a + b*sqrt(5) by the rational block action

    [[A, 5B], [B, A]],

so FLINT can perform the rank and pivot computations over QQ. Polynomial
compatibility equations are then reconstructed over Q(sqrt(5)); no
restriction is placed on the complex-valued carried parameters.
"""

from __future__ import annotations

import sympy as sp
from flint import fmpq, fmpq_mpoly_ctx
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix


def _field_parts(value, field) -> tuple[sp.Rational, sp.Rational]:
    """Return (rational part, sqrt(5) part) of a quadratic-field element."""
    algebraic = field.from_sympy(sp.sympify(value))
    coefficients = list(algebraic.to_list())
    if len(coefficients) > 2:
        raise ValueError("expected an element of a quadratic field")
    coefficients = [field.dom.zero] * (2 - len(coefficients)) + coefficients
    return sp.Rational(coefficients[-1]), sp.Rational(coefficients[-2])


def _flatten_matrix(matrix: sp.Matrix, field):
    rows, columns = matrix.shape
    rational = sp.zeros(rows, columns)
    radical = sp.zeros(rows, columns)
    for i in range(rows):
        for j in range(columns):
            rational[i, j], radical[i, j] = _field_parts(matrix[i, j], field)
    flattened = rational.row_join(5 * radical).col_join(
        radical.row_join(rational)
    )
    return rational, radical, flattened


class _RationalRowBasis:
    """Small exact echelon basis used to select whole quadratic column pairs."""

    def __init__(self, width: int):
        self.width = width
        self.rows: dict[int, list[sp.Rational]] = {}

    def insert(self, vector: list[sp.Rational]) -> bool:
        reduced = list(vector)
        for pivot, basis_row in sorted(self.rows.items()):
            coefficient = reduced[pivot]
            if coefficient:
                reduced = [x - coefficient * y
                           for x, y in zip(reduced, basis_row)]
        pivot = next((i for i, value in enumerate(reduced) if value), None)
        if pivot is None:
            return False
        scale = reduced[pivot]
        reduced = [value / scale for value in reduced]
        for old_pivot, basis_row in list(self.rows.items()):
            coefficient = basis_row[pivot]
            if coefficient:
                self.rows[old_pivot] = [x - coefficient * y
                                        for x, y in zip(basis_row, reduced)]
        self.rows[pivot] = reduced
        return True

    def __len__(self):
        return len(self.rows)


def _quadratic_independent_indices(flattened: sp.Matrix, original_size: int,
                                   field) -> tuple[list[int], int]:
    """Select a K-basis among pairs of rationalized columns.

    ``flattened`` is the block representation of a K-linear map. Its RREF
    gives coordinates of all rationalized columns in a QQ basis. We greedily
    retain an original K column exactly when its two rational columns add two
    dimensions, preserving the K-vector-space structure.
    """
    reduced, pivot_columns = (
        DomainMatrix.from_Matrix(flattened).convert_to(QQ).rref()
    )
    rank_q = len(pivot_columns)
    if rank_q % 2:
        raise AssertionError("quadratic extension rank must be even")
    rank_k = rank_q // 2
    if rank_k == 0:
        return [], 0
    rref_matrix = reduced.to_Matrix()
    basis = _RationalRowBasis(rank_q)
    selected: list[int] = []
    for j in range(original_size):
        first = [rref_matrix[i, j] for i in range(rank_q)]
        second = [rref_matrix[i, j + original_size]
                  for i in range(rank_q)]
        before = len(basis)
        added_first = basis.insert(first)
        added_second = basis.insert(second)
        increment = len(basis) - before
        if increment == 2 and added_first and added_second:
            selected.append(j)
        elif increment != 0 or added_first or added_second:
            raise AssertionError("a quadratic column pair added an odd rank")
        if len(basis) == rank_q:
            break
    if len(selected) != rank_k or len(basis) != rank_q:
        raise AssertionError("failed to extract a full quadratic pivot basis")
    return selected, rank_k


def _quadratic_inverse(matrix: sp.Matrix, field) -> sp.Matrix:
    """Invert a nonsingular K-matrix through its 2x2 rational block form."""
    size = matrix.rows
    if matrix.cols != size:
        raise ValueError("pivot block must be square")
    rational, radical, flattened = _flatten_matrix(matrix, field)
    domain = DomainMatrix.from_Matrix(flattened).convert_to(QQ)
    numerator, denominator = domain.inv_den()
    inverse_flat = numerator.to_Matrix() / field.dom.to_sympy(denominator)
    c = inverse_flat[:size, :size]
    d = inverse_flat[size:, :size]
    if inverse_flat[:size, size:] != 5 * d:
        raise AssertionError("inverse failed quadratic block symmetry (upper)")
    if inverse_flat[size:, size:] != c:
        raise AssertionError("inverse failed quadratic block symmetry (lower)")
    generator = field.ext.as_expr()
    return sp.Matrix(size, size, lambda i, j:
                     sp.Rational(c[i, j]) + generator * sp.Rational(d[i, j]))


def _fmpq(value) -> fmpq:
    value = sp.Rational(value)
    return fmpq(int(value.p), int(value.q))


def _polynomial_pair(expression, parameters, field, polynomial_context):
    # Converting an expanded coefficient such as N*sqrt(5) through
    # AlgebraicField.from_sympy can route through PSLQ/integers.factorint.
    # For the huge exact coefficients produced by later bands that path can
    # overflow even though the coefficient is visibly in Q(sqrt(5)). Work in
    # Q[s, parameters] instead, then reduce powers of s modulo the quadratic
    # minimal polynomial. This keeps the conversion exact and avoids a
    # number-field recognition heuristic.
    alpha = field.ext.as_expr()
    alpha_symbol = sp.Symbol("_quadratic_flat_alpha")
    expression_in_alpha = sp.expand(expression.xreplace({alpha: alpha_symbol}))
    polynomial = sp.Poly(expression_in_alpha, alpha_symbol, *parameters, domain=QQ)
    minimal_polynomial = sp.Poly(
        field.ext.minpoly.as_expr().subs(field.ext.minpoly.gen, alpha_symbol),
        alpha_symbol,
        domain=QQ,
    )
    if minimal_polynomial.degree() != 2:
        raise ValueError("only quadratic number fields are supported")

    power_remainders = {
        0: sp.Poly(1, alpha_symbol, domain=QQ),
        1: sp.Poly(alpha_symbol, alpha_symbol, domain=QQ),
    }

    def reduce_alpha_power(exponent: int):
        remainder = power_remainders.get(exponent)
        if remainder is None:
            remainder = sp.Poly(alpha_symbol**exponent, alpha_symbol, domain=QQ).rem(
                minimal_polynomial
            )
            power_remainders[exponent] = remainder
        return remainder

    rational_terms: dict[tuple[int, ...], fmpq] = {}
    radical_terms: dict[tuple[int, ...], fmpq] = {}
    for powers, coefficient in polynomial.terms():
        alpha_power, parameter_powers = powers[0], powers[1:]
        remainder = reduce_alpha_power(alpha_power)
        rational = coefficient * remainder.coeff_monomial(1)
        radical = coefficient * remainder.coeff_monomial(alpha_symbol)
        if rational:
            rational_terms[parameter_powers] = (
                rational_terms.get(parameter_powers, fmpq(0)) + _fmpq(rational)
            )
        if radical:
            radical_terms[parameter_powers] = (
                radical_terms.get(parameter_powers, fmpq(0)) + _fmpq(radical)
            )
    return (polynomial_context.from_dict(rational_terms),
            polynomial_context.from_dict(radical_terms))


def _polynomial_expression(poly, parameters) -> sp.Expr:
    terms = []
    for powers, coefficient in poly.terms():
        monomial = sp.Integer(1)
        for variable, exponent in zip(parameters, powers):
            if exponent:
                monomial *= variable ** int(exponent)
        rational = sp.Rational(int(coefficient.p), int(coefficient.q))
        terms.append(rational * monomial)
    return sp.Add(*terms) if terms else sp.Integer(0)


def _quadratic_residuals(nonpivot_rows, pivot_rows, pivot_coefficients,
                         forcing, parameters, field) -> list[sp.Expr]:
    """Compute exact K-polynomial row compatibility with sparse FLINT polys."""
    polynomial_context = fmpq_mpoly_ctx.get(tuple(str(p) for p in parameters))
    forcing_pairs = [
        _polynomial_pair(value, parameters, field, polynomial_context)
        for value in forcing
    ]
    coefficient_cache = {}
    residuals = []
    generator = field.ext.as_expr()
    zero = polynomial_context.constant(0)
    for local_row, original_row in enumerate(nonpivot_rows):
        rational = forcing_pairs[original_row][0]
        radical = forcing_pairs[original_row][1]
        for pivot_position, pivot_row in enumerate(pivot_rows):
            coefficient = pivot_coefficients[local_row, pivot_position]
            if coefficient == 0:
                continue
            key = sp.srepr(coefficient)
            parts = coefficient_cache.get(key)
            if parts is None:
                parts = _field_parts(coefficient, field)
                coefficient_cache[key] = parts
            c0, c1 = _fmpq(parts[0]), _fmpq(parts[1])
            b0, b1 = forcing_pairs[pivot_row]
            rational -= c0 * b0 + (5 * c1) * b1
            radical -= c0 * b1 + c1 * b0
        if rational != zero or radical != zero:
            residual = _polynomial_expression(rational, parameters)
            radical_expression = _polynomial_expression(radical, parameters)
            residual = sp.expand(residual + generator * radical_expression)
            residuals.append(residual)
    print(
        f"quadratic_flat_compatibility rows={len(forcing)} "
        f"pivot_rows={len(pivot_rows)} constraints={len(residuals)}",
        flush=True,
    )
    return residuals


def quadratic_affine_solution(M: sp.Matrix, forcing: list[sp.Expr], field,
                              variables: tuple[sp.Symbol, ...],
                              parameters: tuple[sp.Symbol, ...], tag: str):
    """Solve an affine system over Q(sqrt(5)) without number-field RREF.

    Returns the same ``(solution, free, rank, constraints)`` tuple as the
    power-relation lift's existing affine solver.
    """
    if field.ext.minpoly.degree() != 2:
        raise ValueError("only quadratic number fields are supported")
    rational, radical, flattened = _flatten_matrix(M, field)
    columns, rank_columns = _quadratic_independent_indices(
        flattened, M.cols, field
    )
    rows, rank_rows = _quadratic_independent_indices(
        _flatten_matrix(M.T, field)[2], M.rows, field
    )
    if rank_columns != rank_rows:
        raise AssertionError(f"rank mismatch for {tag}")
    rank = rank_columns
    if rank == 0:
        raise NotImplementedError("zero-rank affine block is not expected")

    pivot_block = M.extract(rows, columns)
    inverse = _quadratic_inverse(pivot_block, field)
    nonpivot_columns = [j for j in range(M.cols) if j not in columns]
    nonpivot_rows = [i for i in range(M.rows) if i not in rows]
    free = sp.symbols(f"d0:{len(nonpivot_columns)}")

    pivot_rhs = sp.Matrix([forcing[i] for i in rows])
    if free:
        pivot_rhs += M.extract(rows, nonpivot_columns) * sp.Matrix(free)
    pivot_values = -(inverse * pivot_rhs)
    solution = [sp.Integer(0)] * M.cols
    for column, symbol in zip(nonpivot_columns, free):
        solution[column] = symbol
    for column, value in zip(columns, pivot_values):
        solution[column] = sp.expand(value)

    pivot_coefficients = M.extract(nonpivot_rows, columns) * inverse
    constraints = _quadratic_residuals(
        nonpivot_rows, rows, pivot_coefficients, forcing, parameters, field
    )
    print(
        f"{tag}_affine_profile rank={rank} fibre={len(free)} "
        f"compatibility_constraints={len(constraints)}",
        flush=True,
    )
    return solution, free, rank, constraints
