"""Extend the corrected F2 power-relation calculation through weights 45/55.

This is a support-relaxed necessary system.  It first lifts each selected
carrier through the already-audited 48/58, 47/57, and 46/56 affine blocks,
then computes exact left-kernel compatibility for the next pair.  At this
level the correction -3*kappa*P^3*Q is essential: it removes the resonant
weight-70 face of Q^3-P^5 and contributes at weights 45 and 55.

Even a branch exclusion here would remain conditional on the theorem-facing
F2 support bridge and would not, by itself, exclude the degree pair.
"""

from __future__ import annotations

import hashlib
import pickle
import sys
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import QQ

import f2_75_125_power_relation_fourthbands as fourthbands_module
from f2_75_125_power_relation_fourthbands import lift_first_two_blocks
from f2_75_125_power_relation_pq_nextbands import coefficient_system
from f2_75_125_power_relation_pq_thirdbands import (
    exact_affine_solution,
    nullspace_from_rref,
    residual_polynomials,
)
from f2_75_125_quadratic_flat_solver import quadratic_affine_solution


def trim(poly):
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def add(*polys):
    n = max((len(p) for p in polys), default=1)
    out = [sp.Integer(0)] * n
    for poly in polys:
        for i, coefficient in enumerate(poly):
            out[i] += coefficient
    return trim([sp.expand(value) for value in out])


def scale(coefficient, poly):
    return trim([sp.expand(coefficient * value) for value in poly])


def neg(poly):
    return scale(-1, poly)


def sub(left, right):
    return add(left, neg(right))


def mul(left, right):
    out = [sp.Integer(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        if a == 0:
            continue
        for j, b in enumerate(right):
            if b != 0:
                out[i + j] += a * b
    return trim([sp.expand(value) for value in out])


def deriv(poly):
    if len(poly) <= 1:
        return [sp.Integer(0)]
    return trim([sp.expand(i * poly[i]) for i in range(1, len(poly))])


def shift(poly, amount):
    return [sp.Integer(0)] * amount + list(poly)


def pow_poly(poly, exponent):
    out = [sp.Integer(1)]
    for _ in range(exponent):
        out = mul(out, poly)
    return out


def from_expr(expression, t):
    polynomial = sp.Poly(expression, t)
    return trim([polynomial.nth(i) for i in range(polynomial.degree() + 1)])


def coefficient_rows(poly, variables):
    rows, constants = [], []
    zero = {variable: 0 for variable in variables}
    for coefficient in poly:
        row = [sp.expand(sp.diff(coefficient, variable)) for variable in variables]
        constant = sp.expand(coefficient.subs(zero))
        reconstruction = constant + sum(row[i] * variables[i]
                                        for i in range(len(variables)))
        if sp.expand(coefficient - reconstruction) != 0:
            raise AssertionError("coefficient is not affine-linear in the new slices")
        if any(sp.sympify(value).free_symbols.intersection(variables) for value in row):
            raise AssertionError("new-slice matrix depends on new variables")
        if constant != 0 or any(value != 0 for value in row):
            rows.append(row)
            constants.append(constant)
    return rows, constants


def direct_45_55_system(data, kappa, u5_count=12, w34_degrees=range(7, 67)):
    """Assemble the 45/55 affine block columnwise, without symbolic slices.

    The new coefficients enter linearly. Building one column at a time avoids
    expanding the full expressions in all 72 new slice symbols and then
    differentiating those expressions coefficient-by-coefficient.
    """
    U, h = data["U"], data["h"]
    S, B, U3, U4 = data["S"], data["B"], data["U3"], data["U4"]
    A, A2, A3, A4 = data["A"], data["A2"], data["A3"], data["A4"]
    q24, q23, q22, q21 = (data[name] for name in ("q24", "q23", "q22", "q21"))
    g = data["g"]
    W38, W37, W36, W35 = (data[name] for name in ("W38", "W37", "W36", "W35"))

    q20_binomial_base = add(
        scale(sp.Rational(10, 9), add(mul(mul(h, S), U4), mul(B, U3))),
        scale(sp.Rational(-5, 27), add(mul(pow_poly(S, 2), U3),
                                        mul(S, pow_poly(B, 2)))),
        scale(sp.Rational(20, 243), mul(pow_poly(S, 3), B)),
        scale(sp.Rational(-7, 729), pow_poly(S, 5)),
    )
    q20_base = add(shift(q20_binomial_base, 5),
                   scale(kappa, pow_poly(U, 4)))
    p4_base = add(
        scale(12, mul(pow_poly(U, 6), add(mul(A, A4), mul(A2, A3)))),
        scale(12, mul(pow_poly(U, 3), add(mul(pow_poly(A, 2), A3),
                                          mul(A, pow_poly(A2, 2))))),
        scale(4, mul(pow_poly(A, 3), A2)),
    )

    p_base = add(
        scale(14, mul(A, deriv(W35))),
        scale(-35, mul(deriv(A), W35)),
        scale(13, mul(A2, deriv(W36))),
        scale(-36, mul(deriv(A2), W36)),
        scale(12, mul(A3, deriv(W37))),
        scale(-37, mul(deriv(A3), W37)),
        scale(11, mul(A4, deriv(W38))),
        scale(-38, mul(deriv(A4), W38)),
        scale(-6, mul(pow_poly(U, 5), q20_base)),
        scale(-6, mul(q24, q21)),
        scale(-6, mul(q23, q22)),
        scale(3 * kappa, pow_poly(U, 9)),
    )
    q_base = add(
        scale(24, mul(q24, deriv(W35))),
        scale(-35, mul(deriv(q24), W35)),
        scale(23, mul(q23, deriv(W36))),
        scale(-36, mul(deriv(q23), W36)),
        scale(22, mul(q22, deriv(W37))),
        scale(-37, mul(deriv(q22), W37)),
        scale(21, mul(q21, deriv(W38))),
        scale(-38, mul(deriv(q21), W38)),
        scale(20, mul(q20_base, deriv(g))),
        scale(-39, mul(deriv(q20_base), g)),
        scale(-5, p4_base),
        scale(-9 * kappa, pow_poly(U, 11)),
    )

    p_columns, q_columns = [], []
    for degree in range(u5_count):
        u5_basis = [sp.Integer(0)] * degree + [sp.Integer(1)]
        a5_basis = shift(u5_basis, 3)
        q20_basis = shift(scale(
            sp.Rational(5, 3), mul(pow_poly(h, 2), u5_basis)
        ), 5)
        p_columns.append(add(
            scale(10, mul(a5_basis, deriv(g))),
            scale(-39, mul(deriv(a5_basis), g)),
            scale(-6, mul(pow_poly(U, 5), q20_basis)),
        ))
        q_columns.append(add(
            scale(20, mul(q20_basis, deriv(g))),
            scale(-39, mul(deriv(q20_basis), g)),
            scale(-20, mul(pow_poly(U, 9), a5_basis)),
        ))

    for degree in w34_degrees:
        w34_basis = [sp.Integer(0)] * degree + [sp.Integer(1)]
        p_columns.append(add(
            scale(15, mul(pow_poly(U, 3), deriv(w34_basis))),
            scale(-102, mul(pow_poly(U, 2), mul(deriv(U), w34_basis))),
        ))
        q_columns.append(add(
            scale(25, mul(pow_poly(U, 5), deriv(w34_basis))),
            scale(-170, mul(pow_poly(U, 4), mul(deriv(U), w34_basis))),
        ))

    def rows_for(forcing_poly, column_polys):
        count = max([len(forcing_poly), *(len(poly) for poly in column_polys)])
        rows, constants = [], []
        for degree in range(count):
            constant = sp.expand(forcing_poly[degree]
                                 if degree < len(forcing_poly) else 0)
            row = [sp.expand(poly[degree] if degree < len(poly) else 0)
                   for poly in column_polys]
            if constant != 0 or any(value != 0 for value in row):
                rows.append(row)
                constants.append(constant)
        return rows, constants

    p_rows, p_forcing = rows_for(p_base, p_columns)
    q_rows, q_forcing = rows_for(q_base, q_columns)
    matrix = sp.Matrix([*p_rows, *q_rows])
    forcing = [*p_forcing, *q_forcing]
    return matrix, forcing, len(p_rows)


def exact_affine_solution_with_constraints(
    M: sp.Matrix,
    forcing: list[sp.Expr],
    field,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
    tag: str,
):
    """Parametrize pivot rows and retain, rather than discard, compatibility rows."""
    if field != QQ:
        return quadratic_affine_solution(
            M, forcing, field, variables, parameters, tag
        )
    domain = DomainMatrix.from_Matrix(M).convert_to(field)
    rref, pivot_columns = domain.rref()
    transposed = DomainMatrix.from_Matrix(M.T).convert_to(field)
    transposed_rref, pivot_rows = transposed.rref()
    rank = len(pivot_columns)
    if len(pivot_rows) != rank:
        raise AssertionError(f"inconsistent exact rank profile for {tag}")
    left_kernel = nullspace_from_rref(transposed_rref, pivot_rows, M.rows)
    nonpivot_columns = [j for j in range(M.cols) if j not in pivot_columns]
    free = sp.symbols(f"d0:{len(nonpivot_columns)}")
    pivot_block = M.extract(pivot_rows, pivot_columns)
    free_block = M.extract(pivot_rows, nonpivot_columns)
    forcing_block = sp.Matrix([forcing[i] for i in pivot_rows])
    pivot_domain = DomainMatrix.from_Matrix(pivot_block).convert_to(field)
    inverse_numerator, inverse_denominator = pivot_domain.inv_den()
    denominator = field.to_sympy(inverse_denominator)
    pivot_rhs = forcing_block
    if free:
        pivot_rhs = pivot_rhs + free_block * sp.Matrix(free)
    pivot_values = -(inverse_numerator.to_Matrix() * pivot_rhs) / denominator
    solution = [sp.Integer(0)] * M.cols
    for column, symbol in zip(nonpivot_columns, free):
        solution[column] = symbol
    for column, value in zip(pivot_columns, pivot_values):
        solution[column] = sp.expand(value)
    constraints = residual_polynomials(
        M,
        forcing,
        field,
        parameters,
        f"compatibility {tag}",
        left_kernel=left_kernel,
        rank=rank,
    )
    print(
        f"{tag}_affine_profile rank={rank} fibre={len(free)} "
        f"compatibility_constraints={len(constraints)}",
        flush=True,
    )
    return solution, free, rank, constraints


def first_two_blocks_signature(tag, h_expr, V_expr, field):
    """Fingerprint the exact inputs that determine the cached 48/58,47/57 lift."""
    digest = hashlib.sha256()
    digest.update(tag.encode("utf-8"))
    digest.update(sp.srepr(h_expr).encode("utf-8"))
    digest.update(sp.srepr(V_expr).encode("utf-8"))
    digest.update(str(field).encode("utf-8"))
    module_names = (
        "f2_75_125_power_relation_fourthbands",
        "f2_75_125_power_relation_pq_nextbands",
        "f2_75_125_power_relation_pq_thirdbands",
    )
    for module_name in module_names:
        module = sys.modules[module_name]
        source_path = Path(module.__file__).resolve()
        digest.update(source_path.name.encode("utf-8"))
        digest.update(source_path.read_bytes())
    return digest.hexdigest()


def direct_46_56_system(
    U, h, S, B, U3, A, A2, A3, q24, q23, q22, g, f38, f37, ell,
    W35_coeffs, U4_coeffs,
):
    """Assemble the 46/56 matrix columnwise instead of differentiating a
    fully expanded expression in all 73 new variables.

    Only the forcing depends on carried parameters.  The new-slice columns
    are generated directly from the exact band convolution.
    """
    U2, U3p, U4p, U5p, U6p, U9p = (
        pow_poly(U, exponent) for exponent in (2, 3, 4, 5, 6, 9)
    )
    Up = deriv(U)
    h2 = pow_poly(h, 2)
    q21_0_inner = add(
        scale(54, mul(S, U3)),
        scale(27, pow_poly(B, 2)),
        scale(-9, mul(pow_poly(S, 2), B)),
        pow_poly(S, 4),
    )
    q21_0 = scale(sp.Rational(5, 243), shift(mul(h, q21_0_inner), 5))

    p_columns, q_columns = [], []
    U2Up = mul(U2, Up)
    U4Up = mul(U4p, Up)
    for i in range(7, 68):
        p_columns.append(add(
            scale(15 * i, shift(U3p, i - 1)),
            scale(-105, shift(U2Up, i)),
        ))
        q_columns.append(add(
            scale(25 * i, shift(U5p, i - 1)),
            scale(-175, shift(U4Up, i)),
        ))

    h2_gprime = mul(h2, deriv(g))
    h2_g = mul(h2, g)
    h2prime_g = mul(deriv(h2), g)
    for j in range(12):
        p_columns.append(add(
            scale(11, shift(deriv(g), j + 3)),
            scale(-39 * (j + 3), shift(g, j + 2)),
            scale(-10, shift(mul(U5p, h2), j + 5)),
        ))
        q_columns.append(add(
            scale(35, shift(h2_gprime, j + 5)),
            scale(-65, shift(h2prime_g, j + 5)),
            scale(-65 * (j + 5), shift(h2_g, j + 4)),
            scale(-20, shift(U9p, j + 3)),
        ))

    p_forcing = add(
        scale(14, mul(A, deriv(ell))),
        scale(-36, mul(deriv(A), ell)),
        scale(13, mul(A2, deriv(f37))),
        scale(-37, mul(deriv(A2), f37)),
        scale(12, mul(A3, deriv(f38))),
        scale(-38, mul(deriv(A3), f38)),
        scale(-6, mul(U5p, q21_0)),
        scale(-6, mul(q24, q22)),
        scale(-3, pow_poly(q23, 2)),
    )
    q_forcing = add(
        scale(24, mul(q24, deriv(ell))),
        scale(-36, mul(deriv(q24), ell)),
        scale(23, mul(q23, deriv(f37))),
        scale(-37, mul(deriv(q23), f37)),
        scale(22, mul(q22, deriv(f38))),
        scale(-38, mul(deriv(q22), f38)),
        scale(21, mul(q21_0, deriv(g))),
        scale(-39, mul(deriv(q21_0), g)),
        scale(-60, mul(U6p, mul(A, A3))),
        scale(-30, mul(U6p, pow_poly(A2, 2))),
        scale(-60, mul(U3p, mul(pow_poly(A, 2), A2))),
        scale(-5, pow_poly(A, 4)),
    )

    def block_rows(columns, forcing):
        height = max([len(forcing), *(len(column) for column in columns)])
        rows, constants = [], []
        for degree in range(height):
            row = [column[degree] if degree < len(column) else sp.Integer(0)
                   for column in columns]
            constant = forcing[degree] if degree < len(forcing) else sp.Integer(0)
            if any(value != 0 for value in row) or constant != 0:
                rows.append(row)
                constants.append(constant)
        return rows, constants

    p_rows, p_rhs = block_rows(p_columns, p_forcing)
    q_rows, q_rhs = block_rows(q_columns, q_forcing)
    return (
        sp.Matrix([*p_rows, *q_rows]),
        [*p_rhs, *q_rhs],
        (*W35_coeffs, *U4_coeffs),
        len(p_rows),
    )


def verify_new_band_formulas() -> None:
    """Check the q20 binomial term and corrected 45/55 band convolutions."""
    t, eps = sp.symbols("t eps")
    U, A, A2, A3, A4, A5 = sp.symbols(
        "U A A2 A3 A4 A5", cls=sp.Function
    )
    q24, q23, q22, q21, q20 = sp.symbols(
        "q24 q23 q22 q21 q20", cls=sp.Function
    )
    W39, W38, W37, W36, W35, W34 = sp.symbols(
        "W39 W38 W37 W36 W35 W34", cls=sp.Function
    )
    kappa = sp.symbols("kappa")
    U, A, A2, A3, A4, A5 = (f(t) for f in (U, A, A2, A3, A4, A5))
    q24, q23, q22, q21, q20 = (f(t) for f in (q24, q23, q22, q21, q20))
    W39, W38, W37, W36, W35, W34 = (
        f(t) for f in (W39, W38, W37, W36, W35, W34)
    )

    def bracket_band(left, right, output_weight):
        terms = []
        for left_weight, p in left.items():
            for right_weight, q in right.items():
                if left_weight + right_weight - 4 == output_weight:
                    terms.append(left_weight * p * sp.diff(q, t)
                                 - right_weight * sp.diff(p, t) * q)
        return sum(terms, sp.Integer(0))

    p_bands = {15: U**3, 14: A, 13: A2, 12: A3, 11: A4, 10: A5}
    q_bands = {25: U**5, 24: q24, 23: q23, 22: q22, 21: q21, 20: q20}
    w_bands = {39: W39, 38: W38, 37: W37, 36: W36, 35: W35, 34: W34}
    p_series = sum(p_bands[15 - i] * eps**i for i in range(6))
    q_series = sum(q_bands[25 - i] * eps**i for i in range(6))

    p4_drop5 = (
        4 * U**9 * A5
        + 12 * U**6 * (A * A4 + A2 * A3)
        + 12 * U**3 * (A**2 * A3 + A * A2**2)
        + 4 * A**3 * A2
    )
    direct_p = (
        bracket_band(p_bands, w_bands, 45)
        - 3 * sp.expand(q_series**2).coeff(eps, 5)
        + 3 * kappa * U**9
    )
    expected_p = (
        15 * U**3 * sp.diff(W34, t) - 102 * U**2 * sp.diff(U, t) * W34
        + 14 * A * sp.diff(W35, t) - 35 * sp.diff(A, t) * W35
        + 13 * A2 * sp.diff(W36, t) - 36 * sp.diff(A2, t) * W36
        + 12 * A3 * sp.diff(W37, t) - 37 * sp.diff(A3, t) * W37
        + 11 * A4 * sp.diff(W38, t) - 38 * sp.diff(A4, t) * W38
        + 10 * A5 * sp.diff(W39, t) - 39 * sp.diff(A5, t) * W39
        - 6 * U**5 * q20 - 6 * q24 * q21 - 6 * q23 * q22
        + 3 * kappa * U**9
    )
    direct_q = (
        bracket_band(q_bands, w_bands, 55)
        - 5 * sp.expand(p_series**4).coeff(eps, 5)
        - 9 * kappa * U**11
    )
    expected_q = (
        25 * U**5 * sp.diff(W34, t) - 170 * U**4 * sp.diff(U, t) * W34
        + 24 * q24 * sp.diff(W35, t) - 35 * sp.diff(q24, t) * W35
        + 23 * q23 * sp.diff(W36, t) - 36 * sp.diff(q23, t) * W36
        + 22 * q22 * sp.diff(W37, t) - 37 * sp.diff(q22, t) * W37
        + 21 * q21 * sp.diff(W38, t) - 38 * sp.diff(q21, t) * W38
        + 20 * q20 * sp.diff(W39, t) - 39 * sp.diff(q20, t) * W39
        - 5 * p4_drop5 - 9 * kappa * U**11
    )
    assert sp.simplify(direct_p - expected_p) == 0
    assert sp.simplify(direct_q - expected_q) == 0

    h, S, B, U3, U4, U5 = sp.symbols("h S B U3 U4 U5")
    eps2 = sp.symbols("eps2")
    delta = (
        S / h * eps2 + B / h**2 * eps2**2 + U3 / h**3 * eps2**3
        + U4 / h**3 * eps2**4 + U5 / h**3 * eps2**5
    )
    formal = h**5 * sum(
        sp.binomial(sp.Rational(5, 3), k) * delta**k
        for k in range(1, 6)
    )
    q20_binomial = sp.expand(formal).coeff(eps2, 5)
    q20_formula = (
        sp.Rational(5, 3) * h**2 * U5
        + sp.Rational(10, 9) * (h * S * U4 + B * U3)
        - sp.Rational(5, 27) * (S**2 * U3 + S * B**2)
        + sp.Rational(20, 243) * S**3 * B
        - sp.Rational(7, 729) * S**5
    )
    assert sp.cancel(q20_binomial - q20_formula) == 0
    print("fifth_band_formula_audit=PASS (q20, corrected weights 45/55)", flush=True)


def lift_through_fourth_block(tag: str, h_expr, V_expr, field,
                              checkpoint_path: Path | None = None):
    """Carry the full affine family through weights 46/56."""
    t = sp.symbols("t")
    if checkpoint_path is None:
        data = lift_first_two_blocks(tag, h_expr, V_expr, field)
    else:
        checkpoint_path = Path(checkpoint_path)
        signature = first_two_blocks_signature(tag, h_expr, V_expr, field)
        if checkpoint_path.exists():
            with checkpoint_path.open("rb") as handle:
                checkpoint = pickle.load(handle)
            if checkpoint.get("signature") != signature:
                raise ValueError(f"stale first-two-block checkpoint: {checkpoint_path}")
            data = checkpoint["data"]
            print(f"first_two_blocks_checkpoint=loaded path={checkpoint_path}", flush=True)
        else:
            data = lift_first_two_blocks(tag, h_expr, V_expr, field)
            checkpoint = {"signature": signature, "tag": tag, "data": data}
            temporary_path = checkpoint_path.with_suffix(checkpoint_path.suffix + ".tmp")
            with temporary_path.open("wb") as handle:
                pickle.dump(checkpoint, handle, protocol=pickle.HIGHEST_PROTOCOL)
                handle.flush()
            temporary_path.replace(checkpoint_path)
            print(
                f"first_two_blocks_checkpoint=saved path={checkpoint_path} "
                f"bytes={checkpoint_path.stat().st_size}",
                flush=True,
            )
    U, h = from_expr(data["U"], t), from_expr(h_expr, t)
    S, B = from_expr(data["S"], t), from_expr(data["B"], t)
    U2, U3 = from_expr(data["U2"], t), from_expr(data["U3"], t)
    A, A2, A3 = from_expr(data["A"], t), from_expr(data["A2"], t), from_expr(data["A3"], t)
    q24, q23, q22 = (from_expr(data[name], t) for name in ("q24", "q23", "q22"))
    g, f38 = from_expr(data["g"], t), from_expr(data["f38"], t)
    f37, ell = from_expr(data["f37"], t), from_expr(data["ell"], t)
    parameters = data["parameters"]

    W35_coeffs = sp.symbols("W35_7:68")
    U4_coeffs = sp.symbols("U4_0:12")
    W35 = [sp.Integer(0)] * 7 + list(W35_coeffs)
    U4 = list(U4_coeffs)
    A4 = shift(U4, 3)
    q21 = add(
        scale(sp.Rational(5, 3), shift(mul(pow_poly(h, 2), U4), 5)),
        scale(sp.Rational(5, 243), shift(mul(h, add(
            scale(54, mul(S, U3)),
            scale(27, pow_poly(B, 2)),
            scale(-9, mul(pow_poly(S, 2), B)),
            pow_poly(S, 4),
        )), 5)),
    )

    new_variables = (*W35_coeffs, *U4_coeffs)
    matrix, forcing, returned_variables, _p_row_count = direct_46_56_system(
        U, h, S, B, U3, A, A2, A3, q24, q23, q22, g, f38, f37, ell,
        W35_coeffs, U4_coeffs,
    )
    if returned_variables != new_variables:
        raise AssertionError("direct 46/56 assembler changed the variable order")
    solution, raw_fibre, rank, compatibility_constraints = exact_affine_solution_with_constraints(
        matrix, forcing, field, new_variables, parameters, f"46/56 {tag}"
    )

    fibre = sp.symbols(f"H0:{len(raw_fibre)}")
    renaming = dict(zip(raw_fibre, fibre))
    solution = [sp.expand(value.subs(renaming)) for value in solution]
    substitution = dict(zip(new_variables, solution))
    W35_solved = [sp.expand(value.subs(substitution)) for value in W35]
    U4_solved = [sp.expand(value.subs(substitution)) for value in U4]
    q21_solved = [sp.expand(value.subs(substitution)) for value in q21]
    parameters = (*parameters, *fibre)
    print(
        f"{tag}_46_56_lift=PASS rank={rank} fibre={len(fibre)} "
        f"carried_parameters={len(parameters)}",
        flush=True,
    )
    return {
        "t": t, "U": U, "h": h, "S": S, "B": B,
        "U2": U2, "U3": U3, "U4": U4_solved,
        "A": A, "A2": A2, "A3": A3, "A4": shift(U4_solved, 3),
        "q24": q24, "q23": q23, "q22": q22, "q21": q21_solved,
        "g": g, "W38": f38, "W37": f37, "W36": ell,
        "W35": W35_solved, "parameters": parameters,
        "compatibility_constraints": compatibility_constraints,
    }


def audit_fifth_block(tag: str, h_expr, V_expr, field,
                      checkpoint_path: Path | None = None):
    data = lift_through_fourth_block(tag, h_expr, V_expr, field, checkpoint_path)
    parameters = data["parameters"]

    kappa = sp.Symbol(f"kappa70_{tag}")
    U5_coeffs = sp.symbols(f"U5_{tag}_0:12")
    W34_coeffs = sp.symbols(f"W34_{tag}_7:67")
    new_variables = (*U5_coeffs, *W34_coeffs)
    matrix, forcing, p_row_count = direct_45_55_system(data, kappa)
    if matrix.cols != len(new_variables):
        raise AssertionError("direct 45/55 assembler changed the variable order")
    all_parameters = (*parameters, kappa)
    if any(entry.free_symbols.intersection(all_parameters) for row in matrix.tolist() for entry in row):
        raise AssertionError("45/55 coefficient matrix depends on carried parameters")
    print(
        f"{tag}_45_55_column_system_ready rows={matrix.rows} "
        f"P_rows={p_row_count} variables={matrix.cols}",
        flush=True,
    )

    if field == QQ:
        transposed_domain = DomainMatrix.from_Matrix(matrix.T).convert_to(field)
        left_rref, pivots = transposed_domain.rref()
        rank = len(pivots)
        left_kernel = nullspace_from_rref(left_rref, pivots, matrix.rows)
        fifth_band_constraints = residual_polynomials(
            matrix, forcing, field, all_parameters, f"compatibility 45/55 {tag}",
            left_kernel=left_kernel, rank=rank,
        )
    else:
        _, _, rank, fifth_band_constraints = quadratic_affine_solution(
            matrix, forcing, field, new_variables, all_parameters,
            f"compatibility 45/55 {tag}",
        )
    prior_constraints = data["compatibility_constraints"]
    residuals = [*prior_constraints, *fifth_band_constraints]
    print(
        f"branch={tag} 45/55_rows={matrix.rows} variables={matrix.cols} "
        f"rank={rank} fibre={matrix.cols-rank} carried_parameters={len(all_parameters)} "
        f"new_45_55_constraints={len(fifth_band_constraints)} "
        f"prior_46_56_constraints={len(prior_constraints)} "
        f"combined_constraints={len(residuals)}",
        flush=True,
    )
    for i, residual in enumerate(residuals[:12]):
        poly = sp.Poly(residual, *all_parameters, domain=sp.EX)
        print(
            f"residual_{i} total_degree={poly.total_degree()} "
            f"kappa_degree={sp.degree(residual, kappa)} "
            f"factor={sp.factor_terms(residual)}",
            flush=True,
        )
    if len(residuals) > 12:
        print(f"residuals_printed=12 total={len(residuals)}", flush=True)
    return residuals


def main() -> None:
    verify_new_band_formulas()
    t = sp.symbols("t")
    r = t**2 - 3 * t + 3
    h211 = sp.expand((t - 1) ** 2 * r)
    v211 = sp.expand(-sp.Rational(1, 3) * t * (t - 1) * r)
    checkpoint = Path(__file__).with_name(
        "f2_75_125_211_first_two_blocks_46_56.pkl"
    )
    audit_fifth_block("2+1+1", h211, v211, QQ, checkpoint)
    print("scope=exact next necessary bands on one quartic edge branch; not a pair exclusion")


if __name__ == "__main__":
    main()
