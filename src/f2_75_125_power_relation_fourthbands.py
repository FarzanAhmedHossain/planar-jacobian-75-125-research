"""Couple the next [P,W*] and [Q,W*] bands exactly.

This extends the exact affine-lift calculation through the weight-46 row of
[P,W*] and the weight-56 row of [Q,W*]. It is conditional on the F2 outer
support envelope and the edge power-relation construction; the allowed W*
slice is deliberately support-relaxed. An empty compatibility ideal would
therefore exclude the edge branch, while a nonempty result is only a failure
of this particular rigidity attempt.
"""

from __future__ import annotations

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import QQ

from f2_75_125_power_relation_pq_nextbands import coefficient_system, weight49_parameters
from f2_75_125_power_relation_pq_thirdbands import (
    exact_affine_solution,
    residual_polynomials,
)


def verify_band_formulae() -> None:
    """Independently compare the displayed band rows with convolution rules."""
    t, eps = sp.symbols("t eps")
    U, A, A2, A3, A4 = sp.symbols("U A A2 A3 A4", cls=sp.Function)
    q24, q23, q22, q21 = sp.symbols("q24 q23 q22 q21", cls=sp.Function)
    g, f, k, ell, m = sp.symbols("g f k ell m", cls=sp.Function)
    U, A, A2, A3, A4 = (x(t) for x in (U, A, A2, A3, A4))
    q24, q23, q22, q21 = (x(t) for x in (q24, q23, q22, q21))
    g, f, k, ell, m = (x(t) for x in (g, f, k, ell, m))

    def bracket_band(left, right, weight):
        terms = []
        for p_weight, p in left.items():
            for q_weight, q in right.items():
                if p_weight + q_weight - 4 == weight:
                    terms.append(p_weight * p * sp.diff(q, t)
                                 - q_weight * sp.diff(p, t) * q)
        return sum(terms, sp.Integer(0))

    p_bands = {15: U**3, 14: A, 13: A2, 12: A3, 11: A4}
    q_bands = {25: U**5, 24: q24, 23: q23, 22: q22, 21: q21}
    w_bands = {39: g, 38: f, 37: k, 36: ell, 35: m}

    q_series = sum(q_bands[25 - i] * eps**i for i in range(5))
    p_series = sum(p_bands[15 - i] * eps**i for i in range(5))
    direct_p = bracket_band(p_bands, w_bands, 46) - 3 * sp.expand(q_series**2).coeff(eps, 4)
    expected_p = (
        15 * U**3 * sp.diff(m, t) - 105 * U**2 * sp.diff(U, t) * m
        + 14 * A * sp.diff(ell, t) - 36 * sp.diff(A, t) * ell
        + 13 * A2 * sp.diff(k, t) - 37 * sp.diff(A2, t) * k
        + 12 * A3 * sp.diff(f, t) - 38 * sp.diff(A3, t) * f
        + 11 * A4 * sp.diff(g, t) - 39 * sp.diff(A4, t) * g
        - 6 * U**5 * q21 - 6 * q24 * q22 - 3 * q23**2
    )
    direct_q = bracket_band(q_bands, w_bands, 56) - 5 * sp.expand(p_series**4).coeff(eps, 4)
    expected_q = (
        5 * U**4 * (5 * U * sp.diff(m, t) - 35 * sp.diff(U, t) * m)
        + 24 * q24 * sp.diff(ell, t) - 36 * sp.diff(q24, t) * ell
        + 23 * q23 * sp.diff(k, t) - 37 * sp.diff(q23, t) * k
        + 22 * q22 * sp.diff(f, t) - 38 * sp.diff(q22, t) * f
        + 21 * q21 * sp.diff(g, t) - 39 * sp.diff(q21, t) * g
        - 20 * U**9 * A4 - 60 * U**6 * A * A3 - 30 * U**6 * A2**2
        - 60 * U**3 * A**2 * A2 - 5 * A**4
    )
    assert sp.simplify(direct_p - expected_p) == 0
    assert sp.simplify(direct_q - expected_q) == 0

    h, S, B, U3, U4 = sp.symbols("h S B U3 U4")
    u1, u2 = h**2 * S, h * B
    binomial_v4 = (
        sp.Rational(5, 3) * h**2 * U4
        + sp.Rational(5, 9) * (2 * u1 * U3 + u2**2) / h
        - sp.Rational(5, 81) * (3 * u1**2 * u2) / h**4
        + sp.Rational(5, 243) * u1**4 / h**7
    )
    reduced_v4 = (
        sp.Rational(5, 3) * h**2 * U4
        + sp.Rational(5, 243) * h * (54 * S * U3 + 27 * B**2 - 9 * S**2 * B + S**4)
    )
    assert sp.cancel(binomial_v4 - reduced_v4) == 0
    print("formula_audit=PASS (46/56 band convolution; binomial V4)", flush=True)


def lift_first_two_blocks(tag: str, h: sp.Expr, V: sp.Expr, field):
    """Return solved slices through W*_36 and U_3, plus free parameters."""
    t = sp.symbols("t")
    U = sp.expand(t * h)
    old_parameters, S, f38, _ = weight49_parameters(h, V, field, tag)

    B_coeffs = sp.symbols("B0:8")
    f37_coeffs = sp.symbols("F8:68")
    B = sum(B_coeffs[i] * t**i for i in range(8))
    f37 = sum(f37_coeffs[i - 8] * t**i for i in range(8, 68))
    U1 = sp.expand(h**2 * S)
    U2 = sp.expand(h * B)
    A = sp.expand(t**3 * U1)
    A2 = sp.expand(t**3 * U2)
    q24 = sp.expand(sp.Rational(5, 3) * t**5 * h**2 * U1)
    q23 = sp.expand(
        sp.Rational(5, 3) * t**5 * h**2 * U2
        + sp.Rational(5, 9) * t**5 * h**3 * S**2
    )
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
    Mp, bp = coefficient_system(p48, t, old_parameters, old_variables, f"P48 {tag}")
    Mq, bq = coefficient_system(q58, t, old_parameters, old_variables, f"Q58 {tag}")
    M_old = sp.Matrix.vstack(Mp, Mq)
    b_old = [*bp, *bq]
    old_solution, raw_fibre1, rank_old = exact_affine_solution(
        M_old, b_old, field, old_variables, old_parameters, f"48/58 {tag}"
    )
    fibre1 = sp.symbols(f"D0:{len(raw_fibre1)}")
    rename1 = dict(zip(raw_fibre1, fibre1))
    old_solution = [sp.expand(value.subs(rename1)) for value in old_solution]
    parameters1 = (*old_parameters, *fibre1)
    substitution1 = dict(zip(old_variables, old_solution))
    B_solved = sp.expand(B.subs(substitution1))
    f37_solved = sp.expand(f37.subs(substitution1))
    U2_solved = sp.expand(h * B_solved)
    A2_solved = sp.expand(t**3 * U2_solved)
    q23_solved = sp.expand(
        sp.Rational(5, 3) * t**5 * h**2 * U2_solved
        + sp.Rational(5, 9) * t**5 * h**3 * S**2
    )
    U3_coeffs = sp.symbols("C0:12")
    U3 = sum(U3_coeffs[i] * t**i for i in range(12))
    q22 = sp.expand(
        t**5 * h**2 * (
            sp.Rational(5, 3) * U3
            + sp.Rational(10, 9) * S * B_solved
            - sp.Rational(5, 81) * S**3
        )
    )

    A3 = sp.expand(t**3 * U3)
    L_coeffs = sp.symbols("L8:68")
    ell = sum(L_coeffs[i - 8] * t**i for i in range(8, 68))

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

    variables2 = (*L_coeffs, *U3_coeffs)
    M1, b1 = coefficient_system(p47, t, parameters1, variables2, f"P47 {tag}")
    M2, b2 = coefficient_system(q57, t, parameters1, variables2, f"Q57 {tag}")
    M_next = sp.Matrix.vstack(M1, M2)
    b_next = [*b1, *b2]
    solution2, raw_fibre2, rank_next = exact_affine_solution(
        M_next, b_next, field, variables2, parameters1, f"47/57 {tag}"
    )

    # Give the second block's free variables fresh names; the affine solver's
    # default d0,d1,... names would otherwise collide with the first block.
    fibre2 = sp.symbols(f"E0:{len(raw_fibre2)}")
    rename = dict(zip(raw_fibre2, fibre2))
    solution2 = [sp.expand(value.subs(rename)) for value in solution2]
    parameters2 = (*parameters1, *fibre2)
    substitution2 = dict(zip(variables2, solution2))
    ell_solved = sp.expand(ell.subs(substitution2))
    U3_solved = sp.expand(U3.subs(substitution2))
    return {
        "t": t,
        "U": U,
        "S": S,
        "U1": U1,
        "U2": U2_solved,
        "U3": U3_solved,
        "A": A,
        "A2": A2_solved,
        "A3": sp.expand(t**3 * U3_solved),
        "q24": q24,
        "q23": q23_solved,
        "q22": sp.expand(q22.subs(substitution2)),
        "g": g,
        "f38": f38,
        "f37": f37_solved,
        "ell": ell_solved,
        "B": B_solved,
        "parameters": parameters2,
        "rank48": rank_old,
        "rank47": rank_next,
    }


def audit_fourth_block(tag: str, h: sp.Expr, V: sp.Expr, field) -> None:
    print(f"branch_start={tag}", flush=True)
    data = lift_first_two_blocks(tag, h, V, field)
    t, U = data["t"], data["U"]
    S, U2, U3 = data["S"], data["U2"], data["U3"]
    A, A2, A3 = data["A"], data["A2"], data["A3"]
    q24, q23, q22 = data["q24"], data["q23"], data["q22"]
    g, f38, f37, ell = data["g"], data["f38"], data["f37"], data["ell"]
    parameters = data["parameters"]

    D_coeffs = sp.symbols("U4_0:12")
    U4 = sum(D_coeffs[i] * t**i for i in range(12))
    A4 = sp.expand(t**3 * U4)
    q21 = sp.expand(
        t**5 * (
            sp.Rational(5, 3) * h**2 * U4
            + sp.Rational(5, 243) * h * (
                54 * S * U3 + 27 * data["B"]**2 - 9 * S**2 * data["B"] + S**4
            )
        )
    )

    M_coeffs = sp.symbols("W35_7:68")
    W35 = sum(M_coeffs[i - 7] * t**i for i in range(7, 68))

    p46 = (
        15 * U**3 * sp.diff(W35, t)
        - 105 * U**2 * sp.diff(U, t) * W35
        + 14 * A * sp.diff(ell, t)
        - 36 * sp.diff(A, t) * ell
        + 13 * A2 * sp.diff(f37, t)
        - 37 * sp.diff(A2, t) * f37
        + 12 * A3 * sp.diff(f38, t)
        - 38 * sp.diff(A3, t) * f38
        + 11 * A4 * sp.diff(g, t)
        - 39 * sp.diff(A4, t) * g
        - 6 * U**5 * q21
        - 6 * q24 * q22
        - 3 * q23**2
    )

    q56 = (
        5 * U**4 * (5 * U * sp.diff(W35, t) - 35 * sp.diff(U, t) * W35)
        + 24 * q24 * sp.diff(ell, t)
        - 36 * sp.diff(q24, t) * ell
        + 23 * q23 * sp.diff(f37, t)
        - 37 * sp.diff(q23, t) * f37
        + 22 * q22 * sp.diff(f38, t)
        - 38 * sp.diff(q22, t) * f38
        + 21 * q21 * sp.diff(g, t)
        - 39 * sp.diff(q21, t) * g
        - 20 * U**9 * A4
        - 60 * U**6 * A * A3
        - 30 * U**6 * A2**2
        - 60 * U**3 * A**2 * A2
        - 5 * A**4
    )

    new_variables = (*M_coeffs, *D_coeffs)
    if set(new_variables).intersection(parameters):
        raise AssertionError("new slice variables collide with carried parameters")
    Mp, bp = coefficient_system(p46, t, parameters, new_variables, f"P46 {tag}")
    Mq, bq = coefficient_system(q56, t, parameters, new_variables, f"Q56 {tag}")
    M = sp.Matrix.vstack(Mp, Mq)
    forcing = [*bp, *bq]
    residuals = residual_polynomials(M, forcing, field, parameters, f"46/56 {tag}")
    rank = DomainMatrix.from_Matrix(M).convert_to(field).rank()
    print(
        f"branch={tag} rows={M.rows} variables={M.cols} rank={rank} "
        f"parameter_count={len(parameters)} new_fibre_dimension={M.cols-rank}",
        flush=True,
    )
    if residuals:
        for i, residual in enumerate(residuals):
            print(f"{tag}_residual_{i}={residual}", flush=True)
        if len(parameters) <= 12 and len(residuals) <= 40:
            basis = sp.groebner(residuals, *parameters, domain=field)
            print(f"{tag}_groebner_basis={list(basis)}", flush=True)
            print(f"{tag}_residual_ideal_unit={basis.contains(1)}", flush=True)
        else:
            print("groebner_skipped=parameter_or_generator_budget", flush=True)
    else:
        print(f"{tag}_residual_ideal=zero", flush=True)


def main() -> None:
    verify_band_formulae()
    t = sp.symbols("t")
    sqrt5 = sp.sqrt(5)
    alpha = (3 + sqrt5) / 2
    h22 = sp.expand((t - 1) ** 2 * (t - alpha) ** 2)
    m, n = (4 * alpha - 11) / 3, 3 - alpha
    v22 = sp.expand(t * (t - 1) * (t - alpha) * (m * t + n))
    audit_fourth_block("2+2", h22, v22, QQ.algebraic_field(sqrt5))

    r = t**2 - 3 * t + 3
    h211 = sp.expand((t - 1) ** 2 * r)
    v211 = sp.expand(-sp.Rational(1, 3) * t * (t - 1) * r)
    audit_fourth_block("2+1+1", h211, v211, QQ)
    print("scope=weight-46/56 necessary bands in the relaxed F2 edge model; not a pair exclusion")


if __name__ == "__main__":
    main()
