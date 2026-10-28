"""Couple the F2 terminal ordinary-multiplicity-six selector to edge bands.

This is an exact, support-relaxed necessary-condition calculation for the
2+1+1 quartic edge.  The source-facing implication (F2 chain -> terminal
half-space -> ordinary multiplicity six at the distinguished root) is stated
separately from this finite algebra.  The calculation retains theta=0 and
does not assume the fifth slice U5 is globally polynomial.

Only U1,...,U4 and the corrected-power-relation bands through weights
48/58, 47/57, and 46/56 are used.  A nonunit residual ideal is not a branch
exclusion, and this file does not prove the source chain attachment.
"""

from __future__ import annotations

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import QQ

from f2_75_125_power_relation_pq_nextbands import coefficient_system
from f2_75_125_power_relation_pq_thirdbands import residual_polynomials
from f2_75_125_power_relation_second_band import build_matrix


def taylor_coefficient(poly: sp.Expr, t: sp.Symbol, root: sp.Expr, order: int):
    return sp.factor(sp.diff(poly, t, order).subs(t, root) / sp.factorial(order))


def exact_rank(matrix: sp.Matrix) -> int:
    return DomainMatrix.from_Matrix(matrix).convert_to(QQ).rank()


def solve_block(matrix: sp.Matrix, forcing: list[sp.Expr], parameters,
                variables, tag: str, fibre_prefix: str):
    """Return a pivot-row affine solution and all exact left-kernel residuals.

    Residual equations are retained as parameter constraints; they are not
    silently assumed to vanish.  The pivot-row expression is valid on their
    common zero locus.
    """
    dm = DomainMatrix.from_Matrix(matrix).convert_to(QQ)
    rref, pivots = dm.rref()
    rank = len(pivots)
    tdm = DomainMatrix.from_Matrix(matrix.T).convert_to(QQ)
    _, pivot_rows = tdm.rref()
    if len(pivot_rows) != rank:
        raise AssertionError(f"rank profile mismatch in {tag}")

    left_kernel = tdm.nullspace().to_Matrix()
    residuals = residual_polynomials(
        matrix, forcing, QQ, parameters, tag,
        left_kernel=left_kernel, rank=rank,
    )

    pivot_columns = list(pivots)
    nonpivot_columns = [j for j in range(matrix.cols) if j not in pivot_columns]
    raw_fibre = sp.symbols(f"{fibre_prefix}0:{len(nonpivot_columns)}")
    pivot_block = matrix.extract(pivot_rows, pivot_columns)
    free_block = matrix.extract(pivot_rows, nonpivot_columns)
    forcing_block = sp.Matrix([forcing[i] for i in pivot_rows])
    pivot_domain = DomainMatrix.from_Matrix(pivot_block).convert_to(QQ)
    inverse_numerator, inverse_denominator = pivot_domain.inv_den()
    pivot_values = -(
        inverse_numerator.to_Matrix()
        * (forcing_block + free_block * sp.Matrix(raw_fibre))
    ) / QQ.to_sympy(inverse_denominator)

    solution = [sp.Integer(0)] * matrix.cols
    for column, symbol in zip(nonpivot_columns, raw_fibre):
        solution[column] = symbol
    for column, value in zip(pivot_columns, pivot_values):
        solution[column] = sp.expand(value)

    fibre = sp.symbols(f"{fibre_prefix}free0:{len(raw_fibre)}")
    rename = dict(zip(raw_fibre, fibre))
    solution = [sp.expand(value.subs(rename)) for value in solution]
    print(
        f"{tag}_shape={matrix.rows}x{matrix.cols} rank={rank} "
        f"leftkernel={matrix.rows-rank} fibre={len(fibre)} "
        f"nonzero_residuals={len(residuals)}",
        flush=True,
    )
    for i, residual in enumerate(residuals):
        print(f"{tag}_residual_{i}={sp.factor(residual)}", flush=True)
    return solution, fibre, residuals, rank


def append_equation(matrix_parts, forcing_parts, expression, t, parameters,
                    variables, tag: str):
    M, b = coefficient_system(expression, t, parameters, variables, tag)
    matrix_parts.append(M)
    forcing_parts.extend(b)


def selector_and_first_block():
    t = sp.symbols("t")
    r = t**2 - 3*t + 3
    h = sp.expand((t - 1)**2 * r)
    V = sp.expand(-sp.Rational(1, 3) * t * (t - 1) * r)
    U = sp.expand(t*h)
    theta = sp.symbols("theta211")
    S = sp.expand(theta*(t - 1)*r)
    U1 = sp.expand(h**2*S)
    # The weight-49 two-dimensional selector has this exact v=0 kernel line.
    f38 = sp.expand(-sp.Rational(5, 6)*theta*t**8*(t - 1)**14*r**8)
    M49 = build_matrix(h, V)
    selector_vector = sp.Matrix([
        *[sp.expand(S).coeff(t, j) for j in range(4)],
        *[sp.expand(f38).coeff(t, j) for j in range(8, 68)],
    ])
    if M49.shape != (72, 64) or exact_rank(M49) != 62:
        raise AssertionError("unexpected weight-49 selector rank")
    if M49 * selector_vector != sp.zeros(M49.rows, 1):
        raise AssertionError("v=0 selector is not in the exact weight-49 kernel")

    C0 = taylor_coefficient(h**3, t, 1, 6)
    A1 = taylor_coefficient(U1, t, 1, 5)
    if sp.simplify(C0 - 1) != 0 or sp.simplify(A1 - theta) != 0:
        raise AssertionError("unexpected ordinary sextic normalization")

    # Write U2=h^2*T.  The tangent sextic C(v+A w/(6C))^6 gives
    # [v^4]U2=5*A^2/(12*C), which is T(1)=5*theta^2/12 here.
    Tc = sp.symbols("T48_0:4")
    T = sum(Tc[i]*t**i for i in range(4))
    B = sp.expand(h*T)
    U2 = sp.expand(h*B)

    Fc = sp.symbols("F37_8:68")
    f37 = sum(Fc[i-8]*t**i for i in range(8, 68))
    A = sp.expand(t**3*U1)
    A2 = sp.expand(t**3*U2)
    q24 = sp.expand(sp.Rational(5, 3)*t**5*h**2*U1)
    q23 = sp.expand(
        sp.Rational(5, 3)*t**5*h**2*U2
        + sp.Rational(5, 9)*t**5*h**3*S**2
    )
    g = sp.expand(U**7*V)

    p48 = (
        15*U**3*sp.diff(f37,t) - 111*U**2*sp.diff(U,t)*f37
        + 14*A*sp.diff(f38,t) - 38*sp.diff(A,t)*f38
        + 13*A2*sp.diff(g,t) - 39*sp.diff(A2,t)*g
        - 6*U**5*q23 - 3*q24**2
    )
    q58 = (
        5*U**4*(5*U*sp.diff(f37,t) - 37*sp.diff(U,t)*f37)
        + 24*q24*sp.diff(f38,t) - 38*sp.diff(q24,t)*f38
        + 23*q23*sp.diff(g,t) - 39*sp.diff(q23,t)*g
        - 20*U**9*A2 - 30*U**6*A**2
    )
    variables = (*Fc, *Tc)
    params = (theta,)
    matrices, forcings = [], []
    append_equation(matrices, forcings, p48, t, params, variables, "P48 tangent")
    append_equation(matrices, forcings, q58, t, params, variables, "Q58 tangent")
    append_equation(
        matrices, forcings, T.subs(t, 1)-sp.Rational(5,12)*theta**2,
        t, params, variables, "tangent U2 diagonal",
    )
    matrix = sp.Matrix.vstack(*matrices)
    solution, fibre, residuals, rank = solve_block(
        matrix, forcings, params, variables, "edge211 tangent 48/58", "F48_"
    )
    sub = dict(zip(variables, solution))
    result = {
        "t": t, "r": r, "h": h, "V": V, "U": U, "theta": theta,
        "S": S, "U1": U1, "f38": f38,
        "f37": sp.expand(f37.subs(sub)), "T": sp.expand(T.subs(sub)),
        "B": sp.expand(B.subs(sub)), "U2": sp.expand(U2.subs(sub)),
        "A": A, "A2": sp.expand(A2.subs(sub)), "q24": q24,
        "q23": sp.expand(q23.subs(sub)), "g": g,
        "parameters": (*params, *fibre), "constraints": list(residuals),
        "ranks": [rank],
    }
    assert sp.rem(sp.Poly(result["U1"], t), sp.Poly((t-1)**5, t)).is_zero
    assert sp.rem(sp.Poly(result["U2"], t), sp.Poly((t-1)**4, t)).is_zero
    print(
        f"selector=PASS M49=72x64 rank=62 kernel=2 S=theta*(t-1)*r "
        f"C={C0} A={A1}", flush=True
    )
    print("terminal_tangent_orders_so_far=ord(U1)>=5, ord(U2)>=4", flush=True)
    return result


def add_second_block(data):
    t, r, h = data["t"], data["r"], data["h"]
    U, S, theta = data["U"], data["S"], data["theta"]
    A, A2 = data["A"], data["A2"]
    q24, q23, g = data["q24"], data["q23"], data["g"]
    f38, f37 = data["f38"], data["f37"]
    B = data["B"]
    params = data["parameters"]

    # The terminal tangent selector forces ord_1(U3)>=3.  The simple-root
    # K8 contact already forces r|U3, so use U3=(t-1)^3*r*C3.
    Cc = sp.symbols("C3_0:7")
    C3 = sum(Cc[i]*t**i for i in range(7))
    U3 = sp.expand((t-1)**3*r*C3)
    A3 = sp.expand(t**3*U3)
    q22 = sp.expand(t**5*h**2*(
        sp.Rational(5,3)*U3 + sp.Rational(10,9)*S*B
        - sp.Rational(5,81)*S**3
    ))
    Lc = sp.symbols("W36_8:68")
    ell = sum(Lc[i-8]*t**i for i in range(8,68))

    p47 = (
        15*U**3*sp.diff(ell,t) - 108*U**2*sp.diff(U,t)*ell
        + 14*A*sp.diff(data["f37"],t) - 37*sp.diff(A,t)*data["f37"]
        + 13*A2*sp.diff(f38,t) - 38*sp.diff(A2,t)*f38
        + 12*A3*sp.diff(g,t) - 39*sp.diff(A3,t)*g
        - 6*U**5*q22 - 6*q24*q23
    )
    q57 = (
        5*U**4*(5*U*sp.diff(ell,t) - 36*sp.diff(U,t)*ell)
        + 24*q24*sp.diff(data["f37"],t) - 37*sp.diff(q24,t)*data["f37"]
        + 23*q23*sp.diff(f38,t) - 38*sp.diff(q23,t)*f38
        + 22*q22*sp.diff(g,t) - 39*sp.diff(q22,t)*g
        - 20*U**9*A3 - 60*U**6*A*A2 - 20*U**3*A**3
    )
    variables = (*Lc, *Cc)
    matrices, forcings = [], []
    append_equation(matrices, forcings, p47, t, params, variables, "P47 tangent")
    append_equation(matrices, forcings, q57, t, params, variables, "Q57 tangent")
    append_equation(
        matrices, forcings, C3.subs(t,1)-sp.Rational(5,54)*theta**3,
        t, params, variables, "tangent U3 diagonal",
    )
    matrix = sp.Matrix.vstack(*matrices)
    solution, fibre, residuals, rank = solve_block(
        matrix, forcings, params, variables, "edge211 tangent 47/57", "F47_"
    )
    sub = dict(zip(variables, solution))
    out = dict(data)
    out.update({
        "ell": sp.expand(ell.subs(sub)),
        "C3": sp.expand(C3.subs(sub)),
        "U3": sp.expand(U3.subs(sub)),
        "A3": sp.expand(A3.subs(sub)),
        "q22": sp.expand(q22.subs(sub)),
        "parameters": (*params, *fibre),
        "constraints": [*data["constraints"], *residuals],
        "ranks": [*data["ranks"], rank],
    })
    if not sp.rem(sp.Poly(out["U3"],t),sp.Poly((t-1)**3,t)).is_zero:
        raise AssertionError("terminal U3 order was lost")
    return out


def add_third_block(data):
    t, r, h = data["t"], data["r"], data["h"]
    U, S, theta = data["U"], data["S"], data["theta"]
    A, A2, A3 = data["A"], data["A2"], data["A3"]
    q24, q23, q22, g = data["q24"], data["q23"], data["q22"], data["g"]
    f38, f37, ell = data["f38"], data["f37"], data["ell"]
    U3, B = data["U3"], data["B"]
    params = data["parameters"]

    # At the simple roots K8 forces r|U4; the terminal tangent gives
    # ord_1(U4)>=2, hence h|U4.  This parameterization preserves both.
    D4c = sp.symbols("U4factor_0:8")
    D4 = sum(D4c[i]*t**i for i in range(8))
    U4 = sp.expand(h*D4)
    A4 = sp.expand(t**3*U4)
    q21 = sp.expand(t**5*(
        sp.Rational(5,3)*h**2*U4
        + sp.Rational(5,243)*h*(
            54*S*U3 + 27*B**2 - 9*S**2*B + S**4
        )
    ))
    Wc = sp.symbols("W35_7:68")
    W35 = sum(Wc[i-7]*t**i for i in range(7,68))

    p46 = (
        15*U**3*sp.diff(W35,t) - 105*U**2*sp.diff(U,t)*W35
        + 14*A*sp.diff(ell,t) - 36*sp.diff(A,t)*ell
        + 13*A2*sp.diff(f37,t) - 37*sp.diff(A2,t)*f37
        + 12*A3*sp.diff(f38,t) - 38*sp.diff(A3,t)*f38
        + 11*A4*sp.diff(g,t) - 39*sp.diff(A4,t)*g
        - 6*U**5*q21 - 6*q24*q22 - 3*q23**2
    )
    q56 = (
        5*U**4*(5*U*sp.diff(W35,t)-35*sp.diff(U,t)*W35)
        + 24*q24*sp.diff(ell,t) - 36*sp.diff(q24,t)*ell
        + 23*q23*sp.diff(f37,t) - 37*sp.diff(q23,t)*f37
        + 22*q22*sp.diff(f38,t) - 38*sp.diff(q22,t)*f38
        + 21*q21*sp.diff(g,t) - 39*sp.diff(q21,t)*g
        - 20*U**9*A4 - 60*U**6*A*A3 - 30*U**6*A2**2
        - 60*U**3*A**2*A2 - 5*A**4
    )
    variables = (*Wc, *D4c)
    matrices, forcings = [], []
    append_equation(matrices, forcings, p46, t, params, variables, "P46 tangent")
    append_equation(matrices, forcings, q56, t, params, variables, "Q56 tangent")
    append_equation(
        matrices, forcings, D4.subs(t,1)-sp.Rational(5,432)*theta**4,
        t, params, variables, "tangent U4 diagonal",
    )
    matrix = sp.Matrix.vstack(*matrices)
    solution, fibre, residuals, rank = solve_block(
        matrix, forcings, params, variables, "edge211 tangent 46/56", "F46_"
    )
    sub = dict(zip(variables, solution))
    out = dict(data)
    out.update({
        "W35": sp.expand(W35.subs(sub)),
        "D4": sp.expand(D4.subs(sub)),
        "U4": sp.expand(U4.subs(sub)),
        "A4": sp.expand(A4.subs(sub)),
        "q21": sp.expand(q21.subs(sub)),
        "parameters": (*params, *fibre),
        "constraints": [*data["constraints"], *residuals],
        "ranks": [*data["ranks"], rank],
    })
    if not sp.rem(sp.Poly(out["U4"],t),sp.Poly((t-1)**2,t)).is_zero:
        raise AssertionError("terminal U4 order was lost")

    Djet = sp.expand(3*B-S**2)
    Ejet = sp.expand(27*U3-S**3)
    Delta = sp.expand(Ejet-3*S*Djet)
    Delta1, rem = sp.div(sp.Poly(Delta,t),sp.Poly(h,t))
    if not rem.is_zero:
        raise AssertionError("h does not divide the cube deviation")
    K8 = sp.rem(sp.Poly(S*Delta1-81*U4,t),sp.Poly((t-1)*r,t))
    if not K8.is_zero:
        raise AssertionError("the imposed deeper factors did not satisfy the K8 root selector")
    print("K8_selector=identically_satisfied_under_v0_tangent_factors", flush=True)
    return out


def summarize_parameter_ideal(data):
    constraints = [sp.factor(p) for p in data["constraints"] if sp.expand(p) != 0]
    params = data["parameters"]
    theta = data["theta"]
    print(
        f"total_carried_parameters={len(params)} residual_equations={len(constraints)}",
        flush=True,
    )
    if not constraints:
        print("global_residual_ideal=zero; no coefficient obstruction", flush=True)
        return

    # Exact elimination checks over QQ. The zero-theta specialization is kept
    # separate, and saturation tests whether the nonzero-theta stratum is
    # impossible.  Neither calculation assumes theta != 0.
    gb = sp.groebner(constraints, *params, order="grevlex", domain=QQ)
    unit = any(poly.as_expr() == 1 for poly in gb.polys)
    print(f"global_residual_groebner_basis_size={len(gb.polys)} unit={unit}", flush=True)

    zero_constraints = [sp.expand(p.subs(theta,0)) for p in constraints]
    zero_constraints = [p for p in zero_constraints if p != 0]
    other_params = tuple(p for p in params if p != theta)
    if not zero_constraints:
        zero_unit = False
    elif other_params:
        gb0 = sp.groebner(zero_constraints, *other_params, order="grevlex", domain=QQ)
        zero_unit = any(poly.as_expr() == 1 for poly in gb0.polys)
    else:
        zero_unit = sp.Rational(1) in [sp.factor(p) for p in zero_constraints]
    print(f"theta_zero_specialization_unit={zero_unit}", flush=True)

    z = sp.symbols("theta_inverse")
    nonzero_ideal = [*constraints, z*theta-1]
    gbnz = sp.groebner(nonzero_ideal, z, *params, order="grevlex", domain=QQ)
    nonzero_unit = any(poly.as_expr() == 1 for poly in gbnz.polys)
    print(f"theta_nonzero_saturation_unit={nonzero_unit}", flush=True)
    for i, p in enumerate(constraints):
        print(f"global_residual_{i}={sp.factor(p)}", flush=True)


def main():
    print("source_scope=GGHV F2 terminal attachment is a theorem input, not proved by this arithmetic file", flush=True)
    data = selector_and_first_block()
    data = add_second_block(data)
    data = add_third_block(data)
    summarize_parameter_ideal(data)
    print("scope=necessary exact coupled coefficient conditions only; zero theta retained; no U5 assumption", flush=True)


if __name__ == "__main__":
    main()
