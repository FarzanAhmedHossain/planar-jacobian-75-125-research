"""Exact compatibility audit for the fixed F2 quartic pattern 2+1+1.

This is a new, standalone replay of (i) the previously unsaved K7-divisibility
augmentations to the weight 48/58 and 47/57 power-relation bands, (ii) the
corresponding deeper-contact subbranch, and (iii) the next 46/56 block with
the tuned-root K8 cube-deviation selector.

Everything here is conditional on the normalized F2 edge/support model and
uses necessary coefficient equations only. A compatible relaxed block is not
a Keller pair; a unit-ideal certificate from a relaxation would exclude the
branch. The script never divides by a parameter. The tuned/deep split is
explicit, and the deeper B=h*T computation is labelled as the q>=4
subbranch (or U1=0 after its separate K12 bound), not as all deeper cases.
"""

from __future__ import annotations

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import QQ

from f2_75_125_power_relation_pq_nextbands import (
    coefficient_system,
    weight49_parameters,
)
from f2_75_125_power_relation_pq_thirdbands import (
    exact_affine_solution,
    residual_polynomials,
)


def fixed_carrier():
    t = sp.symbols("t")
    r = t**2 - 3*t + 3
    h = sp.expand((t - 1)**2 * r)
    v = sp.expand(-sp.Rational(1, 3) * t * (t - 1) * r)
    return t, r, h, v


def rank(matrix: sp.Matrix) -> int:
    return DomainMatrix.from_Matrix(matrix).convert_to(QQ).rank()


def rename_solution(solution, free, prefix: str):
    fresh = sp.symbols(f"{prefix}0:{len(free)}")
    rename = dict(zip(free, fresh))
    return [sp.expand(value.subs(rename)) for value in solution], fresh


def tuned_48_58_prefix():
    """Solve 48/58 plus h | D, where D=3B-S^2 (all roots tuned)."""
    t, r, h, v = fixed_carrier()
    U = sp.expand(t*h)
    params, S, f38, _ = weight49_parameters(h, v, QQ, "211_tuned")

    Bc = sp.symbols("B0:8")
    Fc = sp.symbols("F8:68")
    B = sum(Bc[i]*t**i for i in range(8))
    f37 = sum(Fc[i-8]*t**i for i in range(8, 68))
    U1, U2 = sp.expand(h**2*S), sp.expand(h*B)
    A, A2 = sp.expand(t**3*U1), sp.expand(t**3*U2)
    q24 = sp.expand(sp.Rational(5, 3)*t**5*h**2*U1)
    q23 = sp.expand(
        sp.Rational(5, 3)*t**5*h**2*U2
        + sp.Rational(5, 9)*t**5*h**3*S**2
    )
    g = sp.expand(U**7*v)

    p48 = (
        15*U**3*sp.diff(f37, t) - 111*U**2*sp.diff(U, t)*f37
        + 14*A*sp.diff(f38, t) - 38*sp.diff(A, t)*f38
        + 13*A2*sp.diff(g, t) - 39*sp.diff(A2, t)*g
        - 6*U**5*q23 - 3*q24**2
    )
    q58 = (
        5*U**4*(5*U*sp.diff(f37, t) - 37*sp.diff(U, t)*f37)
        + 24*q24*sp.diff(f38, t) - 38*sp.diff(q24, t)*f38
        + 23*q23*sp.diff(g, t) - 39*sp.diff(q23, t)*g
        - 20*U**9*A2 - 30*U**6*A**2
    )
    variables = (*Fc, *Bc)
    Mp, bp = coefficient_system(p48, t, params, variables, "tuned P48")
    Mq, bq = coefficient_system(q58, t, params, variables, "tuned Q58")
    # K7 tuned-root jet: h divides D=3B-S^2 (double-root order two too).
    drem = sp.rem(sp.Poly(3*B-S**2, t), sp.Poly(h, t)).as_expr()
    Md, bd = coefficient_system(drem, t, params, variables, "tuned D mod h")
    M = sp.Matrix.vstack(Mp, Mq, Md)
    forcing = [*bp, *bq, *bd]
    print(f"tuned_48_58_plus_Dmodh shape={M.shape} rank={rank(M)}", flush=True)
    solution, free, rnk = exact_affine_solution(
        M, forcing, QQ, variables, params, "211 tuned 48/58+Dmodh"
    )
    solution, free = rename_solution(solution, free, "Tfree")
    p1 = (*params, *free)
    sub = dict(zip(variables, solution))
    return {
        "t": t, "r": r, "h": h, "v": v, "U": U, "S": S,
        "f38": f38, "f37": sp.expand(f37.subs(sub)),
        "B": sp.expand(B.subs(sub)), "U1": U1,
        "U2": sp.expand(U2.subs(sub)),
        "A": A, "A2": sp.expand(A2.subs(sub)),
        "q24": q24, "q23": sp.expand(q23.subs(sub)), "g": g,
        "parameters": p1, "rank": rnk, "free": free,
    }


def tuned_47_57_with_Emodh(data):
    """Solve 47/57 plus h | E=27U3-S^3, retaining all old free parameters."""
    t, h, v = data["t"], data["h"], data["v"]
    U, S, params = data["U"], data["S"], data["parameters"]
    A, A2, q24, q23, g = (
        data["A"], data["A2"], data["q24"], data["q23"], data["g"]
    )
    f38, f37 = data["f38"], data["f37"]
    U2 = data["U2"]
    q23 = sp.expand(q23)

    Cc = sp.symbols("C0:12")
    U3 = sum(Cc[i]*t**i for i in range(12))
    A3 = sp.expand(t**3*U3)
    q22 = sp.expand(t**5*h**2*(
        sp.Rational(5, 3)*U3 + sp.Rational(10, 9)*S*data["B"]
        - sp.Rational(5, 81)*S**3
    ))
    Lc = sp.symbols("L8:68")
    ell = sum(Lc[i-8]*t**i for i in range(8, 68))

    p47 = (
        15*U**3*sp.diff(ell, t) - 108*U**2*sp.diff(U, t)*ell
        + 14*A*sp.diff(f37, t) - 37*sp.diff(A, t)*f37
        + 13*A2*sp.diff(f38, t) - 38*sp.diff(A2, t)*f38
        + 12*A3*sp.diff(g, t) - 39*sp.diff(A3, t)*g
        - 6*U**5*q22 - 6*q24*q23
    )
    q57 = (
        5*U**4*(5*U*sp.diff(ell, t) - 36*sp.diff(U, t)*ell)
        + 24*q24*sp.diff(f37, t) - 37*sp.diff(q24, t)*f37
        + 23*q23*sp.diff(f38, t) - 38*sp.diff(q23, t)*f38
        + 22*q22*sp.diff(g, t) - 39*sp.diff(q22, t)*g
        - 20*U**9*A3 - 60*U**6*A*A2 - 20*U**3*A**3
    )

    variables = (*Lc, *Cc)
    M1, b1 = coefficient_system(p47, t, params, variables, "tuned P47")
    M2, b2 = coefficient_system(q57, t, params, variables, "tuned Q57")
    erem = sp.rem(sp.Poly(27*U3-S**3, t), sp.Poly(h, t)).as_expr()
    Me, be = coefficient_system(erem, t, params, variables, "tuned E mod h")
    M = sp.Matrix.vstack(M1, M2, Me)
    forcing = [*b1, *b2, *be]
    print(f"tuned_47_57_plus_Emodh shape={M.shape} rank={rank(M)}", flush=True)
    solution, free, rnk = exact_affine_solution(
        M, forcing, QQ, variables, params, "211 tuned 47/57+Emodh"
    )
    solution, free = rename_solution(solution, free, "Ufree")
    sub = dict(zip(variables, solution))
    result = dict(data)
    result.update({
        "ell": sp.expand(ell.subs(sub)),
        "U3": sp.expand(U3.subs(sub)),
        "A3": sp.expand(A3.subs(sub)),
        "q22": sp.expand(q22.subs(sub)),
        "parameters": (*params, *free),
        "rank47": rnk,
        "free47": free,
    })
    return result


def tuned_k8_46_56_attack(data):
    """Add the K8 root-jet selector to the next 46/56 paired bands.

    In the tuned case K7 gives h|D and h|E, so Delta=(E-3SD) is h*Delta1.
    K8's requirements ord(J8)>=3 at t=1 and >=2 at the two simple roots,
    with J8=S*Delta-81*h*U4, are equivalent to
        (t-1)r | S*Delta1 - 81*U4.
    This is the first explicit test here of deviation from the local cube.
    """
    t, r, h, U, S = (data[k] for k in ("t", "r", "h", "U", "S"))
    v, g = data["v"], data["g"]
    f38, f37, ell = data["f38"], data["f37"], data["ell"]
    B, U2, U3 = data["B"], data["U2"], data["U3"]
    A, A2, A3 = data["A"], data["A2"], data["A3"]
    q24, q23, q22 = data["q24"], data["q23"], data["q22"]
    params = data["parameters"]

    D = sp.expand(3*B-S**2)
    E = sp.expand(27*U3-S**3)
    Delta = sp.expand(E-3*S*D)
    quotient, remainder = sp.div(sp.Poly(Delta, t), sp.Poly(h, t))
    if not remainder.is_zero:
        raise AssertionError("the carried K7 equations did not imply h | Delta")
    Delta1 = quotient.as_expr()

    U4c = sp.symbols("U4_0:12")
    U4 = sum(U4c[i]*t**i for i in range(12))
    A4 = sp.expand(t**3*U4)
    q21 = sp.expand(t**5*(
        sp.Rational(5, 3)*h**2*U4
        + sp.Rational(5, 243)*h*(
            54*S*U3 + 27*B**2 - 9*S**2*B + S**4
        )
    ))
    Wc = sp.symbols("W35_7:68")
    W35 = sum(Wc[i-7]*t**i for i in range(7, 68))
    p46 = (
        15*U**3*sp.diff(W35, t) - 105*U**2*sp.diff(U, t)*W35
        + 14*A*sp.diff(ell, t) - 36*sp.diff(A, t)*ell
        + 13*A2*sp.diff(f37, t) - 37*sp.diff(A2, t)*f37
        + 12*A3*sp.diff(f38, t) - 38*sp.diff(A3, t)*f38
        + 11*A4*sp.diff(g, t) - 39*sp.diff(A4, t)*g
        - 6*U**5*q21 - 6*q24*q22 - 3*q23**2
    )
    q56 = (
        5*U**4*(5*U*sp.diff(W35, t) - 35*sp.diff(U, t)*W35)
        + 24*q24*sp.diff(ell, t) - 36*sp.diff(q24, t)*ell
        + 23*q23*sp.diff(f37, t) - 37*sp.diff(q23, t)*f37
        + 22*q22*sp.diff(f38, t) - 38*sp.diff(q22, t)*f38
        + 21*q21*sp.diff(g, t) - 39*sp.diff(q21, t)*g
        - 20*U**9*A4 - 60*U**6*A*A3 - 30*U**6*A2**2
        - 60*U**3*A**2*A2 - 5*A**4
    )
    variables = (*Wc, *U4c)
    M1, b1 = coefficient_system(p46, t, params, variables, "tuned P46")
    M2, b2 = coefficient_system(q56, t, params, variables, "tuned Q56")
    # K8: after removing the common h factor in J8, rad(h) must divide it.
    rad = sp.expand((t-1)*r)
    k8rem = sp.rem(
        sp.Poly(S*Delta1-81*U4, t), sp.Poly(rad, t)
    ).as_expr()
    Mk, bk = coefficient_system(k8rem, t, params, variables, "tuned K8 root jet")
    M = sp.Matrix.vstack(M1, M2, Mk)
    forcing = [*b1, *b2, *bk]
    rnk = rank(M)
    print(
        f"tuned_46_56_plus_K8root shape={M.shape} rank={rnk} "
        f"leftkernel={M.rows-rnk} variables={M.cols}", flush=True
    )
    residuals = residual_polynomials(
        M, forcing, QQ, params, "211 tuned 46/56+K8root", rank=rnk
    )
    for i, residual in enumerate(residuals):
        print(f"tuned_K8_residual_{i}={sp.factor(residual)}", flush=True)
    if residuals:
        basis = sp.groebner(residuals, *params, domain=QQ)
        print(f"tuned_K8_residual_gb={list(basis)}", flush=True)
        print(f"tuned_K8_residual_ideal_unit={basis.contains(1)}", flush=True)
    else:
        print("tuned_K8_residual_ideal=zero", flush=True)
    return residuals


def deeper_q4_prefix():
    """Replay the deeper q>=4 subbranch; B=hT, U3=rad*C.

    This is necessary for U1=0 after the K12 double-root conclusion, and is
    also a subbranch of the nonzero all-deeper stratum. The latter's q=3
    double-root leaf must be checked separately.
    """
    t, r, h, v = fixed_carrier()
    rad = sp.expand((t-1)*r)
    U = sp.expand(t*h)
    theta = sp.symbols("theta")
    S = sp.expand(theta*rad)
    f38 = sp.expand(-sp.Rational(5, 6)*theta*t**8*(t-1)**14*r**8)
    Tc = sp.symbols("T0:4")
    T = sum(Tc[i]*t**i for i in range(4))
    B = sp.expand(h*T)
    Fc = sp.symbols("F8:68")
    f37 = sum(Fc[i-8]*t**i for i in range(8, 68))
    U1, U2 = sp.expand(h**2*S), sp.expand(h*B)
    A, A2 = sp.expand(t**3*U1), sp.expand(t**3*U2)
    q24 = sp.expand(sp.Rational(5, 3)*t**5*h**2*U1)
    q23 = sp.expand(
        sp.Rational(5, 3)*t**5*h**2*U2
        + sp.Rational(5, 9)*t**5*h**3*S**2
    )
    g = sp.expand(U**7*v)
    p48 = (
        15*U**3*sp.diff(f37, t) - 111*U**2*sp.diff(U, t)*f37
        + 14*A*sp.diff(f38, t) - 38*sp.diff(A, t)*f38
        + 13*A2*sp.diff(g, t) - 39*sp.diff(A2, t)*g
        - 6*U**5*q23 - 3*q24**2
    )
    q58 = (
        5*U**4*(5*U*sp.diff(f37, t) - 37*sp.diff(U, t)*f37)
        + 24*q24*sp.diff(f38, t) - 38*sp.diff(q24, t)*f38
        + 23*q23*sp.diff(g, t) - 39*sp.diff(q23, t)*g
        - 20*U**9*A2 - 30*U**6*A**2
    )
    oldvars = (*Fc, *Tc)
    params = (theta,)
    M1, b1 = coefficient_system(p48, t, params, oldvars, "deep P48")
    M2, b2 = coefficient_system(q58, t, params, oldvars, "deep Q58")
    Mold = sp.Matrix.vstack(M1, M2)
    bold = [*b1, *b2]
    print(f"deep_q4_48_58 shape={Mold.shape} rank={rank(Mold)}", flush=True)
    sol, free, rnk = exact_affine_solution(
        Mold, bold, QQ, oldvars, params, "211 deeper q4 48/58"
    )
    sol, free = rename_solution(sol, free, "Dfree")
    params1 = (*params, *free)
    sub = dict(zip(oldvars, sol))
    Bsol, f37sol = sp.expand(B.subs(sub)), sp.expand(f37.subs(sub))
    U2sol, A2sol = sp.expand(h*Bsol), sp.expand(t**3*h*Bsol)
    q23sol = sp.expand(
        sp.Rational(5, 3)*t**5*h**2*U2sol
        + sp.Rational(5, 9)*t**5*h**3*S**2
    )
    Cc = sp.symbols("C0:9")
    C = sum(Cc[i]*t**i for i in range(9))
    U3 = sp.expand(rad*C)
    A3 = sp.expand(t**3*U3)
    q22 = sp.expand(t**5*h**2*(
        sp.Rational(5, 3)*U3 + sp.Rational(10, 9)*S*Bsol
        - sp.Rational(5, 81)*S**3
    ))
    Lc = sp.symbols("L8:68")
    ell = sum(Lc[i-8]*t**i for i in range(8, 68))
    p47 = (
        15*U**3*sp.diff(ell, t) - 108*U**2*sp.diff(U, t)*ell
        + 14*A*sp.diff(f37sol, t) - 37*sp.diff(A, t)*f37sol
        + 13*A2sol*sp.diff(f38, t) - 38*sp.diff(A2sol, t)*f38
        + 12*A3*sp.diff(g, t) - 39*sp.diff(A3, t)*g
        - 6*U**5*q22 - 6*q24*q23sol
    )
    q57 = (
        5*U**4*(5*U*sp.diff(ell, t) - 36*sp.diff(U, t)*ell)
        + 24*q24*sp.diff(f37sol, t) - 37*sp.diff(q24, t)*f37sol
        + 23*q23sol*sp.diff(f38, t) - 38*sp.diff(q23sol, t)*f38
        + 22*q22*sp.diff(g, t) - 39*sp.diff(q22, t)*g
        - 20*U**9*A3 - 60*U**6*A*A2sol - 20*U**3*A**3
    )
    nextvars = (*Lc, *Cc)
    N1, n1 = coefficient_system(p47, t, params1, nextvars, "deep P47")
    N2, n2 = coefficient_system(q57, t, params1, nextvars, "deep Q57")
    N = sp.Matrix.vstack(N1, N2)
    nf = [*n1, *n2]
    nrank = rank(N)
    print(
        f"deep_q4_47_57 shape={N.shape} rank={nrank} "
        f"leftkernel={N.rows-nrank} variables={N.cols}", flush=True
    )
    residuals = residual_polynomials(
        N, nf, QQ, params1, "211 deeper q4 47/57", rank=nrank
    )
    for i, residual in enumerate(residuals):
        print(f"deep_q4_residual_{i}={sp.factor(residual)}", flush=True)
    if not residuals:
        print("deep_q4_residual_ideal=zero", flush=True)
    return residuals


def audit_weight49_selector():
    """Replay the two-dimensional S selector and its synchronized root values."""
    t, r, h, v = fixed_carrier()
    M = __import__("f2_75_125_power_relation_second_band").build_matrix(h, v)
    D = DomainMatrix.from_Matrix(M).convert_to(QQ)
    K = D.nullspace().to_Matrix()
    if M.shape != (72, 64) or D.rank() != 62 or K.shape != (2, 64):
        raise AssertionError("unexpected weight-49 matrix/kernel profile")
    Sbasis = [sp.factor(sum(K[i,j]*t**j for j in range(4))) for i in range(2)]
    target1, target2 = (t-1)*r, 3*t**3-112*t**2+318*t-264
    # Exact span equality for the two-dimensional S projection.
    target_matrix = sp.Matrix([[sp.expand(p).coeff(t,j) for j in range(4)]
                               for p in (target1,target2)])
    basis_matrix = sp.Matrix([[sp.expand(p).coeff(t,j) for j in range(4)]
                              for p in Sbasis])
    if basis_matrix.rank() != 2 or target_matrix.rank() != 2:
        raise AssertionError("S projection is not two-dimensional")
    if basis_matrix.col_join(target_matrix).rank() != 2:
        raise AssertionError("claimed S span differs from computed kernel")
    assert target1.subs(t, 1) == 0 and target2.subs(t, 1) == -55
    assert sp.rem(target1, r, t) == 0 and sp.rem(target2, r, t) == 45
    print("weight49_selector=PASS matrix=72x64 rank=62 kernel=2", flush=True)
    print("S_span={(t-1)r, 3t^3-112t^2+318t-264}", flush=True)
    print("root_values: S(1)=-55v; S(alpha)=S(beta)=45v", flush=True)


def main():
    audit_weight49_selector()
    tuned = tuned_48_58_prefix()
    tuned = tuned_47_57_with_Emodh(tuned)
    tuned_k8_46_56_attack(tuned)
    deeper_q4_prefix()
    print("scope=exact necessary coefficient compatibilities only; no pair/branch exclusion unless an explicit unit certificate is reported", flush=True)


if __name__ == "__main__":
    main()
