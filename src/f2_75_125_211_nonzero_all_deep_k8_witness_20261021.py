#!/usr/bin/env python3
"""Exact replay of a nonzero-U1 all-deep K3--K8 local survivor.

This is deliberately a witness against the sufficiency of the current local
divisibility screen.  It is not a solution of the full band recurrence or a
Keller pair.
"""

from __future__ import annotations

import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
OUT = HERE / "f2_75_125_211_nonzero_all_deep_k8_witness_20261021.json"
sys.path.insert(0, str(HERE))
import f2_75_125_k17_nonzero_first_slice as template_module  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def quotient_by_power(numerator: sp.Expr, divisor: sp.Expr, power: int, t: sp.Symbol):
    den = sp.Poly(divisor**power, t, extension=True)
    quotient, remainder = sp.div(sp.Poly(sp.expand(numerator), t, extension=True), den)
    if remainder:
        raise AssertionError(f"not divisible by h^{power}: remainder={remainder.as_expr()}")
    return quotient.as_expr()


def order_at(poly: sp.Expr, point: sp.Expr, t: sp.Symbol) -> int | str:
    p = sp.Poly(sp.expand(poly), t, domain=sp.QQ)
    if p.is_zero:
        return "infinity"
    u = sp.Symbol("u")
    shifted = sp.Poly(sp.expand(p.as_expr().subs(t, point + u)), u, domain=sp.QQ)
    return min(monomial[0] for monomial, coefficient in shifted.terms() if coefficient)


def order_profile_screen(profile: tuple[int, ...], root_order: int):
    """Replay the existing K9--K24 unique-pole profile screen without NumPy."""
    infinity = 99

    def packed(total: int, alpha: Fraction, exponent_base: int):
        rows = []
        for multiplicities, coefficient, _old_base in template_module.templates(
            total, alpha, tuple(range(1, total + 1)), exponent_base
        ):
            if not coefficient:
                continue
            first = tuple(multiplicities[:5]) + (0,) * max(0, 5 - len(multiplicities))
            base = root_order * (exponent_base - 3 * sum(multiplicities))
            supported = not any(multiplicities[5:])
            rows.append((first, base, coefficient, supported))
        return rows

    def least(terms):
        evaluated = []
        for multiplicities, base, coefficient, supported in terms:
            if any(m and profile[i] >= infinity
                   for i, m in enumerate(multiplicities)):
                continue
            order = base + sum(m * profile[i] for i, m in enumerate(multiplicities))
            evaluated.append((order, multiplicities, coefficient, supported))
        if not evaluated:
            return None, [], False
        minimum = min(row[0] for row in evaluated)
        tied = [row for row in evaluated if row[0] == minimum]
        return minimum, tied, all(row[3] for row in tied)

    by_band = {}
    for n in range(9, 25):
        main = packed(n, Fraction(5, 3), 5)
        corrections = [
            packed(n - 5 * resonance,
                   Fraction(5 - resonance, 3), 5 - resonance)
            for resonance in range(1, min(4, n // 5) + 1)
        ]
        main_min, tied, _ = least(main)
        correction_min = min(
            (least(terms)[0] if least(terms)[0] is not None else 10**9)
            for terms in corrections
        ) if corrections else 10**9
        obstructed = (
            main_min is not None and main_min < 0 and len(tied) == 1
            and tied[0][3] and main_min < correction_min
        )
        by_band[str(n)] = {
            "main_minimum_order": main_min,
            "main_minimum_tie_count": len(tied),
            "main_minimum_supported_on_U1_to_U5": bool(tied) and all(row[3] for row in tied),
            "correction_minimum_order": correction_min,
            "uniquely_obstructed": obstructed,
        }
        if obstructed:
            return {"survives": False, "first_excluded_band": n, "by_band": by_band}
    return {"survives": True, "first_excluded_band": None, "by_band": by_band}


def solve_k8_double_root_slice(t: sp.Symbol, R: sp.Expr, h: sp.Expr, rho: sp.Expr):
    """Solve the first two nontrivial K8 coefficients on a K9--K24 survivor profile."""
    u = sp.symbols("u")
    d0, d1 = sp.symbols("d0 d1")
    # These degree-capped polynomials realize the valuation witness
    # [5,3,2,0,0] at t=1 and [3,2,2,1,1] at each simple root of R.
    A, B = rho, rho
    U1, U2 = sp.expand(h**2 * A), sp.expand(h * B)
    U3 = sp.expand((t - 1) ** 2 * R**2)
    U4 = sp.expand(R * (d0 + d1 * (t - 1)))
    # The band variable is U5=t*(U5/t); t is a unit at every root of h,
    # but the support condition requires retaining this factor globally.
    U5 = t * R
    D = sp.expand(3 * B - A**2)
    E = sp.expand(27 * U3 - A**3)
    Delta = sp.expand(E - 3 * A * D)
    J8 = sp.expand(A * Delta - 81 * U4 * h)
    N8 = sp.Rational(5, 59049) * sp.expand(
        J8**2 - D * Delta**2 + 6 * A * D**2 * Delta
        - 3 * A**2 * D**3 + 3 * D**4 - 243 * D**2 * U4 * h
        + 486 * Delta * U5 * h**2
    )
    at_u = sp.Poly(sp.expand(N8.subs(t, u + 1)), u, domain=sp.QQ[d0, d1])
    coeff4 = sp.factor(at_u.nth(4))
    coeff5 = sp.factor(at_u.nth(5))
    solutions = sp.solve((coeff4, coeff5), (d0, d1), dict=True)
    N3 = sp.expand(18 * h**3 * U1 * U2 - U1**3)
    N4 = sp.expand(
        54 * h**6 * U1 * U3
        + 27 * h**6 * U2**2
        - 9 * h**3 * U1**2 * U2
        + U1**4
    )
    N5 = sp.expand(
        810 * h**9 * (U1 * U4 + U2 * U3)
        - 135 * h**6 * (U1**2 * U3 + U1 * U2**2)
        + 60 * h**3 * U1**3 * U2
        - 7 * U1**5
    )
    W6 = sp.expand(Delta**2 - 3 * D**3)
    W7 = sp.expand(A * W6 + 3 * D**2 * Delta - 162 * U4 * h * Delta)
    screens = {
        "double_root": order_profile_screen((5, 3, 2, 0, 0), 2),
        "simple_root": order_profile_screen((3, 2, 2, 1, 1), 1),
    }
    checked = []
    for solution in solutions:
        specialized = sp.factor(N8.subs(solution))
        quotient = quotient_by_power(specialized, h, 3, t)
        exact_quotients = {
            "K3": quotient_by_power(N3, h, 4, t),
            "K4": quotient_by_power(N4, h, 7, t),
            "K5": quotient_by_power(N5.subs(solution), h, 10, t),
            "K6": quotient_by_power(W6, h, 1, t),
            "K7": quotient_by_power(W7.subs(solution), h, 2, t),
            "K8": quotient,
        }
        checked.append({
            "d0": str(sp.simplify(solution[d0])),
            "d1": str(sp.simplify(solution[d1])),
            "K8_quotient": str(sp.factor(quotient)),
            "all_K3_to_K8_divisibilities_verified": True,
            "K3_to_K7_quotient_degrees": {
                name: int(sp.Poly(expr, t, extension=True).degree())
                for name, expr in exact_quotients.items() if name != "K8"
            },
            "U4_at_double_root_nonzero": sp.simplify(solution[d0]) != 0,
        })
    if not checked:
        raise AssertionError("the K8 order-four/order-five equations had no solution")
    return {
        "profile_at_double_root_U1_to_U5_over_t": [5, 3, 2, 0, 0],
        "profile_at_each_simple_root_U1_to_U5_over_t": [3, 2, 2, 1, 1],
        "assignment_family": {
            "U1": "(t-1)^5 R^3",
            "U2": "(t-1)^3 R^2",
            "U3": "(t-1)^2 R^2",
            "U4": "R*(d0+d1*(t-1))",
            "U5_over_t": "R",
        },
        "K8_order4_coefficient": str(coeff4),
        "K8_order5_coefficient": str(coeff5),
        "solutions_with_exact_global_h3_divisibility": checked,
        "K9_K24_unique_pole_screens": {
            name: {"survives": value["survives"],
                   "first_excluded_band": value["first_excluded_band"]}
            for name, value in screens.items()
        },
        "status": (
            "EXACT_K3_K8_COMPATIBLE_AND_K9_K24_ORDER_SCREEN_SURVIVOR"
            if all(value["survives"] for value in screens.values())
            else "K3_K8_COMPATIBLE_BUT_REMOVED_BY_K9_K24_ORDER_SCREEN"
        ),
    }


def dependent_power_lift(t: sp.Symbol, R: sp.Expr, h: sp.Expr, rho: sp.Expr):
    """Construct a full dependent pair invisible to all homogeneous edge bands."""
    w, x, y = sp.symbols("w x y")
    S = h + rho * w / 3
    F = sp.expand(S**3)
    G = sp.expand(S**5)
    Phi = sp.expand(
        w * (-3 * F * sp.diff(G, w) + 5 * G * sp.diff(F, w))
        + t * (25 * sp.diff(F, t) * G - 15 * F * sp.diff(G, t))
        + t * w * (-sp.diff(F, t) * sp.diff(G, w)
                   + sp.diff(F, w) * sp.diff(G, t))
    )
    if Phi != 0:
        raise AssertionError("the power lift failed the full homogeneous near-edge bracket identity")

    U = [sp.expand(F.coeff(w, i)) for i in range(6)]
    V = [sp.expand(G.coeff(w, i)) for i in range(6)]
    expected_U = [h**3, h**2 * rho, h * rho**2 / 3, rho**3 / 27, 0, 0]
    if any(sp.expand(U[i] - expected_U[i]) != 0 for i in range(6)):
        raise AssertionError("unexpected P slices in the exact cube lift")
    if any(sp.Poly(V[i], t).degree() > 19 for i in range(1, 5) if V[i] != 0):
        raise AssertionError("an exact Q slice exceeded the recorded top-band degree cap")

    # Lift the edge series back to actual x,y polynomials.  This pair is
    # intentionally dependent: it witnesses why [P,Q]=1 must enter the proof.
    tx = x * y**5
    common = sp.expand(x * (h.subs(t, tx) + rho.subs(t, tx) * y / 3))
    P, Q = sp.expand(common**3), sp.expand(common**5)
    deg_common = max(i + j for (i, j), coefficient in sp.Poly(common, x, y).terms() if coefficient)
    deg_P = max(i + j for (i, j), coefficient in sp.Poly(P, x, y).terms() if coefficient)
    deg_Q = max(i + j for (i, j), coefficient in sp.Poly(Q, x, y).terms() if coefficient)
    support_C = {monomial for monomial, coefficient in sp.Poly(common, x, y).terms() if coefficient}
    support_P = {monomial for monomial, coefficient in sp.Poly(P, x, y).terms() if coefficient}
    support_Q = {monomial for monomial, coefficient in sp.Poly(Q, x, y).terms() if coefficient}
    jacobian = sp.expand(sp.diff(P, x) * sp.diff(Q, y) - sp.diff(P, y) * sp.diff(Q, x))
    if (deg_common, deg_P, deg_Q, jacobian) != (25, 75, 125, 0):
        raise AssertionError("the dependent polynomial lift failed its exact degree/Jacobian check")
    if (
        not all(0 <= i <= 5 and 0 <= j <= 20 for i, j in support_C)
        or not all(0 <= i <= 15 and 0 <= j <= 60 for i, j in support_P)
        or not all(0 <= i <= 25 and 0 <= j <= 100 for i, j in support_Q)
        or (5, 20) not in support_C or (15, 60) not in support_P or (25, 100) not in support_Q
    ):
        raise AssertionError("the dependent lift missed the scaled Newton support rectangles/corners")

    return {
        "definition": "C=x*(h(x*y^5)+(x*y^5-1)*R(x*y^5)*y/3), P=C^3, Q=C^5",
        "slice_identity": "F=(h+rho*w/3)^3, G=(h+rho*w/3)^5=F^(5/3)",
        "homogeneous_near_edge_bracket_Phi": "0 identically",
        "P_slices_U1_to_U5": [str(sp.factor(U[i])) for i in range(1, 6)],
        "Q_slice_degrees_V0_to_V5": [
            int(sp.Poly(value, t).degree()) if value != 0 else None for value in V
        ],
        "total_degrees": {"C": deg_common, "P": deg_P, "Q": deg_Q},
        "support_box_verified": {
            "C": "0<=i<=5, 0<=j<=20, top corner (5,20) present",
            "P": "0<=i<=15, 0<=j<=60, top corner (15,60) present",
            "Q": "0<=i<=25, 0<=j<=100, top corner (25,100) present",
        },
        "Jacobian": "0 identically",
        "root_order_profiles": {
            "double_root_U1_to_U5_over_t": [5, 4, 3, "infinity", "infinity"],
            "each_simple_root_U1_to_U5_over_t": [3, 3, 3, "infinity", "infinity"],
        },
        "order_profile_screen": {
            "double_root": order_profile_screen((5, 4, 3, 99, 99), 2),
            "simple_root": order_profile_screen((3, 3, 3, 99, 99), 1),
        },
        "scope": (
            "This is an exact dependent polynomial pair of degrees (75,125) with zero Jacobian, not a Keller pair. It satisfies the homogeneous near-edge bracket identically and all corresponding local root divisibilities. It exhibits the degenerate power component that local homogeneous equations cannot distinguish from a Keller branch; the nonzero constant-Jacobian target must be used."
        ),
    }


def main() -> None:
    t = sp.symbols("t")
    R = t**2 - 3 * t + 3
    rho = (t - 1) * R
    h = (t - 1) ** 2 * R

    # Choose A=rho, B=h, U3=h, U4=rho, U5=0, and recover the original
    # coefficient slices U1=h^2*A and U2=h*B.
    A, B, U3, U4, U5 = rho, h, h, rho, sp.S.Zero
    U1, U2 = sp.expand(h**2 * A), sp.expand(h * B)

    # Exact binomial-power polynomiality conditions K3--K5.
    N3 = sp.expand(18 * h**3 * U1 * U2 - U1**3)
    q3 = quotient_by_power(N3, h, 4, t)
    N4 = sp.expand(
        54 * h**6 * U1 * U3
        + 27 * h**6 * U2**2
        - 9 * h**3 * U1**2 * U2
        + U1**4
    )
    q4 = quotient_by_power(N4, h, 7, t)
    N5 = sp.expand(
        810 * h**9 * (U1 * U4 + U2 * U3)
        - 135 * h**6 * (U1**2 * U3 + U1 * U2**2)
        + 60 * h**3 * U1**3 * U2
        - 7 * U1**5
    )
    q5 = quotient_by_power(N5, h, 10, t)

    # K6/K7 pole numerators after U1=h^2*A and U2=h*B.
    D = sp.expand(3 * B - A**2)
    E = sp.expand(27 * U3 - A**3)
    Delta = sp.expand(E - 3 * A * D)
    W6 = sp.expand(Delta**2 - 3 * D**3)
    W7 = sp.expand(A * W6 + 3 * D**2 * Delta - 162 * U4 * h * Delta)
    q6 = quotient_by_power(W6, h, 1, t)
    q7 = quotient_by_power(W7, h, 2, t)

    # K8 exact pole numerator from the compressed identity.
    J8 = sp.expand(A * Delta - 81 * U4 * h)
    N8 = sp.Rational(5, 59049) * sp.expand(
        J8**2
        - D * Delta**2
        + 6 * A * D**2 * Delta
        - 3 * A**2 * D**3
        + 3 * D**4
        - 243 * D**2 * U4 * h
        + 486 * Delta * U5 * h**2
    )
    q8 = quotient_by_power(N8, h, 3, t)

    pieces = {"R": R, "h": h, "rho": rho, "U1": U1, "U2": U2,
              "U3": U3, "U4": U4, "U5_over_t": U5,
              "K3_quotient": q3, "K4_quotient": q4, "K5_quotient": q5,
              "K6_quotient": q6, "K7_quotient": q7, "K8_quotient": q8}
    encoded = {name: str(sp.factor(expr)) for name, expr in pieces.items()}
    degree_bounds = {
        name: (None if expr == 0 else int(sp.Poly(expr, t).degree()))
        for name, expr in {"U1": U1, "U2": U2, "U3": U3,
                           "U4": U4, "U5_over_t": U5}.items()
    }
    root_orders = {
        "t=1_double_root": {
            name: order_at(expr, sp.Integer(1), t)
            for name, expr in {"h": h, "U1": U1, "U2": U2,
                               "U3": U3, "U4": U4}.items()
        },
        "each_simple_root_of_R": {
            "h": 1, "U1": 3, "U2": 2, "U3": 1, "U4": 1
        },
    }
    if degree_bounds != {"U1": 11, "U2": 8, "U3": 4, "U4": 3, "U5_over_t": None}:
        raise AssertionError(f"unexpected degree profile {degree_bounds}")
    if root_orders["t=1_double_root"] != {
        "h": 2, "U1": 5, "U2": 4, "U3": 2, "U4": 1
    }:
        raise AssertionError(f"unexpected double-root profile {root_orders}")
    profile_screens = {
        "double_root_profile_U1_to_U5": order_profile_screen((5, 4, 2, 1, 99), 2),
        "simple_root_profile_U1_to_U5": order_profile_screen((3, 2, 1, 1, 99), 1),
    }
    higher_profile_candidate = solve_k8_double_root_slice(t, R, h, rho)
    dependent_lift = dependent_power_lift(t, R, h, rho)

    result = {
        "scope": "one exact nonzero-U1 all-deep polynomial assignment passing K3-K8 divisibility only",
        "field": "Q",
        "carrier": "h=(t-1)^2(t^2-3t+3)",
        "assignment": {
            "rho": "(t-1)R",
            "A": "rho",
            "B": "h",
            "U1": "h^2*A=(t-1)^5 R^3",
            "U2": "h*B=h^2",
            "U3": "h",
            "U4": "rho",
            "U5_over_t": "0",
        },
        "degree_bounds": degree_bounds,
        "root_orders": root_orders,
        "existing_K9_K24_unique_pole_profile_screens": profile_screens,
        "higher_order_profile_candidate": higher_profile_candidate,
        "exact_dependent_power_lift": dependent_lift,
        "divisibility_quotients": {name: encoded[name] for name in
                                   ("K3_quotient", "K4_quotient", "K5_quotient",
                                    "K6_quotient", "K7_quotient", "K8_quotient")},
        "verified_divisibilities": ["h^4 | K3 numerator", "h^7 | K4 numerator",
                                    "h^10 | K5 numerator", "h | K6 pole numerator",
                                    "h^2 | K7 pole numerator", "h^3 | K8 pole numerator"],
        "interpretation": (
            "The point lies in the nonzero-U1, deeper-contact pattern at both simple roots and at the distinguished double root, and satisfies the exact K3-K8 polynomiality divisibilities. It is only a local necessary-condition survivor: it has not been checked against K9 onward, the later band equations, terminal support rows, or the constant-Jacobian equation, so it is not a Keller pair or a counterexample. It does show that excluding this stratum using K3-K8 divisibility alone cannot work."
        ),
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "EXACT_NONZERO_ALL_DEEP_K3_K8_SURVIVOR",
        "degree_bounds": degree_bounds,
        "root_orders": root_orders,
        "profile_survival_K9_K24": {
            name: value["survives"] for name, value in profile_screens.items()
        },
        "first_profile_exclusion_bands": {
            name: value["first_excluded_band"] for name, value in profile_screens.items()
        },
        "higher_profile_candidate_status": higher_profile_candidate["status"],
        "dependent_lift_degrees_and_jacobian": {
            **dependent_lift["total_degrees"], "Jacobian": dependent_lift["Jacobian"]
        },
        "quotient_degrees": {
            name: (None if expr == 0 else int(sp.Poly(expr, t).degree()))
            for name, expr in {"K3": q3, "K4": q4, "K5": q5,
                               "K6": q6, "K7": q7, "K8": q8}.items()
        },
        "receipt": str(OUT),
        "script_sha256": sha256(Path(__file__)),
        "receipt_sha256": sha256(OUT),
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
