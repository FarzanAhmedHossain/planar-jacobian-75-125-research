#!/usr/bin/env python3
"""Restrict the exact coupled weight-48/58 compatibility ideal to S(1)=0.

This runs only the active canonical 2+1+1 chart. It reconstructs the relaxed
weight-48 P-band and weight-58 Q-band equations from the pinned source,
eliminates their shared new coefficients by exact rational left-kernel
compatibility, and restricts the resulting parameter ideal to the terminal
line c1=-2*c0 forced by ord_1(U1)>=5. It does not run the paused 2+2 chart.
"""

from __future__ import annotations

import hashlib
import json
from functools import reduce
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import QQ

from f2_75_125_power_relation_pq_nextbands import (
    coefficient_system,
    weight49_parameters,
)


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "f2_75_125_211_weight49_coupled_band_terminal_20261021.json"
WEIGHT49_RECEIPT = HERE / "f2_75_125_211_weight49_root_profile_reduction_20261021.json"
ROOTWISE_RECEIPT = HERE / "f2_75_125_211_rootwise_divisor_reduction_20261021.json"
K10_RECEIPT = HERE / "f2_75_125_211_k10_double_root_profile_exclusion_20261021.json"
DEPENDENCIES = [
    HERE / "f2_75_125_power_relation_pq_nextbands.py",
    HERE / "f2_75_125_power_relation_second_band.py",
    WEIGHT49_RECEIPT,
    ROOTWISE_RECEIPT,
    K10_RECEIPT,
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_audit() -> dict:
    t = sp.symbols("t")
    R = t**2 - 3 * t + 3
    h = sp.expand((t - 1) ** 2 * R)
    V = sp.expand(-sp.Rational(1, 3) * t * (t - 1) * R)
    c, S, f38, _kernel = weight49_parameters(h, V, QQ, "2+1+1")
    if len(c) != 2:
        raise AssertionError(f"expected two weight-49 parameters, found {len(c)}")
    c0, c1 = c
    terminal_S = sp.factor(S.subs(c1, -2 * c0))
    if sp.expand(terminal_S + c0 * (t - 1) * R / 3) != 0:
        raise AssertionError("the terminal line disagrees with the exact weight-49 kernel")

    B_coeffs = sp.symbols("B0:8")
    f37_coeffs = sp.symbols("F8:68")
    B = sum(B_coeffs[i] * t**i for i in range(8))
    f37 = sum(f37_coeffs[i - 8] * t**i for i in range(8, 68))
    U = sp.expand(t * h)
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

    new_variables = (*f37_coeffs, *B_coeffs)
    Mp, bp = coefficient_system(p48, t, c, new_variables, "P48 2+1+1")
    Mq, bq = coefficient_system(q58, t, c, new_variables, "Q58 2+1+1")
    M = sp.Matrix.vstack(Mp, Mq)
    forcing = [*bp, *bq]
    M_domain = DomainMatrix.from_Matrix(M).convert_to(QQ)
    rank = M_domain.rank()
    left_kernel = DomainMatrix.from_Matrix(M.T).convert_to(QQ).nullspace().to_Matrix()
    residuals = []
    for i in range(left_kernel.rows):
        residual = sp.expand(sum(left_kernel[i, j] * forcing[j] for j in range(M.rows)))
        if not sp.Poly(residual, *c, domain=QQ).is_zero:
            residuals.append(sp.factor(residual))

    if not residuals:
        line_polynomials = []
        line_gcd = sp.Integer(0)
        groebner_basis = []
    else:
        line_polynomials = [sp.factor(r.subs(c1, -2 * c0)) for r in residuals]
        nonzero_line = [sp.Poly(p, c0, domain=QQ) for p in line_polynomials if p != 0]
        if not nonzero_line:
            line_gcd = sp.Integer(0)
        else:
            gcd_poly = reduce(lambda left, right: left.gcd(right), nonzero_line)
            line_gcd = sp.factor(gcd_poly.as_expr())
        groebner_basis = [str(p) for p in sp.groebner(residuals, *c, domain=QQ)]

    if line_gcd == 0:
        terminal_slice = "all c0 remain compatible with these two bands"
    elif sp.Poly(line_gcd, c0, domain=QQ).degree() == 0:
        terminal_slice = "no point on c1=-2*c0 satisfies the coupled-band residuals"
    elif sp.factor(line_gcd) == c0 or sp.factor(line_gcd) == c0**2:
        terminal_slice = "only c0=0 survives; the nonzero-U1 slice is excluded by these bands"
    elif sp.factor(line_gcd).has(c0):
        terminal_slice = "the nonzero-U1 slice is restricted to roots of the displayed univariate gcd"
    else:
        terminal_slice = "terminal slice requires inspection"

    return {
        "status": "EXACT_COUPLED_BAND_TERMINAL_RESTRICTION",
        "branch": "canonical normalized F2 2+1+1 only",
        "scope": "support-relaxed exact compatibility of weight-48 P-band and weight-58 Q-band after eliminating shared new coefficients; necessary conditions, not a full branch exclusion",
        "carrier": {"h": "(t-1)^2*(t^2-3t+3)", "V": "-(1/3)*t*(t-1)*(t^2-3t+3)"},
        "kernel": {"dimension": len(c), "terminal_equation": "c1=-2*c0", "S_on_terminal_line": str(terminal_S)},
        "linear_system": {
            "P_rows": Mp.rows,
            "Q_rows": Mq.rows,
            "new_variables": M.cols,
            "rank_over_Q": rank,
            "left_kernel_dimension": left_kernel.rows,
            "nonzero_compatibility_residual_count": len(residuals),
        },
        "residuals_in_c0_c1": [str(r) for r in residuals],
        "terminal_line_residuals": [str(r) for r in line_polynomials],
        "terminal_line_gcd": str(line_gcd),
        "terminal_slice_conclusion": terminal_slice,
        "groebner_basis": groebner_basis,
        "whole_2plus1plus1_branch_excluded": False,
        "input_sha256": {path.name: sha256(path) for path in DEPENDENCIES},
    }


if __name__ == "__main__":
    result = exact_audit()
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "linear_system": result["linear_system"],
        "terminal_line_gcd": result["terminal_line_gcd"],
        "terminal_slice_conclusion": result["terminal_slice_conclusion"],
        "groebner_basis": result["groebner_basis"],
        "whole_branch_excluded": result["whole_2plus1plus1_branch_excluded"],
        "receipt": str(OUTPUT),
    }, indent=2), flush=True)
