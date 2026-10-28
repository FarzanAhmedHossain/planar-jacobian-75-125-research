#!/usr/bin/env python3
"""Intersect the 2+1+1 root-order split with its exact weight-49 kernel.

For the canonical carrier h=(t-1)^2(t^2-3t+3) and the exact edge slice
V=-(1/3)t(t-1)(t^2-3t+3), the weight-49 kernel expresses U1=h^2 S in two
parameters. The terminal condition ord_1(U1)>=5 imposes S(1)=0. This checker
replays the exact kernel, proves that nonzero U1 is tuned at both simple roots,
and that the zero case is deep at both. It eliminates three of the four
root-order strata in this normalized model; it does not exclude the full
2+1+1 branch.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ

from f2_75_125_power_relation_pq_nextbands import weight49_parameters


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "f2_75_125_211_weight49_root_profile_reduction_20261021.json"
ROOTWISE_RECEIPT = HERE / "f2_75_125_211_rootwise_divisor_reduction_20261021.json"
DEPENDENCIES = [
    HERE / "f2_75_125_power_relation_pq_nextbands.py",
    HERE / "f2_75_125_power_relation_second_band.py",
    HERE / "f2_75_125_sparse_ring_lift.py",
    ROOTWISE_RECEIPT,
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    t = sp.symbols("t")
    R = t**2 - 3 * t + 3
    h = sp.expand((t - 1) ** 2 * R)
    V = sp.expand(-sp.Rational(1, 3) * t * (t - 1) * R)
    parameters, S, _f38, _kernel_matrix = weight49_parameters(
        h, V, QQ, "2+1+1"
    )
    if len(parameters) != 2:
        raise AssertionError(f"expected a 2-dimensional weight-49 kernel, got {len(parameters)}")
    c0, c1 = parameters
    expected_S = (
        20 * c0 * t**3 - 40 * c0 * t**2 + 42 * c0
        + 17 * c1 * t**3 - 48 * c1 * t**2 + 42 * c1 * t
    ) / 42
    if sp.expand(S - expected_S) != 0:
        raise AssertionError(
            "weight-49 kernel projection to S changed: "
            f"S={sp.factor(S)}, expected={sp.factor(expected_S)}, "
            f"difference={sp.factor(S - expected_S)}"
        )

    terminal_value = sp.factor(S.subs(t, 1))
    if sp.factor(terminal_value - sp.Rational(11, 42) * (2 * c0 + c1)) != 0:
        raise AssertionError("wrong exact terminal S(1) row")
    terminal_substitution = {c1: -2 * c0}
    reduced_S = sp.factor(S.subs(terminal_substitution))
    expected_reduced_S = -c0 * (t - 1) * R / 3
    if sp.expand(reduced_S - expected_reduced_S) != 0:
        raise AssertionError("terminal specialization did not produce the claimed factorization")

    numerator = sp.together(reduced_S).as_numer_denom()[0]
    quotient, remainder = sp.div(sp.Poly(numerator, t), sp.Poly(R, t))
    if not remainder.is_zero or sp.factor(quotient.as_expr() + c0 * (t - 1)) != 0:
        raise AssertionError("the terminal U1 slice is not divisible by R as claimed")
    quotient_at_double_root = sp.factor(reduced_S / (t - 1)).subs(t, 1)
    if sp.factor(quotient_at_double_root + sp.Rational(1, 3) * c0) != 0:
        raise AssertionError("the factor S/(t-1) has the wrong value at t=1")

    rootwise = json.loads(ROOTWISE_RECEIPT.read_text(encoding="utf-8"))
    profiles = {
        tuple(row["profile_names"]): row for row in rootwise["profiles"]
    }
    tt = profiles[("TT",)]
    mixed = profiles[("DT", "TD")]
    dd = profiles[("DD",)]
    if tt["total_U1_to_U4_coefficients"] != 32:
        raise AssertionError("unexpected TT rootwise coefficient count")
    if mixed["total_U1_to_U4_coefficients"] != 28:
        raise AssertionError("unexpected mixed rootwise coefficient count")
    if dd["total_U1_to_U4_coefficients"] != 24:
        raise AssertionError("unexpected DD rootwise coefficient count")
    if tt["free_quotient_coefficient_counts"]["U1"] != 3:
        raise AssertionError("unexpected TT U1 quotient dimension")
    if dd["free_quotient_coefficient_counts"]["U1"] != 1:
        raise AssertionError("unexpected DD U1 quotient dimension")
    if sp.factor((reduced_S / (t - 1)).subs(t, 1) + sp.Rational(1, 3) * c0) != 0:
        raise AssertionError("ord_1(S) is not exactly one on the nonzero slice")
    dd_after = dd["total_U1_to_U4_coefficients"]

    result = {
        "status": "EXACT_WEIGHT49_U1_ROOT_PROFILE_REDUCTION",
        "branch": "2+1+1",
        "scope": "canonical normalized carrier h=(t-1)^2(t^2-3t+3), V=-(1/3)t(t-1)(t^2-3t+3); necessary reduction only",
        "carrier": {
            "R": "t^2-3t+3",
            "h": "(t-1)^2 R",
            "V": "-(1/3)t(t-1)R",
        },
        "weight49_kernel_dimension": len(parameters),
        "U1_definition": "U1=h^2*S",
        "S_before_terminal": "(20*c0*t^3-40*c0*t^2+42*c0+17*c1*t^3-48*c1*t^2+42*c1*t)/42",
        "terminal_condition": "ord_1(U1)>=5 iff S(1)=0 iff c1=-2*c0",
        "S_after_terminal": "-(c0/3)*(t-1)*R(t)",
        "numerator_S_over_R": "-c0*(t-1)",
        "root_profile_conclusion": {
            "c0_nonzero": "ord_1(U1)=5 and ord_alpha(U1)=ord_beta(U1)=3; K6-K8 force DD",
            "c0_zero": "c1=0, S=0, U1=0; K6-K8 force DD",
            "TT_DT_TD": "excluded by the exact weight-49 kernel plus S(1)=0",
        },
        "rootwise_quotient_coefficients_before_weight49": {
            "TT": 32,
            "DT": 28,
            "TD": 28,
            "DD": 24,
        },
        "rootwise_quotient_coefficients_after_weight49": {
            "DD_all_U1_cases": dd_after,
            "DD_zero_U1_slice": dd_after - dd["free_quotient_coefficient_counts"]["U1"],
            "TT_DT_TD": 0,
        },
        "discrete_root_profile_strata_excluded": 3,
        "discrete_root_profile_strata_total": 4,
        "discrete_stratum_fraction_excluded": "3/4 (75%; not an algebraic-branch percentage)",
        "whole_2plus1plus1_branch_excluded": False,
        "whole_2plus1plus1_branch_exclusion_percent": 0.0,
        "dependency_sha256": {path.name: sha256(path) for path in DEPENDENCIES},
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "S_after_terminal": result["S_after_terminal"],
        "numerator_S_over_R": result["numerator_S_over_R"],
        "surviving_strata": result["root_profile_conclusion"],
        "coefficient_totals": result["rootwise_quotient_coefficients_after_weight49"],
        "discrete_profile_strata_excluded": "3/4; not a branch percentage",
        "whole_branch_coverage_percent": result["whole_2plus1plus1_branch_exclusion_percent"],
        "receipt": OUTPUT.name,
    }, indent=2))


if __name__ == "__main__":
    main()
