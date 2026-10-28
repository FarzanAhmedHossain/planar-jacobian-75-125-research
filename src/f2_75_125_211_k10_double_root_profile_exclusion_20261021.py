#!/usr/bin/env python3
"""Exact leading-jet exclusion of one double-root contact profile.

For ord(F1,...,F5)=(5,3,2,0,0) at a nonzero double root of h, the K8
order-four equation forces a quadratic relation on the leading coefficients.
The K10 u^45 coefficient then cannot vanish.  This is a local profile
exclusion, not an exclusion of the full 2+1+1 branch.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from functools import lru_cache
from math import factorial
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
K9_RECEIPT = HERE / "f2_75_125_211_exact_k9_tied_residue_check_20261021.json"
PROFILE_RECEIPT = HERE / "f2_75_125_k17_k24_nonzero_first_slice.json"
JOINT_PROFILE_RECEIPT = HERE / "f2_75_125_211_joint_root_profile_screen.json"
OUT = HERE / "f2_75_125_211_k10_double_root_profile_exclusion_20261021.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


@lru_cache(maxsize=None)
def partitions(total: int, max_part: int) -> tuple[tuple[int, ...], ...]:
    out: list[tuple[int, ...]] = []

    def rec(weight: int, left: int, current: list[int]) -> None:
        if weight > max_part:
            if left == 0:
                out.append(tuple(current))
            return
        for multiplicity in range(left // weight + 1):
            current.append(multiplicity)
            rec(weight + 1, left - multiplicity * weight, current)
            current.pop()

    rec(1, total, [])
    return tuple(out)


def generalized_binomial(alpha: Fraction, n: int) -> Fraction:
    value = Fraction(1)
    for i in range(n):
        value *= alpha - i
    return value / factorial(n)


def partition_coefficient(alpha: Fraction, mult: tuple[int, ...]) -> Fraction:
    count = sum(mult)
    value = generalized_binomial(alpha, count) * factorial(count)
    for multiplicity in mult:
        value /= factorial(multiplicity)
    return value


def label(mult: tuple[int, ...]) -> str:
    return "+".join(str(weight)
                     for weight, multiplicity in enumerate(mult, 1)
                     for _ in range(multiplicity))


def main() -> None:
    k9_receipt = json.loads(K9_RECEIPT.read_text(encoding="utf-8"))
    profile_receipt = json.loads(PROFILE_RECEIPT.read_text(encoding="utf-8"))
    joint_receipt = json.loads(JOINT_PROFILE_RECEIPT.read_text(encoding="utf-8"))
    known_orders = {1: 5, 2: 3, 3: 2, 4: 0, 5: 0}
    h_order = 2

    c, b, alpha, xi = sp.symbols("c b alpha xi", nonzero=True)
    a, e, f, u = sp.symbols("a e f u")
    d = sp.symbols("d")  # U4(alpha) may vanish; that case is ruled out by K8.

    # Reconstruct K8 from the exact polynomial numerator using only the
    # leading local jets. Higher jets have order >4 and cannot affect [u^4].
    h_local = c * u**2
    U1_local = a * u**5  # a=0 also covers ord(U1)>5 or U1=0 locally.
    U2_local = b * u**3
    U3_local = e * u**2
    U4_local = d
    U5_band_local = alpha * f  # t(alpha)=alpha is a unit; only order matters.
    A = U1_local / h_local**2
    B = U2_local / h_local
    D = 3 * B - A**2
    Delta = 27 * U3_local - A**3 - 3 * A * D
    J8 = A * Delta - 81 * U4_local * h_local
    N8 = sp.Rational(5, 59049) * (
        J8**2 - D * Delta**2 + 6 * A * D**2 * Delta
        - 3 * A**2 * D**3 + 3 * D**4 - 243 * D**2 * U4_local * h_local
        + 486 * Delta * U5_band_local * h_local**2
    )
    k8_u4 = sp.factor(sp.Poly(sp.expand(N8), u).nth(4))
    B1 = b / c  # leading coefficient of U2/h
    expected_k8_u4 = sp.factor(sp.Rational(5, 59049) * (
        6561 * d**2 * c**2 + 243 * B1**4 - 2187 * B1**2 * d * c
    ))
    if sp.simplify(k8_u4 - expected_k8_u4) != 0:
        raise AssertionError("the reconstructed K8 u^4 coefficient is incorrect")
    xi_relation = sp.factor(27 * xi**2 - 9 * xi + 1)
    expected_k8 = sp.factor(sp.Rational(5, 243) * b**4 / c**4
                             * xi_relation.subs(xi, c**3 * d / b**2))
    if sp.simplify(k8_u4 - expected_k8) != 0:
        raise AssertionError("the K8 leading coefficient did not reduce to the xi relation")
    k8_u4_if_d_zero = sp.simplify(k8_u4.subs(d, 0))
    if sp.simplify(k8_u4_if_d_zero - sp.Rational(5, 243) * b**4 / c**4) != 0:
        raise AssertionError("the d=0 K8 leading coefficient should be nonzero")

    # Enumerate every K10 main partition by valuation.  The numerator is
    # multiplied by t*h^25; t(alpha)=alpha is a unit because alpha != 0.
    main_rows: list[dict[str, object]] = []
    main_at_45 = []
    for mult in partitions(10, 10):
        if sum(weight * m for weight, m in enumerate(mult, 1)) != 10:
            continue
        count = sum(mult)
        cleared_order = h_order * (25 + 5 - 3 * count)
        cleared_order += sum(m * known_orders[weight]
                             for weight, m in enumerate(mult, 1) if weight <= 5)
        row: dict[str, object] = {
            "partition": label(mult),
            "cleared_numerator_order_lower_bound": cleared_order,
        }
        if cleared_order == 45:
            coeff = partition_coefficient(Fraction(5, 3), mult)
            row["generalized_binomial_partition_coefficient"] = str(coeff)
            main_at_45.append((label(mult), coeff))
        main_rows.append(row)
    expected_partitions = ["2+4+4", "2+2+2+4", "2+2+2+2+2"]
    if [row[0] for row in main_at_45] != expected_partitions:
        raise AssertionError(f"unexpected K10 leading partitions: {main_at_45}")

    # Reconstruct [u^45] N10 from the three surviving leading partitions.
    leading_slices = {1: a, 2: b, 3: e, 4: d, 5: f}
    k10_u45_terms = []
    for name, coeff in main_at_45:
        mult = next(row for row in partitions(10, 10)
                    if sum(weight * m for weight, m in enumerate(row, 1)) == 10
                    and "+".join(str(weight)
                                 for weight, multiplicity in enumerate(row, 1)
                                 for _ in range(multiplicity)) == name)
        count = sum(mult)
        term = alpha * coeff * c ** (25 + 5 - 3 * count)
        for weight, multiplicity in enumerate(mult, 1):
            if multiplicity:
                if weight not in leading_slices:
                    raise AssertionError("an unknown slice entered the u^45 leading coefficient")
                term *= leading_slices[weight] ** multiplicity
        k10_u45_terms.append(term)
    k10_u45 = sp.factor(sum(k10_u45_terms))
    expected_k10 = sp.factor(alpha * b**5 * c**15 / 729
                              * (-135 * xi**2 + 60 * xi - 7)
                              .subs(xi, c**3 * d / b**2))
    if sp.simplify(k10_u45 - expected_k10) != 0:
        raise AssertionError("the K10 leading coefficient formula failed")
    reduced_factor = sp.rem(-135 * xi**2 + 60 * xi - 7,
                            xi_relation, xi).factor()
    if sp.simplify(reduced_factor - (15 * xi - 2)) != 0:
        raise AssertionError("the K8 relation did not reduce the K10 factor to 15*xi-2")

    # Unknown-term valuation checks. U6 may start at order 0; U7--U10 are
    # checked directly from all weighted partitions, not inferred by sampling.
    u6_rows: list[dict[str, object]] = []
    for mult in partitions(10, 10):
        if not (mult[5] and not any(mult[6:])):
            continue
        count = sum(mult)
        order = h_order * (25 + 5 - 3 * count)
        order += sum(m * known_orders[weight]
                     for weight, m in enumerate(mult, 1) if weight <= 5)
        u6_rows.append({"partition": label(mult), "cleared_order_lower_bound": order})
    if len(u6_rows) != 5 or min(row["cleared_order_lower_bound"] for row in u6_rows) != 48:
        raise AssertionError("unexpected U6 partition valuation")

    u7_rows: list[dict[str, object]] = []
    for mult in partitions(10, 10):
        if not any(mult[6:]):
            continue
        count = sum(mult)
        original_order = h_order * (5 - 3 * count)
        original_order += sum(m * known_orders[weight]
                              for weight, m in enumerate(mult, 1) if weight <= 5)
        u7_rows.append({"partition": label(mult),
                        "preclearing_order_lower_bound": original_order,
                        "cleared_order_lower_bound": original_order + 50})
    if len(u7_rows) != 7 or min(row["preclearing_order_lower_bound"]
                                for row in u7_rows) != 0:
        raise AssertionError("unexpected U7--U10 partition valuation")

    kappa1_rows: list[dict[str, object]] = []
    for mult in partitions(5, 5):
        if sum(weight * m for weight, m in enumerate(mult, 1)) != 5:
            continue
        count = sum(mult)
        order = h_order * (25 + 4 - 3 * count)
        order += sum(m * known_orders[weight]
                     for weight, m in enumerate(mult, 1))
        kappa1_rows.append({"partition": label(mult), "cleared_order_lower_bound": order})
    if min(row["cleared_order_lower_bound"] for row in kappa1_rows) != 51:
        raise AssertionError("unexpected kappa_1 valuation")

    xi_roots = [sp.Rational(1, 6) + sp.sqrt(-3) / 18,
                sp.Rational(1, 6) - sp.sqrt(-3) / 18]
    residuals = [sp.simplify(15 * root - 2) for root in xi_roots]
    if any(value == 0 for value in residuals):
        raise AssertionError("a K8 root unexpectedly annihilates the K10 factor")

    if joint_receipt["input_receipts"][PROFILE_RECEIPT.name] != sha256(PROFILE_RECEIPT):
        raise AssertionError("double-root order grid does not match the joint-screen input hash")
    survivor_profiles = [
        tuple(99 if value == "infinity" or int(value) >= 13 else int(value)
              for value in row)
        for row in profile_receipt["survivor_profiles"]
    ]
    profile_grid_total = len(survivor_profiles)
    excluded_profile_grid = [
        row for row in survivor_profiles
        if row[0] >= 5 and row[1] == 3 and row[2] >= 2
    ]
    if profile_grid_total != joint_receipt["double_root_survivor_counts"]["nonzero_U1_after_K24"]:
        raise AssertionError("double-root profile grid count disagrees with its joint-screen receipt")
    if len(excluded_profile_grid) != 462:
        raise AssertionError("unexpected number of profile-grid rows covered by the local theorem")

    result = {
        "status": "EXACT_DOUBLE_ROOT_PROFILE_EXCLUDED_AT_K10",
        "scope": "all local coefficient data with ord_u(F1)>=5, ord_u(F2)=3, ord_u(F3)>=2 at a nonzero double root; F4 and F5 are arbitrary regular slices; this is not the full 2+1+1 branch",
        "source_K10_valuation_receipt_sha256": sha256(K9_RECEIPT),
        "local_setup": "u=t-alpha; ord(h)=2; h=c*u^2+O(u^3), c!=0; alpha!=0; b=lc(U2,u^3)!=0; d=U4(alpha) is arbitrary",
        "K8_divisibility_consequence": {
            "N8_u4_coefficient": str(k8_u4),
            "if_d_zero": str(k8_u4_if_d_zero),
            "d_zero_case": "if d=0, this is nonzero because b and c are nonzero; K8 already fails",
            "xi_definition": "xi=c^3*d/b^2",
            "forced_relation_if_d_nonzero": "27*xi^2-9*xi+1=0",
            "roots": [str(root) for root in xi_roots],
        },
        "K10_main_partition_valuation_census": main_rows,
        "partitions_at_order_45": [
            {"partition": name, "partition_coefficient": str(coeff)}
            for name, coeff in main_at_45
        ],
        "N10_u45_before_K8_reduction": str(k10_u45),
        "N10_u45_after_K8_reduction_if_d_nonzero": "alpha*b^5*c^15/729*(15*xi-2)",
        "15xi_minus_2_at_K8_roots": [str(value) for value in residuals],
        "U6_partition_census": u6_rows,
        "kappa1_partition_census": kappa1_rows,
        "kappa2_cleared_order": 56,
        "U7_to_U10_partition_census": u7_rows,
        "polynomiality_requirement": "ord_u(N10)>=50 since h^25 has order 50",
        "double_root_profile_grid_audit": {
            "source_receipt": PROFILE_RECEIPT.name,
            "source_receipt_sha256": sha256(PROFILE_RECEIPT),
            "surviving_nonzero_U1_double_root_profiles_after_K24": profile_grid_total,
            "profiles_excluded_by_ord_F1_ge_5_F2_eq_3_F3_ge_2": len(excluded_profile_grid),
            "excluded_percent_of_this_finite_double_root_profile_grid": round(
                100 * len(excluded_profile_grid) / profile_grid_total, 6
            ),
            "interpretation": "finite valuation-grid fraction only, not a percentage of the algebraic branch or coefficient variety",
        },
        "conclusion": "if d=0 K8 fails; if d!=0 K8 forces xi to one of two nonreal roots, while the coefficient of u^45 in N10 is nonzero. Thus the entire stated local contact-profile stratum has no supported K10 completion",
        "branch_coverage": "whole 2+1+1 branch is not excluded; branch-wide exclusion status remains 0% until the other local profiles/components are covered",
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({**result, "receipt": str(OUT),
                      "script_sha256": sha256(Path(__file__)),
                      "receipt_sha256": sha256(OUT)}, indent=2), flush=True)


if __name__ == "__main__":
    main()
