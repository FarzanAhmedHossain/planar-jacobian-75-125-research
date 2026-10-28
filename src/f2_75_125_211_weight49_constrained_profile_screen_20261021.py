#!/usr/bin/env python3
"""Re-screen the canonical 2+1+1 root-order grid after the exact weight-49 slice.

The exact weight-49 kernel and terminal order condition force either
U1=0 at all three roots (the DD root profile) or
U1=-(c0/3)(t-1)^5(t^2-3t+3)^3 with c0 nonzero. In the latter case the
orders are 5 at the double root and 3 at each simple root. This checker
replays that kernel intersection, applies the existing exact K9--K24
single-negative-pole screens to the forced simple-root order profiles, and
joins survivors by polynomial degree caps. It is a necessary local-profile
screen, not a whole-branch exclusion unless a necessary profile set is empty.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ

import f2_75_125_211_joint_root_profile_screen as joint
from f2_75_125_power_relation_pq_nextbands import weight49_parameters


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "f2_75_125_211_weight49_constrained_profile_screen_20261021.json"
WEIGHT49_RECEIPT = HERE / "f2_75_125_211_weight49_root_profile_reduction_20261021.json"
ROOTWISE_RECEIPT = HERE / "f2_75_125_211_rootwise_divisor_reduction_20261021.json"
K10_CERTIFICATE = HERE / "f2_75_125_211_k10_double_root_profile_exclusion_20261021.py"
K10_RECEIPT = HERE / "f2_75_125_211_k10_double_root_profile_exclusion_20261021.json"

DEPENDENCIES = [
    HERE / "f2_75_125_power_relation_pq_nextbands.py",
    HERE / "f2_75_125_power_relation_second_band.py",
    HERE / "f2_75_125_211_joint_root_profile_screen.py",
    HERE / "f2_75_125_k17_k24_nonzero_first_slice.json",
    HERE / "f2_75_125_k15_k24_zero_u1_double_root.json",
    HERE / "f2_75_125_k6_local_jet_certificate.py",
    HERE / "f2_75_125_k7_double_root_certificate.py",
    HERE / "f2_75_125_k8_local_jet_certificate.py",
    K10_CERTIFICATE,
    K10_RECEIPT,
    WEIGHT49_RECEIPT,
    ROOTWISE_RECEIPT,
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def check_exact_weight49_slice() -> dict:
    t = sp.symbols("t")
    c0, c1 = sp.symbols("c0 c1")
    R = t**2 - 3 * t + 3
    h = sp.expand((t - 1) ** 2 * R)
    V = sp.expand(-sp.Rational(1, 3) * t * (t - 1) * R)
    parameters, S, _f38, _kernel = weight49_parameters(h, V, QQ, "2+1+1")
    expected = (
        20 * c0 * t**3 - 40 * c0 * t**2 + 42 * c0
        + 17 * c1 * t**3 - 48 * c1 * t**2 + 42 * c1 * t
    ) / 42
    if len(parameters) != 2 or tuple(parameters) != (c0, c1):
        raise AssertionError(f"unexpected kernel parameters: {parameters}")
    if sp.expand(S - expected) != 0:
        raise AssertionError(f"weight-49 projection changed: {sp.factor(S)}")
    if sp.factor(S.subs(t, 1) - sp.Rational(11, 42) * (2 * c0 + c1)) != 0:
        raise AssertionError("terminal value S(1) changed")
    reduced = sp.factor(S.subs(c1, -2 * c0))
    expected_reduced = -c0 * (t - 1) * R / 3
    if sp.expand(reduced - expected_reduced) != 0:
        raise AssertionError("terminal specialization no longer gives -(c0/3)(t-1)R")
    U1 = sp.factor(h**2 * reduced)
    expected_U1 = -c0 * (t - 1) ** 5 * R**3 / 3
    if sp.expand(U1 - expected_U1) != 0:
        raise AssertionError("U1 factorization does not match the exact kernel slice")
    return {
        "kernel_dimension": len(parameters),
        "S_before_terminal": str(sp.factor(S)),
        "terminal_equation": "c1=-2*c0",
        "S_after_terminal": str(reduced),
        "U1_after_terminal": str(U1),
        "c0_nonzero_root_orders": {"double_root": 5, "simple_roots": [3, 3]},
        "c0_zero_root_orders": "U1=0 at all three roots",
    }


def run() -> dict:
    weight49 = check_exact_weight49_slice()
    nonzero_data, _zero_data, double_nonzero_all, double_zero_all = joint.decode_double_receipts()
    rootwise = json.loads(ROOTWISE_RECEIPT.read_text(encoding="utf-8"))
    k10_receipt = json.loads(K10_RECEIPT.read_text(encoding="utf-8"))
    for name, expected_hash in rootwise.get("root_jet_dependencies", {}).items():
        dependency = HERE / name
        if not dependency.is_file() or sha256(dependency) != expected_hash:
            raise AssertionError(f"rootwise certificate dependency changed: {name}")
    if k10_receipt.get("status") != "EXACT_DOUBLE_ROOT_PROFILE_EXCLUDED_AT_K10":
        raise AssertionError("the pinned exact K10 double-root certificate is unavailable")
    if rootwise.get("K10_profile_exclusion_receipt", {}).get("sha256") != sha256(K10_RECEIPT):
        raise AssertionError("rootwise divisor receipt does not pin the current K10 certificate")
    required_root_condition = (
        "double root: ord(U1)>=5, ord(U2)>=4, U3(1)=0; "
        "K8 and K10 jointly exclude every ord(U2)=3 profile"
    )
    if required_root_condition not in rootwise.get("local_conditions_used", []):
        raise AssertionError("K6-K10 rootwise divisor conditions changed")

    # The polynomial degree cap deg(U1)<=11 and root orders (3,3) at the
    # simple roots force exactly order 5 at the double root; filter explicitly.
    double_fixed_order = [profile for profile in double_nonzero_all if profile[0] == 5]
    # K6--K10's exact double-root divisor receipt forces ord(U2)>=4 and
    # U3(1)=0 in every surviving tuned/deep profile. In particular, all
    # ord(U2)=3 data are gone (K8 handles ord(U3)=1; K10 handles >=2).
    double_fixed = [
        profile for profile in double_fixed_order
        if profile[1] >= 4 and profile[2] >= 1
    ]
    double_zero = [
        profile for profile in double_zero_all
        if profile[1] >= 4 and profile[2] >= 1
    ]
    simple_bands = joint.build_bands(root_order=1)

    simple_deeper_all = joint.simple_profiles("deeper")
    simple_fixed_initial = [profile for profile in simple_deeper_all if profile[0] == 3]
    simple_fixed, nonzero_band_rows = joint.screen_profiles(
        simple_fixed_initial, simple_bands, "weight49 nonzero-U1 simple root (ord U1=3)"
    )

    simple_zero_initial = joint.simple_profiles("zero_u1")
    simple_zero, zero_band_rows = joint.screen_profiles(
        simple_zero_initial, simple_bands, "weight49 zero-U1 simple root"
    )

    nonzero_joint = joint.query_pair(
        double_fixed, simple_fixed, simple_fixed,
        "weight49 nonzero-U1: double order 5, simple orders 3 and 3",
    )
    zero_joint = joint.query_pair(
        double_zero, simple_zero, simple_zero,
        "weight49 U1=0: DD at both simple roots",
    )

    profiles = {tuple(row["profile_names"]): row for row in rootwise["profiles"]}
    if set(profiles) != {("TT",), ("DT", "TD"), ("DD",)}:
        raise AssertionError("rootwise receipt profile partition changed")

    result = {
        "status": "EXACT_WEIGHT49_CONSTRAINED_ROOT_PROFILE_SCREEN",
        "branch": "canonical normalized F2 2+1+1 carrier only",
        "coverage_warning": "Finite exact necessary root-order screen; counts are not algebraic-branch percentages.",
        "weight49_slice": weight49,
        "rootwise_profile_reduction": {
            "TT_DT_TD_excluded": True,
            "surviving_root_profile": "DD",
            "discrete_root_profile_strata_excluded": 3,
            "discrete_root_profile_strata_total": 4,
            "fraction": "3/4; not branch coverage",
            "pre_kernel_quotient_dimensions_U1_to_U4": {
                str(key): profiles[key]["total_U1_to_U4_coefficients"]
                for key in profiles
            },
            "post_kernel_surviving_DD_quotient_coefficients": profiles[("DD",)]["total_U1_to_U4_coefficients"],
            "double_root_jet_conditions": "ord(U1)>=5, ord(U2)>=4, ord(U3)>=1; exact K6-K10 receipt",
        },
        "degree_caps": list(joint.DEGREE_CAPS),
        "nonzero_U1_slice": {
            "double_root_profiles_before_order5_filter": len(double_nonzero_all),
            "double_root_profiles_after_order5_filter": len(double_fixed_order),
            "double_root_profiles_after_K6_K10_divisibility": len(double_fixed),
            "double_root_profiles_with_order5": len(double_fixed),
            "simple_root_profiles_before_K9_K24": len(simple_fixed_initial),
            "simple_root_profiles_after_K9_K24": len(simple_fixed),
            "simple_root_profile_excluded_percent": round(
                100 * (len(simple_fixed_initial) - len(simple_fixed)) / len(simple_fixed_initial), 6
            ) if simple_fixed_initial else 0.0,
            "new_exclusions_by_band": nonzero_band_rows,
            "degree_sum_join": nonzero_joint,
        },
        "zero_U1_slice": {
            "double_root_profiles_after_K24_before_K6_K10_filter": len(double_zero_all),
            "double_root_profiles_after_K6_K10_divisibility": len(double_zero),
            "simple_root_profiles_before_K9_K24": len(simple_zero_initial),
            "simple_root_profiles_after_K9_K24": len(simple_zero),
            "simple_root_profile_excluded_percent": round(
                100 * (len(simple_zero_initial) - len(simple_zero)) / len(simple_zero_initial), 6
            ) if simple_zero_initial else 0.0,
            "new_exclusions_by_band": zero_band_rows,
            "degree_sum_join": zero_joint,
        },
        "whole_2plus1plus1_branch_excluded": bool(
            not nonzero_joint["compatible"] and not zero_joint["compatible"]
        ),
        "whole_branch_coverage_percent": 100.0 if (
            not nonzero_joint["compatible"] and not zero_joint["compatible"]
        ) else 0.0,
        "input_sha256": {path.name: sha256(path) for path in DEPENDENCIES},
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    output = run()
    summary = {
        "status": output["status"],
        "root_profile_split": output["rootwise_profile_reduction"],
        "nonzero_U1": {
            "simple_initial": output["nonzero_U1_slice"]["simple_root_profiles_before_K9_K24"],
            "simple_after": output["nonzero_U1_slice"]["simple_root_profiles_after_K9_K24"],
            "degree_compatible": output["nonzero_U1_slice"]["degree_sum_join"]["compatible"],
            "witness": output["nonzero_U1_slice"]["degree_sum_join"]["witness"],
        },
        "zero_U1": {
            "simple_initial": output["zero_U1_slice"]["simple_root_profiles_before_K9_K24"],
            "simple_after": output["zero_U1_slice"]["simple_root_profiles_after_K9_K24"],
            "degree_compatible": output["zero_U1_slice"]["degree_sum_join"]["compatible"],
            "witness": output["zero_U1_slice"]["degree_sum_join"]["witness"],
        },
        "whole_branch_excluded": output["whole_2plus1plus1_branch_excluded"],
        "whole_branch_coverage_percent": output["whole_branch_coverage_percent"],
        "receipt": str(OUTPUT),
    }
    print(json.dumps(summary, indent=2), flush=True)
