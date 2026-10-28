#!/usr/bin/env python3
"""Intersect exact U1=0 K10 local obstructions with saved F2 profile grids.

This combines the already-derived K9--K24 necessary order profiles with the
exact single-negative-pole cases in f2_75_125_k10_zero_first_slice.py. It
then reapplies the global degree-sum screen from the 2+1+1 joint-root script.
It is a finite order-profile screen only; tied-pole cases are retained and no
coefficient compatibility between distinct roots is asserted.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import f2_75_125_211_joint_root_profile_screen as joint  # noqa: E402


ZERO_DOUBLE = HERE / "f2_75_125_k15_k24_zero_u1_double_root.json"
JOINT_RECEIPT = HERE / "f2_75_125_211_joint_root_profile_screen.json"
JOINT_SCRIPT = HERE / "f2_75_125_211_joint_root_profile_screen.py"
K10_SCRIPT = HERE / "f2_75_125_k10_zero_first_slice.py"
OUT = HERE / "f2_75_125_211_zero_u1_k10_profile_coverage_20261021.json"
INF = 99


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def single_simple_pole(profile: tuple[int, ...]) -> tuple[str, ...]:
    """Return uniquely leading K10 negative-pole terms at a simple root."""
    _, q, r, s, v = profile
    orders = {
        "2+2+3+3": (-7 + 2 * q + 2 * r, (q, r)),
        "2+3+5": (-4 + q + r + v, (q, r, v)),
        "3+3+4": (-4 + 2 * r + s, (r, s)),
        "5+5": (-1 + 2 * v, (v,)),
    }
    negative = tuple(sorted(
        term for term, (order, variables) in orders.items()
        if order == -1 and all(value < INF for value in variables)
    ))
    return negative if len(negative) == 1 else ()


def double_root_k10_exclusion(profile: tuple[int, ...]) -> str | None:
    """Exact unique-pole cases established by the U1=0 K10 census."""
    _, q, _r, s, v = profile
    if q == 3 and s >= 1:
        return "q=3,s>=1: unique 5+5 negative pole"
    if q >= 4 and v == 0 and s >= 2:
        return "q>=4,v=0,s>=2: unique 5+5 negative pole"
    if q == 4 and s <= 1:
        return "q=4,s<=1: unique 2+4+4 negative pole"
    return None


def main() -> None:
    joint_receipt = json.loads(JOINT_RECEIPT.read_text(encoding="utf-8"))
    zero_receipt = json.loads(ZERO_DOUBLE.read_text(encoding="utf-8"))
    _nonzero_data, _zero_data, _double_nonzero, double_zero = joint.decode_double_receipts()

    simple_bands = joint.build_bands(root_order=1)
    initial_simple = joint.simple_profiles("zero_u1")
    simple_after_k24, _bands = joint.screen_profiles(
        initial_simple, simple_bands, "zero-U1 simple-root K9-K24"
    )

    expected_simple_count = joint_receipt["simple_root_local_screens"]["zero_u1"][
        "survivors_after_K24"
    ]
    if len(simple_after_k24) != expected_simple_count:
        raise AssertionError("reconstructed simple-root profile count differs from receipt")
    if len(double_zero) != joint_receipt["double_root_survivor_counts"]["U1_zero_after_K24"]:
        raise AssertionError("reconstructed double-root profile count differs from receipt")
    if len(double_zero) != len(zero_receipt["survivor_profiles"]):
        raise AssertionError("joint and zero-U1 receipts disagree on the profile count")

    baseline = joint.query_pair(
        double_zero, simple_after_k24, simple_after_k24, "U1=0 baseline K9-K24"
    )
    recorded_baseline = joint_receipt["joint_degree_compatibility"][-1]
    if baseline["degree_compatible_double_root_profiles"] != recorded_baseline[
        "degree_compatible_double_root_profiles"
    ]:
        raise AssertionError("baseline degree-sum replay differs from the saved joint receipt")

    simple_excluded = [profile for profile in simple_after_k24 if single_simple_pole(profile)]
    simple_retained = [profile for profile in simple_after_k24 if not single_simple_pole(profile)]
    double_reasons = {
        profile: double_root_k10_exclusion(profile) for profile in double_zero
    }
    double_excluded = [profile for profile, reason in double_reasons.items() if reason]
    double_retained = [profile for profile, reason in double_reasons.items() if not reason]

    after_k10 = joint.query_pair(
        double_retained, simple_retained, simple_retained,
        "U1=0 after exact K10 unique-pole screens",
    )
    result = {
        "status": (
            "U1_ZERO_ORDER_PROFILE_BRANCH_EXCLUDED"
            if not after_k10["compatible"] else "U1_ZERO_ORDER_PROFILES_SURVIVE"
        ),
        "scope": (
            "Finite K9-K24 local order profiles for the 2+1+1 U1=0 stratum, "
            "filtered by exact K10 unique-negative-pole conditions at the two simple roots "
            "and the distinguished double root, then by global degree sums. Tied-pole cases "
            "are retained; this does not impose coefficient compatibility between roots."
        ),
        "input_sha256": {
            ZERO_DOUBLE.name: sha256(ZERO_DOUBLE),
            JOINT_RECEIPT.name: sha256(JOINT_RECEIPT),
            JOINT_SCRIPT.name: sha256(JOINT_SCRIPT),
            K10_SCRIPT.name: sha256(K10_SCRIPT),
        },
        "finite_profile_counts": {
            "double_root_after_K24": len(double_zero),
            "double_root_K10_unique_pole_excluded": len(double_excluded),
            "double_root_retained": len(double_retained),
            "double_root_exclusion_percent": round(100 * len(double_excluded) / len(double_zero), 6),
            "simple_root_after_K24": len(simple_after_k24),
            "simple_root_K10_unique_pole_excluded": len(simple_excluded),
            "simple_root_retained": len(simple_retained),
            "simple_root_exclusion_percent": round(100 * len(simple_excluded) / len(simple_after_k24), 6),
        },
        "degree_sum_screen": {
            "baseline": baseline,
            "after_K10_unique_poles": after_k10,
        },
        "interpretation": (
            "These percentages refer only to finite order-profile grids. Even if the grid were "
            "empty, it would certify only this U1=0 valuation stratum under the recorded local "
            "screens, not the complete 2+1+1 branch or the (75,125) degree pair."
        ),
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({**result, "receipt": str(OUT), "receipt_sha256": sha256(OUT)}, indent=2))


if __name__ == "__main__":
    main()
