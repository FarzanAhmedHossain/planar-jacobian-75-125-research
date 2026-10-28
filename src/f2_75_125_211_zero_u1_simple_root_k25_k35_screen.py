#!/usr/bin/env python3
"""Extend the U1=0 simple-root valuation census from K24 through K35.

This screens the 18,876 exact order profiles surviving the recorded K9--K24
simple-root pass.  A profile is removed only when a unique negative main pole
supported on U1,...,U5 lies strictly below every first-integral resonance
correction.  Higher U_i orders are assigned their worst-case value zero.
This is a necessary local valuation screen, not a coefficient or branch
exclusion.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

import f2_75_125_211_joint_root_profile_screen as joint
import f2_75_125_k17_nonzero_first_slice as template_module


HERE = Path(__file__).resolve().parent
JOINT_RECEIPT = HERE / "f2_75_125_211_joint_root_profile_screen.json"
OUT = HERE / "f2_75_125_211_zero_u1_simple_root_k25_k35_screen.json"
INF = joint.INF
CHUNK = 128
START_BAND = 25
END_BAND = 35


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def bands_for(n: int):
    main = joint.pack_terms(n, Fraction(5, 3), 5, root_order=1)
    corrections = []
    for resonance in range(1, n // 5 + 1):
        total = n - 5 * resonance
        exponent_base = 5 - resonance
        corrections.append(
            joint.pack_terms(
                total,
                Fraction(exponent_base, 3),
                exponent_base,
                root_order=1,
            )
        )
    return main, corrections


def screen_band(profiles, main, corrections):
    survivors = []
    excluded = 0
    examples = []
    for start in range(0, len(profiles), CHUNK):
        batch_profiles = profiles[start : start + CHUNK]
        batch = np.asarray(batch_profiles, dtype=np.int16)
        main_min, main_ties, main_support = joint.least_orders(main, batch)
        correction_min = np.full(batch.shape[0], joint.INF_ORDER, dtype=np.int32)
        for correction in corrections:
            if correction[0].shape[0] == 0:
                continue
            correction_order, _ties, _support = joint.least_orders(correction, batch)
            correction_min = np.minimum(correction_min, correction_order)
        hit = (
            (main_min < 0)
            & (main_ties == 1)
            & main_support
            & (main_min < correction_min)
        )
        excluded += int(hit.sum())
        for profile, remove in zip(batch_profiles, hit.tolist()):
            if remove:
                if len(examples) < 12:
                    examples.append(list(profile))
            else:
                survivors.append(profile)
    return survivors, excluded, examples


def main() -> None:
    joint_receipt = json.loads(JOINT_RECEIPT.read_text(encoding="utf-8"))
    initial_profiles = joint.simple_profiles("zero_u1")
    if len(initial_profiles) != 22_176:
        raise AssertionError(f"unexpected initial U1=0 simple-root grid: {len(initial_profiles)}")
    profiles, k9_k24_counts = joint.screen_profiles(
        initial_profiles,
        joint.build_bands(root_order=1),
        "U1=0 simple-root exact replay",
    )
    expected = joint_receipt["simple_root_local_screens"]["zero_u1"]
    if len(profiles) != expected["survivors_after_K24"]:
        raise AssertionError("replayed K9--K24 survivors disagree with the hash-bound receipt")
    if len(profiles) != 18_876:
        raise AssertionError(f"unexpected K24 survivor count: {len(profiles)}")

    by_band = {}
    examples = {}
    survivors = profiles
    for n in range(START_BAND, END_BAND + 1):
        main, corrections = bands_for(n)
        survivors, excluded, sample = screen_band(survivors, main, corrections)
        by_band[str(n)] = {"newly_excluded": excluded, "survivors": len(survivors)}
        examples[str(n)] = sample
        print(f"U1=0 simple-root K{n}: excluded={excluded} survivors={len(survivors)}", flush=True)

    k9_k24_excluded = len(initial_profiles) - len(profiles)
    result = {
        "status": "finite local order-profile screen only; no whole-branch exclusion",
        "scope": "2+1+1 U1=0 simple-root profiles, exact first-integral valuation screen K25-K35",
        "field": "exact rational generalized-binomial coefficients and integer valuations",
        "root_order": 1,
        "higher_unknown_orders": 0,
        "only_unique_main_poles_supported_on_U1_to_U5_are_counted": True,
        "resonances_included": "all r=1..floor(K/5), including kappa_5 at K25 and later resonances",
        "initial_profile_count": len(initial_profiles),
        "K9_K24_excluded": k9_k24_excluded,
        "K9_K24_excluded_percent": round(100 * k9_k24_excluded / len(initial_profiles), 6),
        "K9_K24_by_band": k9_k24_counts,
        "profiles_entering_K25": len(profiles),
        "new_exclusions_by_band": by_band,
        "K9_K35_excluded": len(initial_profiles) - len(survivors),
        "K9_K35_excluded_percent": round(100 * (len(initial_profiles) - len(survivors)) / len(initial_profiles), 6),
        "profiles_surviving_K35": len(survivors),
        "survivor_profiles": [list(profile) for profile in survivors],
        "examples_by_band": examples,
        "input_sha256": {
            JOINT_RECEIPT.name: sha256(JOINT_RECEIPT),
            Path(joint.__file__).name: sha256(Path(joint.__file__)),
            Path(template_module.__file__).name: sha256(Path(template_module.__file__)),
        },
        "interpretation": (
            "Order profiles are local necessary data, not polynomial coefficient points. "
            "Even if this finite grid became empty, separate proof that the local charts "
            "cover the entire U1=0 algebraic locus would still be required."
        ),
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in result.items() if key != "survivor_profiles"}, indent=2), flush=True)
    print(f"receipt={OUT} sha256={sha256(OUT)}", flush=True)


if __name__ == "__main__":
    main()
