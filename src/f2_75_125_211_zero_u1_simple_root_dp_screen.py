#!/usr/bin/env python3
"""Extend the U1=0 simple-root valuation screen using exact min-plus DP.

The earlier K25--K40 screen enumerates integer partitions of each band
index.  For a profile (ord(U1),...,ord(U5)), a partition with multiplicities
m_i has valuation

    exponent_base + sum_i m_i * (ord(U_i) - 3),

at a simple root.  Orders above U5 are conservatively set to zero.  The
coefficient is nonzero unless alpha is a nonnegative integer smaller than
the number of parts.  An unbounded-coin dynamic program therefore computes
the same minimum, tie count (capped at two), and support flag without listing
all partitions.  This only screens a finite order-profile grid; it is not a
coefficient-system or branch-wide exclusion.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

import f2_75_125_211_joint_root_profile_screen as joint
import f2_75_125_k17_nonzero_first_slice as template_module


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_zero_u1_simple_root_k36_k40_screen.json"
OUT = HERE / "f2_75_125_211_zero_u1_simple_root_k41_k125_dp_screen.json"
BIG = 1 << 28


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def optimal_tables(profiles: list[tuple[int, ...]], through: int):
    """Return min cost, capped tie count, and all-minimizers-supported tables."""
    profile_array = np.asarray(profiles, dtype=np.int16)
    if profile_array.ndim != 2 or profile_array.shape[1] != 5:
        raise AssertionError("expected five-order profiles")
    count = profile_array.shape[0]
    best = np.full((through + 1, count), BIG, dtype=np.int32)
    ways = np.zeros((through + 1, count), dtype=np.uint8)
    supported = np.zeros((through + 1, count), dtype=np.bool_)
    best[0] = 0
    ways[0] = 1
    supported[0] = True

    # Standard unbounded-coin recurrence in increasing part and total order;
    # this counts each integer partition once, independent of part ordering.
    for part in range(1, through + 1):
        if part <= 5:
            value = profile_array[:, part - 1].astype(np.int32)
            available = value < joint.INF
            part_cost = value - 3
            part_supported = True
        else:
            available = np.ones(count, dtype=np.bool_)
            part_cost = -3  # higher unknown orders are assigned their worst case, 0
            part_supported = False

        for total in range(part, through + 1):
            prior_ways = ways[total - part]
            reachable = available & (prior_ways != 0)
            candidate = best[total - part] + part_cost
            better = reachable & (candidate < best[total])
            tied = reachable & (candidate == best[total])

            if better.any():
                best[total, better] = candidate[better]
                ways[total, better] = prior_ways[better]
                supported[total, better] = (
                    supported[total - part, better] & part_supported
                )
            if tied.any():
                ways[total, tied] = np.minimum(
                    2, ways[total, tied] + prior_ways[tied]
                )
                supported[total, tied] &= (
                    supported[total - part, tied] & part_supported
                )
    return profile_array, best, ways, supported


def correction_minimum(
    band: int,
    profile_array: np.ndarray,
    best: np.ndarray,
) -> np.ndarray:
    """Minimum valuation among all nonzero resonance-correction terms."""
    count = profile_array.shape[0]
    correction_min = np.full(count, joint.INF_ORDER, dtype=np.int32)
    for resonance in range(1, band // 5 + 1):
        total = band - 5 * resonance
        exponent_base = 5 - resonance
        alpha = Fraction(exponent_base, 3)

        # For alpha a nonnegative integer, generalized-binomial terms vanish
        # when the part count exceeds alpha.  In this family the only such
        # cases are alpha=1 (resonance 2) and alpha=0 (resonance 5).
        if alpha.denominator == 1 and alpha.numerator >= 0:
            integer_alpha = alpha.numerator
            if total == 0:
                value = np.full(count, exponent_base, dtype=np.int32)
            elif integer_alpha == 0:
                continue
            elif integer_alpha == 1:
                if total <= 5:
                    orders = profile_array[:, total - 1].astype(np.int32)
                    present = orders < joint.INF
                    value = np.where(
                        present, exponent_base + orders - 3, joint.INF_ORDER
                    ).astype(np.int32)
                else:
                    value = np.full(count, exponent_base - 3, dtype=np.int32)
            else:
                raise AssertionError("unexpected nonnegative integral alpha")
        else:
            value = exponent_base + best[total]
        correction_min = np.minimum(correction_min, value)
    return correction_min


def exact_enumeration_audit(profiles: list[tuple[int, ...]]) -> dict:
    """Compare the DP to the original exact partition enumerator on controls."""
    if not profiles:
        raise AssertionError("missing profiles for enumeration audit")
    indices = sorted({0, min(37, len(profiles) - 1), min(777, len(profiles) - 1), len(profiles) - 1})
    sample = [profiles[index] for index in indices]
    arr, best, ways, supported = optimal_tables(sample, 40)
    bands_checked = [9, 10, 12, 25, 40]
    for band in bands_checked:
        main = joint.pack_terms(band, Fraction(5, 3), 5, root_order=1)
        brute_min, brute_ties, brute_support = joint.least_orders(
            main, np.asarray(sample, dtype=np.int16)
        )
        dp_min = 5 + best[band]
        if not np.array_equal(dp_min, brute_min):
            raise AssertionError(f"main minimum mismatch at K{band}")
        if not np.array_equal(np.minimum(2, brute_ties), ways[band]):
            raise AssertionError(f"main tie-count mismatch at K{band}")
        unique = brute_ties == 1
        if not np.array_equal(supported[band, unique], brute_support[unique]):
            raise AssertionError(f"main support mismatch at K{band}")

        brute_correction = np.full(len(sample), joint.INF_ORDER, dtype=np.int32)
        for resonance in range(1, band // 5 + 1):
            total = band - 5 * resonance
            exponent_base = 5 - resonance
            packed = joint.pack_terms(
                total, Fraction(exponent_base, 3), exponent_base, root_order=1
            )
            if packed[0].shape[0]:
                values, _ties, _support = joint.least_orders(
                    packed, np.asarray(sample, dtype=np.int16)
                )
                brute_correction = np.minimum(brute_correction, values)
        dp_correction = correction_minimum(band, arr, best)
        if not np.array_equal(dp_correction, brute_correction):
            raise AssertionError(f"correction minimum mismatch at K{band}")
    return {"profile_indices": indices, "bands": bands_checked, "exact_match": True}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--through", type=int, default=125)
    parser.add_argument("--no-enumeration-audit", action="store_true")
    args = parser.parse_args()
    if args.through < 41:
        parser.error("--through must be at least 41")

    prior = json.loads(SOURCE.read_text(encoding="utf-8"))
    profiles = [tuple(map(int, row)) for row in prior["survivor_profiles"]]
    if (
        len(profiles) != prior["profiles_surviving_K40"]
        or len(profiles) != 18_876
        or prior["initial_profile_count"] != 22_176
        or prior["K9_K40_excluded"] != 3_300
    ):
        raise AssertionError("K40 input receipt/counts are inconsistent")

    audit = None
    if not args.no_enumeration_audit:
        audit = exact_enumeration_audit(profiles)

    profile_array, best, ways, supported = optimal_tables(profiles, args.through)
    alive = np.ones(len(profiles), dtype=np.bool_)
    by_band = {}
    examples = {}
    for band in range(41, args.through + 1):
        main_order = 5 + best[band]
        corrections = correction_minimum(band, profile_array, best)
        hit = (
            (main_order < 0)
            & (ways[band] == 1)
            & supported[band]
            & (main_order < corrections)
        )
        newly_excluded = int(np.count_nonzero(hit & alive))
        examples[str(band)] = [
            list(profiles[i])
            for i in np.flatnonzero(hit & alive)[:12]
        ]
        alive &= ~hit
        by_band[str(band)] = {
            "newly_excluded": newly_excluded,
            "survivors": int(np.count_nonzero(alive)),
        }
        if newly_excluded:
            print(
                f"U1=0 simple-root K{band}: excluded={newly_excluded} "
                f"survivors={int(np.count_nonzero(alive))}",
                flush=True,
            )

    surviving = [list(profile) for profile, keep in zip(profiles, alive) if keep]
    total_excluded = 22_176 - len(surviving)
    result = {
        "status": "finite local order-profile screen only; no whole-branch exclusion",
        "scope": f"2+1+1 U1=0 simple-root valuation profiles, exact DP screen K41-K{args.through}",
        "field": "integer min-plus recurrence; exact generalized-binomial zero criterion",
        "input_receipt": SOURCE.name,
        "input_sha256": sha256(SOURCE),
        "support_rule": "orders above U5 assigned 0; only a unique negative main pole supported on U1..U5 below every correction is pruned",
        "tie_count_cap": 2,
        "initial_profile_count": 22_176,
        "profiles_entering_K41": len(profiles),
        "K9_K40_excluded": 3_300,
        "new_exclusions_by_band": by_band,
        "K9_through_excluded": total_excluded,
        "K9_through_excluded_percent": round(100 * total_excluded / 22_176, 6),
        "profiles_surviving_through": len(surviving),
        "survivor_profiles": surviving,
        "examples_by_band": examples,
        "enumeration_audit": audit,
        "source_sha256": {
            Path(__file__).name: sha256(Path(__file__)),
            Path(joint.__file__).name: sha256(Path(joint.__file__)),
            Path(template_module.__file__).name: sha256(Path(template_module.__file__)),
        },
        "interpretation": (
            "Profiles are valuation patterns, not polynomial coefficient solutions. "
            "This screen is necessary and finite; it does not exclude the complete U1=0 stratum."
        ),
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in result.items() if key != "survivor_profiles"}, indent=2), flush=True)
    print(f"receipt={OUT} sha256={sha256(OUT)}", flush=True)


if __name__ == "__main__":
    main()
