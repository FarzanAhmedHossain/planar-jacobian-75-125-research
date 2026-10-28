#!/usr/bin/env python3
"""Extend the U1=0 deeper-double-root valuation screen with exact min-plus DP.

The K15--K24 screen enumerates partitions.  For a double root, each part of
size j contributes ord(U_j)-6 to the root valuation; orders above U5 are
conservatively set to zero.  The additive valuation permits an unbounded-coin
dynamic program through K125 without enumerating integer partitions.  This
screens only the finite local order-profile grid, not the coefficient ideal or
the whole 2+1+1 branch.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

import f2_75_125_k15_k24_zero_u1_double_root as prior
import f2_75_125_k17_nonzero_first_slice as template_module


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_k15_k24_zero_u1_double_root.json"
OUT = HERE / "f2_75_125_211_zero_u1_double_root_dp_k25_k125.json"
INFINITY = 13
INFINITE_ORDER = 1 << 28


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def optimal_tables(profiles: list[tuple[int, ...]], through: int):
    """Min valuation, capped minimizer count, and all-minimizers support flag."""
    profile_array = np.asarray(profiles, dtype=np.int16)
    if profile_array.ndim != 2 or profile_array.shape[1] != 4:
        raise AssertionError("expected (ord U2,...,ord U5) profiles")
    count = profile_array.shape[0]
    best = np.full((through + 1, count), INFINITE_ORDER, dtype=np.int32)
    ways = np.zeros((through + 1, count), dtype=np.uint8)
    supported = np.zeros((through + 1, count), dtype=np.bool_)
    best[0] = 0
    ways[0] = 1
    supported[0] = True

    # Unbounded-coin recurrence counts each integer partition once.
    for part in range(2, through + 1):
        if part <= 5:
            order = profile_array[:, part - 2].astype(np.int32)
            available = order < INFINITY
            part_cost = order - 6
            part_supported = True
        else:
            available = np.ones(count, dtype=np.bool_)
            part_cost = -6
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
    """Minimum valuation over all potentially nonzero resonance terms."""
    count = profile_array.shape[0]
    correction_min = np.full(count, INFINITE_ORDER, dtype=np.int32)
    for resonance in range(1, band // 5 + 1):
        total = band - 5 * resonance
        exponent_base = 5 - resonance
        alpha = Fraction(exponent_base, 3)

        # The only nonnegative integral alpha values are 1 (r=2) and 0 (r=5).
        # Their generalized-binomial series has only the corresponding finite
        # number of parts, so handle these exactly rather than using the DP.
        if alpha.denominator == 1 and alpha.numerator >= 0:
            integer_alpha = alpha.numerator
            if total == 0:
                value = np.full(
                    count, 2 * exponent_base, dtype=np.int32
                )
            elif integer_alpha == 0:
                continue
            elif integer_alpha == 1:
                if total == 1:
                    continue  # no part of size one because U1 is identically zero
                base = 2 * exponent_base - 6
                if total <= 5:
                    order = profile_array[:, total - 2].astype(np.int32)
                    present = order < INFINITY
                    value = np.where(
                        present, base + order, INFINITE_ORDER
                    ).astype(np.int32)
                else:
                    # A higher U_j is unknown; order zero is the conservative
                    # lower bound, and this term is unsupported by U2..U5.
                    value = np.full(count, base, dtype=np.int32)
            else:
                raise AssertionError("unexpected nonnegative integral alpha")
        else:
            value = 2 * exponent_base + best[total]
        correction_min = np.minimum(correction_min, value)
    return correction_min


def exact_enumeration_audit(profiles: list[tuple[int, ...]]) -> dict:
    """Compare DP minima, ties, support, and corrections against K15--K24 code."""
    if not profiles:
        raise AssertionError("missing survivor profiles for audit")
    indices = sorted({
        0,
        min(37, len(profiles) - 1),
        min(777, len(profiles) - 1),
        len(profiles) - 1,
    })
    sample = [profiles[index] for index in indices]
    arr, best, ways, supported = optimal_tables(sample, 24)
    bands_checked = [15, 16, 20, 24]
    for band in bands_checked:
        weights = tuple(range(2, band + 1))
        main_terms = template_module.templates(
            band, Fraction(5, 3), weights, 5
        )
        main = prior.pack_terms(main_terms, weights)
        for i, profile in enumerate(sample):
            brute_min, brute_ties, brute_support, _ = prior.least_order(
                main, np.asarray(profile, dtype=np.int16)
            )
            dp_min = 10 + int(best[band, i])
            if dp_min != brute_min:
                raise AssertionError(f"main minimum mismatch at K{band}")
            if min(2, len(brute_ties)) != int(ways[band, i]):
                raise AssertionError(f"main tie-count mismatch at K{band}")
            if len(brute_ties) == 1 and bool(supported[band, i]) != bool(
                brute_support[int(brute_ties[0])]
            ):
                raise AssertionError(f"main support mismatch at K{band}")

        brute_correction = np.full(len(sample), INFINITE_ORDER, dtype=np.int32)
        for resonance in range(1, min(4, band // 5) + 1):
            total = band - 5 * resonance
            weights = tuple(range(2, total + 1))
            terms = template_module.templates(
                total,
                Fraction(5 - resonance, 3),
                weights,
                5 - resonance,
            )
            if not terms:
                continue  # no partition when the residual weight is 1
            packed = prior.pack_terms(terms, weights)
            if len(packed[0]) == 0:
                continue
            for i, profile in enumerate(sample):
                order, _ties, _support, _raw = prior.least_order(
                    packed, np.asarray(profile, dtype=np.int16)
                )
                if order is not None:
                    brute_correction[i] = min(brute_correction[i], order)
        dp_correction = correction_minimum(band, arr, best)
        if not np.array_equal(dp_correction, brute_correction):
            raise AssertionError(f"correction minimum mismatch at K{band}")
    return {"profile_indices": indices, "bands": bands_checked, "exact_match": True}


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--through", type=int, default=125)
    args = parser.parse_args()
    if args.through < 25:
        parser.error("--through must be at least 25")

    prior_receipt = json.loads(SOURCE.read_text(encoding="utf-8"))
    profiles = [tuple(map(int, row)) for row in prior_receipt["survivor_profiles"]]
    if (
        prior_receipt["profiles_after_prior_filters"] != 12_848
        or prior_receipt["survivor_profile_count_after_K24"] != 12_128
        or len(profiles) != 12_128
        or prior_receipt["new_profiles_by_first_excluding_band"].get("15") != 720
        or sum(prior_receipt["new_profiles_by_first_excluding_band"].values()) != 720
    ):
        raise AssertionError("K15--K24 input counts are inconsistent")

    audit = exact_enumeration_audit(profiles)
    profile_array, best, ways, supported = optimal_tables(profiles, args.through)
    alive = np.ones(len(profiles), dtype=np.bool_)
    by_band: dict[str, dict[str, int]] = {}
    examples: dict[str, list[list[int]]] = {}
    for band in range(25, args.through + 1):
        main_order = 10 + best[band]
        corrections = correction_minimum(band, profile_array, best)
        hit = (
            (main_order < 0)
            & (ways[band] == 1)
            & supported[band]
            & (main_order < corrections)
        )
        newly_excluded = int(np.count_nonzero(hit & alive))
        examples[str(band)] = [
            list(profiles[i]) for i in np.flatnonzero(hit & alive)[:12]
        ]
        alive &= ~hit
        by_band[str(band)] = {
            "newly_excluded": newly_excluded,
            "survivors": int(np.count_nonzero(alive)),
        }
    surviving = [
        list(profile) for profile, keep in zip(profiles, alive) if keep
    ]
    new_excluded = len(profiles) - len(surviving)
    total_excluded = 720 + new_excluded
    raw_lattice = 12_870
    result = {
        "status": "finite local order-profile screen only; no whole-branch exclusion",
        "scope": f"2+1+1 U1=0 double-root profiles, exact DP screen K25-K{args.through}",
        "field": "integer min-plus recurrence; exact generalized-binomial zero criterion",
        "input_receipt": SOURCE.name,
        "input_sha256": sha256(SOURCE),
        "support_rule": "orders above U5 assigned 0; only a unique negative main pole supported on U2..U5 below every correction is pruned",
        "profile_counts": {
            "raw_bounded_lattice": raw_lattice,
            "prior_exact_residue_exclusions": prior_receipt["profiles_prefiltered_by_prior_exact_residues"],
            "profiles_after_prior_exact_filters": 12_848,
            "K15_K24_screen_exclusions": 720,
            "profiles_entering_K25": len(profiles),
            "new_K25_through_screen_exclusions": new_excluded,
            "profiles_surviving_through": len(surviving),
            "total_profile_exclusions_after_K15_through": total_excluded,
            "profile_grid_percent_after_K15_through": round(
                100 * total_excluded / 12_848, 6
            ),
            "raw_bounded_lattice_percent_including_prior_exact_strata": round(
                100 * (22 + total_excluded) / raw_lattice, 6
            ),
        },
        "new_exclusions_by_band": by_band,
        "examples_by_band": examples,
        "survivor_profiles": surviving,
        "enumeration_audit": audit,
        "source_sha256": {
            Path(__file__).name: sha256(Path(__file__)),
            Path(prior.__file__).name: sha256(Path(prior.__file__)),
            Path(template_module.__file__).name: sha256(Path(template_module.__file__)),
        },
        "interpretation": (
            "These are valuation profiles, not polynomial coefficient solutions. "
            "The conservative screen does not exclude the full U1=0 or 2+1+1 branch."
        ),
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    concise = {key: value for key, value in result.items() if key != "survivor_profiles"}
    print(json.dumps(concise, indent=2), flush=True)
    print(f"receipt={OUT} sha256={sha256(OUT)}", flush=True)


if __name__ == "__main__":
    main()
