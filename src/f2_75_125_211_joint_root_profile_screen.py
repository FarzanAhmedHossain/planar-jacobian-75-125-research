#!/usr/bin/env python3
"""Joint three-root valuation screen for the F2 2+1+1 carrier.

At h=(t-1)^2(t^2-3t+3), the two simple roots and the distinguished double
root must be roots of the same bounded-degree band polynomials U_i.  This
script screens K9--K24 independently at each simple root using exact
rational binomial coefficients and only counts a unique negative main pole
strictly below every resonance correction.  It then asks whether any
surviving local order profiles at t=1, alpha, beta can satisfy the degree
bounds deg(U1..U4)<=11 and deg(U5/t)<=12.

The local screens are deliberately over-approximations: higher unknown slices
are assigned order zero, and coefficient compatibility between distinct
roots is ignored.  Therefore an empty joint profile set would be a valid
local exclusion; a nonempty set is only a necessary-profile survivor.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

import f2_75_125_k17_nonzero_first_slice as template_module


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
NONZERO_DOUBLE = HERE / "f2_75_125_k17_k24_nonzero_first_slice.json"
ZERO_DOUBLE = HERE / "f2_75_125_k15_k24_zero_u1_double_root.json"
OUTPUT = HERE / "f2_75_125_211_joint_root_profile_screen.json"

INF = 99
INF_ORDER = 32767
DEGREE_CAPS = (11, 11, 11, 11, 12)
CHUNK = 512


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_order(value: object) -> int:
    return INF if value == "infinity" or int(value) >= 13 else int(value)


def pack_terms(total: int, alpha: Fraction, exponent_base: int, root_order: int):
    raw = template_module.templates(
        total, alpha, tuple(range(1, total + 1)), exponent_base
    )
    rows = []
    bases = []
    first_five_only = []
    for multiplicities, coefficient, _old_base in raw:
        if not coefficient:
            continue
        first = list(multiplicities[:5]) + [0] * max(0, 5 - len(multiplicities))
        rows.append(first)
        bases.append(root_order * (exponent_base - 3 * sum(multiplicities)))
        first_five_only.append(not any(multiplicities[5:]))
    return (
        np.asarray(rows, dtype=np.int16),
        np.asarray(bases, dtype=np.int16),
        np.asarray(first_five_only, dtype=np.bool_),
    )


def least_orders(packed, profiles: np.ndarray):
    exponents, bases, _support = packed
    finite_values = np.where(profiles < INF, profiles, 0).astype(np.int32)
    values = bases.astype(np.int32)[:, None] + exponents.astype(np.int32) @ finite_values.T
    infinite_mask = (profiles >= INF).astype(np.int32)
    absent = (exponents.astype(np.int32) @ infinite_mask.T) > 0
    values[absent] = INF_ORDER
    minima = values.min(axis=0)
    tie_counts = (values == minima[None, :]).sum(axis=0)
    argmins = values.argmin(axis=0)
    unique_support = packed[2][argmins]
    return minima, tie_counts, unique_support


def build_bands(root_order: int):
    bands = {}
    for n in range(9, 25):
        main = pack_terms(n, Fraction(5, 3), 5, root_order)
        corrections = []
        for resonance in range(1, min(4, n // 5) + 1):
            total = n - 5 * resonance
            exponent_base = 5 - resonance
            corrections.append(
                pack_terms(
                    total,
                    Fraction(exponent_base, 3),
                    exponent_base,
                    root_order,
                )
            )
        bands[n] = (main, corrections)
    return bands


def profile_batches(profiles: list[tuple[int, ...]]):
    for start in range(0, len(profiles), CHUNK):
        yield start, np.asarray(profiles[start : start + CHUNK], dtype=np.int16)


def screen_profiles(profiles: list[tuple[int, ...]], bands, label: str):
    survivors = profiles
    by_band = {}
    for n in range(9, 25):
        main, corrections = bands[n]
        kept = []
        excluded = 0
        for start, batch in profile_batches(survivors):
            main_min, main_ties, main_supported = least_orders(main, batch)
            correction_min = np.full(batch.shape[0], INF_ORDER, dtype=np.int32)
            for correction in corrections:
                corr_min, _corr_ties, _corr_support = least_orders(correction, batch)
                correction_min = np.minimum(correction_min, corr_min)
            hit = (
                (main_min < 0)
                & (main_ties == 1)
                & main_supported
                & (main_min < correction_min)
            )
            excluded += int(hit.sum())
            if not np.any(hit):
                kept.extend(survivors[start : start + batch.shape[0]])
            else:
                kept.extend(
                    profile
                    for profile, remove in zip(
                        survivors[start : start + batch.shape[0]], hit.tolist()
                    )
                    if not remove
                )
        survivors = kept
        by_band[str(n)] = {"newly_excluded": excluded, "survivors": len(survivors)}
        print(
            f"{label} K{n}: excluded={excluded} survivors={len(survivors)}",
            flush=True,
        )
    return survivors, by_band


def simple_profiles(category: str) -> list[tuple[int, ...]]:
    values = [*range(0, 12), INF]
    u5_values = [*range(0, 13), INF]
    profiles = []
    if category == "zero_u1":
        # U1=0: K6/K7/K8 give ord(U2)>=2, ord(U3)>=1,
        # ord(U4)>=1 at a simple root.
        for u2 in [*range(2, 12), INF]:
            for u3 in [*range(1, 12), INF]:
                for u4 in [*range(1, 12), INF]:
                    for u5 in u5_values:
                        profiles.append((INF, u2, u3, u4, u5))
    elif category == "tuned":
        # The tuned simple-root alternative has (ord U1, ord U2, ord U3)=(2,1,0).
        for u4 in values:
            for u5 in u5_values:
                profiles.append((2, 1, 0, u4, u5))
    elif category == "deeper":
        # The deeper alternative has ord(U1)>=3, ord(U2)>=2, U3=0;
        # K8 also forces U4 to vanish at this simple root.
        for u1 in (3, 4):
            for u2 in [*range(2, 12), INF]:
                for u3 in [*range(1, 12), INF]:
                    for u4 in [*range(1, 12), INF]:
                        for u5 in u5_values:
                            profiles.append((u1, u2, u3, u4, u5))
    else:
        raise ValueError(category)
    return profiles


def decode_double_receipts():
    nonzero_data = json.loads(NONZERO_DOUBLE.read_text(encoding="utf-8"))
    zero_data = json.loads(ZERO_DOUBLE.read_text(encoding="utf-8"))
    nonzero = [
        tuple(parse_order(x) for x in row)
        for row in nonzero_data["survivor_profiles"]
    ]
    zero = [
        (INF, *(parse_order(x) for x in row))
        for row in zero_data["survivor_profiles"]
    ]
    return nonzero_data, zero_data, nonzero, zero


def groups(profiles: list[tuple[int, ...]]):
    result = {}
    by_mask = {}
    for profile in profiles:
        mask = sum(1 << i for i, value in enumerate(profile) if value >= INF)
        by_mask.setdefault(mask, set()).add(profile)
    for mask, full_points in by_mask.items():
        active = tuple(i for i in range(5) if not (mask & (1 << i)))
        projected = sorted({tuple(point[i] for i in active) for point in full_points})
        if active:
            shape = tuple(DEGREE_CAPS[i] + 1 for i in active)
            grid = np.zeros(shape, dtype=np.uint32)
            for point in projected:
                grid[point] = 1
            prefix = grid
            for axis in range(len(active)):
                prefix = np.cumsum(prefix, axis=axis, dtype=np.uint32)
            frontier = [point for point in projected if int(prefix[point]) == 1]
        else:
            prefix = np.asarray(1, dtype=np.uint32)
            frontier = [()]
        result[mask] = {
            "active": active,
            "points": projected,
            "frontier": frontier,
            "prefix": prefix,
        }
    return result


def full_profile(mask: int, point: tuple[int, ...]) -> tuple[int, ...]:
    active_values = iter(point)
    return tuple(
        "infinity" if mask & (1 << i) else next(active_values)
        for i in range(5)
    )


def query_pair(double_profiles, simple_a_profiles, simple_b_profiles, label: str):
    dg = groups(double_profiles)
    ag = groups(simple_a_profiles)
    bg = groups(simple_b_profiles)
    matching_masks = sorted(set(dg) & set(ag) & set(bg))
    compatible_count = 0
    first_witness = None
    compatible_by_mask = {}
    for mask in matching_masks:
        dgroup, agroup, bgroup = dg[mask], ag[mask], bg[mask]
        active = dgroup["active"]
        caps = tuple(DEGREE_CAPS[i] for i in active)
        compatible_here = 0
        for dpoint in dgroup["points"]:
            found = False
            for apoint in agroup["frontier"]:
                limits = tuple(c - d - a for c, d, a in zip(caps, dpoint, apoint))
                if any(limit < 0 for limit in limits):
                    continue
                if not active:
                    bpoint = ()
                else:
                    if any(limit > DEGREE_CAPS[i] for i, limit in zip(active, limits)):
                        limits = tuple(
                            min(limit, DEGREE_CAPS[i])
                            for i, limit in zip(active, limits)
                        )
                    if int(bgroup["prefix"][limits]) == 0:
                        continue
                    bpoint = next(
                        point
                        for point in bgroup["points"]
                        if all(value <= limit for value, limit in zip(point, limits))
                    )
                compatible_here += 1
                compatible_count += 1
                found = True
                if first_witness is None:
                    first_witness = {
                        "zero_mask_variables": [f"U{i+1}" for i in range(5) if mask & (1 << i)],
                        "double_root_orders": list(full_profile(mask, dpoint)),
                        "simple_root_alpha_orders": list(full_profile(mask, apoint)),
                        "simple_root_beta_orders": list(full_profile(mask, bpoint)),
                        "degree_sums": [
                            "zero polynomial" if mask & (1 << i)
                            else dpoint[j] + apoint[j] + bpoint[j]
                            for j, i in enumerate(active)
                        ],
                        "active_variables": [f"U{i+1}" for i in active],
                    }
                break
            if not found:
                continue
        compatible_by_mask[str(mask)] = {
            "double_root_profiles": len(dgroup["points"]),
            "degree_compatible_profiles": compatible_here,
        }
    return {
        "label": label,
        "compatible": compatible_count > 0,
        "matching_zero_masks": len(matching_masks),
        "double_root_profiles": len(double_profiles),
        "degree_compatible_double_root_profiles": compatible_count,
        "incompatible_double_root_profiles": len(double_profiles) - compatible_count,
        "excluded_percent_of_double_root_profile_grid": round(
            100 * (len(double_profiles) - compatible_count) / len(double_profiles), 6
        ) if double_profiles else 0.0,
        "compatible_profiles_by_zero_mask": compatible_by_mask,
        "witness": first_witness,
    }


def main() -> None:
    nonzero_data, zero_data, double_nonzero, double_zero = decode_double_receipts()
    simple_bands = build_bands(root_order=1)
    simple_local = {}
    simple_receipts = {}
    for label in ("tuned", "deeper", "zero_u1"):
        initial = simple_profiles(label)
        print(f"simple root {label}: initial_profiles={len(initial)}", flush=True)
        surviving, by_band = screen_profiles(initial, simple_bands, label)
        simple_local[label] = surviving
        simple_receipts[label] = {
            "initial_profiles": len(initial),
            "survivors_after_K24": len(surviving),
            "excluded_percent_of_simple_root_profile_grid": round(
                100 * (len(initial) - len(surviving)) / len(initial), 6
            ),
            "new_exclusions_by_band": by_band,
            "survivor_examples": [list(row) for row in surviving[:20]],
        }

    joint = [
        query_pair(double_nonzero, simple_local["tuned"], simple_local["tuned"], "nonzero-U1 / tuned-tuned"),
        query_pair(double_nonzero, simple_local["tuned"], simple_local["deeper"], "nonzero-U1 / tuned-deeper"),
        query_pair(double_nonzero, simple_local["deeper"], simple_local["deeper"], "nonzero-U1 / deeper-deeper"),
        query_pair(double_zero, simple_local["zero_u1"], simple_local["zero_u1"], "U1=0"),
    ]
    result = {
        "scope": "K9-K24 exact local valuation screen at both simple roots, joined to prior double-root survivors by global degree sums",
        "field": "exact Q generalized-binomial coefficients; integer root orders",
        "script_sha256": sha256(Path(__file__)),
        "carrier": "h=(t-1)^2(t^2-3t+3)",
        "degree_caps_U1_to_U5_over_t": list(DEGREE_CAPS),
        "chunk_size": CHUNK,
        "input_receipts": {
            NONZERO_DOUBLE.name: sha256(NONZERO_DOUBLE),
            ZERO_DOUBLE.name: sha256(ZERO_DOUBLE),
        },
        "double_root_survivor_counts": {
            "nonzero_U1_after_K24": len(double_nonzero),
            "U1_zero_after_K24": len(double_zero),
        },
        "simple_root_local_screens": simple_receipts,
        "joint_degree_compatibility": joint,
        "interpretation": (
            "Every actual coefficient solution must lie among the returned local-profile survivors and satisfy these degree sums. "
            "An incompatible category is excluded by necessary local conditions. A compatible order witness only shows that this valuation screen does not exclude it; it is not a polynomial or Keller-pair witness."
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "simple_root_local_screens": {
            k: {
                "initial_profiles": v["initial_profiles"],
                "survivors_after_K24": v["survivors_after_K24"],
                "excluded_percent_of_simple_root_profile_grid": v["excluded_percent_of_simple_root_profile_grid"],
                "new_exclusions_by_band": v["new_exclusions_by_band"],
            }
            for k, v in simple_receipts.items()
        },
        "joint_degree_compatibility": joint,
        "receipt": str(OUTPUT.relative_to(ROOT)),
        "receipt_sha256": sha256(OUTPUT),
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
