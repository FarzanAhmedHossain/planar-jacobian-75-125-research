#!/usr/bin/env python3
"""Exact K25--K29 local pole screen for the F2 2+1+1 branch.

This continues the finite order-profile audit from K24.  It uses the exact
first-integral expansion

    sum_r kappa_r w^(5r) F^((5-r)/3)

and retains the r=5 constant resonance at K25 as a possible correction.  The
r=6 resonance begins at K30 and is outside this screen.  A profile is pruned
only when the unique main negative-order term is supported on U1,...,U5 and
lies strictly below every retained correction.  The result is local profile
pruning, not an exclusion of the whole branch.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

from f2_75_125_k17_nonzero_first_slice import INFINITY, templates
from f2_75_125_k17_k24_nonzero_first_slice import packed_terms, term_orders


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_k17_k24_nonzero_first_slice.json"
OUT = HERE / "f2_75_125_k25_k29_nonzero_first_slice.json"
PROFILE_COUNT = 51_870


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    prior = json.loads(SOURCE.read_text(encoding="utf-8"))
    profiles = prior["survivor_profiles"]
    if len(profiles) != prior["survivor_profiles_after_K24"]:
        raise AssertionError("K24 survivor receipt count is inconsistent")
    if prior["K9_to_K24_union_profiles"] != 22_830:
        raise AssertionError("the upstream finite-grid union changed")
    if prior["profiles_entering_K17"] + prior["K9_to_K24_union_profiles"] != PROFILE_COUNT:
        raise AssertionError("finite profile census does not close")

    bands = {}
    for n in range(25, 30):
        main_terms = templates(n, Fraction(5, 3), tuple(range(1, n + 1)), 5)
        corrections = []
        for r in range(1, min(5, n // 5) + 1):
            total = n - 5 * r
            alpha = Fraction(5 - r, 3)
            terms = templates(total, alpha, tuple(range(1, total + 1)), 5 - r)
            corrections.append((r, packed_terms(terms)))
        bands[n] = (packed_terms(main_terms), corrections)

    survivors = []
    first_excluded: dict[str, int] = {str(n): 0 for n in range(25, 30)}
    examples: dict[str, list[dict[str, object]]] = {str(n): [] for n in range(25, 30)}
    for raw_profile in profiles:
        if len(raw_profile) != 5:
            raise AssertionError("malformed K24 survivor profile")
        profile = tuple(INFINITY if value == "infinity" else int(value)
                        for value in raw_profile)
        rejected = False
        for n in range(25, 30):
            main_order, main_ties, _finite, main_support, main_raw = term_orders(
                bands[n][0], profile
            )
            if main_order is None or main_order >= 0 or len(main_ties) != 1:
                continue
            main_idx = int(main_ties[0])
            if not main_support[main_idx]:
                continue

            correction_orders = []
            correction_detail = {}
            for r, packed in bands[n][1]:
                if packed[0].size == 0:
                    correction_detail[str(r)] = None
                    continue
                order, _ties, _mask, _support, _raw = term_orders(packed, profile)
                correction_detail[str(r)] = order
                if order is not None:
                    correction_orders.append(order)
            if not correction_orders or main_order >= min(correction_orders):
                continue

            first_excluded[str(n)] += 1
            if len(examples[str(n)]) < 12:
                multiplicities, coefficient, _base = main_raw[main_idx]
                examples[str(n)].append({
                    "orders_U1_to_U5": list(raw_profile),
                    "main_pole_order": -int(main_order),
                    "main_term_multiplicities": list(multiplicities),
                    "main_coefficient": str(coefficient),
                    "correction_orders_by_resonance": correction_detail,
                })
            rejected = True
            break
        if not rejected:
            survivors.append(raw_profile)

    added = sum(first_excluded.values())
    union = prior["K9_to_K24_union_profiles"] + added
    result = {
        "scope": "K25-K29 exact local valuation screen at the distinguished double root of the F2 2+1+1 nonzero-U1 branch",
        "field": "exact rational binomial coefficients and integer valuations",
        "input_receipt": str(SOURCE.relative_to(HERE.parent)),
        "input_receipt_sha256": sha256(SOURCE),
        "profiles_entering_K25": len(profiles),
        "finite_profile_grid_size": PROFILE_COUNT,
        "new_profiles_by_first_excluding_band": first_excluded,
        "K9_to_K29_union_profiles": union,
        "K9_to_K29_percent_of_finite_profile_grid": round(100 * union / PROFILE_COUNT, 6),
        "survivor_profiles_after_K29": len(survivors),
        "survivor_profiles": survivors,
        "resonance_terms_per_band": {
            str(n): [r for r, _ in bands[n][1]] for n in range(25, 30)
        },
        "examples": examples,
        "scope_note": "Kappa_5 is retained conservatively as a possible order-zero K25 correction; Kappa_6 first appears at K30. This screen is not a whole-branch exclusion.",
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    concise = {key: value for key, value in result.items() if key != "survivor_profiles"}
    print(json.dumps(concise, indent=2), flush=True)
    print(f"receipt={OUT} sha256={sha256(OUT)}", flush=True)


if __name__ == "__main__":
    main()
