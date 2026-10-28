"""Exact K26--K35 extension of the F2 2+1+1 double-root profile screen.

The source first-integral identity is used through w^35, with resonance
corrections r=1,...,7. In particular the negative-exponent r=6 and r=7
corrections enter at K30 and K35. Higher slices are assigned order zero;
profiles are removed only by a unique negative main pole strictly below all
possible correction orders. This remains finite local-profile pruning.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

from f2_75_125_k17_k24_nonzero_first_slice import (
    INFINITY,
    packed_terms,
    term_orders,
    templates,
)


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_k17_k24_nonzero_first_slice.json"
OUT = HERE / "f2_75_125_211_k26_k35_nonzero_first_slice_20261021.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    prior = json.loads(SOURCE.read_text(encoding="utf-8"))
    profiles = prior["survivor_profiles"]
    if len(profiles) != prior["survivor_profiles_after_K24"]:
        raise AssertionError("K24 survivor receipt count is inconsistent")

    bands = {}
    for n in range(26, 36):
        main_terms = templates(n, Fraction(5, 3), tuple(range(1, n + 1)), 5)
        corrections = []
        for r in range(1, min(7, n // 5) + 1):
            total = n - 5 * r
            alpha = Fraction(5 - r, 3)
            terms = templates(total, alpha, tuple(range(1, total + 1)), 5 - r)
            corrections.append((r, packed_terms(terms)))
        bands[n] = (packed_terms(main_terms), corrections)

    excluded_by_band = {n: 0 for n in range(26, 36)}
    excluded_by_u2 = {n: {} for n in range(26, 36)}
    examples = {n: [] for n in range(26, 36)}
    survivors = []
    for raw_profile in profiles:
        if len(raw_profile) != 5:
            raise AssertionError("malformed K24 survivor profile")
        profile = tuple(
            INFINITY if value == "infinity" else int(value)
            for value in raw_profile
        )
        removed = False
        for n in range(26, 36):
            main, corrections = bands[n]
            main_order, main_ties, _finite, main_support, main_raw = term_orders(main, profile)
            if main_order is None or main_order >= 0 or len(main_ties) != 1:
                continue
            main_idx = int(main_ties[0])
            if not main_support[main_idx]:
                continue

            correction_orders = {}
            finite_orders = []
            for r, correction in corrections:
                order, _ties, _mask, _support, _terms = term_orders(correction, profile)
                correction_orders[str(r)] = order
                if order is not None:
                    finite_orders.append(order)
            if not finite_orders or main_order >= min(finite_orders):
                continue

            multiplicities, coefficient, _base = main_raw[main_idx]
            excluded_by_band[n] += 1
            u2 = raw_profile[1]
            key = "infinity" if u2 == "infinity" else str(u2)
            excluded_by_u2[n][key] = excluded_by_u2[n].get(key, 0) + 1
            if len(examples[n]) < 8:
                examples[n].append({
                    "orders_U1_to_U5": raw_profile,
                    "main_pole_order": -int(main_order),
                    "main_term_multiplicities": list(multiplicities),
                    "main_coefficient": str(coefficient),
                    "correction_orders_by_r": correction_orders,
                })
            removed = True
            break
        if not removed:
            survivors.append(raw_profile)

    result = {
        "scope": "K26-K35 local valuation screen at the distinguished double root of the F2 2+1+1 nonzero-U1 deeper branch",
        "field": "exact Q generalized-binomial coefficients and integer order arithmetic",
        "input_receipt": SOURCE.name,
        "input_receipt_sha256": digest(SOURCE),
        "profiles_entering_K26": len(profiles),
        "first_excluding_band_counts": {str(n): count for n, count in excluded_by_band.items()},
        "first_excluding_band_by_ord_U2": {str(n): counts for n, counts in excluded_by_u2.items()},
        "nonzero_main_and_correction_term_counts": {
            str(n): {
                "main_F_5_3": len(bands[n][0][3]),
                "corrections_by_r": {
                    str(r): len(packed[3]) for r, packed in bands[n][1]
                },
            }
            for n in range(26, 36)
        },
        "survivor_profiles_after_K35": len(survivors),
        "survivor_profiles": survivors,
        "examples": {str(n): rows for n, rows in examples.items()},
        "status": "finite local-profile pruning only; whole 2+1+1 branch exclusion remains unproved",
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    concise = {key: value for key, value in result.items() if key != "survivor_profiles"}
    print(json.dumps(concise, indent=2), flush=True)
    print(f"receipt={OUT} sha256={digest(OUT)}", flush=True)


if __name__ == "__main__":
    main()
