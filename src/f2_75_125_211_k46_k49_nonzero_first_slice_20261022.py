"""Exact continuation of the finite F2 2+1+1 valuation-profile screen.

This checks K46--K49 for every profile surviving the saved K45 receipt.
The screen is exact over Q and remains finite-profile pruning, not a
whole-branch exclusion.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

from f2_75_125_k17_nonzero_first_slice import INFINITY, templates
from f2_75_125_k17_k24_nonzero_first_slice import packed_terms, term_orders


HERE = Path(__file__).resolve().parent
K45 = HERE / "f2_75_125_211_k36_k45_nonzero_first_slice_20261022.json"
OUT = HERE / "f2_75_125_211_k46_k49_nonzero_first_slice_20261022.json"
FIRST_BAND = 46
LAST_BAND = 49
GRID_COUNT = 51_870


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def save(path: Path, value: dict) -> None:
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    temp.replace(path)


def main() -> None:
    k45 = json.loads(K45.read_text(encoding="utf-8"))
    if k45.get("last_complete_band") != 45:
        raise AssertionError("input receipt is not complete through K45")
    profiles = k45["survivor_profiles"]
    if len(profiles) != k45.get("survivor_profiles_after_last_complete_band"):
        raise AssertionError("K45 survivor count does not match its list")
    if k45.get("finite_profile_grid_count") != GRID_COUNT:
        raise AssertionError("unexpected finite profile grid")

    counts = {n: 0 for n in range(FIRST_BAND, LAST_BAND + 1)}
    by_u2 = {n: {} for n in counts}
    examples = {n: [] for n in counts}
    term_counts = {}
    prior_excluded = int(k45["cumulative_profiles_excluded_K9_K45"])
    k45_hash = digest(K45)
    script_hash = digest(Path(__file__))

    for n in range(FIRST_BAND, LAST_BAND + 1):
        main_terms = packed_terms(
            templates(n, Fraction(5, 3), tuple(range(1, n + 1)), 5)
        )
        corrections = []
        # For K46--K49, r <= 9 is the full correction range (n < 50).
        for r in range(1, min(9, n // 5) + 1):
            total = n - 5 * r
            alpha = Fraction(5 - r, 3)
            corrections.append((r, packed_terms(
                templates(total, alpha, tuple(range(1, total + 1)), 5 - r)
            )))
        term_counts[str(n)] = {
            "main_F_5_3": len(main_terms[3]),
            "corrections_by_r": {str(r): len(packed[3]) for r, packed in corrections},
        }

        survivors = []
        for raw in profiles:
            if len(raw) != 5:
                raise AssertionError("malformed order profile")
            profile = tuple(
                INFINITY if value == "infinity" else int(value)
                for value in raw
            )
            main_order, main_ties, _finite, main_support, main_raw = term_orders(
                main_terms, profile
            )
            if main_order is None or main_order >= 0 or len(main_ties) != 1:
                survivors.append(raw)
                continue
            main_idx = int(main_ties[0])
            if not main_support[main_idx]:
                survivors.append(raw)
                continue

            correction_orders = {}
            finite_orders = []
            for r, correction in corrections:
                order, _ties, _mask, _support, _terms = term_orders(correction, profile)
                correction_orders[str(r)] = order
                if order is not None:
                    finite_orders.append(order)
            if not finite_orders or main_order >= min(finite_orders):
                survivors.append(raw)
                continue

            multiplicities, coefficient, _base = main_raw[main_idx]
            counts[n] += 1
            u2 = raw[1]
            key = "infinity" if u2 == "infinity" else str(u2)
            by_u2[n][key] = by_u2[n].get(key, 0) + 1
            if len(examples[n]) < 8:
                examples[n].append({
                    "orders_U1_to_U5": raw,
                    "main_pole_order": -int(main_order),
                    "main_term_multiplicities": list(multiplicities),
                    "main_coefficient": str(coefficient),
                    "correction_orders_by_r": correction_orders,
                })

        profiles = survivors
        excluded = prior_excluded + sum(counts.values())
        result = {
            "scope": "K46-K49 exact rational local valuation screen at the distinguished double root of the F2 2+1+1 nonzero-U1 deeper branch",
            "field": "exact Q generalized-binomial coefficients and integer order arithmetic",
            "source_receipts": {K45.name: k45_hash},
            "script_sha256": script_hash,
            "finite_profile_grid_count": GRID_COUNT,
            "profiles_entering_K46": k45["survivor_profiles_after_last_complete_band"],
            "last_complete_band": n,
            "first_excluding_band_counts": {str(j): counts[j] for j in counts},
            "first_excluding_band_by_ord_U2": {str(j): by_u2[j] for j in by_u2},
            "term_counts_by_band": term_counts,
            "new_profiles_excluded_K46_K49": sum(counts.values()),
            "cumulative_profiles_excluded_K9_K49": excluded,
            "cumulative_percent_of_finite_profile_grid": round(100 * excluded / GRID_COUNT, 6),
            "survivor_profiles_after_last_complete_band": len(profiles),
            "survivor_profiles": profiles,
            "examples": {str(j): examples[j] for j in examples},
            "status": "finite local-profile pruning only; whole 2+1+1 branch exclusion remains unproved",
        }
        save(OUT, result)
        print(json.dumps({
            "band": n,
            "newly_excluded": counts[n],
            "survivors": len(profiles),
            "cumulative_excluded": excluded,
            "finite_grid_percent": result["cumulative_percent_of_finite_profile_grid"],
            "terms": term_counts[str(n)],
            "receipt": OUT.name,
        }), flush=True)


if __name__ == "__main__":
    main()
