"""Exact finite-profile continuation through K50-K54.

The K50 boundary introduces correction index r=10, so this continuation
includes the full range 1 <= r <= floor(n/5) rather than the r<=9 range that
was complete only through K49. This remains local-profile pruning, not a
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
K49 = HERE / "f2_75_125_211_k46_k49_nonzero_first_slice_20261022.json"
OUT = HERE / "f2_75_125_211_k50_k54_nonzero_first_slice_20261023.json"
FIRST_BAND = 50
LAST_BAND = 54
GRID_COUNT = 51_870


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path: Path, value: dict) -> None:
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    temp.replace(path)


def main() -> None:
    prior = json.loads(K49.read_text(encoding="utf-8"))
    if prior.get("last_complete_band") != 49:
        raise AssertionError("input receipt is not complete through K49")
    profiles = prior["survivor_profiles"]
    if len(profiles) != prior.get("survivor_profiles_after_last_complete_band"):
        raise AssertionError("K49 survivor count does not match its list")
    if prior.get("finite_profile_grid_count") != GRID_COUNT:
        raise AssertionError("unexpected finite profile grid")

    counts = {n: 0 for n in range(FIRST_BAND, LAST_BAND + 1)}
    by_u2 = {n: {} for n in counts}
    examples = {n: [] for n in counts}
    term_counts = {}
    prior_excluded = int(prior["cumulative_profiles_excluded_K9_K49"])
    prior_hash = digest(K49)
    script_hash = digest(Path(__file__))

    for n in range(FIRST_BAND, LAST_BAND + 1):
        main_terms = packed_terms(
            templates(n, Fraction(5, 3), tuple(range(1, n + 1)), 5)
        )
        corrections = []
        # Include the new r=10 resonance beginning at K50.
        for r in range(1, n // 5 + 1):
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
            key = "infinity" if raw[1] == "infinity" else str(raw[1])
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
            "scope": "K50-K54 exact rational local valuation screen at the distinguished double root of the F2 2+1+1 nonzero-U1 branch",
            "field": "exact Q generalized-binomial coefficients and integer order arithmetic",
            "source_receipts": {K49.name: prior_hash},
            "script_sha256": script_hash,
            "finite_profile_grid_count": GRID_COUNT,
            "profiles_entering_K50": prior["survivor_profiles_after_last_complete_band"],
            "last_complete_band": n,
            "first_excluding_band_counts": {str(j): counts[j] for j in counts},
            "first_excluding_band_by_ord_U2": {str(j): by_u2[j] for j in by_u2},
            "term_counts_by_band": term_counts,
            "new_profiles_excluded_K50_K54": sum(counts.values()),
            "cumulative_profiles_excluded_K9_K54": excluded,
            "cumulative_percent_of_finite_profile_grid": round(100 * excluded / GRID_COUNT, 6),
            "survivor_profiles_after_last_complete_band": len(profiles),
            "survivor_profiles": profiles,
            "examples": {str(j): examples[j] for j in examples},
            "status": "finite local-profile pruning only; whole 2+1+1 branch exclusion remains unproved",
        }
        atomic_json(OUT, result)
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
