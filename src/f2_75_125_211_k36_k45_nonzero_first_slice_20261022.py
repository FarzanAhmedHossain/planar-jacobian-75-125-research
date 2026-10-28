"""Exact extension of the finite F2 2+1+1 double-root profile screen.

This continues the K26--K35 valuation screen through K36--K45. It works only
on the finite order-profile grid surviving K35, using exact rational binomial
coefficients and integer valuations. Each completed band is atomically saved
so a long run can be resumed. This is profile pruning, not a whole-branch
exclusion.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

from f2_75_125_k17_nonzero_first_slice import INFINITY, templates
from f2_75_125_k17_k24_nonzero_first_slice import packed_terms, term_orders


HERE = Path(__file__).resolve().parent
K24 = HERE / "f2_75_125_k17_k24_nonzero_first_slice.json"
K35 = HERE / "f2_75_125_211_k26_k35_nonzero_first_slice_20261021.json"
OUT = HERE / "f2_75_125_211_k36_k45_nonzero_first_slice_20261022.json"
FIRST_BAND = 36
LAST_BAND = 45


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def checkpoint(
    *, script_hash: str, k24_hash: str, k35_hash: str,
    last_complete_band: int, profiles: list, counts: dict,
    by_u2: dict, examples: dict, term_counts: dict, grid_count: int,
    prior_union: int, k26_k35_excluded: int,
) -> dict:
    new_excluded = sum(counts.values())
    cumulative_excluded = prior_union + k26_k35_excluded + new_excluded
    return {
        "scope": "K36-K45 exact rational local valuation screen at the distinguished double root of the F2 2+1+1 nonzero-U1 deeper branch",
        "field": "exact Q generalized-binomial coefficients and integer order arithmetic",
        "source_receipts": {
            K24.name: k24_hash,
            K35.name: k35_hash,
        },
        "script_sha256": script_hash,
        "finite_profile_grid_count": grid_count,
        "profiles_entering_K36": grid_count - prior_union - k26_k35_excluded,
        "last_complete_band": last_complete_band,
        "first_excluding_band_counts": {str(n): counts[n] for n in range(FIRST_BAND, LAST_BAND + 1)},
        "first_excluding_band_by_ord_U2": {str(n): by_u2[n] for n in range(FIRST_BAND, LAST_BAND + 1)},
        "term_counts_by_band": {str(n): term_counts.get(n, {}) for n in range(FIRST_BAND, LAST_BAND + 1)},
        "new_profiles_excluded_K36_K45": new_excluded,
        "cumulative_profiles_excluded_K9_K45": cumulative_excluded,
        "cumulative_percent_of_finite_profile_grid": round(100 * cumulative_excluded / grid_count, 6),
        "survivor_profiles_after_last_complete_band": len(profiles),
        "survivor_profiles": profiles,
        "examples": {str(n): examples[n] for n in range(FIRST_BAND, LAST_BAND + 1)},
        "status": "finite local-profile pruning only; whole 2+1+1 branch exclusion remains unproved",
    }


def main() -> None:
    script_hash = digest(Path(__file__))
    k24_hash = digest(K24)
    k35_hash = digest(K35)
    k24 = json.loads(K24.read_text(encoding="utf-8"))
    k35 = json.loads(K35.read_text(encoding="utf-8"))
    if k35["input_receipt_sha256"] != k24_hash:
        raise AssertionError("K35 receipt does not match the K24 input")
    if k35["survivor_profiles_after_K35"] != len(k35["survivor_profiles"]):
        raise AssertionError("K35 survivor count is inconsistent")
    if k35["survivor_profiles"] != k24["survivor_profiles"]:
        raise AssertionError("K26-K35 profile list changed from K24 despite its zero-kill receipt")
    grid_count = 51_870
    prior_union = int(k24["K9_to_K24_union_profiles"])
    k26_k35_excluded = sum(int(v) for v in k35["first_excluding_band_counts"].values())

    counts = {n: 0 for n in range(FIRST_BAND, LAST_BAND + 1)}
    by_u2 = {n: {} for n in range(FIRST_BAND, LAST_BAND + 1)}
    examples = {n: [] for n in range(FIRST_BAND, LAST_BAND + 1)}
    term_counts = {}
    current_profiles = k35["survivor_profiles"]
    start_band = FIRST_BAND

    if OUT.exists():
        previous = json.loads(OUT.read_text(encoding="utf-8"))
        if previous.get("script_sha256") != script_hash:
            raise AssertionError("existing checkpoint belongs to a different script")
        if previous.get("source_receipts") != {K24.name: k24_hash, K35.name: k35_hash}:
            raise AssertionError("existing checkpoint has different input receipts")
        last = int(previous["last_complete_band"])
        if last < FIRST_BAND or last >= LAST_BAND:
            if last == LAST_BAND:
                print(json.dumps({k: v for k, v in previous.items() if k != "survivor_profiles"}, indent=2))
                return
            raise AssertionError("invalid last-complete-band marker")
        current_profiles = previous["survivor_profiles"]
        counts = {n: int(previous["first_excluding_band_counts"][str(n)]) for n in counts}
        by_u2 = {n: previous["first_excluding_band_by_ord_U2"][str(n)] for n in by_u2}
        examples = {n: previous["examples"][str(n)] for n in examples}
        term_counts = {int(n): value for n, value in previous["term_counts_by_band"].items() if value}
        start_band = last + 1

    for n in range(start_band, LAST_BAND + 1):
        main = packed_terms(templates(n, Fraction(5, 3), tuple(range(1, n + 1)), 5))
        corrections = []
        for r in range(1, min(9, n // 5) + 1):
            total = n - 5 * r
            alpha = Fraction(5 - r, 3)
            packed = packed_terms(templates(total, alpha, tuple(range(1, total + 1)), 5 - r))
            corrections.append((r, packed))
        term_counts[n] = {
            "main_F_5_3": len(main[3]),
            "corrections_by_r": {str(r): len(packed[3]) for r, packed in corrections},
        }

        survivors = []
        for raw_profile in current_profiles:
            if len(raw_profile) != 5:
                raise AssertionError("malformed order profile")
            profile = tuple(
                INFINITY if value == "infinity" else int(value)
                for value in raw_profile
            )
            main_order, main_ties, _finite, main_support, main_raw = term_orders(main, profile)
            if main_order is None or main_order >= 0 or len(main_ties) != 1:
                survivors.append(raw_profile)
                continue
            main_idx = int(main_ties[0])
            if not main_support[main_idx]:
                survivors.append(raw_profile)
                continue

            correction_orders = {}
            finite_orders = []
            for r, correction in corrections:
                order, _ties, _mask, _support, _terms = term_orders(correction, profile)
                correction_orders[str(r)] = order
                if order is not None:
                    finite_orders.append(order)
            if not finite_orders or main_order >= min(finite_orders):
                survivors.append(raw_profile)
                continue

            multiplicities, coefficient, _base = main_raw[main_idx]
            counts[n] += 1
            u2 = raw_profile[1]
            key = "infinity" if u2 == "infinity" else str(u2)
            by_u2[n][key] = by_u2[n].get(key, 0) + 1
            if len(examples[n]) < 8:
                examples[n].append({
                    "orders_U1_to_U5": raw_profile,
                    "main_pole_order": -int(main_order),
                    "main_term_multiplicities": list(multiplicities),
                    "main_coefficient": str(coefficient),
                    "correction_orders_by_r": correction_orders,
                })

        current_profiles = survivors
        value = checkpoint(
            script_hash=script_hash, k24_hash=k24_hash, k35_hash=k35_hash,
            last_complete_band=n, profiles=current_profiles, counts=counts,
            by_u2=by_u2, examples=examples, term_counts=term_counts,
            grid_count=grid_count, prior_union=prior_union,
            k26_k35_excluded=k26_k35_excluded,
        )
        atomic_json(OUT, value)
        print(json.dumps({
            "band": n,
            "newly_excluded": counts[n],
            "survivors": len(current_profiles),
            "cumulative_excluded": value["cumulative_profiles_excluded_K9_K45"],
            "finite_grid_percent": value["cumulative_percent_of_finite_profile_grid"],
            "term_counts": term_counts[n],
            "receipt": str(OUT),
        }), flush=True)


if __name__ == "__main__":
    main()
