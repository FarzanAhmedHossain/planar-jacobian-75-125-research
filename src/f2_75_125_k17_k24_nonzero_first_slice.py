"""Vectorized exact-coefficient local screen through K24, F2 2+1+1 branch.

This starts from the 29,040 order profiles surviving K16.  The first-integral
identity is expanded as sum_r kappa_r w^(5r) F^((5-r)/3).  At the centered
root t=1, the powers of t in kappa_r are units, so K_n's r-th correction
has the valuation of [w^(n-5r)]F^((5-r)/3).  We count only a unique negative
main term supported on U1,...,U5 and strictly below every correction term.
This is a finite local-profile screen, not a whole-branch exclusion.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

from f2_75_125_k17_nonzero_first_slice import INFINITY, templates


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_k15_k16_nonzero_first_slice.json"
OUT = HERE / "f2_75_125_k17_k24_nonzero_first_slice.json"
PROFILE_COUNT = 51_870
PRIOR_UNION = 22_830


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def packed_terms(terms):
    """Pack nonzero exact terms as (first-five exponents, base order, support)."""
    rows = []
    base = []
    support = []
    raw_nonzero = []
    for multiplicities, coefficient, base_order in terms:
        if not coefficient:
            continue
        first_five = list(multiplicities[:5]) + [0] * max(0, 5 - len(multiplicities))
        rows.append(first_five)
        base.append(base_order)
        support.append(not any(multiplicities[5:]))
        raw_nonzero.append((multiplicities, coefficient, base_order))
    return (
        np.asarray(rows, dtype=np.int16),
        np.asarray(base, dtype=np.int16),
        np.asarray(support, dtype=np.bool_),
        raw_nonzero,
    )


def term_orders(packed, profile):
    exponents, base, support, raw_terms = packed
    finite = np.ones(len(base), dtype=np.bool_)
    for j, order in enumerate(profile):
        if order >= INFINITY:
            finite &= exponents[:, j] == 0
    if not np.any(finite):
        return None, np.empty(0, dtype=np.int64), finite, support, raw_terms
    orders = base.astype(np.int32) + exponents.astype(np.int32) @ np.asarray(profile, dtype=np.int32)
    least = int(orders[finite].min())
    tied = np.flatnonzero(finite & (orders == least))
    return least, tied, finite, support, raw_terms


def main() -> None:
    prior = json.loads(SOURCE.read_text(encoding="utf-8"))
    profiles = prior["survivor_profiles_after_K16"]
    if len(profiles) != prior["survivor_profile_count_after_K16"]:
        raise AssertionError("K16 survivor receipt count is inconsistent")
    if prior["profile_count"] != PROFILE_COUNT or prior["K9_to_K16_union_profiles"] != PRIOR_UNION:
        raise AssertionError("input finite-grid census changed")

    bands = {}
    counts_by_band: dict[int, int] = {}
    counts_by_u2: dict[int, dict[str, int]] = {}
    examples: dict[int, list[dict[str, object]]] = {}
    for n in range(17, 25):
        main_terms = templates(n, Fraction(5, 3), tuple(range(1, n + 1)), 5)
        main = packed_terms(main_terms)
        corrections = []
        correction_counts = {}
        for r in range(1, min(4, n // 5) + 1):
            total = n - 5 * r
            alpha = Fraction(5 - r, 3)
            terms = templates(total, alpha, tuple(range(1, total + 1)), 5 - r)
            corrections.append((r, packed_terms(terms)))
            correction_counts[str(r)] = sum(bool(entry[1]) for entry in terms)
        bands[n] = (main, corrections)
        counts_by_band[n] = 0
        counts_by_u2[n] = {}
        examples[n] = []

    survivors = []
    for raw_profile in profiles:
        if len(raw_profile) != 5:
            raise AssertionError("malformed K16 survivor profile")
        profile = tuple(int(value) for value in raw_profile)
        excluded = False
        for n in range(17, 25):
            main_order, main_ties, _finite, main_support, main_raw = term_orders(bands[n][0], profile)
            if main_order is None or main_order >= 0 or len(main_ties) != 1:
                continue
            main_idx = int(main_ties[0])
            if not main_support[main_idx]:
                continue

            correction_orders = []
            by_r = {}
            for r, correction in bands[n][1]:
                correction_order, _ties, _mask, _support, _terms = term_orders(correction, profile)
                by_r[str(r)] = correction_order
                if correction_order is not None:
                    correction_orders.append(correction_order)
            if not correction_orders or main_order >= min(correction_orders):
                continue

            multiplicities, coefficient, _base = main_raw[main_idx]
            counts_by_band[n] += 1
            u2 = "infinity" if profile[1] >= INFINITY else str(profile[1])
            counts_by_u2[n][u2] = counts_by_u2[n].get(u2, 0) + 1
            if len(examples[n]) < 12:
                examples[n].append({
                    "orders_U1_to_U5": list(profile),
                    "main_pole_order": -main_order,
                    "main_term_multiplicities": [i for i, power in enumerate(multiplicities, start=1) for _ in range(power)],
                    "main_coefficient": str(coefficient),
                    "correction_orders_by_r": by_r,
                })
            excluded = True
            break
        if not excluded:
            survivors.append(raw_profile)

    added = sum(counts_by_band.values())
    cumulative = PRIOR_UNION + added
    result = {
        "scope": "K17-K24 exact rational-coefficient local valuation screen at the distinguished double root of the F2 2+1+1 nonzero-U1 deeper branch",
        "field": "exact Q binomial coefficients and integer valuations; vectorized integer order arithmetic",
        "input_receipt": str(SOURCE.relative_to(HERE.parent)),
        "input_receipt_sha256": digest(SOURCE),
        "profiles_entering_K17": len(profiles),
        "new_profiles_by_first_excluding_band": {str(k): v for k, v in counts_by_band.items()},
        "new_profiles_by_ord_U2": {str(k): v for k, v in counts_by_u2.items()},
        "K9_to_K24_union_profiles": cumulative,
        "K9_to_K24_percent_of_finite_profile_grid": round(100 * cumulative / PROFILE_COUNT, 6),
        "survivor_profiles_after_K24": len(survivors),
        "survivor_profiles": survivors,
        "per_band_partition_counts": {
            str(n): {
                "main_F_5_3": sum(bool(entry[1]) for entry in bands[n][0][3]),
                "corrections_by_r": {
                    str(r): sum(bool(entry[1]) for entry in packed[3])
                    for r, packed in bands[n][1]
                },
            }
            for n in range(17, 25)
        },
        "examples": {str(k): v for k, v in examples.items()},
        "status": "finite local-profile pruning only; whole 2+1+1 branch exclusion remains unproved",
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    concise = {key: value for key, value in result.items() if key != "survivor_profiles"}
    print(json.dumps(concise, indent=2), flush=True)
    print(f"receipt={OUT} sha256={digest(OUT)}", flush=True)


if __name__ == "__main__":
    main()
