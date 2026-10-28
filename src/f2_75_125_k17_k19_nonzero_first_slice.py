"""Exact K17--K19 local pole screen for the F2 2+1+1 nonzero-U1 stratum.

This extends the K9--K16 finite valuation-profile screen.  For n=17,18,19
the first-integral coefficient has [w^n]F^(5/3), the kappa_1 term
[w^(n-5)]F^(4/3)/t, a kappa_2 term proportional to U_(n-10), and the
kappa_3 term [w^(n-15)]F^(2/3).  Since the chart is centered at t=1, any
power of t in these corrections is a unit.  Only unique negative main poles
on U1,...,U5 strictly below every correction order are counted.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

from f2_75_125_k17_nonzero_first_slice import (
    INFINITY,
    ORDERS,
    SUPPORT_CUTOFF,
    least_terms,
    sha256,
    templates,
)


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_k15_k16_nonzero_first_slice.json"
OUT = HERE / "f2_75_125_k17_k19_nonzero_first_slice.json"


def main() -> None:
    prior = json.loads(SOURCE.read_text(encoding="utf-8"))
    profiles = prior["survivor_profiles_after_K16"]
    if len(profiles) != prior["survivor_profile_count_after_K16"]:
        raise AssertionError("K16 survivor receipt count is inconsistent")

    band_terms = {}
    for n in range(17, 20):
        band_terms[n] = {
            "main": templates(n, Fraction(5, 3), tuple(range(1, n + 1)), 5),
            "kappa1": templates(n - 5, Fraction(4, 3), tuple(range(1, n - 4)), 4),
            "kappa3": templates(n - 15, Fraction(2, 3), tuple(range(1, n - 14)), 2),
        }
    expected = {
        17: (297, 77, 2),
        18: (385, 101, 3),
        19: (490, 135, 5),
    }
    for n, terms in band_terms.items():
        actual = (len(terms["main"]), len(terms["kappa1"]), len(terms["kappa3"]))
        if actual != expected[n]:
            raise AssertionError(f"K{n} partition census changed: {actual}")

    survivors = []
    first_hits = {17: 0, 18: 0, 19: 0}
    first_hits_by_u2 = {17: {}, 18: {}, 19: {}}
    examples = {17: [], 18: [], 19: []}
    for profile in profiles:
        ORDERS.clear()
        ORDERS.update({i: 0 for i in range(1, 20)})
        ORDERS.update({i: int(value) for i, value in enumerate(profile, start=1)})
        excluded = False
        for n in range(17, 20):
            main_order, tied = least_terms(band_terms[n]["main"])
            kappa1_order, _ = least_terms(band_terms[n]["kappa1"])
            kappa3_order, _ = least_terms(band_terms[n]["kappa3"])
            kappa2_order = ORDERS[n - 10]
            correction_candidates = [kappa2_order]
            if kappa1_order is not None:
                correction_candidates.append(kappa1_order)
            if kappa3_order is not None:
                correction_candidates.append(kappa3_order)
            correction_order = min(correction_candidates)
            obstructed = (
                main_order is not None
                and main_order < 0
                and len(tied) == 1
                and main_order < correction_order
                and not any(tied[0][0][SUPPORT_CUTOFF:])
            )
            if not obstructed:
                continue

            excluded = True
            first_hits[n] += 1
            q = "infinity" if profile[1] >= INFINITY else str(profile[1])
            first_hits_by_u2[n][q] = first_hits_by_u2[n].get(q, 0) + 1
            if len(examples[n]) < 12:
                multiplicities, coefficient = tied[0]
                examples[n].append({
                    "orders_U1_to_U5": profile,
                    "main_pole_order": -int(main_order),
                    "multiplicities": [i for i, k in enumerate(multiplicities, start=1) for _ in range(k)],
                    "coefficient": str(coefficient),
                    "kappa1_order": kappa1_order,
                    "kappa2_order": kappa2_order,
                    "kappa3_order": kappa3_order,
                    "full_correction_order": correction_order,
                })
            break
        if not excluded:
            survivors.append(profile)

    prior_union = prior["K9_to_K16_union_profiles"]
    added = sum(first_hits.values())
    cumulative = prior_union + added
    result = {
        "scope": "K17-K19 exact local valuation screen at the distinguished double root of the F2 2+1+1 nonzero-U1 deeper branch",
        "field": "exact Q coefficients and integer valuations",
        "input_receipt": str(SOURCE.relative_to(HERE.parent)),
        "input_receipt_sha256": sha256(SOURCE),
        "profiles_entering_K17": len(profiles),
        "new_profiles_by_first_excluding_band": first_hits,
        "new_profiles_by_ord_U2": first_hits_by_u2,
        "K9_to_K19_union_profiles": cumulative,
        "K9_to_K19_percent_of_finite_profile_grid": round(100 * cumulative / prior["profile_count"], 6),
        "survivor_profiles_after_K19": len(survivors),
        "survivor_profiles": survivors,
        "partition_counts": {str(n): {key: len(value) for key, value in terms.items()}
                              for n, terms in band_terms.items()},
        "examples": examples,
        "status": "finite local-profile pruning only; whole 2+1+1 branch exclusion remains unproved",
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    concise = {key: value for key, value in result.items()
               if key not in {"survivor_profiles"}}
    print(json.dumps(concise, indent=2), flush=True)
    print(f"receipt={OUT} sha256={hashlib.sha256(OUT.read_bytes()).hexdigest()}", flush=True)


if __name__ == "__main__":
    main()
