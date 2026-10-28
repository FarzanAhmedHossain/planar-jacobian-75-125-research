"""Exact local valuation census for the U1=0, 2+1+1 double-root branch.

The K9--K15 note already excludes q=3 and several q=4 local jets.  This
checker enumerates the remaining bounded order profiles with q>=4, s>=2,
v>=1, r>=3, omitting the two strata already eliminated by the exact
K11--K15 residue identities.  It screens K15--K24 using the full first-
integral expansion sum_r kappa_r t^(-r) w^(5r) F^((5-r)/3).  Only unique
negative main poles supported on U2,...,U5 and below all resonance terms
are counted.  The result is local profile pruning, not branch exclusion.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from f2_75_125_k17_nonzero_first_slice import templates


HERE = Path(__file__).resolve().parent
OUT = HERE / "f2_75_125_k15_k24_zero_u1_double_root.json"
INFINITY = 13
H_ORDER = 2


def pack_terms(terms, weights):
    rows = []
    base = []
    supported = []
    raw = []
    for multiplicities, coefficient, base_order in terms:
        if not coefficient:
            continue
        exponents = [0, 0, 0, 0]
        for weight, multiplicity in zip(weights, multiplicities):
            if 2 <= weight <= 5:
                exponents[weight - 2] = multiplicity
        rows.append(exponents)
        base.append(base_order)
        supported.append(not any(
            multiplicity
            for weight, multiplicity in zip(weights, multiplicities)
            if weight >= 6
        ))
        raw.append((weights, multiplicities, coefficient, base_order))
    return (
        np.asarray(rows, dtype=np.int16),
        np.asarray(base, dtype=np.int16),
        np.asarray(supported, dtype=np.bool_),
        raw,
    )


def least_order(packed, profile):
    exponents, base, supported, raw = packed
    finite = np.ones(len(base), dtype=np.bool_)
    for j, order in enumerate(profile):
        if order >= INFINITY:
            finite &= exponents[:, j] == 0
    if not np.any(finite):
        return None, np.empty(0, dtype=np.int64), supported, raw
    values = base.astype(np.int32) + exponents.astype(np.int32) @ np.asarray(profile, dtype=np.int32)
    minimum = int(values[finite].min())
    tied = np.flatnonzero(finite & (values == minimum))
    return minimum, tied, supported, raw


def main() -> None:
    q_values = [*range(4, 12), INFINITY]
    r_values = [*range(3, 12), INFINITY]
    s_values = [*range(2, 12), INFINITY]
    v_values = [*range(1, 13), INFINITY]
    profiles = []
    prefiltered = 0
    for q, r, s, v in product(q_values, r_values, s_values, v_values):
        # Exact K11--K15 residues already exclude q=4,s=2,r=3 for all v,
        # and q=4,s=2,r>=4,v=1.
        if q == 4 and s == 2 and r == 3:
            prefiltered += 1
            continue
        if q == 4 and s == 2 and r >= 4 and v == 1:
            prefiltered += 1
            continue
        profiles.append((q, r, s, v))

    bands = {}
    for n in range(15, 25):
        main_weights = tuple(range(2, n + 1))
        main_terms = templates(n, Fraction(5, 3), main_weights, 5)
        main = pack_terms(main_terms, main_weights)
        corrections = []
        for resonance in range(1, min(4, n // 5) + 1):
            total = n - 5 * resonance
            weights = tuple(range(2, total + 1))
            terms = templates(total, Fraction(5 - resonance, 3), weights, 5 - resonance)
            corrections.append((resonance, pack_terms(terms, weights)))
        bands[n] = (main, corrections)

    by_band = {n: 0 for n in range(15, 25)}
    by_q = {n: {} for n in range(15, 25)}
    examples = {n: [] for n in range(15, 25)}
    survivors = []
    for profile in profiles:
        excluded = False
        for n in range(15, 25):
            main_order, main_ties, main_support, main_raw = least_order(bands[n][0], profile)
            if main_order is None or main_order >= 0 or len(main_ties) != 1:
                continue
            main_index = int(main_ties[0])
            if not main_support[main_index]:
                continue
            correction_orders = {}
            for resonance, correction in bands[n][1]:
                order, _ties, _support, _raw = least_order(correction, profile)
                correction_orders[str(resonance)] = order
            finite_corrections = [x for x in correction_orders.values() if x is not None]
            if not finite_corrections or main_order >= min(finite_corrections):
                continue

            weights, multiplicities, coefficient, _base = main_raw[main_index]
            by_band[n] += 1
            q = "infinity" if profile[0] >= INFINITY else str(profile[0])
            by_q[n][q] = by_q[n].get(q, 0) + 1
            if len(examples[n]) < 12:
                term_label = [
                    weight
                    for weight, multiplicity in zip(weights, multiplicities)
                    for _ in range(multiplicity)
                ]
                examples[n].append({
                    "orders_U2_to_U5": list(profile),
                    "main_pole_order": -main_order,
                    "main_term_partition": term_label,
                    "main_coefficient": str(coefficient),
                    "correction_orders_by_resonance": correction_orders,
                })
            excluded = True
            break
        if not excluded:
            survivors.append(profile)

    result = {
        "scope": "K15-K24 exact local valuation screen on the U1=0 deeper-double-root stratum of the F2 2+1+1 quartic carrier",
        "field": "exact Q generalized-binomial coefficients and integer valuations",
        "order_caps": {"U2": [4, 11, "infinity"], "U3": [3, 11, "infinity"], "U4": [2, 11, "infinity"], "U5_over_t": [1, 12, "infinity"]},
        "prior_exact_strata_excluded": ["q=4,s=2,r=3, all v", "q=4,s=2,r>=4,v=1"],
        "profiles_after_prior_filters": len(profiles),
        "profiles_prefiltered_by_prior_exact_residues": prefiltered,
        "new_profiles_by_first_excluding_band": {str(k): v for k, v in by_band.items()},
        "new_profiles_by_q": {str(k): v for k, v in by_q.items()},
        "survivor_profile_count_after_K24": len(survivors),
        "survivor_profiles": survivors,
        "examples": {str(k): v for k, v in examples.items()},
        "status": "local order-profile screen only; not an exclusion of the U1=0 branch or the full 2+1+1 branch",
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    concise = {key: value for key, value in result.items() if key != "survivor_profiles"}
    print(json.dumps(concise, indent=2), flush=True)
    print(f"receipt={OUT} sha256={hashlib.sha256(OUT.read_bytes()).hexdigest()}", flush=True)


if __name__ == "__main__":
    main()
