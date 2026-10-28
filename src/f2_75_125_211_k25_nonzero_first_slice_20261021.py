"""Exact K25 extension of the F2 2+1+1 double-root profile screen.

Starts from the local order profiles surviving K24. It applies the K25
coefficient of the first-integral identity, including the r=5 correction
which is the regular scalar kappa_5 [w^0] F^0. Higher slices are assigned
order zero, so a profile is removed only by a unique negative main pole
strictly below every possible correction order. This is finite local-profile
pruning, not a whole-branch exclusion.
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
OUT = HERE / "f2_75_125_211_k25_nonzero_first_slice_20261021.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    prior = json.loads(SOURCE.read_text(encoding="utf-8"))
    profiles = prior["survivor_profiles"]
    if len(profiles) != prior["survivor_profiles_after_K24"]:
        raise AssertionError("K24 survivor receipt count is inconsistent")

    n = 25
    main = packed_terms(templates(n, Fraction(5, 3), tuple(range(1, n + 1)), 5))
    corrections = []
    for r in range(1, 6):
        total = n - 5 * r
        alpha = Fraction(5 - r, 3)
        terms = templates(total, alpha, tuple(range(1, total + 1)), 5 - r)
        corrections.append((r, packed_terms(terms)))

    excluded_profiles = []
    survivors = []
    by_order_u2: dict[str, int] = {}
    examples = []
    for raw_profile in profiles:
        if len(raw_profile) != 5:
            raise AssertionError("malformed K24 survivor profile")
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
        excluded_profiles.append(raw_profile)
        u2 = raw_profile[1]
        key = "infinity" if u2 == "infinity" else str(u2)
        by_order_u2[key] = by_order_u2.get(key, 0) + 1
        if len(examples) < 16:
            examples.append({
                "orders_U1_to_U5": raw_profile,
                "main_pole_order": -int(main_order),
                "main_term_multiplicities": list(multiplicities),
                "main_coefficient": str(coefficient),
                "correction_orders_by_r": correction_orders,
            })

    result = {
        "scope": "K25 local valuation screen at the distinguished double root of the F2 2+1+1 nonzero-U1 deeper branch",
        "field": "exact Q generalized-binomial coefficients and integer order arithmetic",
        "input_receipt": SOURCE.name,
        "input_receipt_sha256": digest(SOURCE),
        "profiles_entering_K25": len(profiles),
        "main_F_5_3_partition_count": len(main[3]),
        "correction_partition_counts_by_r": {
            str(r): len(correction[3]) for r, correction in corrections
        },
        "K25_r5_correction": "kappa_5 [w^0]F^0 = kappa_5, regular of order 0 at the nonzero root",
        "new_profiles_excluded": len(excluded_profiles),
        "new_profiles_excluded_by_ord_U2": by_order_u2,
        "profiles_surviving_K25": len(survivors),
        "survivor_profiles_after_K25": survivors,
        "examples": examples,
        "status": "finite local-profile pruning only; whole 2+1+1 branch exclusion remains unproved",
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    concise = {key: value for key, value in result.items() if key != "survivor_profiles_after_K25"}
    print(json.dumps(concise, indent=2), flush=True)
    print(f"receipt={OUT} sha256={digest(OUT)}", flush=True)


if __name__ == "__main__":
    main()
