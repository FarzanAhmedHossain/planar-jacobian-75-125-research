"""Exact K17 local pole screen for the F2 2+1+1 nonzero-U1 stratum.

The input is the exact K9--K16 profile receipt.  At K17 the first-integral
coefficient has main term [w^17]F^(5/3), a kappa_1 correction
[w^12]F^(4/3)/t, a kappa_2 correction proportional to U7/t^a (a is
irrelevant at t=1), and a kappa_3 correction [w^2]F^(2/3).  We count only
unique negative main poles supported on U1,...,U5 and strictly below every
possible correction order.  This prunes a finite local order-profile set; it
does not by itself exclude the full quartic branch.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from functools import lru_cache
from math import factorial
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_k15_k16_nonzero_first_slice.json"
OUT = HERE / "f2_75_125_k17_nonzero_first_slice.json"
H_ORDER = 2
SUPPORT_CUTOFF = 5
INFINITY = 13
ORDERS: dict[int, int] = {}


@lru_cache(maxsize=None)
def partitions(total: int, weights: tuple[int, ...]):
    if not weights:
        return ((),) if total == 0 else ()
    weight, *tail = weights
    out = []
    for multiplicity in range(total // weight + 1):
        for suffix in partitions(total - multiplicity * weight, tuple(tail)):
            out.append((multiplicity,) + suffix)
    return tuple(out)


def generalized_binomial(alpha: Fraction, n: int) -> Fraction:
    value = Fraction(1)
    for i in range(n):
        value *= alpha - i
    return value / factorial(n)


def coefficient(alpha: Fraction, multiplicities: tuple[int, ...]) -> Fraction:
    count = sum(multiplicities)
    value = generalized_binomial(alpha, count) * factorial(count)
    for multiplicity in multiplicities:
        value /= factorial(multiplicity)
    return value


def templates(total: int, alpha: Fraction, weights: tuple[int, ...], exponent_base: int):
    return tuple(
        (
            multiplicities,
            coefficient(alpha, multiplicities),
            H_ORDER * (exponent_base - 3 * sum(multiplicities)),
        )
        for multiplicities in partitions(total, weights)
    )


def least_terms(terms):
    evaluated = []
    for multiplicities, value, base_order in terms:
        if not value:
            continue
        if any(
            multiplicity and ORDERS[weight] >= INFINITY
            for weight, multiplicity in enumerate(multiplicities, start=1)
        ):
            continue
        order = base_order + sum(
            multiplicity * ORDERS[weight]
            for weight, multiplicity in enumerate(multiplicities, start=1)
        )
        evaluated.append((order, multiplicities, value))
    if not evaluated:
        return None, []
    least = min(order for order, _, _ in evaluated)
    return least, [(m, value) for order, m, value in evaluated if order == least]


def label(multiplicities: tuple[int, ...]) -> str:
    return "+".join(
        str(weight)
        for weight, multiplicity in enumerate(multiplicities, start=1)
        for _ in range(multiplicity)
    ) or "1"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    prior = json.loads(SOURCE.read_text(encoding="utf-8"))
    profiles = prior["survivor_profiles_after_K16"]
    if len(profiles) != prior["survivor_profile_count_after_K16"]:
        raise AssertionError("K16 survivor receipt count is inconsistent")

    main_terms = templates(17, Fraction(5, 3), tuple(range(1, 18)), 5)
    kappa1_terms = templates(12, Fraction(4, 3), tuple(range(1, 13)), 4)
    kappa3_terms = templates(2, Fraction(2, 3), (1, 2), 2)
    if len(main_terms) != 297 or len(kappa1_terms) != 77 or len(kappa3_terms) != 2:
        raise AssertionError("K17 partition census changed")

    new_profiles = []
    survivors = []
    examples = []
    by_ord_u2: dict[str, int] = {}
    for raw in profiles:
        if len(raw) != 5:
            raise AssertionError("malformed K16 survivor profile")
        ORDERS.clear()
        ORDERS.update({i: 0 for i in range(1, 18)})
        ORDERS.update({i: (INFINITY if value == "infinity" else int(value))
                       for i, value in enumerate(raw, start=1)})

        main_order, tied = least_terms(main_terms)
        kappa1_order, _ = least_terms(kappa1_terms)
        kappa3_order, _ = least_terms(kappa3_terms)
        # U7 has unconstrained order >= 0 in the K17 chart; t is a unit at t=1.
        correction_order = min(kappa1_order, 0, kappa3_order)
        obstructed = (
            main_order is not None
            and main_order < 0
            and len(tied) == 1
            and main_order < correction_order
            and not any(tied[0][0][SUPPORT_CUTOFF:])
        )
        if not obstructed:
            survivors.append(raw)
            continue

        multiplicities, value = tied[0]

        new_profiles.append(raw)
        q = "infinity" if raw[1] == "infinity" else str(raw[1])
        by_ord_u2[q] = by_ord_u2.get(q, 0) + 1
        if len(examples) < 20:
            examples.append({
                "orders_U1_to_U5": raw,
                "main_pole_order": -int(main_order),
                "unique_main_term": label(multiplicities),
                "coefficient": str(value),
                "kappa1_order": kappa1_order,
                "kappa2_order_lower_bound": 0,
                "kappa3_order": kappa3_order,
                "full_correction_order": correction_order,
            })

    cumulative = prior["K9_to_K16_union_profiles"] + len(new_profiles)
    result = {
        "scope": "K17 local valuation screen at the distinguished double root of the F2 2+1+1 nonzero-U1 deeper branch",
        "field": "exact Q coefficients and integer valuations",
        "input_receipt": str(SOURCE.relative_to(HERE.parent)),
        "input_receipt_sha256": sha256(SOURCE),
        "prior_profile_count": len(profiles),
        "K17_new_profiles": len(new_profiles),
        "K17_new_by_ord_U2": by_ord_u2,
        "K9_to_K17_union_profiles": cumulative,
        "K9_to_K17_percent_of_finite_profile_grid": round(100 * cumulative / prior["profile_count"], 6),
        "survivor_profile_count_after_K17": len(survivors),
        "survivor_profiles_after_K17": survivors,
        "K17_partition_counts": {
            "main_w17_F_5_3": len(main_terms),
            "correction_w12_F_4_3": len(kappa1_terms),
            "correction_w2_F_2_3": len(kappa3_terms),
        },
        "examples": examples,
        "status": "finite local-profile pruning only; whole 2+1+1 branch exclusion remains unproved",
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in result.items()
                      if key not in {"survivor_profiles_after_K17"}}, indent=2), flush=True)
    print(f"receipt={OUT} sha256={sha256(OUT)}", flush=True)


if __name__ == "__main__":
    main()
