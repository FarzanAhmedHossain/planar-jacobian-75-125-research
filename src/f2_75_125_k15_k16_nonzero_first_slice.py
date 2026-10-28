"""Exact K9--K16 local valuation census on the 2+1+1 nonzero-U1 branch.

This replays the K9--K14 order-profile flags, then tests K15 and K16 on the
profiles that survive. The first-integral expansions add only regular terms
besides the kappa_1 correction at these two bands:
  G_15 = [w^15]F^(5/3) + kappa_1/t [w^10]F^(4/3)
         + kappa_2 U5/t + kappa_3 h^2,
  G_16 = [w^16]F^(5/3) + kappa_1/t [w^11]F^(4/3)
         + kappa_2 U6/t + kappa_3 [w]F^(2/3).
The checker counts only a unique negative main pole supported on U1,...,U5
and strictly below every resonance correction. Higher-slice orders use zero
as the worst-case baseline. Results are local profile pruning, not whole-
branch coverage.
"""

from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from math import factorial
from pathlib import Path


DEGREE_CAP = 11
U5_DEGREE_CAP = 12
U1_ORDER_MAX = 7
U2_ORDER_MAX = 9
INFINITY = U5_DEGREE_CAP + 1
H_ORDER = 2
SUPPORT_CUTOFF = 5


@lru_cache(maxsize=None)
def partitions(total: int, weights: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    if not weights:
        return ((),) if total == 0 else ()
    weight, *tail = weights
    result: list[tuple[int, ...]] = []
    for multiplicity in range(total // weight + 1):
        for suffix in partitions(total - multiplicity * weight, tuple(tail)):
            result.append((multiplicity,) + suffix)
    return tuple(result)


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


def evaluated_terms(terms, orders: dict[int, int]):
    out = []
    for multiplicities, value, base_order in terms:
        if value == 0:
            continue
        if any(
            multiplicity and orders[weight] >= INFINITY
            for weight, multiplicity in enumerate(multiplicities, start=1)
        ):
            continue
        order = base_order + sum(
            multiplicity * orders[weight]
            for weight, multiplicity in enumerate(multiplicities, start=1)
        )
        out.append((order, multiplicities, value))
    return out


def least_terms(terms, orders: dict[int, int]):
    evaluated = evaluated_terms(terms, orders)
    if not evaluated:
        return None, []
    least = min(order for order, _, _ in evaluated)
    return least, [(m, value) for order, m, value in evaluated if order == least]


def robust_pole(main_terms, correction_order: int | None):
    main_order, tied = least_terms(main_terms, ORDERS)
    if main_order is None or main_order >= 0 or len(tied) != 1:
        return False, main_order, tied
    multiplicities, _ = tied[0]
    if any(multiplicities[SUPPORT_CUTOFF:]):
        return False, main_order, tied
    if correction_order is not None and main_order >= correction_order:
        return False, main_order, tied
    return True, main_order, tied


def label(multiplicities: tuple[int, ...]) -> str:
    return "+".join(
        str(weight)
        for weight, multiplicity in enumerate(multiplicities, start=1)
        for _ in range(multiplicity)
    ) or "1"


ORDERS: dict[int, int] = {}


def main() -> None:
    main_terms = {
        n: templates(n, Fraction(5, 3), tuple(range(1, n + 1)), 5)
        for n in range(9, 17)
    }
    kappa1_terms = {
        n: templates(n - 5, Fraction(4, 3), tuple(range(1, n - 4)), 4)
        for n in range(9, 17)
    }
    expected_main = {9: 30, 10: 42, 11: 56, 12: 77, 13: 101, 14: 135, 15: 176, 16: 231}
    expected_kappa1 = {9: 5, 10: 7, 11: 11, 12: 15, 13: 22, 14: 30, 15: 42, 16: 56}
    for n in range(9, 17):
        assert len(main_terms[n]) == expected_main[n]
        assert len(kappa1_terms[n]) == expected_kappa1[n]

    lower_receipt = json.loads(
        Path("research/f2_75_125_k10_nonzero_first_slice.json").read_text(encoding="utf-8")
    )
    mid_receipt = json.loads(
        Path("research/f2_75_125_k11_k12_nonzero_first_slice.json").read_text(encoding="utf-8")
    )
    high_receipt = json.loads(
        Path("research/f2_75_125_k13_k14_nonzero_first_slice.json").read_text(encoding="utf-8")
    )
    assert lower_receipt["profile_count"] == mid_receipt["profile_count"] == high_receipt["profile_count"] == 51_870

    q_values = list(range(3, U2_ORDER_MAX + 1)) + [INFINITY]
    r_values = list(range(1, DEGREE_CAP + 1)) + [INFINITY]
    s_values = list(range(DEGREE_CAP + 1)) + [INFINITY]
    v_values = list(range(U5_DEGREE_CAP + 1)) + [INFINITY]

    profile_count = 0
    prior_band_counts = Counter()
    first_band_counts = Counter()
    prior_union = 0
    new_counts = Counter()
    correction_minima = {15: INFINITY, 16: INFINITY}
    new_by_q = {15: Counter(), 16: Counter()}
    examples: dict[str, list[dict[str, object]]] = {"K15": [], "K16": []}
    signatures = {"K15": set(), "K16": set()}
    survivors_after_k16: list[list[int]] = []

    for u1 in range(5, U1_ORDER_MAX + 1):
        for u2 in q_values:
            for u3 in r_values:
                if u2 == 3 and u3 < 2:
                    continue
                for u4 in s_values:
                    for u5 in v_values:
                        profile_count += 1
                        ORDERS.clear()
                        ORDERS.update({1: u1, 2: u2, 3: u3, 4: u4, 5: u5})
                        ORDERS.update({weight: 0 for weight in range(6, 17)})

                        eliminated = False
                        for n in range(9, 15):
                            correction_order, _ = least_terms(kappa1_terms[n], ORDERS)
                            is_robust, _, _ = robust_pole(main_terms[n], correction_order)
                            if is_robust:
                                prior_band_counts[n] += 1
                                if not eliminated:
                                    first_band_counts[n] += 1
                                eliminated = True

                        if eliminated:
                            prior_union += 1

                        if not eliminated:
                            for n in (15, 16):
                                kappa1_order, _ = least_terms(kappa1_terms[n], ORDERS)
                                if n == 15:
                                    kappa2_order = u5 if u5 < INFINITY else None
                                    kappa3_order = 4  # ord(h^2), with ord(h)=2
                                else:
                                    kappa2_order = 0  # ord(U6) lower bound
                                    kappa3_order = u1 - 2  # ord((2/3)U1/h)
                                corrections = [
                                    order for order in (kappa1_order, kappa2_order, kappa3_order)
                                    if order is not None
                                ]
                                correction_order = min(corrections) if corrections else None
                                correction_minima[n] = min(
                                    correction_minima[n],
                                    correction_order if correction_order is not None else INFINITY,
                                )
                                is_robust, main_order, tied = robust_pole(
                                    main_terms[n], correction_order
                                )
                                if not is_robust:
                                    continue
                                band = f"K{n}"
                                new_counts[n] += 1
                                new_by_q[n][str(u2)] += 1
                                multiplicities, value = tied[0]
                                signature = (u1, u2, u3, u4, u5, label(multiplicities))
                                if signature not in signatures[band] and len(examples[band]) < 20:
                                    signatures[band].add(signature)
                                    examples[band].append(
                                        {
                                            "ord_U1": u1,
                                            "ord_U2": "infinity" if u2 >= INFINITY else u2,
                                            "ord_U3": "infinity" if u3 >= INFINITY else u3,
                                            "ord_U4": "infinity" if u4 >= INFINITY else u4,
                                            "ord_U5_over_t": "infinity" if u5 >= INFINITY else u5,
                                            "band": band,
                                            "main_pole_order": -int(main_order),
                                            "unique_main_term": label(multiplicities),
                                            "coefficient": str(value),
                                            "kappa1_correction_order": kappa1_order,
                                            "full_correction_order": correction_order,
                                            "strictly_below_all_corrections": True,
                                        }
                                    )
                                eliminated = True
                                break
                        if not eliminated:
                            survivors_after_k16.append([u1, u2, u3, u4, u5])
        print(f"completed ord(U1)={u1}; profiles={profile_count}", flush=True)

    assert profile_count == 51_870
    assert prior_band_counts[9] == lower_receipt["K9_robust_unique_negative_pole_profiles"]
    assert prior_band_counts[10] == lower_receipt[
        "K10_robust_unique_negative_pole_profiles_including_K9_overlap"
    ]
    assert prior_band_counts[11] == mid_receipt["K11_robust_unique_poles"]
    assert prior_band_counts[12] == mid_receipt["K12_robust_unique_poles"]
    assert first_band_counts[13] == high_receipt["K13_new_after_K9_K12"]
    assert first_band_counts[14] == high_receipt["K14_new_after_K9_K13"]
    assert prior_union == high_receipt["K9_to_K14_union_profiles"]
    cumulative_union = prior_union + new_counts[15] + new_counts[16]
    assert len(survivors_after_k16) == 51_870 - cumulative_union

    result = {
        "scope": "local K9-K16 order screen at the distinguished double root of the F2 2+1+1 nonzero-U1 deeper branch",
        "field": "exact Q coefficients and integer valuations",
        "profile_count": profile_count,
        "K9_to_K14_union_profiles": prior_union,
        "K15_new_after_K9_K14": new_counts[15],
        "K16_new_after_K9_K15": new_counts[16],
        "K9_to_K16_union_profiles": cumulative_union,
        "union_profile_ratio_not_branch_coverage": str(Fraction(cumulative_union, profile_count)),
        "union_profile_percent_not_branch_coverage": round(100 * cumulative_union / profile_count, 4),
        "K15_kappa1_correction_minimum_on_survivors": correction_minima[15],
        "K16_kappa1_correction_minimum_on_survivors": correction_minima[16],
        "K15_regular_resonance_orders": {"kappa_2_U5": "ord(U5/t)>=0", "kappa_3_h2": 4},
        "K16_regular_resonance_orders": {"kappa_2_U6": 0, "kappa_3_w1_F2_3": "ord(U1)-2>=3"},
        "later_slice_baseline_orders": "ord(U6)=...=ord(U16)=0; only unique main minima supported on U1..U5 and below all resonance corrections are counted",
        "new_profile_counts_by_ord_U2": {"K15": dict(new_by_q[15]), "K16": dict(new_by_q[16])},
        "new_examples": examples,
        "survivor_profile_count_after_K16": len(survivors_after_k16),
        "survivor_profiles_after_K16": survivors_after_k16,
        "status": "local pruning only; does not exclude the full 2+1+1 branch",
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    summary_keys = (
        "profile_count",
        "K9_to_K14_union_profiles",
        "K15_new_after_K9_K14",
        "K16_new_after_K9_K15",
        "K9_to_K16_union_profiles",
        "union_profile_percent_not_branch_coverage",
        "survivor_profile_count_after_K16",
        "status",
    )
    print(json.dumps({key: result[key] for key in summary_keys}, indent=2))
    print(f"receipt={output}")


if __name__ == "__main__":
    main()
