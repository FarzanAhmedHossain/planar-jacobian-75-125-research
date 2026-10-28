"""Exact K13/K14 valuation screen on the nonzero-U1 2+1+1 branch.

At bands 13 and 14 the first-integral identities are
  G_13 = [w^13]F^(5/3) + kappa_1/t [w^8]F^(4/3) + kappa_2/t^2 U3,
  G_14 = [w^14]F^(5/3) + kappa_1/t [w^9]F^(4/3) + kappa_2/t^2 U4.
The kappa_2 terms are regular. Count only a unique negative main pole
supported on U1,...,U5 and strictly below the full kappa_1 correction order.
Later slices take order zero as a worst-case baseline. The result is local
profile pruning, not whole-branch coverage.
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


def robust_pole(main_terms, correction_terms, orders: dict[int, int]):
    main_order, tied = least_terms(main_terms, orders)
    correction_order, _ = least_terms(correction_terms, orders)
    if main_order is None or main_order >= 0 or len(tied) != 1:
        return False, main_order, correction_order, tied
    multiplicities, _ = tied[0]
    if any(multiplicities[SUPPORT_CUTOFF:]):
        return False, main_order, correction_order, tied
    if correction_order is not None and main_order >= correction_order:
        return False, main_order, correction_order, tied
    return True, main_order, correction_order, tied


def label(multiplicities: tuple[int, ...]) -> str:
    return "+".join(
        str(weight)
        for weight, multiplicity in enumerate(multiplicities, start=1)
        for _ in range(multiplicity)
    ) or "1"


def main() -> None:
    main_terms = {
        n: templates(n, Fraction(5, 3), tuple(range(1, n + 1)), 5)
        for n in range(9, 15)
    }
    correction_terms = {
        n: templates(n - 5, Fraction(4, 3), tuple(range(1, n - 4)), 4)
        for n in range(9, 15)
    }
    expected_main_counts = {9: 30, 10: 42, 11: 56, 12: 77, 13: 101, 14: 135}
    expected_correction_counts = {9: 5, 10: 7, 11: 11, 12: 15, 13: 22, 14: 30}
    for n in range(9, 15):
        assert len(main_terms[n]) == expected_main_counts[n]
        assert len(correction_terms[n]) == expected_correction_counts[n]

    q_values = list(range(3, U2_ORDER_MAX + 1)) + [INFINITY]
    r_values = list(range(1, DEGREE_CAP + 1)) + [INFINITY]
    s_values = list(range(DEGREE_CAP + 1)) + [INFINITY]
    v_values = list(range(U5_DEGREE_CAP + 1)) + [INFINITY]

    profile_count = 0
    new_counts = Counter()
    prior_counts = Counter()
    prior_union_profiles = 0
    correction_minima = {13: INFINITY, 14: INFINITY}
    new_by_q = {13: Counter(), 14: Counter()}
    examples: dict[str, list[dict[str, object]]] = {"K13": [], "K14": []}
    signatures = {"K13": set(), "K14": set()}

    for u1 in range(5, U1_ORDER_MAX + 1):
        for u2 in q_values:
            for u3 in r_values:
                if u2 == 3 and u3 < 2:
                    continue
                for u4 in s_values:
                    for u5 in v_values:
                        profile_count += 1
                        orders = {1: u1, 2: u2, 3: u3, 4: u4, 5: u5}
                        orders.update({weight: 0 for weight in range(6, 15)})

                        eliminated_before = False
                        for n in range(9, 13):
                            is_robust, _, _, _ = robust_pole(
                                main_terms[n], correction_terms[n], orders
                            )
                            if is_robust:
                                prior_counts[n] += 1
                                eliminated_before = True
                        if eliminated_before:
                            prior_union_profiles += 1

                        if profile_count % 10_000 == 0:
                            print(
                                json.dumps(
                                    {
                                        "stage": "profile-progress",
                                        "profiles": profile_count,
                                        "prior_union": prior_union_profiles,
                                        "prior_K11": prior_counts[11],
                                        "prior_K12": prior_counts[12],
                                        "K13_new": new_counts["K13"],
                                        "K14_new": new_counts["K14"],
                                    }
                                ),
                                flush=True,
                            )

                        if eliminated_before:
                            continue

                        for n in (13, 14):
                            band = f"K{n}"
                            is_robust, main_order, correction_order, tied = robust_pole(
                                main_terms[n], correction_terms[n], orders
                            )
                            if correction_order is not None:
                                correction_minima[n] = min(
                                    correction_minima[n], correction_order
                                )
                            if not is_robust:
                                continue
                            new_counts[band] += 1
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
                                        "kappa1_correction_order": correction_order,
                                        "strictly_below_correction": True,
                                    }
                                )
                            if n == 13:
                                break

    assert profile_count == 51_870, profile_count
    previous = json.loads(
        Path("research/f2_75_125_k11_k12_nonzero_first_slice.json").read_text(encoding="utf-8")
    )
    assert previous["profile_count"] == profile_count
    k9_to_k12 = previous["K9_to_K12_union_profiles"]
    assert prior_union_profiles == k9_to_k12, (prior_union_profiles, k9_to_k12)
    assert prior_counts[11] == previous["K11_robust_unique_poles"], (
        prior_counts[11], previous["K11_robust_unique_poles"]
    )
    assert prior_counts[12] == previous["K12_robust_unique_poles"], (
        prior_counts[12], previous["K12_robust_unique_poles"]
    )
    union = k9_to_k12 + new_counts["K13"] + new_counts["K14"]

    result = {
        "scope": "local K9-K14 order screen at the distinguished double root of the F2 2+1+1 nonzero-U1 deeper branch",
        "field": "exact Q coefficients and integer valuations",
        "profile_count": profile_count,
        "K9_to_K12_union_profiles": k9_to_k12,
        "K13_new_after_K9_K12": new_counts["K13"],
        "K14_new_after_K9_K13": new_counts["K14"],
        "K9_to_K14_union_profiles": union,
        "union_profile_ratio_not_branch_coverage": str(Fraction(union, profile_count)),
        "union_profile_percent_not_branch_coverage": round(100 * union / profile_count, 4),
        "kappa1_correction_minimum_orders_on_surviving_profiles": {
            "K13": correction_minima[13], "K14": correction_minima[14]
        },
        "regular_kappa2_terms": {
            "K13": "kappa_2*U3, order >=1",
            "K14": "kappa_2*U4, order >=0",
        },
        "later_slice_baseline_orders": "ord(U6)=...=ord(U14)=0; unique main minima must use U1..U5 and lie strictly below the full kappa_1 correction",
        "new_profile_counts_by_ord_U2": {
            "K13": dict(new_by_q[13]), "K14": dict(new_by_q[14])
        },
        "new_examples": examples,
        "status": "local pruning only; does not exclude the full 2+1+1 branch",
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"receipt={output}")


if __name__ == "__main__":
    main()
