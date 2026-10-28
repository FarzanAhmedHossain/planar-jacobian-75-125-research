"""Exact K9/K10 valuation screen for the nonzero-U1 2+1+1 stratum.

This enumerates supported root-order profiles at the distinguished double root.
K9 and K10 coefficients use exact rational generalized-binomial coefficients.
For K10, the kappa_1 correction has order at least one and the kappa_2 h^3
term has order six, so a unique negative pole in [w^10]F^(5/3) cannot cancel.
Only unique least poles supported on U1,...,U5 are counted as robust against
all allowed orders of later slices. This is local pruning, not whole-branch
coverage.
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


@lru_cache(maxsize=None)
def partitions(total: int, weights: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    if not weights:
        return ((),) if total == 0 else ()
    weight, *tail = weights
    rest = tuple(tail)
    result: list[tuple[int, ...]] = []
    for multiplicity in range(total // weight + 1):
        for suffix in partitions(total - multiplicity * weight, rest):
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


def least_terms(terms, orders: dict[int, int]):
    evaluated = []
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


def main() -> None:
    terms9 = templates(9, Fraction(5, 3), tuple(range(1, 10)), 5)
    terms10 = templates(10, Fraction(5, 3), tuple(range(1, 11)), 5)
    correction9 = templates(4, Fraction(4, 3), tuple(range(1, 5)), 4)
    correction10 = templates(5, Fraction(4, 3), tuple(range(1, 6)), 4)
    assert len(terms9) == 30
    assert len(terms10) == 42
    assert len(correction9) == 5
    assert len(correction10) == 7

    q_values = list(range(3, U2_ORDER_MAX + 1)) + [INFINITY]
    r_values = list(range(1, DEGREE_CAP + 1)) + [INFINITY]
    s_values = list(range(DEGREE_CAP + 1)) + [INFINITY]
    v_values = list(range(U5_DEGREE_CAP + 1)) + [INFINITY]

    profile_count = 0
    k9_robust = 0
    k10_robust = 0
    k10_new = 0
    correction9_min = INFINITY
    correction10_min = INFINITY
    examples: list[dict[str, object]] = []
    by_signature: set[tuple[int, int, int, int, int, str]] = set()
    new_by_q: Counter[str] = Counter()

    for u1 in range(5, U1_ORDER_MAX + 1):
        for u2 in q_values:
            for u3 in r_values:
                # K8 excludes q=3,r=1 on the deeper double-root stratum.
                if u2 == 3 and u3 < 2:
                    continue
                for u4 in s_values:
                    for u5 in v_values:
                        profile_count += 1
                        orders = {1: u1, 2: u2, 3: u3, 4: u4, 5: u5}
                        # Later slices take the least restrictive finite order.
                        orders.update({weight: 0 for weight in range(6, 11)})

                        correction_order9, _ = least_terms(correction9, orders)
                        correction_order10, _ = least_terms(correction10, orders)
                        if correction_order9 is not None:
                            correction9_min = min(correction9_min, correction_order9)
                            assert correction_order9 >= 2, (orders, correction_order9)
                        if correction_order10 is not None:
                            correction10_min = min(correction10_min, correction_order10)
                            assert correction_order10 >= 1, (orders, correction_order10)

                        k9_order, k9_terms = least_terms(terms9, orders)
                        k9_bad = bool(
                            k9_order is not None
                            and k9_order < 0
                            and len(k9_terms) == 1
                            and all(m == 0 for m in k9_terms[0][0][5:])
                        )
                        if k9_bad:
                            k9_robust += 1

                        k10_order, k10_terms = least_terms(terms10, orders)
                        k10_bad = bool(
                            k10_order is not None
                            and k10_order < 0
                            and len(k10_terms) == 1
                            and all(m == 0 for m in k10_terms[0][0][5:])
                        )
                        if k10_bad:
                            k10_robust += 1
                            if not k9_bad:
                                k10_new += 1
                                new_by_q[str(u2)] += 1
                                multiplicities, value = k10_terms[0]
                                signature = (u1, u2, u3, u4, u5, label(multiplicities))
                                if signature not in by_signature and len(examples) < 24:
                                    by_signature.add(signature)
                                    examples.append(
                                        {
                                            "ord_U1": u1,
                                            "ord_U2": "infinity" if u2 >= INFINITY else u2,
                                            "ord_U3": "infinity" if u3 >= INFINITY else u3,
                                            "ord_U4": "infinity" if u4 >= INFINITY else u4,
                                            "ord_U5_over_t": "infinity" if u5 >= INFINITY else u5,
                                            "pole_order": -k10_order,
                                            "unique_term": label(multiplicities),
                                            "coefficient": str(value),
                                            "k9_robust_obstruction": False,
                                        }
                                    )

    assert profile_count == 51_870, profile_count
    union_count = k9_robust + k10_new
    result = {
        "scope": "local K9/K10 order screen at the distinguished double root of the F2 2+1+1 nonzero-U1 deeper branch",
        "field": "exact Q coefficients and integer valuations",
        "profile_count": profile_count,
        "K9_robust_unique_negative_pole_profiles": k9_robust,
        "K10_robust_unique_negative_pole_profiles_including_K9_overlap": k10_robust,
        "K10_new_profiles_after_K9": k10_new,
        "K9_or_K10_robust_profiles": union_count,
        "union_profile_ratio_not_branch_coverage": str(Fraction(union_count, profile_count)),
        "union_profile_percent_not_branch_coverage": round(100 * union_count / profile_count, 4),
        "K9_kappa1_correction_minimum_order": correction9_min,
        "K10_kappa1_correction_minimum_order": correction10_min,
        "K10_kappa2_h3_order": 6,
        "later_slice_baseline_orders": "ord(U6)=...=ord(U10)=0; count only unique minima supported on U1..U5",
        "K10_new_profile_counts_by_ord_U2": dict(new_by_q),
        "K10_new_obstruction_examples": examples,
        "status": "local pruning only; does not exclude the full 2+1+1 branch",
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"receipt={output}")


if __name__ == "__main__":
    main()
