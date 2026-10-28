"""Exact K9 valuation screen for the nonzero-U1 deeper double-root stratum.

Scope: the distinguished double root in the F2 2+1+1 branch. The existing
K7/K8 conditions give ord(U1)>=5, ord(U2)>=3, ord(U3)>=1, and if
ord(U2)=3 then ord(U3)>=2. We enumerate all supported orders through degree
11 for U1,...,U4 and through degree 12 for U5. In this fixed 2+1+1 chart,
U1=(t-1)^5 R^2 B with deg(B)<=2, so a nonzero U1 has root order 5, 6,
or 7; U2=(t-1)^2 R C with deg(C)<=7, so its surviving root order is
3..9 or infinity. For U6,...,U9 we use order zero as a worst-case baseline;
we certify an obstruction only when the unique least pole uses no U6,...,U9,
so increasing any of those orders cannot remove it.

This is a local necessary-condition reduction, not a whole-branch exclusion.
The K9 first-integral correction is separately checked to have order >=2 on
every enumerated stratum, so negative poles come from [w^9]F^(5/3).
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import factorial
import json
from pathlib import Path


DEGREE_CAP = 11
U5_DEGREE_CAP = 12
U1_ORDER_MAX = 7
U2_ORDER_MAX = 9
INFINITY = U5_DEGREE_CAP + 1
H_ORDER = 2
WEIGHTS_9 = tuple(range(1, 10))
WEIGHTS_4 = tuple(range(1, 5))


@lru_cache(maxsize=None)
def partitions(total: int, weights: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    """Return multiplicity vectors for weighted partitions of ``total``."""
    if not weights:
        return ((),) if total == 0 else ()
    weight, *tail = weights
    tail_tuple = tuple(tail)
    result: list[tuple[int, ...]] = []
    for multiplicity in range(total // weight + 1):
        remainder = total - multiplicity * weight
        for suffix in partitions(remainder, tail_tuple):
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


def term_order(
    multiplicities: tuple[int, ...],
    weights: tuple[int, ...],
    orders: dict[int, int],
    exponent_base: int,
) -> int | None:
    if any(
        multiplicity and orders[weight] >= INFINITY
        for weight, multiplicity in zip(weights, multiplicities)
    ):
        return None
    count = sum(multiplicities)
    return H_ORDER * (exponent_base - 3 * count) + sum(
        multiplicity * orders[weight]
        for weight, multiplicity in zip(weights, multiplicities)
    )


def least_terms(
    total: int,
    alpha: Fraction,
    weights: tuple[int, ...],
    orders: dict[int, int],
    exponent_base: int,
) -> tuple[int | None, list[tuple[tuple[int, ...], Fraction]]]:
    evaluated: list[tuple[int, tuple[int, ...], Fraction]] = []
    for multiplicities in partitions(total, weights):
        order = term_order(multiplicities, weights, orders, exponent_base)
        if order is not None:
            evaluated.append((order, multiplicities, coefficient(alpha, multiplicities)))
    if not evaluated:
        return None, []
    minimum = min(item[0] for item in evaluated)
    return minimum, [
        (multiplicities, value)
        for order, multiplicities, value in evaluated
        if order == minimum
    ]


def term_templates(
    total: int, alpha: Fraction, weights: tuple[int, ...], exponent_base: int
) -> tuple[tuple[tuple[int, ...], int, Fraction, int], ...]:
    """Precompute partition data so the large order census stays lightweight."""
    return tuple(
        (m, sum(m), coefficient(alpha, m), H_ORDER * (exponent_base - 3 * sum(m)))
        for m in partitions(total, weights)
    )


def label(multiplicities: tuple[int, ...], weights: tuple[int, ...]) -> str:
    return "+".join(
        str(weight)
        for weight, multiplicity in zip(weights, multiplicities)
        for _ in range(multiplicity)
    ) or "1"


def main() -> None:
    degree9 = partitions(9, WEIGHTS_9)
    degree4 = partitions(4, WEIGHTS_4)
    terms9 = term_templates(9, Fraction(5, 3), WEIGHTS_9, 5)
    terms4 = term_templates(4, Fraction(4, 3), WEIGHTS_4, 4)
    assert len(degree9) == 30
    assert len(degree4) == 5

    checked = 0
    robust_obstructions = 0
    profile_examples: list[dict[str, object]] = []
    seen_example_signatures: set[tuple[int, int, str]] = set()
    obstruction_by_qr: Counter[tuple[int, int]] = Counter()
    resonance_minimum = INFINITY

    q_values = list(range(3, U2_ORDER_MAX + 1)) + [INFINITY]
    r_values = list(range(1, DEGREE_CAP + 1)) + [INFINITY]
    s_values = list(range(DEGREE_CAP + 1)) + [INFINITY]
    v_values = list(range(U5_DEGREE_CAP + 1)) + [INFINITY]

    for u1_order in range(5, U1_ORDER_MAX + 1):
        for u2_order in q_values:
            for u3_order in r_values:
                # K8 excludes q=3,r=1 on this deeper double-root branch.
                if u2_order == 3 and u3_order < 2:
                    continue
                for u4_order in s_values:
                    for u5_order in v_values:
                        orders = {
                            1: u1_order,
                            2: u2_order,
                            3: u3_order,
                            4: u4_order,
                            5: u5_order,
                            # Baseline orders. Certified obstructions below
                            # are required to be independent of these entries.
                            6: 0,
                            7: 0,
                            8: 0,
                            9: 0,
                        }
                        checked += 1

                        # The K9 correction is (kappa_1/t)[w^4]F^(4/3).
                        # t is a unit at the root; this verifies its order >=2.
                        correction_orders = [
                            base + sum(m[i - 1] * orders[i] for i in range(1, 5))
                            for m, _, _, base in terms4
                            if not any(m[i - 1] and orders[i] >= INFINITY for i in range(1, 5))
                        ]
                        correction_order = min(correction_orders) if correction_orders else None
                        if correction_order is not None:
                            resonance_minimum = min(resonance_minimum, correction_order)
                            assert correction_order >= 2, (orders, correction_order)

                        evaluated = [
                            (
                                base + sum(m[i - 1] * orders[i] for i in range(1, 10)),
                                m,
                                value,
                            )
                            for m, _, value, base in terms9
                            if not any(m[i - 1] and orders[i] >= INFINITY for i in range(1, 10))
                        ]
                        if not evaluated:
                            continue
                        minimum = min(order for order, _, _ in evaluated)
                        if minimum >= 0:
                            continue
                        lowest = [
                            (m, value)
                            for order, m, value in evaluated
                            if order == minimum
                        ]
                        if len(lowest) != 1:
                            continue
                        multiplicities, value = lowest[0]
                        support = [
                            weight
                            for weight, multiplicity in zip(WEIGHTS_9, multiplicities)
                            if multiplicity
                        ]
                        # If the unique least term uses only U1,...,U5, every
                        # term involving U6,...,U9 starts strictly higher at
                        # order zero and can only rise if those orders increase.
                        if max(support) > 5:
                            continue
                        assert value != 0
                        robust_obstructions += 1
                        obstruction_by_qr[(u2_order, u3_order)] += 1
                        term_label = label(multiplicities, WEIGHTS_9)
                        signature = (u2_order, u3_order, term_label)
                        if signature not in seen_example_signatures and len(profile_examples) < 24:
                            seen_example_signatures.add(signature)
                            profile_examples.append(
                                {
                                    "ord_U1": u1_order,
                                    "ord_U2": "infinity" if u2_order >= INFINITY else u2_order,
                                    "ord_U3": "infinity" if u3_order >= INFINITY else u3_order,
                                    "ord_U4": "infinity" if u4_order >= INFINITY else u4_order,
                                    "ord_U5_over_t": (
                                        "infinity" if u5_order >= INFINITY else u5_order
                                    ),
                                    "pole_order": -minimum,
                                    "unique_term": term_label,
                                    "coefficient": str(value),
                                }
                            )

    total_profiles = checked
    result = {
        "scope": "local K9 order screen at the distinguished double root of F2 2+1+1, nonzero U1 deeper branch only",
        "field": "Q exact rational binomial coefficients and integer valuations",
        "degree_cap_for_root_orders": DEGREE_CAP,
        "degree9_weighted_partition_count": len(degree9),
        "degree4_weighted_partition_count": len(degree4),
        "profile_count": total_profiles,
        "profiles_with_robust_unique_negative_pole": robust_obstructions,
        "raw_profile_ratio_not_branch_coverage": f"{Fraction(robust_obstructions, total_profiles)}",
        "raw_profile_ratio_percent_not_branch_coverage": round(
            100 * robust_obstructions / total_profiles, 4
        ),
        "kappa1_correction_minimum_order": resonance_minimum,
        "baseline_orders_for_U6_to_U9": 0,
        "robustness_rule": "count only a unique least negative term whose support uses U1..U5; increasing ord(U6..U9) cannot remove it",
        "order_bounds": {
            "U1": "5..7 from U1=(t-1)^5 R^2 B, deg(B)<=2; U1=0 handled separately",
            "U2": "3..9 from U2=(t-1)^2 R C, deg(C)<=7, or identically zero",
            "U3": "1..11 or identically zero; if ord(U2)=3 then ord(U3)>=2",
            "U4": "0..11 or identically zero",
            "U5_over_t": "0..12 or identically zero",
        },
        "unique_obstruction_examples": profile_examples,
        "status": "local pruning only; does not exclude the whole 2+1+1 branch",
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"receipt={output}")


if __name__ == "__main__":
    main()
