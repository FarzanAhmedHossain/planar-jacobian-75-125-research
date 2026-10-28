"""Exact K9--K12 valuation screen on the nonzero-U1 2+1+1 double-root branch.

K11 and K12 contain the first-integral corrections
  G_11 = [w^11]F^(5/3) + kappa_1/t [w^6]F^(4/3) + kappa_2/t^2 U1,
  G_12 = [w^12]F^(5/3) + kappa_1/t [w^7]F^(4/3) + kappa_2/t^2 U2.
The script counts only unique negative minima of the binomial term supported
on U1,...,U5 and strictly below the entire kappa_1 correction. The kappa_2
terms are regular at the double root. Later coefficients use order zero as a
worst-case baseline, so these certificates survive any increase of their
orders. This is local profile pruning, not whole-branch coverage.
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


def label(multiplicities: tuple[int, ...]) -> str:
    return "+".join(
        str(weight)
        for weight, multiplicity in enumerate(multiplicities, start=1)
        for _ in range(multiplicity)
    ) or "1"


def robust_main_pole(main_terms, orders: dict[int, int], correction_min: int | None):
    main_order, tied = least_terms(main_terms, orders)
    if main_order is None or main_order >= 0 or len(tied) != 1:
        return False, main_order, tied
    multiplicities, _ = tied[0]
    if any(multiplicities[SUPPORT_CUTOFF:]):
        return False, main_order, tied
    if correction_min is not None and main_order >= correction_min:
        return False, main_order, tied
    return True, main_order, tied


def main() -> None:
    main_terms = {
        n: templates(n, Fraction(5, 3), tuple(range(1, n + 1)), 5)
        for n in range(9, 13)
    }
    corrections = {
        n: templates(n - 5, Fraction(4, 3), tuple(range(1, n - 4)), 4)
        for n in range(9, 13)
    }
    expected_partition_counts = {9: 30, 10: 42, 11: 56, 12: 77}
    for n in range(9, 13):
        assert len(main_terms[n]) == expected_partition_counts[n]
    assert [len(corrections[n]) for n in range(9, 13)] == [5, 7, 11, 15]

    q_values = list(range(3, U2_ORDER_MAX + 1)) + [INFINITY]
    r_values = list(range(1, DEGREE_CAP + 1)) + [INFINITY]
    s_values = list(range(DEGREE_CAP + 1)) + [INFINITY]
    v_values = list(range(U5_DEGREE_CAP + 1)) + [INFINITY]

    profile_count = 0
    band_counts = Counter()
    band_new_counts = Counter()
    correction_minima = {n: INFINITY for n in range(9, 13)}
    new_by_q = {n: Counter() for n in range(11, 13)}
    examples: dict[str, list[dict[str, object]]] = {"K11": [], "K12": []}
    example_signatures = {"K11": set(), "K12": set()}

    for u1 in range(5, U1_ORDER_MAX + 1):
        for u2 in q_values:
            for u3 in r_values:
                if u2 == 3 and u3 < 2:
                    continue
                for u4 in s_values:
                    for u5 in v_values:
                        profile_count += 1
                        orders = {1: u1, 2: u2, 3: u3, 4: u4, 5: u5}
                        orders.update({weight: 0 for weight in range(6, 13)})

                        correction_orders: dict[int, int | None] = {}
                        for n in range(9, 13):
                            correction_order, _ = least_terms(corrections[n], orders)
                            correction_orders[n] = correction_order
                            if correction_order is not None:
                                correction_minima[n] = min(
                                    correction_minima[n], correction_order
                                )

                        eliminated_before = False
                        for n in range(9, 13):
                            is_robust, main_order, tied = robust_main_pole(
                                main_terms[n], orders, correction_orders[n]
                            )
                            if not is_robust:
                                continue
                            band = f"K{n}"
                            band_counts[band] += 1
                            if eliminated_before:
                                continue
                            if n in (11, 12):
                                band_new_counts[band] += 1
                                new_by_q[n][str(u2)] += 1
                                multiplicities, value = tied[0]
                                signature = (u1, u2, u3, u4, u5, label(multiplicities))
                                if signature not in example_signatures[band] and len(examples[band]) < 24:
                                    example_signatures[band].add(signature)
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
                                            "kappa1_correction_order": correction_orders[n],
                                            "strictly_below_correction": True,
                                        }
                                    )
                            eliminated_before = True

    assert profile_count == 51_870, profile_count
    # Reconstruct K9/K10 robust counts from the recorded exact census so the
    # cumulative union uses the already independently replayed lower screens.
    previous = json.loads(
        Path("research/f2_75_125_k10_nonzero_first_slice.json").read_text(encoding="utf-8")
    )
    assert previous["profile_count"] == profile_count
    k9_count = previous["K9_robust_unique_negative_pole_profiles"]
    k10_new = previous["K10_new_profiles_after_K9"]
    assert band_counts["K9"] == k9_count, (band_counts["K9"], k9_count)
    assert band_counts["K10"] == previous[
        "K10_robust_unique_negative_pole_profiles_including_K9_overlap"
    ]
    k11_new = band_new_counts["K11"]
    k12_new = band_new_counts["K12"]
    union = k9_count + k10_new + k11_new + k12_new

    result = {
        "scope": "local K9-K12 order screen at the distinguished double root of the F2 2+1+1 nonzero-U1 deeper branch",
        "field": "exact Q coefficients and integer valuations",
        "profile_count": profile_count,
        "K9_robust_unique_negative_pole_profiles": k9_count,
        "K10_new_profiles_after_K9": k10_new,
        "K11_robust_unique_poles": band_counts["K11"],
        "K11_new_profiles_after_K9_K10": k11_new,
        "K12_robust_unique_poles": band_counts["K12"],
        "K12_new_profiles_after_K9_K10_K11": k12_new,
        "K9_to_K12_union_profiles": union,
        "union_profile_ratio_not_branch_coverage": str(Fraction(union, profile_count)),
        "union_profile_percent_not_branch_coverage": round(100 * union / profile_count, 4),
        "kappa1_correction_minimum_orders": {f"K{n}": correction_minima[n] for n in range(9, 13)},
        "regular_kappa2_terms": {"K11": "kappa_2*U1, order >=5", "K12": "kappa_2*U2, order >=3"},
        "later_slice_baseline_orders": "ord(U6)=...=ord(U12)=0; robust terms must be unique and supported on U1..U5",
        "new_profile_counts_by_ord_U2": {f"K{n}": dict(new_by_q[n]) for n in (11, 12)},
        "new_examples": examples,
        "status": "local pruning only; does not exclude the full 2+1+1 branch",
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"receipt={output}")


if __name__ == "__main__":
    main()
