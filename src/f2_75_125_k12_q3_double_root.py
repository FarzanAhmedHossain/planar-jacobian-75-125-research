"""Exact K12 pruning on the double-root branch.

K9--K11 leave only (r,v)=(5,2) or r>=6,v>=3 (allowing zero slices)
when q=ord(U2)=3 at a double root and U1=0. This checker enumerates all
weighted partitions of 12, verifies the unique order -8 residue there, and
reduces it modulo the K10 quadratic. On q>=5, it also shows that s=0,1
are impossible; q=4 already has s>=2 from K10. It does not settle q>=4
or the pair.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import factorial


WEIGHTS_12 = tuple(range(2, 13))
WEIGHTS_7 = tuple(range(2, 8))
INFINITY = 10_000


def partitions(total: int, weights: tuple[int, ...], index: int = 0):
    if index == len(weights):
        if total == 0:
            yield ()
        return
    weight = weights[index]
    for multiplicity in range(total // weight + 1):
        for tail in partitions(total - multiplicity * weight, weights, index + 1):
            yield (multiplicity,) + tail


def binomial(alpha: Fraction, n: int) -> Fraction:
    value = Fraction(1)
    for j in range(n):
        value *= (alpha - j) / (j + 1)
    return value


def partition_coefficient(alpha: Fraction, multiplicities: tuple[int, ...]) -> Fraction:
    count = sum(multiplicities)
    value = binomial(alpha, count) * factorial(count)
    for multiplicity in multiplicities:
        value /= factorial(multiplicity)
    return value


def label(multiplicities: tuple[int, ...], weights: tuple[int, ...]) -> str:
    return "+".join(
        str(weight)
        for weight, multiplicity in zip(weights, multiplicities)
        for _ in range(multiplicity)
    )


def term_order(
    multiplicities: tuple[int, ...],
    weights: tuple[int, ...],
    orders: dict[int, int],
    h_order: int,
    exponent_numerator: int,
) -> int | None:
    if any(
        multiplicity and orders.get(weight, INFINITY) >= INFINITY
        for weight, multiplicity in zip(weights, multiplicities)
    ):
        return None
    count = sum(multiplicities)
    return h_order * (exponent_numerator - 3 * count) + sum(
        multiplicity * orders[weight]
        for weight, multiplicity in zip(weights, multiplicities)
        if multiplicity
    )


def least_terms(
    total: int,
    weights: tuple[int, ...],
    orders: dict[int, int],
    h_order: int,
    exponent_numerator: int,
    alpha: Fraction,
) -> tuple[int | None, list[tuple[str, Fraction]]]:
    evaluated = []
    for multiplicities in partitions(total, weights):
        order = term_order(
            multiplicities, weights, orders, h_order, exponent_numerator
        )
        if order is not None:
            evaluated.append(
                (
                    order,
                    label(multiplicities, weights),
                    partition_coefficient(alpha, multiplicities),
                )
            )
    if not evaluated:
        return None, []
    minimum = min(item[0] for item in evaluated)
    return minimum, [
        (term_label, coefficient)
        for order, term_label, coefficient in evaluated
        if order == minimum
    ]


def main() -> None:
    degree12_partitions = list(partitions(12, WEIGHTS_12))
    assert len(degree12_partitions) == 21

    expected_coefficients = {
        "4+4+4": Fraction(-5, 81),
        "2+2+4+4": Fraction(10, 81),
        "2+2+2+2+4": Fraction(-35, 729),
        "2+2+2+2+2+2": Fraction(35, 6561),
    }
    actual_coefficients = {
        label(multiplicities, WEIGHTS_12): partition_coefficient(
            Fraction(5, 3), multiplicities
        )
        for multiplicities in degree12_partitions
    }
    for term, coefficient in expected_coefficients.items():
        assert actual_coefficients[term] == coefficient

    checked = 0
    case_summary: dict[str, int] = {}
    r_values = list(range(6, 12)) + [INFINITY]
    v_values = list(range(3, 13)) + [INFINITY]
    # The q=3 survivors are r=5,v=2 or r>=6,v>=3. Later F-slices take
    # their worst-case local order zero; any higher actual order only removes
    # terms or raises their order.
    strata = [(5, 2)] + [
        (r, v) for r in r_values for v in v_values
    ]
    expected_labels = sorted(expected_coefficients)
    for r, v in strata:
        orders = {2: 3, 3: r, 4: 0, 5: v}
        orders.update({weight: 0 for weight in range(6, 13)})
        minimum, terms = least_terms(
            12,
            WEIGHTS_12,
            orders,
            h_order=2,
            exponent_numerator=5,
            alpha=Fraction(5, 3),
        )
        assert minimum == -8, (r, v, minimum, terms)
        assert sorted(term for term, _ in terms) == expected_labels

        # At K12, the kappa_1 mode is [w^7]F^(4/3), while the kappa_2 mode
        # is [w^2]F=U2. Neither has a pole as deep as -8.
        correction_orders = {2: 3, 3: r, 4: 0, 5: v}
        correction_orders.update({weight: 0 for weight in range(6, 8)})
        correction_minimum, _ = least_terms(
            7,
            WEIGHTS_7,
            correction_orders,
            h_order=2,
            exponent_numerator=4,
            alpha=Fraction(4, 3),
        )
        assert correction_minimum >= 1
        assert 3 > -8  # order of the kappa_2 U2 term
        pattern = f"r={'infinity' if r == INFINITY else r},v={'infinity' if v == INFINITY else v}"
        case_summary[pattern] = minimum
        checked += 1

    # In the q>=4 branch, K10 gives v>=1 and q=4 => s>=2; K9 gives
    # s+v>=2. If q>=5, only s=0,1 remain to test. Their unique K12
    # lowest term is 4+4+4, with order -8 or -5 respectively.
    qge5_checked = 0
    qge5_by_s = {"s=0": 0, "s=1": 0}
    for q in list(range(5, 12)) + [INFINITY]:
        for r in list(range(3, 12)) + [INFINITY]:
            for s in (0, 1):
                v_minimum = 2 if s == 0 else 1
                for v in list(range(v_minimum, 13)) + [INFINITY]:
                    orders = {2: q, 3: r, 4: s, 5: v}
                    orders.update({weight: 0 for weight in range(6, 13)})
                    minimum, terms = least_terms(
                        12,
                        WEIGHTS_12,
                        orders,
                        h_order=2,
                        exponent_numerator=5,
                        alpha=Fraction(5, 3),
                    )
                    expected_order = -8 if s == 0 else -5
                    assert minimum == expected_order
                    assert terms == [("4+4+4", Fraction(-5, 81))]

                    correction_orders = {2: q, 3: r, 4: s, 5: v}
                    correction_orders.update({weight: 0 for weight in range(6, 8)})
                    correction_minimum, _ = least_terms(
                        7,
                        WEIGHTS_7,
                        correction_orders,
                        h_order=2,
                        exponent_numerator=4,
                        alpha=Fraction(4, 3),
                    )
                    assert correction_minimum >= (-1 if s == 0 else 0)
                    assert q > minimum  # kappa_2 U2 is also higher order
                    qge5_by_s[f"s={s}"] += 1
                    qge5_checked += 1

    # With h=c*u^2, U2=d*u^3, U4=a, and z=c^3*a/d^2, the order -8
    # coefficient is (a^3/c^4)*C(z), where C is the displayed Laurent
    # polynomial. Clear 6561*z^3 and reduce its numerator by the K10 relation.
    numerator = {
        3: Fraction(-405),
        2: Fraction(810),
        1: Fraction(-315),
        0: Fraction(35),
    }
    # More explicitly, -405*z^3+810*z^2-315*z+35 reduces to 7/3-14*z.
    z3_reduction = {1: Fraction(59, 405), 0: Fraction(-28, 1215)}
    z2_reduction = {1: Fraction(4, 9), 0: Fraction(-7, 135)}
    reduced = {
        1: numerator[3] * z3_reduction[1]
        + Fraction(810) * z2_reduction[1]
        + numerator[1],
        0: numerator[3] * z3_reduction[0]
        + Fraction(810) * z2_reduction[0]
        + numerator[0],
    }
    assert reduced == {1: Fraction(-14), 0: Fraction(7, 3)}
    # z=1/6 is not a root of 135*z^2-60*z+7.
    k10_at_sixth = Fraction(135, 36) - Fraction(60, 6) + 7
    assert k10_at_sixth == Fraction(3, 4)

    print(
        json.dumps(
            {
                "degree12_weighted_partition_count": len(degree12_partitions),
                "surviving_q3_order_strata_checked": checked,
                "qge5_small_s_strata_checked": qge5_checked,
                "qge5_small_s_strata_by_s": qge5_by_s,
                "K12_unique_pole_order": -8,
                "K12_lowest_terms": expected_labels,
                "K12_normalized_residue": "(-5/81 + 10/(81*z) - 35/(729*z^2) + 35/(6561*z^3))",
                "cleared_numerator_mod_K10": "(7/3)*(1-6*z)",
                "K10_at_z_equals_1_over_6": str(k10_at_sixth),
                "correction_orders": {
                    "q3_kappa1_mode": "at least 1 on all checked strata",
                    "q3_kappa2_mode": "3",
                    "qge5_s0_kappa1_mode": "at least -1",
                    "qge5_s1_kappa1_mode": "at least 0",
                    "qge5_kappa2_mode": "at least 5",
                },
                "qge5_consequence": "s=0,1 have a unique uncancellable K12 pole; hence s>=2",
                "all_arithmetic_exact": True,
                "scope": "excludes q=3 and forces s>=2 for q>=5 at a double root on U1=0; q=4 and higher-q cases remain",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
