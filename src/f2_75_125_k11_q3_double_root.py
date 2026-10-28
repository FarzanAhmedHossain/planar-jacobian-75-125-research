"""Exact K11 pruning on the q=3 double-root branch.

On the U1=0 branch at a double root of h, K9 leaves (q,r,v)=(3,3,0)
and (3,4,1) subject to 4*d*e=9*c^3*g. K10 pins
z=c^3*a/d^2 by 135*z^2-60*z+7=0. This checker enumerates every weighted
partition in the K11 binomial coefficient. It excludes r=3,4; forces v=2
when r=5 (with one further exact residue equation); and excludes v=2 for
r>=6. It is only a local-pruning result.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import factorial


WEIGHTS_11 = tuple(range(2, 12))
WEIGHTS_6 = tuple(range(2, 7))
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
    degree11_partitions = list(partitions(11, WEIGHTS_11))
    assert len(degree11_partitions) == 14

    expected_k11_coefficients = {
        "3+4+4": Fraction(-5, 27),
        "2+4+5": Fraction(-10, 27),
        "2+2+3+4": Fraction(20, 81),
        "2+2+2+5": Fraction(20, 243),
        "2+2+2+2+3": Fraction(-35, 729),
    }
    actual_k11_coefficients = {
        label(multiplicities, WEIGHTS_11): partition_coefficient(
            Fraction(5, 3), multiplicities
        )
        for multiplicities in degree11_partitions
    }
    for term, coefficient in expected_k11_coefficients.items():
        assert actual_k11_coefficients[term] == coefficient

    expected_leading_terms = sorted(expected_k11_coefficients)
    a_group = sorted(("3+4+4", "2+2+3+4", "2+2+2+2+3"))
    b_group = sorted(("2+4+5", "2+2+2+5"))
    cases = []
    cases_to_check = [(3, 0, -5), (4, 1, -4), (5, 2, -3)]
    cases_to_check.extend((5, v, -3) for v in range(3, 13))
    cases_to_check.append((5, INFINITY, -3))
    cases_to_check.extend((r, 2, -3) for r in range(6, 12))
    cases_to_check.append((INFINITY, 2, -3))
    for r, v, expected_order in cases_to_check:
        orders = {2: 3, 3: r, 4: 0, 5: v}
        orders.update({weight: 0 for weight in range(6, 12)})
        minimum, terms = least_terms(
            11,
            WEIGHTS_11,
            orders,
            h_order=2,
            exponent_numerator=5,
            alpha=Fraction(5, 3),
        )
        assert minimum == expected_order
        leading_labels = sorted(term for term, _ in terms)
        if r in (3, 4) or (r == 5 and v == 2):
            assert leading_labels == expected_leading_terms
        elif r == 5:
            assert leading_labels == a_group
        else:
            assert leading_labels == b_group

        # The only K11 correction on U1=0 is kappa_1/t [w^6]F^(4/3).
        # Its least order is -1 on both strata, so it cannot affect the
        # degree -5 or -4 binomial residue.
        correction_orders = {2: 3, 3: r, 4: 0, 5: v, 6: 0}
        correction_minimum, correction_terms = least_terms(
            6,
            WEIGHTS_6,
            correction_orders,
            h_order=2,
            exponent_numerator=4,
            alpha=Fraction(4, 3),
        )
        if r in (3, 4, 5):
            assert correction_minimum == -1
            assert sorted(term for term, _ in correction_terms) == ["2+2+2", "2+4"]
        cases.append(
            {
                "q": 3,
                "r": r,
                "s": 0,
                "v": v,
                "K11_binomial_pole_order": minimum,
                "K11_leading_terms": leading_labels,
                "kappa1_correction_order": correction_minimum,
            }
        )

    # Normalize h=c*u^2, U2=d*u^3, U3=e*u^r, U4=a, F5=g*u^v.
    # Both K9 ties impose gamma=c^3*g/(d*e)=4/9. The K11 leading
    # coefficient divided by d^4*e/c^10 is A(z)+gamma*B(z), where
    # z=c^3*a/d^2 and the exact rational polynomials are recorded below.
    A = {
        2: Fraction(-5, 27),
        1: Fraction(20, 81),
        0: Fraction(-35, 729),
    }
    B = {1: Fraction(-10, 27), 0: Fraction(20, 243)}
    gamma = Fraction(4, 9)
    residue_polynomial = {
        power: A.get(power, Fraction(0)) + gamma * B.get(power, Fraction(0))
        for power in (2, 1, 0)
    }
    assert residue_polynomial == {
        2: Fraction(-5, 27),
        1: Fraction(20, 243),
        0: Fraction(-25, 2187),
    }

    # K10 gives 135*z^2-60*z+7=0. Hence the numerator of the
    # K11 residue over 2187 is -405*z^2+180*z-25 = -4.
    reduced_constant = Fraction(405) * Fraction(7, 135) + Fraction(-25)
    # Use the K10 relation to replace z^2 by (60*z-7)/135; the z terms cancel.
    reduced_z_coefficient = Fraction(180) - Fraction(405 * 60, 135)
    assert reduced_z_coefficient == 0
    assert reduced_constant == -4
    assert residue_polynomial[2] == Fraction(-405, 2187)
    assert residue_polynomial[1] == Fraction(180, 2187)
    assert residue_polynomial[0] == Fraction(-25, 2187)

    # For r=5, v=2, K9 supplies no gamma relation. K11 requires
    # A(z)+gamma*B(z)=0. At K10 roots both A and B are nonzero, so this
    # determines gamma but does not itself contradict the branch.
    A_after_K10 = {1: Fraction(120, 729), 0: Fraction(-28, 729)}
    B_after_scale = {1: Fraction(-90, 243), 0: Fraction(20, 243)}
    assert Fraction(135) * Fraction(4, 81) - Fraction(60) * Fraction(2, 9) + 7 == Fraction(1, 3)
    assert Fraction(135) * Fraction(49, 900) - Fraction(60) * Fraction(7, 30) + 7 == Fraction(7, 20)
    assert A_after_K10[1] != 0 and B_after_scale[1] != 0
    assert Fraction(120, 729) * Fraction(7, 30) - Fraction(28, 729) == 0
    assert Fraction(-90, 243) * Fraction(2, 9) + Fraction(20, 243) == 0

    print(
        json.dumps(
            {
                "degree11_weighted_partition_count": len(degree11_partitions),
                "checked_local_strata": cases,
                "K9_tied_cases": "r=3,v=0 and r=4,v=1 give gamma=c^3*g/(d*e)=4/9",
                "K10_normalized_relation": "135*z^2 - 60*z + 7 = 0",
                "K11_normalized_residue_after_K9_K10": "-4/2187",
                "nonzero_factor": "d^4*e/c^10",
                "K11_branch_consequences": {
                    "r_3_or_4": "excluded by the nonzero residue -4/2187",
                    "r_5": "v=2 is necessary; then gamma=2*(30*z-7)/(15*(9*z-2))",
                    "r_at_least_6": "v=2 is impossible; v>=3 or U3=0 remains",
                },
                "all_arithmetic_exact": True,
                "scope": "K11 local pruning on q=3 double-root strata only; not a pair exclusion",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
