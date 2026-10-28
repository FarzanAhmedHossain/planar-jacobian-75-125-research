"""Exact K10 local pruning on the U1=0 branch at simple and double roots.

The first integral gives the two resonant corrections at band distance ten.
This checker independently enumerates the weighted partitions in
[w^10]F^(5/3), computes exact rational coefficients and local orders, and
checks the exact pole orders and leading residues, including tied terms. It
does not eliminate the (75,125) pair or the remaining local strata.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import factorial


WEIGHTS = tuple(range(2, 11))
INFINITY = 10_000


def partitions(total: int, weights: tuple[int, ...] = WEIGHTS, index: int = 0):
    """Yield multiplicities of weighted partitions of total."""
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


def order_of_term(
    multiplicities: tuple[int, ...],
    weights: tuple[int, ...],
    orders: dict[int, int],
    h_order: int = 2,
) -> int | None:
    if any(
        multiplicity and orders.get(weight, INFINITY) >= INFINITY
        for weight, multiplicity in zip(weights, multiplicities)
    ):
        return None
    count = sum(multiplicities)
    return h_order * (5 - 3 * count) + sum(
        multiplicity * orders[weight]
        for weight, multiplicity in zip(weights, multiplicities)
        if multiplicity
    )


def least_terms(
    orders: dict[int, int],
    h_order: int = 2,
) -> tuple[int | None, list[tuple[str, Fraction]]]:
    evaluated = []
    for multiplicities in partitions(10):
        term_order = order_of_term(
            multiplicities, WEIGHTS, orders, h_order=h_order
        )
        if term_order is not None:
            evaluated.append(
                (
                    term_order,
                    label(multiplicities, WEIGHTS),
                    partition_coefficient(Fraction(5, 3), multiplicities),
                )
            )
    if not evaluated:
        return None, []
    minimum = min(item[0] for item in evaluated)
    return minimum, [
        (term_label, coefficient)
        for term_order, term_label, coefficient in evaluated
        if term_order == minimum
    ]


def make_orders(q: int, r: int, s: int, v: int) -> dict[int, int]:
    """Known lower-slice orders; later polynomial slices take worst-case 0."""
    return {2: q, 3: r, 4: s, 5: v, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0}


def main() -> None:
    degree10_partitions = list(partitions(10))
    assert len(degree10_partitions) == 12

    # These are the only lowest terms in the q=3, s=0 survivor. The exact
    # coefficients come from the generalized binomial expansion.
    expected = {
        "2+2+2+2+2": Fraction(-7, 729),
        "2+2+2+4": Fraction(20, 243),
        "2+4+4": Fraction(-5, 27),
    }
    actual = {
        label(m, WEIGHTS): partition_coefficient(Fraction(5, 3), m)
        for m in degree10_partitions
    }
    for term, coefficient in expected.items():
        assert actual[term] == coefficient, (term, actual[term], coefficient)

    # First-integral formula at distance ten:
    # G_10 = [w^10]F^(5/3) + kappa_1/t [w^5]F^(4/3)
    #         + kappa_2/t^2 h^3.
    # With F_1=0, the only partitions of 5 in the middle coefficient are 5
    # and 2+3. Verify their exact coefficients independently.
    middle = {
        label(m, (2, 3, 4, 5)): partition_coefficient(Fraction(4, 3), m)
        for m in partitions(5, (2, 3, 4, 5))
    }
    assert middle == {"2+3": Fraction(4, 9), "5": Fraction(4, 3)}

    checked_q3 = checked_qge4_v0 = checked_q4_small_s = 0
    minimum_resonance_order = INFINITY
    q_values = list(range(4, 12)) + [INFINITY]
    r_values = list(range(3, 12)) + [INFINITY]
    s_values = list(range(12)) + [INFINITY]

    # The q=3 K9 survivors are (r,v)=(3,0), (4,1), or r>=5,v>=2.
    for r in r_values:
        if r == 3:
            allowed_v = [0]
        elif r == 4:
            allowed_v = [1]
        else:
            allowed_v = list(range(2, 13)) + [INFINITY]
        for v in allowed_v:
            for s in s_values:
                minimum, terms = least_terms(make_orders(3, r, s, v))
                names = sorted(term for term, _ in terms)
                resonance_order = min(
                    2 + v if v < INFINITY else INFINITY,
                    3 + r - 4 if r < INFINITY else INFINITY,
                )
                minimum_resonance_order = min(
                    minimum_resonance_order, resonance_order
                )
                assert resonance_order >= 2
                if s == 0:
                    assert minimum == -5
                    assert names == sorted(expected)
                else:
                    assert minimum == -5
                    assert names == ["2+2+2+2+2"]
                checked_q3 += 1

    # If q>=4 and v=0, K9 says s>=2. The unique lowest K10 pole is 5+5,
    # involving only the known nonzero leading coefficient of F5.
    for q in q_values:
        for r in r_values:
            for s in list(range(2, 12)) + [INFINITY]:
                minimum, terms = least_terms(make_orders(q, r, s, 0))
                assert minimum == -2
                assert terms == [("5+5", Fraction(5, 9))]
                checked_qge4_v0 += 1

    # For q=4 and s=0 or 1, the K9 condition s+v>=2 applies. In either
    # case the unique least term is 2+4+4, so q=4 requires s>=2.
    for r in r_values:
        for s in (0, 1):
            allowed_vs = range(2, 13) if s == 0 else range(1, 13)
            for v in list(allowed_vs) + [INFINITY]:
                minimum, terms = least_terms(make_orders(4, r, s, v))
                expected_order = -4 if s == 0 else -2
                assert minimum == expected_order
                assert terms == [("2+4+4", Fraction(-5, 27))]
                checked_q4_small_s += 1

    # Leading residue in the q=3,s=0 case. Write h=c*u^2+..., F2=d*u^3+...,
    # F4=a+..., with c,d nonzero. Clearing 729*c^10 gives the exact equation.
    leading_labels = set(expected)
    computed_residue: dict[tuple[int, int, int], Fraction] = {}
    by_label = {label(m, WEIGHTS): m for m in degree10_partitions}
    for term in leading_labels:
        multiplicities = by_label[term]
        h_power = 5 - 3 * sum(multiplicities)
        monomial = (
            10 + h_power,
            multiplicities[0],
            multiplicities[2],
        )  # exponents of (c,d,a) after clearing 729*c^10
        computed_residue[monomial] = (
            729 * actual[term]
        )
    assert computed_residue == {
        (0, 5, 0): Fraction(-7),
        (3, 3, 1): Fraction(60),
        (6, 1, 2): Fraction(-135),
    }
    discriminant = Fraction((-60) ** 2 - 4 * 135 * 7)
    assert discriminant == -180

    # Simple roots occur in the 2+1+1 quartic branch. For U1=0, the
    # preceding K6--K8 results give q>=2, r>=1, s>=1; v is unrestricted.
    # At h-order one, only four degree-ten partitions can have negative order.
    simple_labels = {
        "2+2+3+3": Fraction(10, 81),
        "2+3+5": Fraction(-10, 27),
        "3+3+4": Fraction(-5, 27),
        "5+5": Fraction(5, 9),
    }
    simple_q_values = list(range(2, 12)) + [INFINITY]
    simple_r_values = list(range(1, 12)) + [INFINITY]
    simple_s_values = list(range(1, 12)) + [INFINITY]
    simple_v_values = list(range(13)) + [INFINITY]
    simple_strata_checked = 0
    simple_pattern_counts: dict[str, int] = {}
    for q in simple_q_values:
        for r in simple_r_values:
            for s in simple_s_values:
                for v in simple_v_values:
                    orders = make_orders(q, r, s, v)
                    minimum, terms = least_terms(orders, h_order=1)
                    candidate_orders = {
                        "2+2+3+3": -7 + 2 * q + 2 * r,
                        "2+3+5": -4 + q + r + v,
                        "3+3+4": -4 + 2 * r + s,
                        "5+5": -1 + 2 * v,
                    }
                    candidate_orders = {
                        term: order
                        for term, order in candidate_orders.items()
                        if all(
                            value < INFINITY
                            for value in {
                                "2+2+3+3": (q, r),
                                "2+3+5": (q, r, v),
                                "3+3+4": (r, s),
                                "5+5": (v,),
                            }[term]
                        )
                    }
                    expected_negative = {
                        term
                        for term, order in candidate_orders.items()
                        if order == -1
                    }
                    observed_negative = {term for term, _ in terms} if minimum == -1 else set()
                    assert observed_negative == expected_negative
                    assert minimum is None or minimum >= -1
                    assert expected_negative <= set(simple_labels)
                    if expected_negative:
                        assert minimum == -1
                        for term, coefficient in terms:
                            assert simple_labels[term] == coefficient
                    else:
                        assert minimum is None or minimum >= 0

                    resonance_order = min(
                        1 + v if v < INFINITY else INFINITY,
                        q + r - 2 if q < INFINITY and r < INFINITY else INFINITY,
                    )
                    assert resonance_order >= 1
                    pattern = ",".join(sorted(expected_negative)) or "regular"
                    simple_pattern_counts[pattern] = simple_pattern_counts.get(pattern, 0) + 1
                    simple_strata_checked += 1

    expected_simple_patterns = {
        "regular",
        "2+2+3+3",
        "3+3+4",
        "5+5",
        "2+2+3+3,3+3+4",
        "2+2+3+3,2+3+5,5+5",
        "2+2+3+3,2+3+5,3+3+4,5+5",
        "3+3+4,5+5",
    }
    assert set(simple_pattern_counts) == expected_simple_patterns, simple_pattern_counts

    # After writing h=c*u, F2=d*u^q, F3=e*u^r, F4=f*u^s,
    # F5=g*u^v, the four normalized coefficients under the common scale
    # 81*c^7/5 are 2, -6, -3, and 9 on their respective monomials.
    normalized_simple_coefficients = {
        "2+2+3+3": Fraction(81, 5) * Fraction(10, 81),
        "2+3+5": Fraction(81, 5) * Fraction(-10, 27),
        "3+3+4": Fraction(81, 5) * Fraction(-5, 27),
        "5+5": Fraction(81, 5) * Fraction(5, 9),
    }
    assert normalized_simple_coefficients == {
        "2+2+3+3": Fraction(2),
        "2+3+5": Fraction(-6),
        "3+3+4": Fraction(-3),
        "5+5": Fraction(9),
    }

    print(
        json.dumps(
            {
                "degree10_partition_count": len(degree10_partitions),
                "q3_order_strata_checked": checked_q3,
                "qge4_v0_order_strata_checked": checked_qge4_v0,
                "q4_s0_or_s1_order_strata_checked": checked_q4_small_s,
                "resonance_correction_minimum_order": minimum_resonance_order,
                "simple_root_order_strata_checked": simple_strata_checked,
                "simple_root_negative_residue_patterns": simple_pattern_counts,
                "simple_root_tie_equations": {
                    "A+C": "q=2,r=1,s=1,v>=1: 2*d^2 - 3*c^3*f = 0",
                    "A+B+D": "q=2,r=1,v=0,s>=2: 2*d^2*e^2 - 6*c^3*d*e*g + 9*c^6*g^2 = 0",
                    "A+B+C+D": "q=2,r=1,s=1,v=0: 2*d^2*e^2 - 6*c^3*d*e*g - 3*c^3*e^2*f + 9*c^6*g^2 = 0",
                    "C+D": "q>=3,r=1,s=1,v=0: e^2*f - 3*c^3*g^2 = 0",
                },
                "q3_survivor_condition": {
                    "ord_U4": 0,
                    "normalized_leading_value": "z = c^3*U4(alpha)/d^2",
                    "equation": "135*z^2 - 60*z + 7 = 0",
                    "discriminant": int(discriminant),
                    "roots_over_C": ["(10 + i*sqrt(5))/45", "(10 - i*sqrt(5))/45"],
                },
                "qge4_conditions": {
                    "ord_F5_zero_is_impossible": True,
                    "q_equals_4_requires_ord_U4_at_least_2": True,
                },
                "all_exact_rational_arithmetic": True,
                "scope": "K10 local pruning on the U1=0 branch at simple and double roots; not a pair exclusion",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
