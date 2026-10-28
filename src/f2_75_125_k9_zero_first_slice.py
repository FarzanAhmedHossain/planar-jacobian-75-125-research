"""Exact K9 valuation split for the U1=0 deeper-double-root branch.

The checker enumerates the degree-nine weighted partitions in F^(5/3),
checks the one resonant contribution at degree nine, and replays the local
order classification at a double root of h.  It proves necessary K9
conditions only; it does not eliminate the whole (75,125) system.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import factorial


WEIGHTS = tuple(range(2, 10))
INFINITY = 10_000


def partitions(total: int, weights: tuple[int, ...] = WEIGHTS, index: int = 0):
    """Yield multiplicity tuples for partitions using the given weights."""
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


def multinomial_binomial(alpha: Fraction, multiplicities: tuple[int, ...]) -> Fraction:
    k = sum(multiplicities)
    value = binomial(alpha, k) * factorial(k)
    for multiplicity in multiplicities:
        value /= factorial(multiplicity)
    return value


def order_of_term(
    multiplicities: tuple[int, ...],
    weights: tuple[int, ...],
    orders: dict[int, int],
    h_order: int,
    h_power_base: int,
) -> int | None:
    if any(multiplicity and orders.get(weight, INFINITY) >= INFINITY for weight, multiplicity in zip(weights, multiplicities)):
        return None
    k = sum(multiplicities)
    return h_order * (h_power_base - 3 * k) + sum(
        multiplicity * orders[weight]
        for weight, multiplicity in zip(weights, multiplicities)
        if multiplicity
    )


def leading_terms(
    total: int,
    weights: tuple[int, ...],
    orders: dict[int, int],
    h_order: int,
    h_power_base: int,
) -> tuple[int | None, list[tuple[tuple[int, ...], Fraction]]]:
    evaluated: list[tuple[int, tuple[int, ...], Fraction]] = []
    for multiplicities in partitions(total, weights):
        term_order = order_of_term(
            multiplicities, weights, orders, h_order, h_power_base
        )
        if term_order is not None:
            evaluated.append(
                (
                    term_order,
                    multiplicities,
                    multinomial_binomial(Fraction(5, 3), multiplicities),
                )
            )
    if not evaluated:
        return None, []
    minimum = min(item[0] for item in evaluated)
    return minimum, [
        (multiplicities, coefficient)
        for term_order, multiplicities, coefficient in evaluated
        if term_order == minimum
    ]


def fmt(multiplicities: tuple[int, ...], weights: tuple[int, ...]) -> str:
    return "+".join(
        str(weight)
        for weight, multiplicity in zip(weights, multiplicities)
        for _ in range(multiplicity)
    )


def main() -> None:
    degree9_partitions = list(partitions(9))
    coefficients = {
        fmt(m, WEIGHTS): multinomial_binomial(Fraction(5, 3), m)
        for m in degree9_partitions
    }
    expected_coefficients = {
        "3+3+3": Fraction(-5, 81),
        "2+2+2+3": Fraction(20, 243),
        "2+2+5": Fraction(-5, 27),
        "4+5": Fraction(10, 9),
        "2+3+4": Fraction(-10, 27),
    }
    for key, value in expected_coefficients.items():
        assert coefficients[key] == value, (key, coefficients[key], value)

    # [w^4]F^(4/3) is the only factor multiplying the K5 resonance at K9.
    mode1_coefficients = {
        fmt(m, (2, 3, 4)): multinomial_binomial(Fraction(4, 3), m)
        for m in partitions(4, (2, 3, 4))
    }
    assert mode1_coefficients == {
        "2+2": Fraction(2, 9),
        "4": Fraction(4, 3),
    }
    # At h-order 2, q=ord(U2)>=3 and ord(U4)>=1 make this resonance regular.
    for q in range(3, 12):
        for s4 in range(1, 12):
            mode_orders = {2: q, 3: 99, 4: s4}
            term_orders = [
                order_of_term(m, (2, 3, 4), mode_orders, 2, 4)
                for m in partitions(4, (2, 3, 4))
            ]
            assert min(v for v in term_orders if v is not None) >= 2

    # K7/K8 leave the following U1=0 deeper-double-root hypotheses:
    # q=ord(U2)>=3, r=ord(U3)>=1, s=ord(U4)>=1; if q=3 then r>=2.
    # The following finite census verifies all exact leading-term cases up to
    # the support caps (zero polynomials are represented by INFINITY).
    checked = 0
    unique_obstructions = 0
    tie_cases = 0
    allowed_cases = 0
    q_values = list(range(3, 12)) + [INFINITY]
    r_values = list(range(1, 12)) + [INFINITY]
    s_values = list(range(1, 12)) + [INFINITY]
    v_values = list(range(13)) + [INFINITY]

    for q in q_values:
        for r in r_values:
            if q == 3 and r < 2:
                continue
            for s4 in s_values:
                for v5 in v_values:
                    # Later coefficients are regular at this nonzero root.
                    # Order zero is the worst case for pole cancellation.
                    orders = {2: q, 3: r, 4: s4, 5: v5, 6: 0, 7: 0, 8: 0, 9: 0}
                    minimum, terms = leading_terms(9, WEIGHTS, orders, 2, 5)
                    checked += 1

                    if q == 3 and r == 2:
                        assert minimum == -3
                        assert len(terms) == 1 and fmt(terms[0][0], WEIGHTS) == "2+2+2+3"
                        unique_obstructions += 1
                        continue

                    if q >= 4 and r in (1, 2):
                        expected = -5 if r == 1 else -2
                        assert minimum == expected
                        assert len(terms) == 1 and fmt(terms[0][0], WEIGHTS) == "3+3+3"
                        unique_obstructions += 1
                        continue

                    if q == 3 and r == 3:
                        if v5 == 0:
                            assert minimum == -2
                            assert {fmt(m, WEIGHTS) for m, _ in terms} == {
                                "2+2+2+3",
                                "2+2+5",
                            }
                            tie_cases += 1
                        else:
                            assert minimum == -2
                            assert len(terms) == 1 and fmt(terms[0][0], WEIGHTS) == "2+2+2+3"
                            unique_obstructions += 1
                        continue

                    if q == 3 and r == 4:
                        if v5 == 1:
                            assert minimum == -1
                            assert {fmt(m, WEIGHTS) for m, _ in terms} == {
                                "2+2+2+3",
                                "2+2+5",
                            }
                            tie_cases += 1
                        else:
                            assert minimum < 0 and len(terms) == 1
                            assert fmt(terms[0][0], WEIGHTS) == (
                                "2+2+5" if v5 == 0 else "2+2+2+3"
                            )
                            unique_obstructions += 1
                        continue

                    if q == 3 and r >= 5:
                        if v5 < 2:
                            assert minimum == v5 - 2
                            assert len(terms) == 1 and fmt(terms[0][0], WEIGHTS) == "2+2+5"
                            unique_obstructions += 1
                        else:
                            assert minimum is None or minimum >= 0
                            allowed_cases += 1
                        continue

                    if q >= 4 and r >= 3:
                        if s4 + v5 == 1:
                            assert minimum == -1
                            assert len(terms) == 1 and fmt(terms[0][0], WEIGHTS) == "4+5"
                            unique_obstructions += 1
                        else:
                            assert minimum is None or minimum >= 0
                            allowed_cases += 1
                        continue

                    raise AssertionError((q, r, s4, v5, minimum, terms))

    # Exact cancellation in the two tied cases. Write h=c*u^2,
    # U2=d*u^3, U3=e*u^r, U5=g*u^(r-3), for r=3 or 4.
    coeff_a = Fraction(20, 243)  # (2+2+2+3)
    coeff_b = Fraction(-5, 27)   # (2+2+5)
    # The coefficient is (20/243)c^-7*d^3*e-(5/27)c^-4*g*d^2;
    # after clearing (5/243)c^-7*d^2 it is 4*d*e-9*c^3*g.
    assert coeff_a / Fraction(5, 243) == 4
    assert -coeff_b / Fraction(5, 243) == 9

    print(
        json.dumps(
            {
                "degree9_partition_count": len(degree9_partitions),
                "checked_order_strata": checked,
                "unique_nonzero_pole_obstructions": unique_obstructions,
                "leading_coefficient_tie_cases": tie_cases,
                "nonnegative_or_absent_pole_cases": allowed_cases,
                "resonance_mode_order_minimum": 2,
                "necessary_K9_conditions": {
                    "q_equals_3": (
                        "r=2 is impossible; if r=3 then v=0 and 4*d*e=9*c^3*g; "
                        "if r=4 then v=1 and the same relation holds; "
                        "if r>=5 then v>=2"
                    ),
                    "q_at_least_4": "r>=3 and s+v>=2",
                    "definitions": "q=ord(U2), r=ord(U3), s=ord(U4), v=ord(U5)",
                    "root_normalization": "ord(h)=2; c,d,e,g are the first nonzero local coefficients",
                },
                "all_exact": True,
                "scope": "K9 local necessary conditions on the U1=0 deeper-double-root branch only",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
