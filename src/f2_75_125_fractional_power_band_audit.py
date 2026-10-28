"""Replay the corrected fractional-power law for the F2 edge bands.

It verifies q_s = [z^s](P(z)^(5/3) + kappa*z^5*P(z)^(4/3))
for saved bands and provides an exact polynomial continuation formula.
Every required division by U=t*h is checked with exact polynomial division.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
from math import factorial
from pathlib import Path

import sympy as sp

from f2_75_125_sparse_terminal_q_postcheck import load_sparse


@lru_cache(None)
def partitions(distance: int, smallest: int = 1):
    if distance == 0:
        return ((),)
    result = []
    for i in range(smallest, distance + 1):
        for tail in partitions(distance - i, i):
            result.append((i, *tail))
    return tuple(result)


def fractional_coefficient(series, exponent, distance, U, R, field):
    """Return [z^distance](sum p_i z^i)^exponent, or reject nonpolynomiality."""
    total = R.zero
    base_power = int(3 * exponent)
    if sp.Rational(base_power, 3) != exponent:
        raise ValueError("only the 5/3 and 4/3 exponents are supported")
    for part in partitions(distance):
        counts = {}
        for i in part:
            counts[i] = counts.get(i, 0) + 1
        k = len(part)
        scalar = sp.binomial(exponent, k) * factorial(k)
        for count in counts.values():
            scalar /= factorial(count)
        value = R.one
        for index, count in counts.items():
            value *= series[index] ** count
        power = base_power - 3 * k
        if power >= 0:
            value *= U ** power
        else:
            divisor = U ** (-power)
            quotient, remainder = value.div(divisor)
            if remainder:
                raise AssertionError(
                    f"nonpolynomial {exponent} coefficient at z^{distance}: "
                    f"remainder has {len(remainder)} terms after division by U^{(-power)}")
            value = quotient
        total += value * field.convert(scalar)
    return total


def corrected_q(series, distance, U, kappa, R, field):
    value = fractional_coefficient(series, sp.Rational(5, 3), distance, U, R, field)
    if distance >= 5:
        value += kappa * fractional_coefficient(
            series, sp.Rational(4, 3), distance - 5, U, R, field)
    return value


def product_coefficient(left, right, distance, R):
    return sum((left[i] * right[distance - i]
                for i in range(distance + 1)
                if i < len(left) and distance - i < len(right)), R.zero)


def integer_power_coefficient(series, exponent, distance, R):
    values = [R.one] + [R.zero] * distance
    for _ in range(exponent):
        values = [sum((values[j] * series[i - j]
                       for j in range(i + 1)
                       if i - j < len(series)), R.zero)
                  for i in range(distance + 1)]
    return values[distance]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--through", type=int)
    args = parser.parse_args()
    s = load_sparse(args.checkpoint)
    R, field = s["R"], s["field"]
    t = s["gens"][0]
    U = t * s["h"]
    series = s["p"]
    depth = len(s["q"]) - 1 if args.through is None else args.through
    if depth >= len(series):
        raise ValueError("P bands do not reach the requested distance")
    comparisons = []
    corrected_relation = []
    for distance in range(depth + 1):
        if distance <= 9:
            expected = corrected_q(series, distance, U, s["kappa"], R, field)
            actual = s["q"][distance]
            if expected != actual:
                difference = expected - actual
                raise AssertionError(
                    f"q_{distance} disagrees with corrected power law: "
                    f"{len(difference)} terms")
            comparisons.append({"distance": distance, "terms": len(expected), "equal": True})
        q3 = integer_power_coefficient(s["q"], 3, distance, R)
        p5 = integer_power_coefficient(series, 5, distance, R)
        correction = R.zero
        if distance >= 5:
            p3 = [integer_power_coefficient(series, 3, i, R)
                  for i in range(distance - 4)]
            correction = 3 * s["kappa"] * product_coefficient(
                p3, s["q"], distance - 5, R)
        if distance >= 10:
            target = distance - 10
            q2 = [integer_power_coefficient(s["q"], 2, i, R)
                  for i in range(target + 1)]
            correction += 3*s["kappa2"]*product_coefficient(
                series, q2, target, R)
        if distance >= 15:
            correction += 3*s["kappa3"]*integer_power_coefficient(
                series, 4, distance - 15, R)
        relation_residual = q3 - p5 - correction
        if relation_residual:
            raise AssertionError(
                f"corrected cubic relation fails at distance {distance}: "
                f"{len(relation_residual)} terms")
        corrected_relation.append({"distance": distance, "terms": len(q3), "equal": True})
    print({"branch": s["data"]["branch"], "depth": depth,
           "law": "Q^3=P^5+3*kappa*z^5*P^3*Q+3*kappa2*z^10*P*Q^2+3*kappa3*z^15*P^4",
           "checks": comparisons,
           "corrected_cubic_relation":"Q^3=P^5+3*kappa*z^5*P^3*Q+3*kappa2*z^10*P*Q^2+3*kappa3*z^15*P^4",
           "corrected_cubic_checks":corrected_relation,
           "status": "PASS"})


if __name__ == "__main__":
    main()
