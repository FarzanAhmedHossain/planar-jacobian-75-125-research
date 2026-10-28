"""Exact three-band exclusion for q=4,s=2,v=1,r>=4 at a double root.

At a double root of h on U1=0, this checker enumerates all weighted
partitions for K11, K12, and K13 with q=ord(U2)=4, s=ord(U4)=2,
v=ord(U5/t)=1, and r=ord(U3)>=4 (including U3=0). It verifies that only
the displayed lower-slice terms attain the relevant poles, checks that
resonance corrections have higher order, and reduces the exact residues:
K11 gives phi=(9z-2)/27, K13 gives z=1/3, and K12 then forces g=0,
contradicting v=1. This is a local branch exclusion only.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import factorial

import sympy as sp


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


def coefficient(alpha: Fraction, multiplicities: tuple[int, ...]) -> Fraction:
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


def local_census(
    total: int,
    alpha: Fraction,
    orders: dict[int, int],
    h_order: int = 2,
) -> tuple[int, list[tuple[str, Fraction]]]:
    weights = tuple(range(2, total + 1))
    evaluated = []
    for multiplicities in partitions(total, weights):
        if any(
            multiplicity and orders.get(weight, INFINITY) >= INFINITY
            for weight, multiplicity in zip(weights, multiplicities)
        ):
            continue
        count = sum(multiplicities)
        order = h_order * (3 * alpha - 3 * count) + sum(
            multiplicity * orders[weight]
            for weight, multiplicity in zip(weights, multiplicities)
            if multiplicity
        )
        evaluated.append(
            (
                order,
                label(multiplicities, weights),
                coefficient(alpha, multiplicities),
            )
        )
    minimum = min(item[0] for item in evaluated)
    return minimum, [
        (term, coefficient)
        for order, term, coefficient in evaluated
        if order == minimum
    ]


def main() -> None:
    expected = {
        11: (-1, {"5+6": Fraction(10, 9), "2+4+5": Fraction(-10, 27), "2+2+2+5": Fraction(20, 243)}),
        12: (
            -2,
            {
                "6+6": Fraction(5, 9),
                "4+4+4": Fraction(-5, 81),
                "2+5+5": Fraction(-5, 27),
                "2+4+6": Fraction(-10, 27),
                "2+2+4+4": Fraction(10, 81),
                "2+2+2+6": Fraction(20, 243),
                "2+2+2+2+4": Fraction(-35, 729),
                "2+2+2+2+2+2": Fraction(35, 6561),
            },
        ),
        13: (
            -3,
            {
                "4+4+5": Fraction(-5, 27),
                "2+5+6": Fraction(-10, 27),
                "2+2+4+5": Fraction(20, 81),
                "2+2+2+2+5": Fraction(-35, 729),
            },
        ),
    }

    checked_orders = []
    for r in list(range(4, 12)) + [INFINITY]:
        orders = {2: 4, 3: r, 4: 2, 5: 1}
        orders.update({weight: 0 for weight in range(6, 14)})
        for total, alpha in ((11, Fraction(5, 3)), (12, Fraction(5, 3)), (13, Fraction(5, 3))):
            minimum, terms = local_census(total, alpha, orders)
            expected_order, expected_terms = expected[total]
            assert minimum == expected_order, (r, total, minimum, terms)
            assert dict(terms) == expected_terms, (r, total, terms)

        # Resonant corrections at these three bands have orders at least
        # 2, 1, and 0, respectively; the kappa_2 modes are U1=0, U2, U3.
        corr_orders = {**orders}
        correction_mins = []
        for total, lower_bound in ((6, 2), (7, 1), (8, 0)):
            corr_min, _ = local_census(
                total,
                Fraction(4, 3),
                corr_orders,
            )
            assert corr_min >= lower_bound, (r, total, corr_min)
            correction_mins.append(corr_min)
        checked_orders.append(r)

    z, phi, W = sp.symbols("z phi W")

    # Normalize by the nonzero leading monomials:
    # z=c^3*a/d^2, phi=c^6*f/d^3, W=c^9*g^2/d^5.
    # The exact K11 residue divided by d^3*g/c^7 is this expression.
    k11 = sp.Rational(10, 9) * phi - sp.Rational(10, 27) * z + sp.Rational(20, 243)
    phi_forced = sp.Rational(1, 27) * (9 * z - 2)
    assert sp.factor(sp.solve(k11, phi)[0] - phi_forced) == 0

    # The K13 leading residue divided by d^4*g/c^10.
    k13 = (
        -sp.Rational(5, 27) * z**2
        - sp.Rational(10, 27) * phi
        + sp.Rational(20, 81) * z
        - sp.Rational(35, 729)
    )
    k13_after_k11 = sp.factor(k13.subs(phi, phi_forced))
    assert k13_after_k11 == -sp.Rational(15, 729) * (3 * z - 1) ** 2
    assert sp.solve(k13_after_k11, z) == [sp.Rational(1, 3)]

    # The K12 leading residue divided by d^6/c^13.
    k12 = (
        sp.Rational(5, 9) * phi**2
        - sp.Rational(5, 81) * z**3
        - sp.Rational(5, 27) * W
        - sp.Rational(10, 27) * z * phi
        + sp.Rational(10, 81) * z**2
        + sp.Rational(20, 243) * phi
        - sp.Rational(35, 729) * z
        + sp.Rational(35, 6561)
    )
    k12_after_k11 = sp.factor(k12.subs(phi, phi_forced))
    assert sp.simplify(
        k12_after_k11 + sp.Rational(5, 2187) * (81 * W + (3 * z - 1) ** 3)
    ) == 0
    k12_after_k13 = sp.factor(k12_after_k11.subs(z, sp.Rational(1, 3)))
    assert k12_after_k13 == -sp.Rational(5, 27) * W
    assert sp.solve(k12_after_k13, W) == [0]

    print(
        json.dumps(
            {
                "r_orders_checked": ["4..11", "U3=0"],
                "K11_K12_K13_minimum_orders": {"K11": -1, "K12": -2, "K13": -3},
                "resonance_correction_minima": {
                    "K11_kappa1": "at least 2",
                    "K12_kappa1": "at least 1",
                    "K13_kappa1": "at least 0",
                },
                "K11_relation": "phi=(9*z-2)/27",
                "K13_after_K11": "-15*(3*z-1)^2/729",
                "K12_after_K11": "-5*(81*W+(3*z-1)^3)/2187",
                "forced": {"z": "1/3", "W": "0", "g": "0"},
                "contradiction": "v=1 requires g!=0",
                "all_coefficients_exact": True,
                "scope": "excludes only q=4,s=2,v=1,r>=4 at a double root on U1=0",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
