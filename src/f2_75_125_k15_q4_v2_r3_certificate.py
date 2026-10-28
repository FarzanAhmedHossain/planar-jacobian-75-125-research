"""Exact K11--K15 exclusion of the q=4, s=2, v=2, r=3 jet.

At a double root, the order profile is

    ord(h)=2, ord(U2)=4, ord(U3)=3, ord(U4)=2,
    ord(F5)=2, ord(U6)=0,

with nonzero leading coefficients for U3 and F5.  The lowest K11--K15
residues involve no part 5, so they depend only on

    z = c^3*a/d^2, epsilon = e^2*c^3/d^3, phi = c^6*f/d^3.

The script regenerates those residues from the weighted-partition formula,
checks that all resonance corrections have higher local order, and verifies
the exact certificate epsilon = sum(Q_i*A_i) over QQ.  The branch has
epsilon != 0, so it is impossible.

Use --discover to solve the degree-at-most-two multiplier system and rewrite
the sparse certificate JSON.  The default invocation only replays it.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import sympy as sp


Z, EPSILON, PHI = sp.symbols("z epsilon phi")
VARS = (Z, EPSILON, PHI)
ORDERS = {2: 4, 3: 3, 4: 2, 5: 2, 6: 0}
TOTALS = tuple(range(11, 16))
CERTIFICATE_PATH = Path(__file__).with_name(
    "f2_75_125_k15_q4_v2_r3_certificate.json"
)


@lru_cache(None)
def partitions(total: int, minimum_part: int = 2) -> tuple[tuple[int, ...], ...]:
    if total == 0:
        return ((),)
    return tuple(
        (part,) + rest
        for part in range(minimum_part, total + 1)
        for rest in partitions(total - part, part)
    )


def local_order(partition: Iterable[int]) -> int:
    parts = tuple(partition)
    return sum(ORDERS.get(part, 0) for part in parts) + 10 - 6 * len(parts)


def lowest_residue(total: int) -> tuple[int, tuple[tuple[int, ...], ...], sp.Expr]:
    """Return minimum order, its exact partition list, and normalized sum."""
    census = [(local_order(part), part) for part in partitions(total)]
    minimum = min(order for order, _ in census)
    lowest = tuple(part for order, part in census if order == minimum)

    base = Counter(lowest[0])
    base_m = sum(base.values())
    base_d = base[2] + 2 * base[4] + 3 * base[6]
    base_e = base[3]
    base_c = 5 - 3 * base_m - 3 * base[4] - 6 * base[6]

    normalized = sp.Integer(0)
    for part in lowest:
        counts = Counter(part)
        multiplicity = len(part)
        coefficient = sp.binomial(sp.Rational(5, 3), multiplicity) * sp.factorial(
            multiplicity
        )
        for count in counts.values():
            coefficient /= sp.factorial(count)

        if counts[5]:
            raise AssertionError((total, part, "F5 unexpectedly enters the residue"))
        delta_d = counts[2] + 2 * counts[4] + 3 * counts[6] - base_d
        delta_e = counts[3] - base_e
        delta_c = (
            5 - 3 * multiplicity - 3 * counts[4] - 6 * counts[6] - base_c
        )
        if delta_e % 2:
            raise AssertionError((total, part, "epsilon exponent is not integral"))
        if delta_d + 3 * delta_e // 2 != 0:
            raise AssertionError((total, part, "d scale did not cancel"))
        if delta_c - 3 * delta_e // 2 != 0:
            raise AssertionError((total, part, "c scale did not cancel"))

        normalized += (
            coefficient
            * Z**counts[4]
            * EPSILON ** (delta_e // 2)
            * PHI**counts[6]
        )

    return minimum, lowest, sp.factor(normalized)


def band_equations() -> tuple[sp.Expr, ...]:
    """Polynomial equations, each a nonzero scalar multiple of a residue."""
    a11 = -12 * EPSILON + 54 * PHI + 27 * Z**2 - 36 * Z + 7
    a12 = (
        -27 * EPSILON**2
        + 243 * EPSILON * PHI
        - 324 * EPSILON * Z
        + 126 * EPSILON
        - 729 * PHI**2
        + 486 * PHI * Z
        - 108 * PHI
        + 81 * Z**3
        - 162 * Z**2
        + 63 * Z
        - 7
    )
    a13 = (
        18 * EPSILON * Z
        - 21 * EPSILON
        - 81 * PHI * Z
        + 54 * PHI
        + 54 * Z**2
        - 42 * Z
        + 7
    )
    a14 = (
        -189 * EPSILON**2
        + 972 * EPSILON * PHI
        + 486 * EPSILON * Z**2
        - 1134 * EPSILON * Z
        + 315 * EPSILON
        - 729 * PHI**2
        - 729 * PHI * Z**2
        + 972 * PHI * Z
        - 189 * PHI
        + 324 * Z**3
        - 378 * Z**2
        + 126 * Z
        - 13
    )
    a15 = (
        -189 * EPSILON**2
        + 1620 * EPSILON * PHI
        - 3780 * EPSILON * Z
        + 2100 * EPSILON
        - 3645 * PHI**2
        + 9720 * PHI * Z
        - 3780 * PHI
        + 1620 * Z**3
        - 5670 * Z**2
        + 3150 * Z
        - 455
    )
    return a11, a12, a13, a14, a15


def correction_orders(total: int) -> dict[str, int | None]:
    """Lower bounds for resonance terms in G_total at the double root."""
    kappa1_degree = total - 5
    kappa1 = min(
        8 + sum(ORDERS.get(part, 0) for part in part) - 6 * len(part)
        for part in partitions(kappa1_degree)
    )
    f_index = total - 10
    kappa2 = None if f_index == 1 else ORDERS.get(f_index, 0)
    kappa3 = 4 if total == 15 else None
    return {"kappa1": kappa1, "kappa2": kappa2, "kappa3": kappa3}


def residue_census() -> dict[str, object]:
    equations = band_equations()
    expected_scales = (
        -sp.Rational(5, 729),
        -sp.Rational(5, 6561),
        sp.Rational(10, 2187),
        sp.Rational(5, 19683),
        sp.Rational(1, 19683),
    )
    expected_minima = (-1, -2, -3, -4, -5)
    expected_counts = (5, 11, 7, 13, 11)
    minima: dict[str, int] = {}
    counts: dict[str, int] = {}
    corrections: dict[str, dict[str, int | None]] = {}

    for index, total in enumerate(TOTALS):
        minimum, lowest, expression = lowest_residue(total)
        minima[str(total)] = minimum
        counts[str(total)] = len(lowest)
        if minimum != expected_minima[index] or len(lowest) != expected_counts[index]:
            raise AssertionError((total, minimum, len(lowest)))
        if sp.expand(expression - expected_scales[index] * equations[index]) != 0:
            raise AssertionError(f"K{total} normalized residue mismatch")

        bounds = correction_orders(total)
        corrections[str(total)] = bounds
        for key, bound in bounds.items():
            if bound is not None and bound <= minimum:
                raise AssertionError((total, key, bound, minimum))

    return {
        "minimum_orders_K11_to_K15": minima,
        "lowest_partition_counts_K11_to_K15": counts,
        "resonance_order_lower_bounds": corrections,
    }


def monomials_at_most_degree(max_degree: int) -> list[tuple[int, int, int]]:
    return [
        (i, j, k)
        for total in range(max_degree + 1)
        for i in range(total + 1)
        for j in range(total - i + 1)
        for k in (total - i - j,)
    ]


def sparse_terms(poly: sp.Expr) -> list[dict[str, object]]:
    terms = []
    for exponents, coefficient in sp.Poly(poly, *VARS).terms():
        coefficient = sp.Rational(coefficient)
        terms.append(
            {
                "exponents": list(exponents),
                "numerator": int(coefficient.p),
                "denominator": int(coefficient.q),
            }
        )
    return terms


def discover_cofactors() -> list[dict[str, object]]:
    equations = band_equations()
    degree_bound = 2
    monomials = monomials_at_most_degree(degree_bound)
    unknowns = sp.symbols(f"q0:{len(equations) * len(monomials)}")
    expression = sum(
        unknowns[i * len(monomials) + j]
        * Z**monomial[0]
        * EPSILON**monomial[1]
        * PHI**monomial[2]
        * equations[i]
        for i in range(len(equations))
        for j, monomial in enumerate(monomials)
    ) - EPSILON
    poly = sp.Poly(sp.expand(expression), *VARS)
    matrix, target = sp.linear_eq_to_matrix(
        [poly.coeff_monomial(monomial) for monomial in poly.monoms()], unknowns
    )
    if matrix.rank() != matrix.row_join(target).rank():
        raise AssertionError("degree-two cofactor system has no solution")

    solution = sp.linsolve((matrix, target)).args[0]
    # Set all free parameters to zero to get one deterministic exact witness.
    free_parameters = set().union(*(value.free_symbols for value in solution))
    solution = tuple(value.subs({symbol: 0 for symbol in free_parameters}) for value in solution)
    cofactors = []
    for i in range(len(equations)):
        cofactor = sum(
            solution[i * len(monomials) + j]
            * Z**monomial[0]
            * EPSILON**monomial[1]
            * PHI**monomial[2]
            for j, monomial in enumerate(monomials)
        )
        cofactors.append({"terms": sparse_terms(sp.factor(cofactor))})
    return cofactors


def read_cofactors() -> list[sp.Expr]:
    data = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
    if data.get("target") != "epsilon" or data.get("variables") != [
        "z",
        "epsilon",
        "phi",
    ]:
        raise AssertionError("certificate header mismatch")
    if data.get("equations") != ["A11", "A12", "A13", "A14", "A15"]:
        raise AssertionError("certificate equation labels mismatch")
    if data.get("cofactor_degree_bound") != 2:
        raise AssertionError("certificate degree bound mismatch")

    cofactors = []
    for entry in data["cofactors"]:
        poly = sp.Integer(0)
        for term in entry["terms"]:
            numerator = int(term["numerator"])
            denominator = int(term["denominator"])
            if denominator == 0:
                raise AssertionError("zero denominator in certificate")
            monomial = sp.prod(var**power for var, power in zip(VARS, term["exponents"]))
            poly += sp.Rational(numerator, denominator) * monomial
        cofactors.append(sp.expand(poly))
    return cofactors


def certificate_data(cofactors: list[dict[str, object]]) -> dict[str, object]:
    return {
        "target": "epsilon",
        "variables": ["z", "epsilon", "phi"],
        "equations": ["A11", "A12", "A13", "A14", "A15"],
        "cofactor_degree_bound": 2,
        "cofactors": cofactors,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--discover",
        action="store_true",
        help="solve and persist the degree-two exact cofactor certificate",
    )
    args = parser.parse_args()

    census = residue_census()
    if args.discover:
        CERTIFICATE_PATH.write_text(
            json.dumps(certificate_data(discover_cofactors()), indent=2) + "\n",
            encoding="utf-8",
        )

    cofactors = read_cofactors()
    equations = band_equations()
    identity = sp.expand(sum(q * equation for q, equation in zip(cofactors, equations)))
    if sp.expand(identity - EPSILON) != 0:
        raise AssertionError("cofactor identity failed over QQ")
    if any(sp.Poly(cofactor, *VARS).total_degree() > 2 for cofactor in cofactors):
        raise AssertionError("cofactor exceeds the recorded degree bound")

    print(
        json.dumps(
            {
                "profile": "q=4,s=2,v=2,r=3",
                "identity_verified_over_QQ": True,
                "identity_target": "epsilon",
                "epsilon_nonzero_on_branch": True,
                **census,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
