"""Exact K15 exclusion of the q=4, s=2, v=1, r=3 double-root jet.

The script derives the lowest local residues at bands K11--K15 from the
weighted-partition expansion of F^(5/3), verifies that the first-integral
resonance corrections have strictly higher local order, and checks an exact
ideal-membership certificate

    epsilon**3 = sum(Q_i * E_i, i=11,...,15).

The nonzero branch has epsilon != 0, so this certificate excludes it.

Run normally to verify the checked-in certificate. Use --discover to solve
the bounded degree-3 exact linear system and rewrite the certificate JSON.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import sympy as sp


Z, GAMMA, EPSILON, PHI = sp.symbols("z gamma epsilon phi")
VARS = (Z, GAMMA, EPSILON, PHI)
ORDERS = {2: 4, 3: 3, 4: 2, 5: 1, 6: 0}
CERTIFICATE_PATH = Path(__file__).with_name(
    "f2_75_125_k15_q4_v1_r3_certificate.json"
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
    """Order of a term of F^(5/3), with ord(h)=2 and ord(F_i) as in ORDERS."""
    parts = tuple(partition)
    return sum(ORDERS.get(part, 0) for part in parts) + 10 - 6 * len(parts)


def lowest_residue(total: int) -> tuple[int, tuple[tuple[int, ...], ...], sp.Expr]:
    """Return the minimum order, its exact partition census, and normalized sum.

    Locally write h=c*u^2, U_2=d*u^4, U_3=e*u^3, U_4=a*u^2,
    F_5=g*u, U_6=f. The normalized invariants are

        z=c^3*a/d^2, gamma=c^3*g/(d*e),
        epsilon=e^2*c^3/d^3, phi=c^6*f/d^3.

    The returned sum is the coefficient of the common nonzero local scale
    selected by the first minimum-order partition.
    """
    census = [(local_order(part), part) for part in partitions(total)]
    minimum = min(order for order, _ in census)
    lowest = tuple(part for order, part in census if order == minimum)
    if any(max(part) > 6 for part in lowest):
        raise AssertionError((total, "a new band enters the minimum residue", lowest))

    base = Counter(lowest[0])
    base_m = sum(base.values())
    base_d = base[2] + 2 * base[4] + base[5] + 3 * base[6]
    base_e = base[3] + base[5]
    base_c = 5 - 3 * base_m - 3 * base[4] - 3 * base[5] - 6 * base[6]

    normalized = sp.Integer(0)
    for part in lowest:
        counts = Counter(part)
        m = len(part)
        coefficient = sp.binomial(sp.Rational(5, 3), m) * sp.factorial(m)
        for multiplicity in counts.values():
            coefficient /= sp.factorial(multiplicity)

        d_exp = counts[2] + 2 * counts[4] + counts[5] + 3 * counts[6]
        e_exp = counts[3] + counts[5]
        c_exp = 5 - 3 * m - 3 * counts[4] - 3 * counts[5] - 6 * counts[6]
        delta_d = d_exp - base_d
        delta_e = e_exp - base_e
        delta_c = c_exp - base_c
        if delta_e % 2:
            raise AssertionError((total, part, "odd epsilon exponent"))
        if delta_d + 3 * delta_e // 2 != 0:
            raise AssertionError((total, part, "d-scale did not cancel"))
        if delta_c - 3 * delta_e // 2 != 0:
            raise AssertionError((total, part, "c-scale did not cancel"))

        monomial = (
            Z ** counts[4]
            * GAMMA ** counts[5]
            * PHI ** counts[6]
            * EPSILON ** (delta_e // 2)
        )
        normalized += coefficient * monomial

    return minimum, lowest, sp.factor(normalized)


def band_equations() -> tuple[sp.Expr, ...]:
    """Normalized K11--K15 lowest-residue equations."""
    e11 = (
        27 * Z**2
        + 54 * Z * GAMMA
        - 36 * Z
        - 162 * GAMMA * PHI
        + 54 * PHI
        + 27 * EPSILON * GAMMA
        - 12 * GAMMA
        - 12 * EPSILON
        + 7
    )
    e12 = (
        81 * Z**3
        - 162 * Z**2
        + 486 * Z * PHI
        + 486 * Z * EPSILON * GAMMA
        - 324 * Z * EPSILON
        + 63 * Z
        - 729 * PHI**2
        + 243 * EPSILON * GAMMA**2
        + 243 * EPSILON * PHI
        - 108 * PHI
        - 324 * EPSILON * GAMMA
        - 27 * EPSILON**2
        + 126 * EPSILON
        - 7
    )
    e13 = (
        81 * Z**2 * GAMMA
        - 108 * Z**2
        + 162 * Z * PHI
        - 108 * Z * GAMMA
        - 36 * Z * EPSILON
        + 84 * Z
        + 162 * GAMMA * PHI
        + 81 * EPSILON * GAMMA**2
        - 108 * PHI
        - 108 * EPSILON * GAMMA
        + 21 * GAMMA
        + 42 * EPSILON
        - 14
    )
    e14 = sp.Rational(5, 19683) * (
        324 * EPSILON**2 * GAMMA
        - 189 * EPSILON**2
        - 729 * EPSILON * GAMMA**2 * Z
        + 486 * EPSILON * GAMMA**2
        - 1458 * EPSILON * GAMMA * PHI
        + 1944 * EPSILON * GAMMA * Z
        - 756 * EPSILON * GAMMA
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
    e15 = sp.Rational(1, 19683) * (
        -189 * EPSILON**2
        - 1215 * EPSILON * GAMMA**3
        + 4860 * EPSILON * GAMMA**2
        + 4860 * EPSILON * GAMMA * Z
        - 5670 * EPSILON * GAMMA
        + 1620 * EPSILON * PHI
        - 3780 * EPSILON * Z
        + 2100 * EPSILON
        - 7290 * GAMMA * PHI * Z
        + 4860 * GAMMA * PHI
        + 4860 * GAMMA * Z**2
        - 3780 * GAMMA * Z
        + 630 * GAMMA
        - 3645 * PHI**2
        + 9720 * PHI * Z
        - 3780 * PHI
        + 1620 * Z**3
        - 5670 * Z**2
        + 3150 * Z
        - 455
    )
    return e11, e12, e13, e14, e15


def low_order_census() -> dict[str, object]:
    equations = band_equations()
    expected_scales = (
        -sp.Rational(5, 729),
        -sp.Rational(5, 6561),
        -sp.Rational(5, 2187),
        sp.Integer(1),
        sp.Integer(1),
    )
    minima: dict[int, int] = {}
    partition_counts: dict[int, int] = {}
    for offset, total in enumerate(range(11, 16)):
        minimum, lowest, expression = lowest_residue(total)
        minima[total] = minimum
        partition_counts[total] = len(lowest)
        if sp.expand(expression - expected_scales[offset] * equations[offset]) != 0:
            raise AssertionError(f"K{total} residue formula mismatch")

    # At K15 the new resonant correction and all earlier resonance modes are
    # strictly above the order -5 binomial residue.
    k15_min = minima[15]
    kappa1_min = min(
        sum(ORDERS.get(part, 0) for part in part)
        + 8
        - 6 * len(part)
        for part in partitions(10)
    )
    corrections = {
        "kappa1_over_t_order_lower_bound": kappa1_min,
        "kappa2_over_t2_order": ORDERS[5],
        "kappa3_over_t3_order": 4,  # ord(h^2) at a double root
    }
    if not all(order > k15_min for order in corrections.values()):
        raise AssertionError((k15_min, corrections))

    return {
        "minimum_orders_K11_to_K15": minima,
        "lowest_partition_counts_K11_to_K15": partition_counts,
        "K15_correction_order_lower_bounds": corrections,
        "all_exact": True,
    }


def monomials_at_most_degree(max_degree: int) -> list[tuple[int, ...]]:
    result: list[tuple[int, ...]] = []

    def append_exact(position: int, remaining: int, prefix: list[int]) -> None:
        if position == len(VARS):
            for exponent in range(remaining + 1):
                result.append(tuple(prefix + [exponent]))
            return
        for exponent in range(remaining + 1):
            append_exact(position + 1, remaining - exponent, prefix + [exponent])

    for total in range(max_degree + 1):
        append_exact(0, total, [])
    return sorted(set(result), key=lambda exponents: (sum(exponents), exponents))


def discover_cofactors() -> list[dict[str, object]]:
    """Find degree-at-most-three Q_i with epsilon^3 = sum Q_i E_i over QQ."""
    equations = band_equations()
    mons = monomials_at_most_degree(3)
    equation_polys = [sp.Poly(eq, *VARS, domain="QQ") for eq in equations]
    columns: list[dict[tuple[int, ...], sp.Rational]] = []
    labels: list[tuple[int, tuple[int, ...]]] = []
    row_support: set[tuple[int, ...]] = {(0, 0, 3, 0)}

    for equation_index, equation in enumerate(equation_polys):
        equation_terms = {
            tuple(monomial): sp.Rational(coefficient)
            for monomial, coefficient in equation.terms()
        }
        for monomial in mons:
            column = {
                tuple(tuple_monomial[i] + monomial[i] for i in range(4)): coefficient
                for tuple_monomial, coefficient in equation_terms.items()
            }
            columns.append(column)
            labels.append((equation_index, monomial))
            row_support.update(column)

    rows = sorted(row_support, key=lambda exponents: (sum(exponents), exponents))
    matrix = sp.Matrix(
        [
            [column.get(row, sp.Integer(0)) for column in columns]
            + [sp.Integer(1) if row == (0, 0, 3, 0) else sp.Integer(0)]
            for row in rows
        ]
    )
    reduced, pivots = matrix.rref()
    if any(
        all(reduced[row, col] == 0 for col in range(len(columns)))
        and reduced[row, len(columns)] != 0
        for row in range(reduced.rows)
    ):
        raise AssertionError("degree-three certificate system is inconsistent")

    solution = [sp.Integer(0)] * len(columns)
    for row, pivot in enumerate(pivots):
        if pivot < len(columns):
            solution[pivot] = sp.Rational(reduced[row, len(columns)])

    cofactors = [sp.Integer(0)] * len(equations)
    for coefficient, (equation_index, monomial) in zip(solution, labels):
        if coefficient:
            cofactors[equation_index] += coefficient * sp.prod(
                variable**exponent for variable, exponent in zip(VARS, monomial)
            )

    if sp.expand(
        sum(cofactor * equation for cofactor, equation in zip(cofactors, equations))
        - EPSILON**3
    ) != 0:
        raise AssertionError("discovered cofactor identity failed exact expansion")

    return [serialize_polynomial(cofactor) for cofactor in cofactors]


def serialize_polynomial(expression: sp.Expr) -> dict[str, object]:
    poly = sp.Poly(expression, *VARS, domain="QQ")
    return {
        "terms": [
            {
                "exponents": list(monomial),
                "numerator": int(coefficient.p),
                "denominator": int(coefficient.q),
            }
            for monomial, coefficient in poly.terms()
        ]
    }


def load_polynomial(data: dict[str, object]) -> sp.Expr:
    expression = sp.Integer(0)
    for term in data["terms"]:
        exponents = term["exponents"]
        coefficient = sp.Rational(term["numerator"], term["denominator"])
        expression += coefficient * sp.prod(
            variable**exponent for variable, exponent in zip(VARS, exponents)
        )
    return expression


def verify_certificate(path: Path) -> dict[str, object]:
    census = low_order_census()
    data = json.loads(path.read_text(encoding="utf-8"))
    if data["target"] != "epsilon^3":
        raise AssertionError("wrong certificate target")
    cofactors = [load_polynomial(item) for item in data["cofactors"]]
    if len(cofactors) != 5:
        raise AssertionError("expected five cofactor polynomials")
    equations = band_equations()
    identity = sp.expand(
        sum(cofactor * equation for cofactor, equation in zip(cofactors, equations))
        - EPSILON**3
    )
    if identity != 0:
        raise AssertionError("stored characteristic-zero ideal identity failed")

    term_counts = [len(sp.Poly(cofactor, *VARS).terms()) for cofactor in cofactors]
    degrees = [
        max((sum(monomial) for monomial, _ in sp.Poly(cofactor, *VARS).terms()), default=-1)
        for cofactor in cofactors
    ]
    return {
        **census,
        "certificate_identity": "epsilon^3 = sum(Q_i*E_i, i=11,...,15)",
        "cofactor_term_counts": term_counts,
        "cofactor_degree_bounds": degrees,
        "identity_verified_over_QQ": True,
        "branch_contradiction": "epsilon is nonzero on q=4,s=2,v=1,r=3, but epsilon^3=0",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--discover",
        action="store_true",
        help="solve the bounded exact linear system and write the cofactor certificate",
    )
    args = parser.parse_args()
    if args.discover:
        payload = {
            "target": "epsilon^3",
            "variables": ["z", "gamma", "epsilon", "phi"],
            "equations": ["E11", "E12", "E13", "E14", "E15"],
            "cofactor_degree_bound": 3,
            "cofactors": discover_cofactors(),
        }
        CERTIFICATE_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    receipt = verify_certificate(CERTIFICATE_PATH)
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
