"""Exact K12/K14 selector and post-cube residues for the 2+2 F2 branch.

Scope: the U1=0 double-root stratum q=ord(U2)=4, s=ord(U4)=2,
r=ord(U3)>=4, v=ord(U5/t)>=2, including U3=0 and U5=0. It verifies the
K12/K14 lowest-pole selector, then subtracts the tuned local cube
(h + (B2/3) w^2)^3. K14 forces F7 to vanish at each double root; K14 and
K15 then give further exact necessary equations. This is not an exclusion.

No assumption is made that t divides U5: F5=U5/t may have a pole at zero.
The local computations are at roots tau=1,alpha, where t is a unit.
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from fractions import Fraction
from math import factorial

import sympy as sp


INFINITY = 10_000


def partitions(total: int, minimum_part: int = 2):
    """Yield nondecreasing weighted partitions of total."""
    if total == 0:
        yield ()
        return
    for part in range(minimum_part, total + 1):
        for tail in partitions(total - part, part):
            yield (part,) + tail


def partition_coefficient(alpha: sp.Expr, parts: tuple[int, ...]) -> sp.Expr:
    """Coefficient of the indicated monomial in (1+X)^alpha."""
    counts = Counter(parts)
    return sp.binomial(alpha, len(parts)) * factorial(len(parts)) / sp.prod(
        factorial(multiplicity) for multiplicity in counts.values()
    )


def local_census(
    total: int,
    alpha: sp.Expr,
    orders: dict[int, int],
    h_order: int = 2,
) -> tuple[int, list[tuple[tuple[int, ...], sp.Expr]]]:
    evaluated = []
    for parts in partitions(total):
        count = len(parts)
        order = h_order * (3 * alpha - 3 * count) + sum(
            orders.get(part, 0) for part in parts
        )
        evaluated.append((order, parts, partition_coefficient(alpha, parts)))
    minimum = min(order for order, _, _ in evaluated)
    return minimum, [
        (parts, coefficient)
        for order, parts, coefficient in evaluated
        if order == minimum
    ]


def leading_residue(
    terms: list[tuple[tuple[int, ...], sp.Expr]],
    alpha: sp.Expr,
    leading_coefficients: dict[int, sp.Expr],
    h_lead: sp.Symbol,
) -> sp.Expr:
    result = sp.Integer(0)
    for parts, coefficient in terms:
        term = coefficient * h_lead ** (3 * alpha - 3 * len(parts))
        for part in parts:
            if part not in leading_coefficients:
                raise AssertionError(("unexpected part in lowest residue", part, parts))
            term *= leading_coefficients[part]
        result += term
    return sp.factor(result)


def correction_minimum(
    total: int,
    alpha: sp.Expr,
    orders: dict[int, int],
    h_order: int = 2,
) -> int:
    return local_census(total, alpha, orders, h_order)[0]


def cube_subtracted_residues(
    total: int,
    orders: dict[int, int],
    leading_coefficients: dict[int, sp.Expr],
    h_lead: sp.Symbol,
    A_lead: sp.Symbol,
) -> dict[int, sp.Expr]:
    """Lowest-coefficient Laurent map for F=(h+A*w^2)^3+E.

    For m residual factors with total w-weight L, let j=(total-L)/2.
    The exact term is
      binom(5/3,m) * multinomial * binom(5-3m,j) * h^(5-3m-j) A^j E^m.
    """
    by_order: dict[int, sp.Expr] = defaultdict(lambda: sp.Integer(0))
    for m in range(1, total // 3 + 1):
        base_power = sp.Rational(5) - 3 * m
        for parts in partitions(total, minimum_part=3):
            if len(parts) != m:
                continue
            remainder = total - sum(parts)
            if remainder < 0 or remainder % 2:
                continue
            j = remainder // 2
            if base_power >= 0 and j > base_power:
                continue
            if any(part not in leading_coefficients for part in parts):
                continue
            counts = Counter(parts)
            multinomial = factorial(m) / sp.prod(
                factorial(multiplicity) for multiplicity in counts.values()
            )
            coefficient = (
                sp.binomial(sp.Rational(5, 3), m)
                * multinomial
                * sp.binomial(base_power, j)
            )
            order = 2 * (base_power - j) + sum(orders[part] for part in parts)
            term = coefficient * h_lead ** (base_power - j) * A_lead**j
            for part in parts:
                term *= leading_coefficients[part]
            by_order[order] += term
    return {order: sp.factor(value) for order, value in by_order.items() if value != 0}


def full_k14_jet_coefficients(e7_order: int) -> tuple[sp.Expr, sp.Expr]:
    """Expand K14 through u^-1 with arbitrary first local jet corrections."""
    u = sp.symbols("u")
    c, c1, c2, c3 = sp.symbols("c c1 c2 c3", nonzero=True)
    A, A1, A2 = sp.symbols("A A1 A2", nonzero=True)
    h = c * u**2 + c1 * u**3 + c2 * u**4 + c3 * u**5
    a_series = A + A1 * u + A2 * u**2

    residual_orders = {3: 4, 4: 3, 5: 2, 6: 1, 7: e7_order}
    residuals: dict[int, sp.Expr] = {}
    first_coeffs: dict[int, sp.Symbol] = {}
    for index in range(3, 15):
        order = residual_orders.get(index, 0)
        a0, a1, a2 = sp.symbols(f"e{index}0 e{index}1 e{index}2")
        first_coeffs[index] = a0
        residuals[index] = (
            a0 * u**order + a1 * u ** (order + 1) + a2 * u ** (order + 2)
        )

    total = 14
    expression = sp.Integer(0)
    for m in range(1, total // 3 + 1):
        base_power = sp.Rational(5) - 3 * m
        for parts in partitions(total, minimum_part=3):
            if len(parts) != m:
                continue
            remainder = total - sum(parts)
            if remainder < 0 or remainder % 2:
                continue
            j = remainder // 2
            if base_power >= 0 and j > base_power:
                continue
            counts = Counter(parts)
            multinomial = factorial(m) / sp.prod(
                factorial(multiplicity) for multiplicity in counts.values()
            )
            coefficient = (
                sp.binomial(sp.Rational(5, 3), m)
                * multinomial
                * sp.binomial(base_power, j)
            )
            term = coefficient * h ** (base_power - j) * a_series**j
            for part in parts:
                term *= residuals[part]
            expression += term

    negative_series = sp.expand(sp.series(expression, u, 0, 0).removeO())
    return (
        sp.factor(negative_series.coeff(u, -2)),
        sp.factor(negative_series.coeff(u, -1)),
    )


def main() -> None:
    # The local q=4,s=2,r>=4,v>=2 profile includes U3=0 and F5=0.
    # Orders for parts >=7 are only lower bounds; a zero/higher slice removes
    # its terms and cannot create a more singular partition.
    profile = {2: 4, 3: 4, 4: 2, 5: 2, 6: 0}
    profile.update({index: 0 for index in range(7, 15)})

    expected_lowest = {
        12: {
            (2, 2, 2, 2, 2, 2),
            (2, 2, 2, 2, 4),
            (2, 2, 2, 6),
            (2, 2, 4, 4),
            (2, 4, 6),
            (4, 4, 4),
            (6, 6),
        },
        14: {
            (2, 2, 2, 2, 2, 2, 2),
            (2, 2, 2, 2, 2, 4),
            (2, 2, 2, 2, 6),
            (2, 2, 2, 4, 4),
            (2, 2, 4, 6),
            (2, 4, 4, 4),
            (2, 6, 6),
            (4, 4, 6),
        },
    }
    census = {}
    for total, expected_order in ((12, -2), (14, -4)):
        minimum, terms = local_census(total, sp.Rational(5, 3), profile)
        assert minimum == expected_order, (total, minimum)
        assert {parts for parts, _ in terms} == expected_lowest[total]
        assert all(set(parts) <= {2, 4, 6} for parts, _ in terms)
        census[str(total)] = {
            "minimum_order": int(minimum),
            "lowest_partition_count": len(terms),
            "lowest_partitions": ["+".join(map(str, parts)) for parts, _ in terms],
        }

    # Exact K12/K14 normalized selector equations.
    c, d, a, f = sp.symbols("c d a f", nonzero=True)
    leading = {2: d, 4: a, 6: f}
    low12 = local_census(12, sp.Rational(5, 3), profile)[1]
    low14 = local_census(14, sp.Rational(5, 3), profile)[1]
    raw12 = leading_residue(low12, sp.Rational(5, 3), leading, c)
    raw14 = leading_residue(low14, sp.Rational(5, 3), leading, c)

    z, phi = sp.symbols("z phi")
    e12 = (
        81 * z**3 - 162 * z**2 + 486 * z * phi + 63 * z
        - 729 * phi**2 - 108 * phi - 7
    )
    e14 = (
        324 * z**3 - 729 * z**2 * phi - 378 * z**2 + 972 * z * phi
        + 126 * z - 729 * phi**2 - 189 * phi - 13
    )
    normalized12 = sp.factor(
        raw12.subs({a: z * d**2 / c**3, f: phi * d**3 / c**6})
        * (-6561 * c**13 / (5 * d**6))
    )
    normalized14 = sp.factor(
        raw14.subs({a: z * d**2 / c**3, f: phi * d**3 / c**6})
        * (19683 * c**16 / (5 * d**7))
    )
    assert sp.expand(normalized12 - e12) == 0
    assert sp.expand(normalized14 - e14) == 0
    resultant = sp.factor(sp.resultant(e12, e14, phi))
    assert resultant == -3**15 * (3 * z - 1) ** 7
    assert sp.factor(e12.subs(z, sp.Rational(1, 3))) == -(27 * phi - 1) ** 2
    z_forced = sp.Rational(1, 3)
    phi_forced = sp.Rational(1, 27)

    # K12 and K14 first-integral resonance modes lie strictly above their
    # lowest residues. For K14, [w^9]F^(4/3) is regular and [w^4]F=U4.
    k12_kappa1_min = int(correction_minimum(
        7, sp.Rational(4, 3), profile
    ))
    k14_kappa1_min = int(correction_minimum(
        9, sp.Rational(4, 3), profile
    ))
    assert k12_kappa1_min == 2
    assert k14_kappa1_min == 0
    assert profile[2] == 4 and profile[4] == 2  # kappa_2 modes at K12/K14

    # Global divisibilities at both double roots give a polynomial cube
    # subtraction. Let A=B2/3 and C=h+A*w^2. Then F=C^3+E, with
    # E3=U3, E4=h(3C4-B2^2)/3, E5=U5/t, and E6=F6-B2^3/27.
    # The rootwise orders are ord(E3)>=4, ord(E4)>=3, ord(E5)>=2,
    # ord(E6)>=1. No t|U5 assumption is used.
    t = sp.symbols("t")
    B2, B3, C4, L, T, U5, U6, U7, U8, N = sp.symbols(
        "B2 B3 C4 L T U5 U6 U7 U8 N"
    )
    r0 = sp.symbols("r0")
    h_global = r0**2
    assert sp.expand(
        h_global * (3 * ((B2**2 + r0 * L) / 3) - B2**2) / 3
        - r0**3 * L / 3
    ) == 0
    assert sp.expand(
        ((r0 * N) / (27 * t)) - r0 * N / (27 * t)
    ) == 0
    # The first relation is r0|(3C4-B2^2). The second is the polynomialized
    # form of r0|(27F6-B2^3), with F6=U6/t and gcd(t,r0)=1.
    assert sp.expand(U5 - r0**2 * T).subs(U5, r0**2 * T) == 0

    # Exact transformed K14 Laurent coefficients. The first certificate
    # allows E7(tau) to be nonzero; the second uses its forced vanishing and
    # checks arbitrary h,A,E_i first-order jets, so no hidden u^-1 term is
    # omitted after the u^-2 coefficient is set to zero.
    k14_order_m2, _ = full_k14_jet_coefficients(e7_order=0)
    assert sp.factor(k14_order_m2 - 5 * sp.Symbol("e70") ** 2 / (9 * c)) == 0
    k14_after_e7_m2, k14_after_e7_m1 = full_k14_jet_coefficients(e7_order=1)
    e30, e40, e50, e60, e80 = sp.symbols("e30 e40 e50 e60 e80")
    expected_k14_m1 = 5 * (
        6 * c**3 * e60 * e80
        - 2 * e30 * e50 * e60
        - e40**2 * e60
        - e40 * e50**2
    ) / (27 * c**4)
    assert k14_after_e7_m2 == 0
    assert sp.factor(k14_after_e7_m1 - expected_k14_m1) == 0

    # K14's first-integral corrections are regular at tau, so its u^-2
    # coefficient forces F7(tau)=0. At both roots this gives r0|U7.
    # K14's next u^-1 coefficient is the equation shown below.
    k14_next_equation = sp.expand(
        6 * c**3 * e60 * e80
        - 2 * e30 * e50 * e60
        - e40**2 * e60
        - e40 * e50**2
    )

    # K15 after the same cube subtraction. E7 is now at least order one.
    # The kappa_1 coefficient [w^10]F^(4/3) has order -2, but its exact
    # leading residue vanishes after z=1/3, phi=1/27. The other K15
    # resonance modes have orders >=2 and 4.
    corr10_orders = dict(profile)
    corr10_min = correction_minimum(10, sp.Rational(4, 3), corr10_orders)
    assert corr10_min == -2
    low10_corr = local_census(10, sp.Rational(4, 3), corr10_orders)[1]
    raw10_corr = leading_residue(low10_corr, sp.Rational(4, 3), leading, c)
    normalized10_corr = sp.factor(
        raw10_corr.subs({a: z * d**2 / c**3, f: phi * d**3 / c**6})
    )
    expected10_corr = -4 * d**5 * (3 * z - 1) * (-27 * phi + 9 * z - 2) / (
        729 * c**11
    )
    assert sp.factor(normalized10_corr - expected10_corr) == 0
    assert normalized10_corr.subs(
        {z: z_forced, phi: phi_forced}
    ) == 0

    cube15_orders = {3: 4, 4: 3, 5: 2, 6: 1, 7: 1}
    cube15_orders.update({index: 0 for index in range(8, 16)})
    e3, e4, e5, e6, eta7 = sp.symbols("e3 e4 e5 e6 eta7")
    x8_15, x9_15 = sp.symbols("x8_15 x9_15")
    cube15_coeffs = {3: e3, 4: e4, 5: e5, 6: e6, 7: eta7}
    cube15_coeffs.update({index: sp.Symbol(f"x{index}_15") for index in range(8, 16)})
    cube15 = cube_subtracted_residues(
        15, cube15_orders, cube15_coeffs, c, sp.Symbol("A0", nonzero=True)
    )
    expected_k15_m2 = -5 * (3 * e3 * e6**2 + 6 * e4 * e5 * e6 + e5**3) / (
        81 * c**4
    )
    assert sp.factor(cube15[-2] - expected_k15_m2) == 0
    k15_equation = 3 * e3 * e6**2 + 6 * e4 * e5 * e6 + e5**3
    assert sp.factor(k15_equation.subs(e5, 0)) == 3 * e3 * e6**2
    assert sp.factor(k15_equation.subs(e6, 0)) == e5**3

    # Convert the root equations to polynomial congruences using
    # h=r0^2, c=lambda^2, lambda=r0'(tau), and the global quotients
    # U3=h^2*B3, U5=r0^2*T,
    # 27U6-t*B2^3=r0*N, 3C4-B2^2=r0*L.
    lam, tau = sp.symbols("lam tau", nonzero=True)
    Bt, Lt, Tt, Nt, U8t = sp.symbols("Bt Lt Tt Nt U8t")
    e3_global = lam**4 * Bt
    e4_global = lam**3 * Lt / 3
    e5_global = lam**2 * Tt / tau
    e6_global = lam * Nt / (27 * tau)
    c_global = lam**2
    k15_global = sp.factor(
        k15_equation.subs(
            {e3: e3_global, e4: e4_global, e5: e5_global, e6: e6_global}
        )
        - lam**6
        * (tau * Bt * Nt**2 + 18 * tau * Lt * Tt * Nt + 243 * Tt**3)
        / (243 * tau**3)
    )
    assert k15_global == 0
    k14_global_expr = (
        6 * c**3 * e6 * sp.Symbol("x8")
        - 2 * e3 * e5 * e6
        - e4**2 * e6
        - e4 * e5**2
    )
    k14_global = sp.factor(
        k14_global_expr.subs(
            {
                c: c_global,
                e3: e3_global,
                e4: e4_global,
                e5: e5_global,
                e6: e6_global,
                sp.Symbol("x8"): U8t / tau,
            }
        )
        - lam**7
        * (
            54 * Nt * U8t
            - 18 * Bt * Tt * Nt
            - tau * Lt**2 * Nt
            - 81 * Lt * Tt**2
        )
        / (243 * tau**2)
    )
    assert k14_global == 0

    print(
        json.dumps(
            {
                "selector": {
                    "K12_equation": "81*z^3-162*z^2+486*z*phi+63*z-729*phi^2-108*phi-7=0",
                    "K14_equation": "324*z^3-729*z^2*phi-378*z^2+972*z*phi+126*z-729*phi^2-189*phi-13=0",
                    "resultant_phi": str(resultant),
                    "forced": {"z": "1/3", "phi": "1/27"},
                    "f0_nonzero_is_forced": True,
                    "lowest_partition_census": census,
                    "K12_kappa1_minimum": k12_kappa1_min,
                    "K14_kappa1_minimum": k14_kappa1_min,
                },
                "guards": {
                    "h_leading_c_nonzero": True,
                    "U2_leading_d_nonzero": True,
                    "U4_leading_a_nonzero": True,
                    "U3_zero_or_higher_order_allowed": True,
                    "U5_over_t_zero_or_higher_order_allowed": True,
                    "t_is_a_unit_at_both_roots": True,
                    "no_assumption_t_divides_U5": True,
                    "no_division_by_U3_or_U5_or_U7": True,
                },
                "global_cube_subtraction": {
                    "A": "B2/3",
                    "E4": "h*(3*C4-B2^2)/3 = r0^3*L/3",
                    "E6": "F6-B2^3/27; polynomialized as 27*U6-t*B2^3 = r0*N when F6=U6/t",
                    "root_orders": {"E3": ">=4", "E4": ">=3", "E5=U5/t": ">=2", "E6": ">=1"},
                    "base_divisibilities": [
                        "r0 | (3*C4-B2^2)",
                        "r0 | (27*F6-B2^3) in the localization at t != 0",
                        "r0 | (27*U6-t*B2^3) in K[t]",
                        "U3=h^2*B3",
                        "U5=r0^2*T, without requiring t|U5",
                    ],
                },
                "post_cube": {
                    "K14_order_minus_2": "5*F7(tau)^2/(9*c)",
                    "forced_at_each_double_root": "F7(tau)=0; hence r0|U7",
                    "K14_order_minus_1_after_F7_zero": "5*(6*c^3*e6*F8-2*e3*e5*e6-e4^2*e6-e4*e5^2)/(27*c^4)=0",
                    "K15_order_minus_2": "3*e3*e6^2+6*e4*e5*e6+e5^3=0",
                    "K15_kappa1_order_minus_2": "0 after z=1/3, phi=1/27",
                    "zero_cases": {
                        "e6=0": "forces e5=0; leaves e3 unconstrained at K15",
                        "e5=0": "forces e3*e6^2=0",
                        "e5_nonzero": "forces e6 nonzero; no division by e3 or e4 is valid globally",
                    },
                },
                "globalized_next_constraints": {
                    "quotients": {
                        "L": "(3*C4-B2^2)/r0",
                        "T": "U5/r0^2",
                        "N": "(27*U6-t*B2^3)/r0",
                    },
                    "from_K15": "r0 | (t*B3*N^2 + 18*t*L*T*N + 243*T^3)",
                    "from_K14_next_residue": "r0 | (54*N*U8 - 18*B3*T*N - t*L^2*N - 81*L*T^2)",
                    "degree_pair_excluded": False,
                },
                "scope": "exact necessary conditions for the 2+2 q=4,s=2,r>=4,v>=2 local branch; not an exclusion of the branch or of (75,125)",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
