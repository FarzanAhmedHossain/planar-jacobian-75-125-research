"""Exact K6/K7 double-root jet classification for the F2 (75,125) bands.

The checker uses only Python's standard library. It replays the weighted
partition valuations, verifies the local residue polynomials, and computes the
quadratic resultant that isolates the exceptional (r,q)=(4,2) jet. It checks
these finite arithmetic steps; the band-model derivation and its upstream
Jacobian-pair reduction remain separate proof obligations.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import permutations
from math import factorial


Multiplicity = tuple[int, ...]
Polynomial = dict[int, Fraction]  # exponent -> coefficient, in Y
MultiPolynomial = dict[tuple[int, ...], Fraction]


def weighted_partitions(total: int, index: int = 1) -> list[Multiplicity]:
    """Return all tuples m with sum(j*m_j)=total."""
    if index > total:
        return [()] if total == 0 else []
    result: list[Multiplicity] = []
    for multiplicity in range(total // index + 1):
        for tail in weighted_partitions(total - index * multiplicity, index + 1):
            result.append((multiplicity,) + tail)
    return result


def pad(m: Multiplicity, total: int) -> Multiplicity:
    return m + (0,) * (total - len(m))


def binomial_5_over_3(k: int) -> Fraction:
    value = Fraction(1)
    for j in range(k):
        value *= Fraction(5, 3) - j
    return value / factorial(k)


def partition_coefficient(m: Multiplicity) -> Fraction:
    k = sum(m)
    value = binomial_5_over_3(k) * factorial(k)
    for multiplicity in m:
        value /= factorial(multiplicity)
    return value


def valuation(m: Multiplicity, orders: tuple[int, ...], h_order: int = 2) -> int:
    k = sum(m)
    return 5 * h_order - 3 * h_order * k + sum(a * b for a, b in zip(m, orders))


def least_terms(
    total: int, orders: tuple[int, ...], h_order: int = 2
) -> tuple[int, list[tuple[Multiplicity, Fraction]]]:
    census = [
        (
            valuation(pad(m, total), orders, h_order),
            pad(m, total),
            partition_coefficient(pad(m, total)),
        )
        for m in weighted_partitions(total)
    ]
    lowest = min(v for v, _, _ in census)
    return lowest, [(m, c) for v, m, c in census if v == lowest]


def normalized_residue(
    terms: list[tuple[Multiplicity, Fraction]], scale: Fraction, expected_h_power: int
) -> dict[tuple[int, int, int], Fraction]:
    """Substitute A=X*c^2, B=Y*c, C=Z and remove c^expected_h_power."""
    out: dict[tuple[int, int, int], Fraction] = {}
    for m, coefficient in terms:
        assert not any(m[3:])
        h_power = 5 - 3 * sum(m) + 2 * m[0] + m[1]
        assert h_power == expected_h_power
        monomial = m[:3]
        out[monomial] = out.get(monomial, Fraction(0)) + scale * coefficient
    return {m: c for m, c in out.items() if c}


def p_add(a: Polynomial, b: Polynomial) -> Polynomial:
    out = dict(a)
    for e, c in b.items():
        out[e] = out.get(e, Fraction(0)) + c
    return {e: c for e, c in out.items() if c}


def p_scale(a: Polynomial, c: Fraction) -> Polynomial:
    return {e: value * c for e, value in a.items() if value * c}


def p_mul(a: Polynomial, b: Polynomial) -> Polynomial:
    out: Polynomial = {}
    for e, c in a.items():
        for f, d in b.items():
            out[e + f] = out.get(e + f, Fraction(0)) + c * d
    return {e: c for e, c in out.items() if c}


def p_pow(a: Polynomial, n: int) -> Polynomial:
    out: Polynomial = {0: Fraction(1)}
    for _ in range(n):
        out = p_mul(out, a)
    return out


def determinant4(matrix: list[list[Polynomial]]) -> Polynomial:
    out: Polynomial = {}
    for permutation in permutations(range(4)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(4)
            for j in range(i + 1, 4)
        )
        term: Polynomial = {0: Fraction(-1 if inversions % 2 else 1)}
        for row, column in enumerate(permutation):
            term = p_mul(term, matrix[row][column])
        out = p_add(out, term)
    return out


def specialize_y(coefficients: list[Polynomial], y: Fraction) -> list[Fraction]:
    return [sum(c * y**e for e, c in p.items()) for p in coefficients]


def mp_add(a: MultiPolynomial, b: MultiPolynomial) -> MultiPolynomial:
    out = dict(a)
    for monomial, coefficient in b.items():
        out[monomial] = out.get(monomial, Fraction(0)) + coefficient
    return {m: c for m, c in out.items() if c}


def mp_scale(a: MultiPolynomial, coefficient: Fraction) -> MultiPolynomial:
    return {m: c * coefficient for m, c in a.items() if c * coefficient}


def mp_mul(a: MultiPolynomial, b: MultiPolynomial) -> MultiPolynomial:
    out: MultiPolynomial = {}
    for m, c in a.items():
        for n, d in b.items():
            monomial = tuple(x + y for x, y in zip(m, n))
            out[monomial] = out.get(monomial, Fraction(0)) + c * d
    return {m: c for m, c in out.items() if c}


def mp_pow(a: MultiPolynomial, exponent: int) -> MultiPolynomial:
    dimension = len(next(iter(a))) if a else 6
    out: MultiPolynomial = {(0,) * dimension: Fraction(1)}
    for _ in range(exponent):
        out = mp_mul(out, a)
    return out


def mp_var(index: int, dimension: int = 6) -> MultiPolynomial:
    monomial = [0] * dimension
    monomial[index] = 1
    return {tuple(monomial): Fraction(1)}


def mp_substitute(
    polynomial: MultiPolynomial, substitutions: tuple[MultiPolynomial, ...]
) -> MultiPolynomial:
    out: MultiPolynomial = {}
    for monomial, coefficient in polynomial.items():
        term: MultiPolynomial = {(0,) * len(substitutions): coefficient}
        for variable, exponent in enumerate(monomial):
            term = mp_mul(term, mp_pow(substitutions[variable], exponent))
        out = mp_add(out, term)
    return out


def partition_pole_map(total: int, h_power_target: int) -> MultiPolynomial:
    """Pole numerator terms after U1=h^2*A, U2=h*B; U3=C, U4=U4."""
    out: MultiPolynomial = {}
    for raw in weighted_partitions(total):
        m = pad(raw, total)
        k = sum(m)
        h_power = 5 - 3 * k + 2 * m[0] + m[1]
        if h_power != h_power_target:
            continue
        assert not any(m[4:]), "higher bands must not enter the negative h-pole part"
        monomial = (m[0], m[1], m[2], m[3], 0, 0)
        out = mp_add(out, {monomial: partition_coefficient(m)})
    return out


def main() -> None:
    # In the binomial coefficient [w^s]F^(5/3), a partition m has order
    # 10 - 6*sum(m_j) + sum(m_j*ord(U_j)) at a double root of h.
    orders_42 = (4, 2, 0, 0, 0, 0, 0)
    k6_min, k6_terms = least_terms(6, orders_42)
    k7_min, k7_terms = least_terms(7, orders_42)
    assert k6_min == -2
    assert k7_min == -4
    expected_k6 = {
        (0, 0, 2, 0, 0, 0): Fraction(5, 9),
        (0, 3, 0, 0, 0, 0): Fraction(-5, 81),
        (1, 1, 1, 0, 0, 0): Fraction(-10, 27),
        (2, 2, 0, 0, 0, 0): Fraction(10, 81),
        (3, 0, 1, 0, 0, 0): Fraction(20, 243),
        (4, 1, 0, 0, 0, 0): Fraction(-35, 729),
        (6, 0, 0, 0, 0, 0): Fraction(35, 6561),
    }
    expected_k7 = {
        (0, 2, 1, 0, 0, 0, 0): Fraction(-5, 27),
        (1, 0, 2, 0, 0, 0, 0): Fraction(-5, 27),
        (1, 3, 0, 0, 0, 0, 0): Fraction(20, 243),
        (2, 1, 1, 0, 0, 0, 0): Fraction(20, 81),
        (3, 2, 0, 0, 0, 0, 0): Fraction(-70, 729),
        (4, 0, 1, 0, 0, 0, 0): Fraction(-35, 729),
        (5, 1, 0, 0, 0, 0, 0): Fraction(70, 2187),
        (7, 0, 0, 0, 0, 0, 0): Fraction(-65, 19683),
    }
    assert dict(k6_terms) == expected_k6
    assert dict(k7_terms) == expected_k7

    # The cleared local residue equations for r=4, q=2.
    # X=A/c^2, Y=B/c, Z=U_3(alpha), where
    # h=c*u^2+..., U_1=A*u^4+..., U_2=B*u^2+....
    k6 = {
        (0, 0, 2): 729,
        (0, 3, 0): -81,
        (1, 1, 1): -486,
        (2, 2, 0): 162,
        (3, 0, 1): 108,
        (4, 1, 0): -63,
        (6, 0, 0): 7,
    }
    k7 = {
        (0, 2, 1): 729,
        (1, 0, 2): 729,
        (1, 3, 0): -324,
        (2, 1, 1): -972,
        (3, 2, 0): 378,
        (4, 0, 1): 189,
        (5, 1, 0): -126,
        (7, 0, 0): 13,
    }
    assert normalized_residue(k6_terms, Fraction(6561, 5), -1) == {
        (0, 0, 2): Fraction(729),
        (0, 3, 0): Fraction(-81),
        (1, 1, 1): Fraction(-486),
        (2, 2, 0): Fraction(162),
        (3, 0, 1): Fraction(108),
        (4, 1, 0): Fraction(-63),
        (6, 0, 0): Fraction(7),
    }
    assert normalized_residue(k7_terms, Fraction(-19683, 5), -2) == {
        (0, 2, 1): Fraction(729),
        (1, 0, 2): Fraction(729),
        (1, 3, 0): Fraction(-324),
        (2, 1, 1): Fraction(-972),
        (3, 2, 0): Fraction(378),
        (4, 0, 1): Fraction(189),
        (5, 1, 0): Fraction(-126),
        (7, 0, 0): Fraction(13),
    }
    assert all(sum(e * w for e, w in zip(m, (1, 2, 3))) == 6 for m in k6)
    assert all(sum(e * w for e, w in zip(m, (1, 2, 3))) == 7 for m in k7)

    # Substitution of R7=kappa*(-sigma/3+2*rho^2/9)/S in the w^7
    # coefficient equation cancels these four exact rational coefficients.
    assert Fraction(-1, 3) + Fraction(5, 15) == 0
    assert Fraction(4, 9) - Fraction(20, 45) == 0
    assert Fraction(2, 5) * Fraction(-1, 3) + Fraction(2, 15) == 0
    assert Fraction(2, 5) * Fraction(2, 9) - Fraction(4, 45) == 0
    # The coefficient correction G7-[w^7]F^(5/3) has the stated factors.
    assert Fraction(-1, 3) + Fraction(5, 3) == Fraction(4, 3)
    assert Fraction(2, 9) + Fraction(5, 9) - Fraction(5, 9) == Fraction(2, 9)

    # Exact pole decomposition of the full binomial coefficients after
    # U1=h^2*A and U2=h*B. K6 has only a h^-1 pole; K7 has h^-2 and h^-1.
    assert partition_pole_map(6, -1) == {
        (0, 0, 2, 0, 0, 0): Fraction(5, 9),
        (0, 3, 0, 0, 0, 0): Fraction(-5, 81),
        (1, 1, 1, 0, 0, 0): Fraction(-10, 27),
        (2, 2, 0, 0, 0, 0): Fraction(10, 81),
        (3, 0, 1, 0, 0, 0): Fraction(20, 243),
        (4, 1, 0, 0, 0, 0): Fraction(-35, 729),
        (6, 0, 0, 0, 0, 0): Fraction(35, 6561),
    }
    assert partition_pole_map(7, -2) == {
        (0, 2, 1, 0, 0, 0): Fraction(-5, 27),
        (1, 0, 2, 0, 0, 0): Fraction(-5, 27),
        (1, 3, 0, 0, 0, 0): Fraction(20, 243),
        (2, 1, 1, 0, 0, 0): Fraction(20, 81),
        (3, 2, 0, 0, 0, 0): Fraction(-70, 729),
        (4, 0, 1, 0, 0, 0): Fraction(-35, 729),
        (5, 1, 0, 0, 0, 0): Fraction(70, 2187),
        (7, 0, 0, 0, 0, 0): Fraction(-65, 19683),
    }
    assert partition_pole_map(7, -1) == {
        (0, 0, 1, 1, 0, 0): Fraction(10, 9),
        (1, 1, 0, 1, 0, 0): Fraction(-10, 27),
        (3, 0, 0, 1, 0, 0): Fraction(20, 243),
    }

    # Verify the compressed D,E forms exactly, with variables ordered
    # (A,B,C,U4,D,E).
    A, B, C, U4, D, E = (mp_var(i) for i in range(6))
    raw_w6 = {
        (0, 0, 2, 0, 0, 0): Fraction(5, 9),
        (0, 3, 0, 0, 0, 0): Fraction(-5, 81),
        (1, 1, 1, 0, 0, 0): Fraction(-10, 27),
        (2, 2, 0, 0, 0, 0): Fraction(10, 81),
        (3, 0, 1, 0, 0, 0): Fraction(20, 243),
        (4, 1, 0, 0, 0, 0): Fraction(-35, 729),
        (6, 0, 0, 0, 0, 0): Fraction(35, 6561),
    }
    raw_w7 = {
        (0, 2, 1, 0, 0, 0): Fraction(-5, 27),
        (1, 0, 2, 0, 0, 0): Fraction(-5, 27),
        (1, 3, 0, 0, 0, 0): Fraction(20, 243),
        (2, 1, 1, 0, 0, 0): Fraction(20, 81),
        (3, 2, 0, 0, 0, 0): Fraction(-70, 729),
        (4, 0, 1, 0, 0, 0): Fraction(-35, 729),
        (5, 1, 0, 0, 0, 0): Fraction(70, 2187),
        (7, 0, 0, 0, 0, 0): Fraction(-65, 19683),
    }
    raw_w71 = {
        (0, 0, 1, 1, 0, 0): Fraction(10, 9),
        (1, 1, 0, 1, 0, 0): Fraction(-10, 27),
        (3, 0, 0, 1, 0, 0): Fraction(20, 243),
    }
    sub = (A, mp_scale(mp_add(mp_pow(A, 2), D), Fraction(1, 3)),
           mp_scale(mp_add(mp_pow(A, 3), E), Fraction(1, 27)), U4, D, E)
    hat_w6 = mp_scale(mp_substitute(raw_w6, sub), Fraction(6561, 5))
    hat_w7 = mp_scale(mp_substitute(raw_w7, sub), Fraction(-19683, 5))
    hat_w71 = mp_scale(mp_substitute(raw_w71, sub), Fraction(-19683, 5))
    delta = mp_add(E, mp_scale(mp_mul(A, D), Fraction(-3)))
    expected_w6 = mp_add(
        mp_add(mp_pow(delta, 2), mp_scale(mp_pow(D, 3), Fraction(-3))), {}
    )
    expected_w7 = mp_add(
        mp_mul(A, expected_w6), mp_scale(mp_mul(mp_pow(D, 2), delta), Fraction(3))
    )
    expected_w71 = mp_scale(mp_mul(U4, delta), Fraction(-162))
    assert hat_w6 == expected_w6
    assert hat_w7 == expected_w7
    assert hat_w71 == expected_w71

    # For X != 0, weighted homogeneity permits X=1. Compute the exact
    # Sylvester resultant in Z of the two quadratics over Q[Y].
    f2 = {0: {0: Fraction(7), 1: Fraction(-63), 2: Fraction(162), 3: Fraction(-81)},
          1: {0: Fraction(108), 1: Fraction(-486)},
          2: {0: Fraction(729)}}
    g2 = {0: {0: Fraction(13), 1: Fraction(-126), 2: Fraction(378), 3: Fraction(-324)},
          1: {0: Fraction(189), 1: Fraction(-972), 2: Fraction(729)},
          2: {0: Fraction(729)}}
    zero: Polynomial = {}
    sylvester = [
        [f2[2], f2[1], f2[0], zero],
        [zero, f2[2], f2[1], f2[0]],
        [g2[2], g2[1], g2[0], zero],
        [zero, g2[2], g2[1], g2[0]],
    ]
    resultant = determinant4(sylvester)
    expected_resultant = p_scale(
        p_pow({0: Fraction(-1), 1: Fraction(3)}, 7), Fraction(-(3**15))
    )
    assert resultant == expected_resultant

    # Its only root is Y=1/3; both quadratics then become (27Z-1)^2.
    y0 = Fraction(1, 3)
    f_at_y0 = specialize_y([f2[0], f2[1], f2[2]], y0)
    g_at_y0 = specialize_y([g2[0], g2[1], g2[2]], y0)
    square = [Fraction(1), Fraction(-54), Fraction(729)]
    assert f_at_y0 == square and g_at_y0 == square

    # The other contact strata are disposed of by their exact lowest terms.
    r4_q3 = (4, 3, 0, 0, 0, 0, 0)
    assert least_terms(6, r4_q3) == (
        -2,
        [
            ((0, 0, 2, 0, 0, 0), Fraction(5, 9)),
            ((3, 0, 1, 0, 0, 0), Fraction(20, 243)),
            ((6, 0, 0, 0, 0, 0), Fraction(35, 6561)),
        ],
    )
    assert least_terms(7, r4_q3) == (
        -4,
        [
            ((1, 0, 2, 0, 0, 0, 0), Fraction(-5, 27)),
            ((4, 0, 1, 0, 0, 0, 0), Fraction(-35, 729)),
            ((7, 0, 0, 0, 0, 0, 0), Fraction(-65, 19683)),
        ],
    )
    # After X=1 these give 729Z^2+108Z+7=0 and
    # 729Z^2+189Z+13=0; their difference forces Z=-2/27,
    # which leaves 3 in the first equation.
    assert Fraction(729) * Fraction(-2, 27) ** 2 + 108 * Fraction(-2, 27) + 7 == 3

    r5_q2 = (5, 2, 0, 0, 0, 0, 0)
    assert least_terms(6, r5_q2) == (
        -2,
        [
            ((0, 0, 2, 0, 0, 0), Fraction(5, 9)),
            ((0, 3, 0, 0, 0, 0), Fraction(-5, 81)),
        ],
    )
    assert least_terms(7, r5_q2) == (
        -4,
        [((0, 2, 1, 0, 0, 0, 0), Fraction(-5, 27))],
    )
    r5_q3 = (5, 3, 0, 0, 0, 0, 0)
    assert least_terms(6, r5_q3) == (
        -2,
        [((0, 0, 2, 0, 0, 0), Fraction(5, 9))],
    )

    # At a simple root e=1 the minimal orders are (r,q)=(2,1). The
    # valuations are exactly the double-root (4,2) valuations divided by 2,
    # so they produce the same normalized K6/K7 residue equations.
    r2_q1 = (2, 1, 0, 0, 0, 0, 0)
    assert least_terms(6, r2_q1, h_order=1) == (-1, k6_terms)
    assert least_terms(7, r2_q1, h_order=1) == (-2, k7_terms)
    assert normalized_residue(k6_terms, Fraction(6561, 5), -1) == k6
    assert normalized_residue(k7_terms, Fraction(-19683, 5), -2) == k7

    # The simple-root mixed contacts (r=2,q>=2) and (r>=3,q=1)
    # inherit the same incompatible residue pairs as the double-root cases.
    r2_q2 = (2, 2, 0, 0, 0, 0, 0)
    r3_q1 = (3, 1, 0, 0, 0, 0, 0)
    r3_q2 = (3, 2, 0, 0, 0, 0, 0)
    assert least_terms(6, r2_q2, h_order=1) == (
        -1,
        [
            ((0, 0, 2, 0, 0, 0), Fraction(5, 9)),
            ((3, 0, 1, 0, 0, 0), Fraction(20, 243)),
            ((6, 0, 0, 0, 0, 0), Fraction(35, 6561)),
        ],
    )
    assert least_terms(7, r2_q2, h_order=1) == (
        -2,
        [
            ((1, 0, 2, 0, 0, 0, 0), Fraction(-5, 27)),
            ((4, 0, 1, 0, 0, 0, 0), Fraction(-35, 729)),
            ((7, 0, 0, 0, 0, 0, 0), Fraction(-65, 19683)),
        ],
    )
    assert least_terms(6, r3_q1, h_order=1) == (
        -1,
        [
            ((0, 0, 2, 0, 0, 0), Fraction(5, 9)),
            ((0, 3, 0, 0, 0, 0), Fraction(-5, 81)),
        ],
    )
    assert least_terms(7, r3_q1, h_order=1) == (
        -2,
        [((0, 2, 1, 0, 0, 0, 0), Fraction(-5, 27))],
    )
    assert least_terms(6, r3_q2, h_order=1) == (
        -1,
        [((0, 0, 2, 0, 0, 0), Fraction(5, 9))],
    )

    print("K6/K7 exact partition censuses: PASS")
    print("r=4, q>=3 is impossible; r>=5, q=2 is impossible")
    print("r=4, q=2 forces Y=X^2/3 and U3(alpha)=X^3/27")
    print("r>=5, q>=3 forces U3(alpha)=0")
    print("the same exhaustive split holds at simple roots with r>=2, q>=1")
    print("full K6/K7 pole numerators and their D,E compression: PASS")


if __name__ == "__main__":
    main()
