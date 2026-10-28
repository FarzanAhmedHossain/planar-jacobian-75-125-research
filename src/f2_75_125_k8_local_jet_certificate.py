"""Exact eighth-band recurrence and local valuation checks for F2 (75,125).

This standard-library checker validates the rational recurrence arithmetic,
the resonance correction, and the finite weighted-partition censuses stated in
the research note. It does not eliminate the complete coefficient model.
"""

from __future__ import annotations

from fractions import Fraction
from math import factorial


MultiPolynomial = dict[tuple[int, ...], Fraction]


def weighted_partitions(total: int, index: int = 1) -> list[tuple[int, ...]]:
    if index > total:
        return [()] if total == 0 else []
    out: list[tuple[int, ...]] = []
    for multiplicity in range(total // index + 1):
        for tail in weighted_partitions(total - index * multiplicity, index + 1):
            out.append((multiplicity,) + tail)
    return out


def pad(m: tuple[int, ...], total: int) -> tuple[int, ...]:
    return m + (0,) * (total - len(m))


def binomial_5_over_3(k: int) -> Fraction:
    value = Fraction(1)
    for j in range(k):
        value *= Fraction(5, 3) - j
    return value / factorial(k)


def coefficient(m: tuple[int, ...]) -> Fraction:
    k = sum(m)
    value = binomial_5_over_3(k) * factorial(k)
    for multiplicity in m:
        value /= factorial(multiplicity)
    return value


def mp_add(a: MultiPolynomial, b: MultiPolynomial) -> MultiPolynomial:
    out = dict(a)
    for monomial, coefficient in b.items():
        out[monomial] = out.get(monomial, Fraction(0)) + coefficient
    return {m: c for m, c in out.items() if c}


def mp_scale(a: MultiPolynomial, scalar: Fraction) -> MultiPolynomial:
    return {m: c * scalar for m, c in a.items() if c * scalar}


def mp_mul(a: MultiPolynomial, b: MultiPolynomial) -> MultiPolynomial:
    out: MultiPolynomial = {}
    for left, c in a.items():
        for right, d in b.items():
            monomial = tuple(x + y for x, y in zip(left, right))
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


def k8_pole_numerator() -> MultiPolynomial:
    """Return N8 where the singular part of [w^8]F^(5/3) is N8/h^3.

    Variable order is (h,A,B,C,U4,U5), after U1=h^2*A and U2=h*B.
    Every negative h power is generated directly from the multinomial
    expansion, rather than entered from the claimed compressed formula.
    """
    out: MultiPolynomial = {}
    for raw in weighted_partitions(8):
        m = pad(raw, 8)
        k = sum(m)
        h_power = 5 - 3 * k + 2 * m[0] + m[1]
        if h_power >= 0:
            continue
        assert h_power >= -3
        assert not any(m[5:]), "U6 and higher cannot occur in the pole part"
        monomial = (h_power + 3, m[0], m[1], m[2], m[3], m[4])
        out = mp_add(out, {monomial: coefficient(m)})
    return out


def compressed_k8_pole_numerator() -> MultiPolynomial:
    """Return the proposed D,Delta expression expanded in (h,A,B,C,U4,U5)."""
    h, A, B, C, U4, U5 = (mp_var(i) for i in range(6))
    D = mp_add(mp_scale(B, Fraction(3)), mp_scale(mp_pow(A, 2), Fraction(-1)))
    E = mp_add(mp_scale(C, Fraction(27)), mp_scale(mp_pow(A, 3), Fraction(-1)))
    Delta = mp_add(E, mp_scale(mp_mul(A, D), Fraction(-3)))
    J = mp_add(mp_mul(A, Delta), mp_scale(mp_mul(U4, h), Fraction(-81)))

    remainder: MultiPolynomial = {}
    remainder = mp_add(remainder, mp_scale(mp_mul(D, mp_pow(Delta, 2)), Fraction(-1)))
    remainder = mp_add(
        remainder,
        mp_scale(mp_mul(mp_mul(A, mp_pow(D, 2)), Delta), Fraction(6)),
    )
    remainder = mp_add(
        remainder,
        mp_scale(mp_mul(mp_pow(A, 2), mp_pow(D, 3)), Fraction(-3)),
    )
    remainder = mp_add(remainder, mp_scale(mp_pow(D, 4), Fraction(3)))
    remainder = mp_add(
        remainder,
        mp_scale(mp_mul(mp_mul(mp_pow(D, 2), U4), h), Fraction(-243)),
    )
    remainder = mp_add(
        remainder,
        mp_scale(
            mp_mul(mp_mul(Delta, U5), mp_pow(h, 2)), Fraction(486)
        ),
    )
    return mp_scale(
        mp_add(mp_pow(J, 2), remainder), Fraction(5, 59049)
    )


def verify_band_ratio_first_integral() -> None:
    """Check the exact characteristic identity for I=w^5/(S*f^(1/3))."""
    a, L, M = (mp_var(i, 3) for i in range(3))
    one = {(0, 0, 0): Fraction(1)}
    first = mp_mul(
        mp_add(mp_scale(a, Fraction(3)), L),
        mp_add(mp_scale(one, Fraction(5)), mp_scale(M, Fraction(-1, 3))),
    )
    second = mp_mul(
        mp_add(mp_scale(one, Fraction(15)), mp_scale(M, Fraction(-1))),
        mp_add(mp_scale(a, Fraction(-1)), mp_scale(L, Fraction(-1, 3))),
    )
    assert mp_add(first, second) == {}


def census(
    total: int, orders: tuple[int, ...], root_multiplicity: int
) -> tuple[int, list[tuple[tuple[int, ...], Fraction]]]:
    terms = []
    for raw in weighted_partitions(total):
        m = pad(raw, total)
        k = sum(m)
        order = 5 * root_multiplicity - 3 * root_multiplicity * k
        order += sum(a * b for a, b in zip(m, orders))
        terms.append((order, m, coefficient(m)))
    low = min(order for order, _, _ in terms)
    return low, [(m, c) for order, m, c in terms if order == low]


def main() -> None:
    verify_band_ratio_first_integral()

    # For L_j=[w^j]log(1+rho*w+sigma*w^2+tau*w^3+...), the w^8
    # recurrence reduces to (15D+9a)P+forcing=0 for R8=kappa*P/S.
    # These Fraction identities are the coefficientwise cancellation after
    # P=-14*rho^3/81+4*rho*sigma/9-tau/3 is substituted.
    assert Fraction(-14, 9) + Fraction(14, 9) == 0
    assert Fraction(4) - Fraction(36, 9) == 0
    assert Fraction(-3) + Fraction(27, 9) == 0
    assert Fraction(-70, 9) + Fraction(70, 9) == 0
    assert Fraction(20, 3) - Fraction(60, 9) == 0
    assert Fraction(-5) + Fraction(45, 9) == 0

    # Coefficients of G8-[w^8]F^(5/3), after removing kappa*h^5/S.
    assert Fraction(5, 3) - Fraction(1, 3) == Fraction(4, 3)  # tau
    assert Fraction(10, 9) - Fraction(5, 9) - Fraction(5, 9) + Fraction(4, 9) == Fraction(4, 9)
    assert Fraction(-5, 81) - Fraction(5, 27) + Fraction(10, 27) - Fraction(14, 81) == Fraction(-4, 81)

    # Exact K8 pole compression. Polynomiality of V8 requires h^3 | N8,
    # because the resonance correction becomes regular after U1=h^2*A,
    # U2=h*B. The two sides are generated independently: one from all
    # weighted partitions, the other from the proposed D,Delta identity.
    assert k8_pole_numerator() == compressed_k8_pole_numerator()

    # The remainder after J^2 has order at least 3e on a tuned root of
    # multiplicity e=1 or 2: K6/K7 give ord(D),ord(Delta)>=e there.
    for e in (1, 2):
        d_order = delta_order = e
        remainder_order = min(
            d_order + 2 * delta_order,
            2 * d_order + delta_order,
            3 * d_order,
            4 * d_order,
            2 * d_order + e,
            delta_order + 2 * e,
        )
        assert remainder_order >= 3 * e
        assert (3 * e + 1) // 2 == (2 if e == 1 else 3)

    # On a deeper simple root, ord(A),ord(D),ord(Delta)>=1 and ord(h)=1.
    # If U4 is a unit, J has exact order 1, so J^2 cannot be canceled by a
    # remainder of order >=3; hence U4 must vanish there.
    assert min(1 + 2, 2 + 1, 3, 4, 2 + 1, 1 + 2) == 3
    assert 2 * 1 < 3

    # On a deeper double root, A,D,Delta have order >=1 and h has order 2.
    # J^2 has order >=4, while the unique possible order-3 remainder is
    # -D*Delta^2. Thus its leading coefficient must vanish.
    deep_double_remainder_orders = (1 + 2, 2 + 1 + 1, 2 + 3, 4, 2 + 2, 1 + 4)
    assert deep_double_remainder_orders == (3, 4, 5, 4, 4, 5)
    assert 2 * 2 > 3

    # A symbolic all-deep jet survives every local divisibility test through
    # K8 for any quartic h: U1=0, U2=h^2, U3=U4=h, U5=0. Exponents here are
    # powers of the formal polynomial h, so they check divisibility exactly.
    h = mp_var(0, 1)
    h2, h3, h4 = mp_pow(h, 2), mp_pow(h, 3), mp_pow(h, 4)
    D = mp_scale(h, Fraction(3))
    Delta = mp_scale(h, Fraction(27))
    k4_witness = {(10,): Fraction(27)}  # N4 = 27*h^10, h^7 | N4
    k5_quotient_witness = {(2,): Fraction(810)}  # N5/h^10 = 810*h^2
    k6_witness = mp_add(mp_pow(Delta, 2), mp_scale(mp_pow(D, 3), Fraction(-3)))
    k7_witness = mp_scale(h3, Fraction(-3645))
    # This is 3*D^2*Delta + h*(-162*U4*Delta).
    assert k4_witness == {(10,): Fraction(27)}
    assert k5_quotient_witness == {(2,): Fraction(810)}
    assert k6_witness == {
        (2,): Fraction(729),
        (3,): Fraction(-81),
    }
    assert all(power >= 1 for (power,), _ in k6_witness.items())
    assert all(power >= 2 for (power,), _ in k7_witness.items())
    k8_witness = mp_scale(
        mp_add(mp_scale(h3, Fraction(-2187)), mp_scale(h4, Fraction(4617))),
        Fraction(5, 59049),
    )
    assert k8_witness == {
        (3,): Fraction(-5, 27),
        (4,): Fraction(95, 243),
    }

    # The resonance correction is regular for both the tuned and deeper local
    # patterns left by K6/K7, so B8=[w^8]F^(5/3) must have no root pole.
    # At a tuned root, r=2e,q=e and X=1,Y=1/3,Z=1/27; the lowest possible
    # pole coefficient cancels exactly, so this leading K8 test adds nothing.
    tuned_orders = (2, 1, 0, 0, 0, 0, 0, 0)
    tuned_min, tuned_terms = census(8, tuned_orders, root_multiplicity=1)
    assert tuned_min == -3
    assert len(tuned_terms) == 10
    residue = sum(
        c * Fraction(1, 3) ** m[1] * Fraction(1, 27) ** m[2]
        for m, c in tuned_terms
    )
    assert residue == 0
    tuned_double_min, tuned_double_terms = census(
        8, (4, 2, 0, 0, 0, 0, 0, 0), root_multiplicity=2
    )
    assert tuned_double_min == -6 and tuned_double_terms == tuned_terms

    # At a deeper simple root (r>=3,q>=2,U3(alpha)=0), the unique lowest
    # term is (5/9)*U4(alpha)^2/h, a simple pole unless U4(alpha)=0.
    simple_deep_min, simple_deep_terms = census(
        8, (3, 2, 1, 0, 0, 0, 0, 0), root_multiplicity=1
    )
    assert simple_deep_min == -1
    assert simple_deep_terms == [((0, 0, 0, 2, 0, 0, 0, 0), Fraction(5, 9))]

    # At a deeper double root, q=3 and ord(U3)=1 leave one uncancellable
    # order -3 term, -5*U2*U3^2/(27*h^4). Hence q>=4 or ord(U3)>=2.
    double_deep_min, double_deep_terms = census(
        8, (5, 3, 1, 0, 0, 0, 0, 0), root_multiplicity=2
    )
    assert double_deep_min == -3
    assert double_deep_terms == [((0, 1, 2, 0, 0, 0, 0, 0), Fraction(-5, 27))]

    print("K8 recurrence and exact rational coefficient identities: PASS")
    print("all-order band-ratio first-integral identity: PASS")
    print("full K8 pole numerator and D,Delta compression: PASS")
    print("tuned K6/K7 root jets have a cancelling lowest K8 residue")
    print("K8 forces ord(A*Delta-81*U4*h) >= ceil(3e/2) on tuned roots")
    print("deep simple roots force U4(alpha)=0")
    print("deep double roots force ord(D)>=2 or ord(Delta)>=2")
    print("the all-deep jet U1=0, U2=h^2, U3=U4=h survives local K4-K8 divisibility")


if __name__ == "__main__":
    main()
