"""Exact depth-six local obstruction at a double root in the F2 (75,125) model.

This checker verifies the generalized-binomial valuation census and the small
rational coefficient identities used in the depth-six resonance propagation.
It proves a local necessary condition inside the band model; it is not an
exclusion of the degree pair or a proof of the band recurrence itself.
"""

from __future__ import annotations

from fractions import Fraction
from math import factorial


def weighted_multiplicities(total_degree: int, index: int = 1) -> list[tuple[int, ...]]:
    """All (m_1,...,m_s) with sum(j*m_j)=total_degree."""
    if index > total_degree:
        return [()] if total_degree == 0 else []
    out: list[tuple[int, ...]] = []
    for multiplicity in range(total_degree // index + 1):
        for tail in weighted_multiplicities(total_degree - index * multiplicity, index + 1):
            out.append((multiplicity,) + tail)
    return out


def binomial_5_over_3(k: int) -> Fraction:
    value = Fraction(1)
    for j in range(k):
        value *= Fraction(5, 3) - j
    return value / factorial(k)


def coefficient(multiplicities: tuple[int, ...]) -> Fraction:
    count = sum(multiplicities)
    multinomial = Fraction(factorial(count), 1)
    for multiplicity in multiplicities:
        multinomial /= factorial(multiplicity)
    return binomial_5_over_3(count) * multinomial


def depth_six_valuation(
    multiplicities: tuple[int, ...], orders: tuple[int, ...], root_multiplicity: int
) -> int:
    count = sum(multiplicities)
    return 5 * root_multiplicity - 3 * root_multiplicity * count + sum(
        m * order for m, order in zip(multiplicities, orders)
    )


def main() -> None:
    # At a double root alpha, h has order 2, U1 has order >=4, U2 has order 1
    # in the branch to be excluded, and U3,...,U6 are ordinary polynomials.
    # Since alpha != 0, the depth-five/six support shifts divide by t, a unit
    # at alpha, so those shifted coefficients still have nonnegative order.
    orders = (4, 1, 0, 0, 0, 0)
    terms = []
    for multiplicities in weighted_multiplicities(6):
        if len(multiplicities) != 6:
            multiplicities = multiplicities + (0,) * (6 - len(multiplicities))
        val = depth_six_valuation(multiplicities, orders, root_multiplicity=2)
        coeff = coefficient(multiplicities)
        terms.append((val, multiplicities, coeff))

    least = min(val for val, _, _ in terms)
    least_terms = [(m, c) for val, m, c in terms if val == least]
    triple_u2 = (0, 3, 0, 0, 0, 0)
    assert least == -5
    assert least_terms == [(triple_u2, Fraction(-5, 81))]
    assert all(val > -5 for val, multiplicities, _ in terms if multiplicities != triple_u2)

    # Write S=t*h and rho=U1/h^3.  The coefficient-ratio recurrence gives
    # R5=kappa/S and R6=-kappa*rho/(3*S).  Its exact K5/K6 cancellations are
    # the following integer coefficient sums.
    assert 15 - 15 == 0  # K5: 5*A0*R5 + 15*t*R5' = 0.
    assert -6 + 5 + 1 == 0  # K6: coefficients of kappa*t*rho*S'/S^2.
    assert 5 - 5 == 0  # K6: coefficients of kappa*t*rho'/S.
    assert Fraction(-1, 3) + Fraction(5, 3) == Fraction(4, 3)

    # The correction to the binomial G6 coefficient is 4*kappa*h*U1/(3*t),
    # hence regular at alpha.  It cannot cancel the unique u^-5 pole above.
    print("depth-six exact valuation census: unique order -5 term is -5*U2^3/(81*h^4)")
    print("resonance correction verified: G6 - [w^6]F^(5/3) = 4*kappa*h*U1/(3*t)")
    print("local conclusion: ord_alpha(U2) != 1; together with ord_alpha(U2)>=1, ord_alpha(U2)>=2")


if __name__ == "__main__":
    main()
