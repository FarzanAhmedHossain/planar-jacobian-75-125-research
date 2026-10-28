"""Exact support and coefficient audit linking the F2 quartic edge to Magnus R.

The accompanying note proves the degree-gap step.  This checker verifies the
integer support crosswalk and the two explicit quartic coefficient vectors.
It is a necessary-condition bridge only, not an exclusion certificate.
"""

from __future__ import annotations

import json
from fractions import Fraction
from itertools import product


Alg = tuple[Fraction, Fraction]  # a + b*sqrt(5)
ZERO: Alg = (Fraction(0), Fraction(0))
ONE: Alg = (Fraction(1), Fraction(0))


def add(a: Alg, b: Alg) -> Alg:
    return (a[0] + b[0], a[1] + b[1])


def neg(a: Alg) -> Alg:
    return (-a[0], -a[1])


def sub(a: Alg, b: Alg) -> Alg:
    return add(a, neg(b))


def mul(a: Alg, b: Alg) -> Alg:
    return (a[0] * b[0] + 5 * a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def scale(q: int, a: Alg) -> Alg:
    return (q * a[0], q * a[1])


def poly_mul(a: list[Alg], b: list[Alg]) -> list[Alg]:
    out = [ZERO for _ in range(len(a) + len(b) - 1)]
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] = add(out[i + j], mul(ai, bj))
    while len(out) > 1 and out[-1] == ZERO:
        out.pop()
    return out


def fmt(a: Alg) -> str:
    return f"({a[0]})+({a[1]})*sqrt(5)"


def main() -> None:
    # The F2 edge polynomial is P_15(t)=a*t^3*h(t)^3.  With h monic,
    # its two allowed quartic coefficient vectors are as follows.
    alpha_plus = (Fraction(3, 2), Fraction(1, 2))
    alpha_minus = (Fraction(3, 2), Fraction(-1, 2))
    one_minus = [neg(ONE), ONE]

    h22_vectors = []
    for alpha in (alpha_plus, alpha_minus):
        t_minus_alpha = [neg(alpha), ONE]
        h22 = poly_mul(poly_mul(one_minus, one_minus),
                       poly_mul(t_minus_alpha, t_minus_alpha))
        expected = [
            sub(scale(3, alpha), ONE),
            sub((Fraction(2), Fraction(0)), scale(8, alpha)),
            scale(7, alpha),
            scale(-2, add(alpha, ONE)),
            ONE,
        ]
        if h22 != expected:
            raise AssertionError("2+2 coefficient vector mismatch")
        h22_vectors.append([fmt(c) for c in h22])

    t2_minus_3t_plus_3 = [(Fraction(3), Fraction(0)),
                           (Fraction(-3), Fraction(0)), ONE]
    h211 = poly_mul(poly_mul(one_minus, one_minus), t2_minus_3t_plus_3)
    expected_211 = [(Fraction(c), Fraction(0)) for c in (3, -9, 10, -5, 1)]
    if h211 != expected_211:
        raise AssertionError("2+1+1 coefficient vector mismatch")

    # In t=xy^5, z=y^-1 coordinates, ell=5*i-j.  The Magnus remainder
    # rectangle [0,9]x[0,36] intersects ell=15 only at t exponents 3..9.
    z15_exponents = [
        i for i in range(10)
        if 0 <= 5 * i - 15 <= 36
    ]
    if z15_exponents != list(range(3, 10)):
        raise AssertionError("unexpected Z edge support")

    # The Magnus pre-generator has 0<=i<=5, 0<=j<=20, ell<=5.
    # Its ell=5 edge is precisely the five monomials t^i z^5, i=1..5.
    r_edge = [
        (i, j) for i in range(6) for j in range(21)
        if 5 * i - j == 5
    ]
    if r_edge != [(i, 5 * i - 5) for i in range(1, 6)]:
        raise AssertionError("unexpected Magnus generator edge support")

    # Since every Magnus-generator monomial has ell<=5, total ell=15 in R^3
    # can only arise from a product of three ell=5 edge monomials.
    possible_lambdas = {5 * i - j for i in range(6) for j in range(21)
                        if 5 * i - j <= 5 and 0 <= 5 * i - j}
    triples_at_15 = [triple for triple in product(possible_lambdas, repeat=3)
                     if sum(triple) == 15]
    if triples_at_15 != [(5, 5, 5)]:
        raise AssertionError("R^3 edge is not the cube of the R edge")

    print(json.dumps({
        "status": "PASS",
        "F2_edge": "P_15=a*t^3*h^3, deg(h)=4",
        "Magnus_Z_edge_t_exponents": z15_exponents,
        "Magnus_R_edge_support_xy": r_edge,
        "Magnus_R_edge_t_exponents": list(range(1, 6)),
        "degree_gap": {
            "deg(h_P^3-h_R^3)": "<=6 after matching leading coefficient",
            "distinct monic quartics imply degree": ">=8",
            "conclusion": "h_R is proportional to h_P",
        },
        "2+2_coefficients_ascending": h22_vectors,
        "2+1+1_coefficients_ascending": [str(c[0]) for c in h211],
        "scope": "exact necessary-condition bridge; not an exclusion",
    }, indent=2))


if __name__ == "__main__":
    main()
