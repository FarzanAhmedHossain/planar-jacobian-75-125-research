"""Exact arithmetic audit of the F2 reduction's total-top-form consequence.

This checks only finite integer/support arithmetic, conditional on the cited
standard-pair reduction and the F2 corner theorem. It does not prove those
theorem inputs and it does not exclude a (75,125) Keller pair.
"""

from __future__ import annotations

import json

import sympy as sp


def main() -> None:
    x, y, t = sp.symbols("x y t")

    degree_p, degree_q = 75, 125
    degree_gcd = int(sp.gcd(degree_p, degree_q))
    ratio = (degree_p // degree_gcd, degree_q // degree_gcd)
    assert degree_gcd == 25
    assert ratio == (3, 5)

    a_bar, b_bar = 5, 20  # the F2, j=1 corner A0
    transformed_gcd = a_bar + b_bar
    transports = [
        (M, N)
        for M in range(1, degree_gcd // a_bar + 1)
        for N in range(1, degree_gcd // b_bar + 1)
        if M * a_bar + N * b_bar == degree_gcd
    ]
    assert transformed_gcd == degree_gcd == 25
    assert transports == [(1, 1)]

    # The standardized support lies in these rectangles and contains the
    # forced corner monomials. Only the northeast vertex has maximal total
    # degree in either rectangle.
    p_box = [(i, j) for i in range(16) for j in range(61)]
    q_box = [(i, j) for i in range(26) for j in range(101)]
    p_top = [(i, j) for i, j in p_box if i + j == degree_p]
    q_top = [(i, j) for i, j in q_box if i + j == degree_q]
    assert p_top == [(15, 60)]
    assert q_top == [(25, 100)]
    assert p_top[0] == (3 * 5, 3 * 20)
    assert q_top[0] == (5 * 5, 5 * 20)

    # The common degree-25 base in standard coordinates has exactly two
    # distinct factors, with multiplicities 5 and 20.
    common_base = x**5 * y**20
    assert sp.Poly(common_base, x, y).total_degree() == 25
    base_factors = {
        str(factor): multiplicity
        for factor, multiplicity in sp.factor_list(common_base, x, y)[1]
    }
    assert base_factors == {"x": 5, "y": 20}
    assert common_base == (x * y**4) ** 5

    # x^25+y^25 is squarefree as a binary form: in the x-chart, its roots
    # are those of t^25+1, whose derivative is relatively prime to it.
    fermat_univariate = t**25 + 1
    fermat_derivative = sp.diff(fermat_univariate, t)
    assert sp.gcd(fermat_univariate, fermat_derivative) == 1
    fermat_distinct_linear_factors = 25
    common_base_distinct_linear_factors = 2
    assert fermat_distinct_linear_factors != common_base_distinct_linear_factors

    print(
        json.dumps(
            {
                "raw_degree_pair": [degree_p, degree_q],
                "reduced_ratio": list(ratio),
                "raw_gcd": degree_gcd,
                "F2_corner_A0": [a_bar, b_bar],
                "positive_degree_transport_solutions": [list(pair) for pair in transports],
                "standardized_top_support_P": [list(point) for point in p_top],
                "standardized_top_support_Q": [list(point) for point in q_top],
                "standardized_common_degree_25_base": "x^5*y^20",
                "standardized_base_factor_multiplicities": [5, 20],
                "standardized_base_is_also_a_fifth_power": "x^5*y^20 = (x*y^4)^5",
                "fermat_form_squarefree": True,
                "fermat_distinct_linear_factor_count": fermat_distinct_linear_factors,
                "factor_count_invariant_under_invertible_linear_change": True,
                "conditional_consequence": (
                    "the F2-normalized top-form branch is two-factor type 5+20; "
                    "H=x^25+y^25 cannot be that carrier"
                ),
                "scope": (
                    "exact arithmetic conditional on the standardization/F2 inputs; "
                    "not a proof of the inputs and not an exclusion of (75,125)"
                ),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
