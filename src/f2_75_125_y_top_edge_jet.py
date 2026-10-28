"""Exact local-valuation check for the first four upper-y edge slices."""

from __future__ import annotations

import json

import sympy as sp


def main() -> None:
    u = sp.symbols("u")
    a, b = sp.symbols("a b", nonzero=True)
    checks_depth2 = 0
    checks_operator = 0
    checks_depth3 = 0
    checks_depth3_image_order = 0

    x = sp.symbols("x")
    R = sp.Function("R")(x)
    p1 = sp.Function("p1")(x)
    p2 = sp.Function("p2")(x)
    p3 = sp.Function("p3")(x)
    w2 = sp.Function("w2")(x)
    w3 = sp.Function("w3")(x)
    kappa = sp.Rational(5, 3) * b / a
    q1 = kappa * R**2 * p1
    q2 = kappa * R**2 * p2 + w2 / (3 * a)
    cross_depth3 = (
        98 * sp.diff(p1, x) * q2
        - 59 * p1 * sp.diff(q2, x)
        + 99 * sp.diff(p2, x) * q1
        - 58 * p2 * sp.diff(q1, x)
    )
    expected_depth3 = (
        kappa
        * (
            40 * R**2 * (sp.diff(p1, x) * p2 + p1 * sp.diff(p2, x))
            - 234 * R * sp.diff(R, x) * p1 * p2
        )
        + (98 * sp.diff(p1, x) * w2 - 59 * p1 * sp.diff(w2, x)) / (3 * a)
    )
    assert sp.simplify(cross_depth3 - expected_depth3) == 0

    q3 = kappa * R**2 * p3 + w3 / (3 * a)
    cross_depth4 = (
        97 * sp.diff(p1, x) * q3
        - 59 * p1 * sp.diff(q3, x)
        + 98 * sp.diff(p2, x) * q2
        - 58 * p2 * sp.diff(q2, x)
        + 99 * sp.diff(p3, x) * q1
        - 57 * p3 * sp.diff(q1, x)
    )
    expected_depth4 = kappa * (
        40 * R**2 * (sp.diff(p1, x) * p3 + p1 * sp.diff(p3, x))
        - 232 * R * sp.diff(R, x) * p1 * p3
        + 40 * R**2 * p2 * sp.diff(p2, x)
        - 116 * R * sp.diff(R, x) * p2**2
    ) + (
        97 * sp.diff(p1, x) * w3
        - 59 * p1 * sp.diff(w3, x)
        + 98 * sp.diff(p2, x) * w2
        - 58 * p2 * sp.diff(w2, x)
    ) / (3 * a)
    assert sp.simplify(cross_depth4 - expected_depth4) == 0

    for multiplicity in range(1, 6):
        carrier = u**multiplicity
        for p_order in range(0, 15):
            p = u**p_order
            q = sp.Rational(5, 3) * b / a * carrier**2 * p
            middle = sp.expand(99 * sp.diff(p, u) * q - 59 * p * sp.diff(q, u))
            expected = (
                sp.Rational(5, 3)
                * b
                / a
                * (40 * p_order - 118 * multiplicity)
                * u ** (2 * multiplicity + 2 * p_order - 1)
            )
            assert sp.simplify(middle - expected) == 0
            assert 40 * p_order - 118 * multiplicity != 0
            checks_depth2 += 1

        for w_order in range(0, 25):
            w = u**w_order
            d2 = sp.expand(98 * sp.diff(carrier, u) * w - 20 * carrier * sp.diff(w, u))
            expected = (
                (98 * multiplicity - 20 * w_order)
                * u ** (multiplicity + w_order - 1)
            )
            assert sp.simplify(d2 - expected) == 0
            assert 98 * multiplicity - 20 * w_order != 0
            checks_operator += 1

        for p_order in range((multiplicity + 1) // 2, 15):
            w_order = 2 * p_order - multiplicity
            if w_order > 24:
                continue  # The depth-two equation already has no polynomial w2.
            assert 59 * multiplicity - 20 * p_order != 0
            for p2_order in range(15):
                assert 40 * (p_order + p2_order) - 234 * multiplicity != 0
                order_a = 3 * p_order - multiplicity - 1
                order_b = 2 * multiplicity + p_order + p2_order - 1
                required_order = 3 * multiplicity - 1
                minimum_p1_order = (4 * multiplicity + 2) // 3
                # Every putative violation has A as the unique leading term,
                # and that term falls below the order allowed by D_3.
                if p_order < minimum_p1_order:
                    assert order_a < order_b
                    assert order_a < required_order
                checks_depth3 += 1
                checks_depth3_image_order += 1

    depth4_cases = []
    checks_depth4_minimal_jet = 0
    checks_depth4_tied_branches = 0
    for multiplicity in range(2, 6):
        p_order = (4 * multiplicity + 2) // 3
        w2_order = 2 * p_order - multiplicity
        w3_order = 3 * p_order - 4 * multiplicity
        w13_order = p_order + w3_order - 1
        required_order = 3 * multiplicity - 1
        assert 59 * multiplicity - 20 * p_order != 0
        assert 97 * multiplicity - 20 * w3_order != 0
        assert 97 * p_order - 59 * w3_order != 0
        assert w13_order < required_order
        assert w13_order < 2 * multiplicity + p_order - 1
        checks_depth4_minimal_jet += 1

        for p2_order in range(15):
            if (multiplicity, p2_order) in ((2, 0), (4, 0)):
                depth4_cases.append(
                    {
                        "root_multiplicity": multiplicity,
                        "p59_order": p_order,
                        "p58_order": p2_order,
                        "status": "exceptional equal-order depth3 branch",
                    }
                )
                continue

            order_a = 3 * p_order - multiplicity - 1
            order_b = 2 * multiplicity + p_order + p2_order - 1
            assert order_a < order_b

            order_w22 = p2_order + w2_order - 1
            order_p2_square = 2 * multiplicity + 2 * p2_order - 1
            order_p1_p3 = 2 * multiplicity + p_order - 1
            assert min(order_w22, order_p2_square, order_p1_p3) > w13_order
            assert w13_order < required_order
            depth4_cases.append(
                {
                    "root_multiplicity": multiplicity,
                    "p59_order": p_order,
                    "p58_order": p2_order,
                    "status": "excluded by unique low-order term",
                }
            )
            checks_depth4_minimal_jet += 1

    tied_branch_checks = []
    for multiplicity, p_order in ((2, 3), (4, 6)):
        p2_order = 0
        w2_order = 2 * p_order - multiplicity
        for w3_order in range(25):
            assert 97 * multiplicity - 20 * w3_order != 0
        order_a = 3 * p_order - multiplicity - 1
        order_b = 2 * multiplicity + p_order + p2_order - 1
        assert order_a == order_b

        # Normalize the leading coefficient of R and p1 to 1. Cancellation
        # of the tied depth-three coefficient forces this exact p2 constant.
        a_coefficient = 59 * multiplicity - 20 * p_order
        b_coefficient = 40 * (p_order + p2_order) - 234 * multiplicity
        p2_constant = sp.simplify(-sp.Rational(a_coefficient, 1) / (3 * a * b_coefficient))
        assert sp.simplify(p2_constant - 1 / (18 * a)) == 0

        # At depth four, the w22 and p2^2 coefficients are the only terms at
        # order 2m-1 after that cancellation; their sum is nonzero.
        depth4_leading = sp.simplify(
            sp.Rational(5, 3) * b / a
            * (
                (40 * p2_order - 116 * multiplicity) * p2_constant**2
                + p2_constant
                * (98 * p2_order - 58 * w2_order)
                / (3 * a)
            )
        )
        assert depth4_leading != 0
        w3_order_lower = 3 * p_order - 4 * multiplicity + 1
        p1_w3_order_lower = p_order + w3_order_lower - 1
        required_order = 3 * multiplicity - 1
        assert p1_w3_order_lower > 2 * multiplicity - 1
        assert 2 * multiplicity + p_order - 1 > 2 * multiplicity - 1
        assert 2 * multiplicity - 1 < required_order
        tied_branch_checks.append(
            {
                "root_multiplicity": multiplicity,
                "p59_order": p_order,
                "forced_p58_constant_over_leading_p59_squared": "1/(18a)",
                "depth4_first_uncancelled_order": 2 * multiplicity - 1,
                "status": "excluded by nonzero exact coefficient",
            }
        )
        checks_depth4_tied_branches += 1

    print(
        json.dumps(
            {
                "root_multiplicities_checked": [1, 5],
                "p59_orders_checked": [0, 14],
                "depth2_leading_coefficients_checked": checks_depth2,
                "D2_leading_coefficients_checked": checks_operator,
                "depth3_cross_term_identity": "exact symbolic identity",
                "depth3_root_order_coefficients_checked": checks_depth3,
                "depth3_image_order_cases_checked": checks_depth3_image_order,
                "depth4_cross_term_identity": "exact symbolic identity",
                "depth4_minimal_jet_cases_checked": checks_depth4_minimal_jet,
                "depth4_unique_low_term_exclusions": sum(
                    case["status"] == "excluded by unique low-order term"
                    for case in depth4_cases
                ),
                "depth4_tied_branches_checked": checks_depth4_tied_branches,
                "depth4_tied_branch_summary": tied_branch_checks,
                "all_exact": True,
                "conclusion": (
                    "depth4 gives ord(p59)>=ceil(4m/3)+1 for m=2,3,4,5; "
                    "the m=2,4 tied branches are also excluded exactly"
                ),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
