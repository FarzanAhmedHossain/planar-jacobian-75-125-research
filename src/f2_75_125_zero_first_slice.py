"""Exact first constraints on the p59=0 upper-y branch."""

from __future__ import annotations

import json

import sympy as sp


def main() -> None:
    x = sp.symbols("x")
    a, b = sp.symbols("a b", nonzero=True)
    R = sp.Function("R")(x)
    p2 = sp.Function("p2")(x)
    p3 = sp.Function("p3")(x)
    p4 = sp.Function("p4")(x)
    p5 = sp.Function("p5")(x)
    pk = sp.Function("pk")(x)
    w4 = sp.Function("w4")(x)
    w5 = sp.Function("w5")(x)
    kappa = sp.Rational(5, 3) * b / a
    q2 = kappa * R**2 * p2
    q3 = kappa * R**2 * p3

    c4 = 98 * sp.diff(p2, x) * q2 - 58 * p2 * sp.diff(q2, x)
    expected_c4 = kappa * (
        40 * R**2 * p2 * sp.diff(p2, x)
        - 116 * R * sp.diff(R, x) * p2**2
    )
    assert sp.simplify(c4 - expected_c4) == 0

    c5 = (
        97 * sp.diff(p2, x) * q3
        - 58 * p2 * sp.diff(q3, x)
        + 98 * sp.diff(p3, x) * q2
        - 57 * p3 * sp.diff(q2, x)
    )
    expected_c5 = kappa * (
        40 * R**2 * (sp.diff(p2, x) * p3 + p2 * sp.diff(p3, x))
        - 230 * R * sp.diff(R, x) * p2 * p3
    )
    assert sp.simplify(c5 - expected_c5) == 0

    q4 = kappa * R**2 * p4 + w4 / (3 * a)
    c6 = (
        96 * sp.diff(p2, x) * q4
        - 58 * p2 * sp.diff(q4, x)
        + 97 * sp.diff(p3, x) * q3
        - 57 * p3 * sp.diff(q3, x)
        + 98 * sp.diff(p4, x) * q2
        - 56 * p4 * sp.diff(q2, x)
    )
    expected_c6 = kappa * (
        40 * R**2 * (sp.diff(p2, x) * p4 + p2 * sp.diff(p4, x))
        - 228 * R * sp.diff(R, x) * p2 * p4
        + 40 * R**2 * p3 * sp.diff(p3, x)
        - 114 * R * sp.diff(R, x) * p3**2
    ) + (96 * sp.diff(p2, x) * w4 - 58 * p2 * sp.diff(w4, x)) / (3 * a)
    assert sp.simplify(c6 - expected_c6) == 0

    q5 = kappa * R**2 * p5 + w5 / (3 * a)
    c7 = (
        95 * sp.diff(p2, x) * q5
        - 58 * p2 * sp.diff(q5, x)
        + 96 * sp.diff(p3, x) * q4
        - 57 * p3 * sp.diff(q4, x)
        + 97 * sp.diff(p4, x) * q3
        - 56 * p4 * sp.diff(q3, x)
        + 98 * sp.diff(p5, x) * q2
        - 55 * p5 * sp.diff(q2, x)
    )
    expected_c7 = kappa * (
        40 * R**2 * (sp.diff(p2 * p5, x) + sp.diff(p3 * p4, x))
        - 226 * R * sp.diff(R, x) * (p2 * p5 + p3 * p4)
    ) + (
        95 * sp.diff(p2, x) * w5
        - 58 * p2 * sp.diff(w5, x)
        + 96 * sp.diff(p3, x) * w4
        - 57 * p3 * sp.diff(w4, x)
    ) / (3 * a)
    assert sp.simplify(c7 - expected_c7) == 0

    c8 = (
        95 * sp.diff(p3, x) * q5
        - 57 * p3 * sp.diff(q5, x)
        + 97 * sp.diff(p5, x) * q3
        - 55 * p5 * sp.diff(q3, x)
    )
    expected_c8 = kappa * (
        40 * R**2 * sp.diff(p3 * p5, x) - 224 * R * sp.diff(R, x) * p3 * p5
    ) + (95 * sp.diff(p3, x) * w5 - 57 * p3 * sp.diff(w5, x)) / (3 * a)
    assert sp.simplify(c8 - expected_c8) == 0

    qk = kappa * R**2 * pk
    k_symbol = sp.Symbol("k")
    c2k = (100 - k_symbol) * sp.diff(pk, x) * qk - (
        60 - k_symbol
    ) * pk * sp.diff(qk, x)
    expected_c2k = kappa * (
        20 * R**2 * sp.diff(pk**2, x)
        - (120 - 2 * k_symbol) * R * sp.diff(R, x) * pk**2
    )
    assert sp.simplify(c2k - expected_c2k) == 0

    depth4_cases = 0
    depth5_cases = 0
    depth6_cases = 0
    for s in (1, 2, 3):
        for m in range(1, 6):
            assert ((100 - s) * m) % 20 != 0

    for m in range(1, 6):
        for t in range(15):
            coeff4 = 40 * t - 116 * m
            assert coeff4 != 0  # 10t=29m has no solution for 1<=m<=5.
            if 2 * t < m:
                assert 2 * m + 2 * t - 1 < 3 * m - 1
            depth4_cases += 1

            for z in range(15):
                total_order = t + z
                coeff5 = 40 * total_order - 230 * m
                if total_order < m:
                    assert coeff5 != 0
                    assert 2 * m + total_order - 1 < 3 * m - 1
                depth5_cases += 1

        for z in range(15):
            coeff6 = 40 * z - 114 * m
            assert coeff6 != 0  # 20z=57m has no solution for 1<=m<=5.
            if 2 * z < m:
                assert 2 * m + 2 * z - 1 < 3 * m - 1
            depth6_cases += 1

    multiplicity_partitions = [
        [5],
        [4, 1],
        [3, 2],
        [3, 1, 1],
        [2, 2, 1],
        [2, 1, 1, 1],
        [1, 1, 1, 1, 1],
    ]
    d4_kernel_partitions = [
        parts for parts in multiplicity_partitions if all(m % 5 == 0 for m in parts)
    ]
    assert d4_kernel_partitions == [[5]]
    d5_kernel_partitions = [
        parts for parts in multiplicity_partitions if all(m % 4 == 0 for m in parts)
    ]
    assert d5_kernel_partitions == []

    depth7_mixed_cross_cases = 0
    depth8_mixed_cross_cases = 0
    for m in range(1, 5):  # a mixed degree-five carrier has root multiplicity at most four
        for product_order in range(29):  # each slice degree is at most fourteen
            if product_order < m:
                assert 40 * product_order - 226 * m != 0
                assert 40 * product_order - 224 * m != 0
            depth7_mixed_cross_cases += 1
            depth8_mixed_cross_cases += 1

    mixed_first_slice_cases = 0
    mixed_half_order_exclusions = 0
    mixed_leading_resonances = []
    for k in range(2, 20):
        degree_cap = 3 + (60 - k) // 5
        for m in range(1, 6):
            for t in range(degree_cap + 1):
                leading_coefficient = 40 * t - (120 - 2 * k) * m
                if leading_coefficient == 0:
                    mixed_leading_resonances.append(
                        {"first_nonzero_depth": k, "root_multiplicity": m, "slice_order": t}
                    )
                elif 2 * t < m:
                    assert 2 * m + 2 * t - 1 < 3 * m - 1
                    mixed_half_order_exclusions += 1
                mixed_first_slice_cases += 1

    assert mixed_leading_resonances == [
        {"first_nonzero_depth": 4, "root_multiplicity": 5, "slice_order": 14},
        {"first_nonzero_depth": 5, "root_multiplicity": 4, "slice_order": 11},
        {"first_nonzero_depth": 8, "root_multiplicity": 5, "slice_order": 13},
        {"first_nonzero_depth": 10, "root_multiplicity": 2, "slice_order": 5},
        {"first_nonzero_depth": 10, "root_multiplicity": 4, "slice_order": 10},
        {"first_nonzero_depth": 12, "root_multiplicity": 5, "slice_order": 12},
        {"first_nonzero_depth": 15, "root_multiplicity": 4, "slice_order": 9},
        {"first_nonzero_depth": 16, "root_multiplicity": 5, "slice_order": 11},
    ]

    pure_first_slice_resonances = []
    for k in (2, 3):
        degree_cap = 3 + (60 - k) // 5
        for m in range(1, 6):
            for t in range(degree_cap + 1):
                coefficient = 40 * t - (120 - 2 * k) * m
                if coefficient == 0:
                    pure_first_slice_resonances.append((k, m, t))
                if coefficient != 0 and 2 * t < m:
                    assert 2 * m + 2 * t - 1 < 3 * m - 1
    assert pure_first_slice_resonances == []

    p58_nonzero_cases = 0
    p57_nonzero_forced_t_cases = 0
    p57_zero_forced_t_cases = 0
    for m in range(1, 6):
        for t in range((m + 1) // 2, 15):
            ell4 = 2 * t - m
            assert 40 * t - 116 * m != 0
            assert 96 * m - 20 * ell4 != 0
            w_term_coefficient = 58 * m - 20 * t
            assert w_term_coefficient != 0
            p58_nonzero_cases += 1

            for z in range(15):
                if t + z < m:
                    order_w_term = 3 * t - m - 1
                    order_pair_24_min = 2 * m + t - 1
                    order_p3_square = 2 * m + 2 * z - 1
                    assert order_w_term < order_pair_24_min
                    assert order_w_term < order_p3_square
                    assert order_w_term < 3 * m - 1
                    p57_nonzero_forced_t_cases += 1

            if t < (4 * m + 2) // 3:
                order_w_term = 3 * t - m - 1
                order_pair_24_min = 2 * m + t - 1
                assert order_w_term < order_pair_24_min
                assert order_w_term < 3 * m - 1
                p57_zero_forced_t_cases += 1

    # At depth six, compare the three local orders when p58, p57, and p56
    # are nonzero.  A low minimum can survive only as a tie; the bounded
    # census isolates the single tied stratum t=m,z=0.  The p56 term cannot
    # participate in that tie because its order is at least 3m-1 there.
    depth6_low_order_survivors = []
    depth6_order_cases = 0
    for m in range(1, 6):
        required = 3 * m - 1
        for t in range((m + 1) // 2, 15):
            for z in range(15):
                if t + z < m:  # the depth-five condition
                    continue
                order_w = 3 * t - m - 1
                order_p58_p56 = 2 * m + t - 1  # v=ord(p56)>=0
                order_p57_square = 2 * m + 2 * z - 1
                assert 40 * (t + 0) - 228 * m != 0
                assert 40 * z - 114 * m != 0
                for v in [None, *range(15)]:
                    if v is None:  # p56=0, so its term is absent
                        orders = [order_w, order_p57_square]
                    else:
                        order_p58_p56_v = order_p58_p56 + v
                        assert 40 * (t + v) - 228 * m != 0
                        orders = [order_w, order_p58_p56_v, order_p57_square]
                    low = min(orders)
                    tied = sum(order == low for order in orders)
                    depth6_order_cases += 1
                    if low < required and tied == 1:
                        continue
                    if low < required:
                        depth6_low_order_survivors.append((m, t, z, v, low, tied))

    assert all((m, t, z) == (m, m, 0) for m, t, z, *_ in depth6_low_order_survivors)
    assert len(depth6_low_order_survivors) == 5 * 16

    # In the exceptional stratum, write R=c*u^m, p58=e*u^m+..., and
    # p57=d+..., with c,e,d nonzero.  The depth-four equation gives the
    # leading coefficient F=kappa*e^2/c of w4=F*u^m+..., and cancellation
    # at order 2m-1 in depth six is exactly e^3=9*a*c^3*d^2.
    c, e, d = sp.symbols("c e d", nonzero=True)
    F = kappa * e**2 / c
    low_order_coefficient = e * F / (3 * a) - 3 * kappa * c**2 * d**2
    assert sp.simplify(
        low_order_coefficient - kappa * (e**3 - 9 * a * c**3 * d**2) / (3 * a * c)
    ) == 0

    # The next layer removes that exceptional tie.  Depth four and five
    # determine F and G; the only depth-seven coefficient at order m-1 is
    # the displayed w5/w4 contribution, and it is nonzero in characteristic 0.
    m = sp.symbols("m", positive=True, integer=True)
    G = 2 * kappa * e * d / c
    depth4_c = -76 * m * kappa * c**2 * e**2
    depth4_new = 76 * m * c**3 * F
    depth5_c = -190 * m * kappa * c**2 * e * d
    depth5_new = 95 * m * c**3 * G
    assert sp.simplify(depth4_new + depth4_c) == 0
    assert sp.simplify(depth5_new + depth5_c) == 0
    depth7_low = m * (95 * e * G - 57 * d * F) / (3 * a)
    expected_depth7_low = 133 * m * kappa * e**2 * d / (3 * a * c)
    assert sp.simplify(depth7_low - expected_depth7_low) == 0
    for root_multiplicity in range(1, 6):
        for p4_order in range(15):
            for p5_order in range(15):
                assert 2 * root_multiplicity + p4_order - 1 > root_multiplicity - 1
                assert 3 * root_multiplicity + p5_order - 1 > root_multiplicity - 1

    print(
        json.dumps(
            {
                "c4_identity": "exact symbolic identity",
                "c5_identity": "exact symbolic identity",
                "c6_identity": "exact symbolic identity",
                "c7_identity": "exact symbolic identity",
                "c8_cross_identity": "exact symbolic identity",
                "general_first_nonzero_slice_identity": "exact symbolic identity",
                "D1_D2_D3_root_kernel_cases_checked": 15,
                "depth4_local_cases_checked": depth4_cases,
                "depth5_local_cases_checked": depth5_cases,
                "depth6_local_cases_checked": depth6_cases,
                "D4_homogeneous_kernel_root_pattern": d4_kernel_partitions,
                "D5_homogeneous_kernel_root_patterns": d5_kernel_partitions,
                "mixed_p58_zero_p57_nonzero_cross_order_cases": {
                    "depth7_coefficient_40s_minus_226m": depth7_mixed_cross_cases,
                    "depth8_coefficient_40s_minus_224m": depth8_mixed_cross_cases,
                    "depth7_resonance_for_m_le_4": False,
                    "depth8_resonance_for_m_le_4": False,
                },
                "mixed_first_nonzero_slice_cases_checked": mixed_first_slice_cases,
                "mixed_half_order_exclusions": mixed_half_order_exclusions,
                "mixed_leading_resonances": mixed_leading_resonances,
                "pure_first_nonzero_depth_2_or_3_resonances": pure_first_slice_resonances,
                "p58_nonzero_local_cases_checked": p58_nonzero_cases,
                "p57_nonzero_unique_term_cases": p57_nonzero_forced_t_cases,
                "p57_zero_unique_term_cases": p57_zero_forced_t_cases,
                "depth6_mixed_order_cases_checked": depth6_order_cases,
                "depth6_low_order_tie_stratum": {
                    "root_multiplicities": sorted(
                        {item[0] for item in depth6_low_order_survivors}
                    ),
                    "ord_p58_equals_root_multiplicity": True,
                    "ord_p57": 0,
                    "ord_p56": "zero or any value from 0 through 14",
                    "forced_leading_relation": "e^3 = 9*a*c^3*d^2",
                },
                "depth7_tie_exclusion": {
                    "forced_ord_w5": 0,
                    "forced_w5_leading_coefficient": "G = 2*kappa*e*d/c",
                    "depth7_first_uncancelled_order": "m-1",
                    "depth7_leading_coefficient": "133*m*kappa*e^2*d/(3*a*c)",
                    "new_layer_image_minimum_order": "3*m-1",
                    "status": "exceptional depth-six tie is impossible",
                },
                "all_exact": True,
                "conclusion": (
                    "if p59=0 then p58=0 or ord_root(p58)>=ceil(m/2); "
                    "if p58,p57 are both nonzero then ord(p58)>=ceil(4m/3) and "
                    "ord(p57)>=ceil(m/2) (the only depth-six tie is excluded at depth seven); "
                    "if p58!=0 and p57=0 then ord_root(p58)>=ceil(4m/3); "
                    "if p58=0 and p57!=0 then ord_root(p57)>=ceil(m/2)"
                ),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
