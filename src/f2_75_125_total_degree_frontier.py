"""Exact top-total-degree recurrence for the (75,125) outer support model.

The script expands the homogeneous coefficient equations of [P,Q]=1 down
through total-degree deficit 15 from the forced top monomials.  Only the
positive-degree rows are used, so nonzero rescaling of P and Q lets us set
their top coefficients to 1.  It checks the complete support slices, the
linear coefficient of each new layer, every resonance compatibility, and
every output index outside the current support layer.

This is a flexibility audit, not an exclusion: the rows through deficit 15
have a free-parameter solution and do not yet constrain the degree pair.
"""

from __future__ import annotations

import sympy as sp


def layer_cap(deficit: int) -> int:
    """Largest k in the degree-(75-d) or degree-(125-d) support slice."""
    return (5 * deficit) // 6


def main() -> None:
    # P_d,k multiplies x^(15-d+k)y^(60-k); Q_d,k multiplies
    # x^(25-d+k)y^(100-k).  These are exactly the allowed support monomials
    # in the indicated total-degree slices for 0 <= d <= 15.
    p_layers: dict[int, dict[int, sp.Expr]] = {0: {0: sp.Integer(1)}}
    q_layers: dict[int, dict[int, sp.Expr]] = {0: {0: sp.Integer(1)}}
    free_symbols: list[sp.Symbol] = []
    resonances: list[tuple[int, int]] = []
    outside_support_rows = 0
    checked_principal_coefficients = 0

    # In the full F2 rectangles the maximum total-degree layer is a single
    # monomial in each component.  Their common degree-25 power base is thus
    # H_tot=x^5*y^20 in the standardized coordinates, not a generic binary
    # form such as x^25+y^25.
    p_top_support = [
        (i, j)
        for i in range(16)
        for j in range(61)
        if 5 * i - j <= 15 and i + j == 75
    ]
    q_top_support = [
        (i, j)
        for i in range(26)
        for j in range(101)
        if 5 * i - j <= 25 and i + j == 125
    ]
    assert p_top_support == [(15, 60)]
    assert q_top_support == [(25, 100)]
    x, y = sp.symbols("x y")
    H_tot = x**5 * y**20
    assert sp.expand(H_tot**3 - x**15 * y**60) == 0
    assert sp.expand(H_tot**5 - x**25 * y**100) == 0

    for deficit in range(16):
        expected_p = [
            (15 - deficit + k, 60 - k)
            for k in range(layer_cap(deficit) + 1)
        ]
        actual_p = [
            (i, j)
            for i in range(16)
            for j in range(61)
            if i + j == 75 - deficit and 5 * i - j <= 15
        ]
        expected_q = [
            (25 - deficit + k, 100 - k)
            for k in range(layer_cap(deficit) + 1)
        ]
        actual_q = [
            (i, j)
            for i in range(26)
            for j in range(101)
            if i + j == 125 - deficit and 5 * i - j <= 25
        ]
        assert actual_p == expected_p
        assert actual_q == expected_q

    def fixed_residual(deficit: int, output_index: int) -> sp.Expr:
        """Coefficient from already chosen layers, excluding P_0,Q_0."""
        result = sp.Integer(0)
        for p_deficit in range(1, deficit):
            q_deficit = deficit - p_deficit
            for p_index, p_coeff in p_layers[p_deficit].items():
                for q_index, q_coeff in q_layers[q_deficit].items():
                    if p_index + q_index != output_index:
                        continue
                    p_x = 15 - p_deficit + p_index
                    p_y = 60 - p_index
                    q_x = 25 - q_deficit + q_index
                    q_y = 100 - q_index
                    result += (p_x * q_y - p_y * q_x) * p_coeff * q_coeff
        return sp.expand(result)

    for deficit in range(1, 16):
        p_layers[deficit] = {}
        q_layers[deficit] = {}
        cap = layer_cap(deficit)

        # No already-fixed pair may contribute beyond this output cap.
        for output_index in range(cap + 1, deficit + 1):
            assert fixed_residual(deficit, output_index) == 0
            outside_support_rows += 1

        for index in range(cap + 1):
            residual = fixed_residual(deficit, index)
            principal_factor = 4 * deficit - 5 * index

            # Directly replay the principal coefficient from [P_0,Q_d] and
            # [P_d,Q_0]: 5(4d-5k)(3 Q_d,k - 5 P_d,k).
            p0_qd = 15 * (100 - index) - 60 * (25 - deficit + index)
            pd_q0 = (15 - deficit + index) * 100 - (60 - index) * 25
            assert p0_qd == 15 * principal_factor
            assert pd_q0 == -25 * principal_factor
            checked_principal_coefficients += 1

            if principal_factor == 0:
                # The new P and Q coefficients disappear from this row.
                # The old-layer residual must therefore vanish on this
                # branch; keep both new coefficients free afterward.
                assert residual == 0
                p_coeff = sp.Symbol(f"p_{deficit}_{index}")
                q_coeff = sp.Symbol(f"q_{deficit}_{index}")
                p_layers[deficit][index] = p_coeff
                q_layers[deficit][index] = q_coeff
                free_symbols.extend((p_coeff, q_coeff))
                resonances.append((deficit, index))
            else:
                p_coeff = sp.Symbol(f"p_{deficit}_{index}")
                p_layers[deficit][index] = p_coeff
                free_symbols.append(p_coeff)
                q_layers[deficit][index] = sp.expand(
                    (5 * p_coeff - residual / (5 * principal_factor)) / 3
                )

    expected_free_dimension = sum(layer_cap(d) + 1 for d in range(1, 16)) + len(
        resonances
    )
    assert expected_free_dimension == 111
    assert len(free_symbols) == expected_free_dimension
    assert resonances == [(5, 4), (10, 8), (15, 12)]

    # Replay every constructed coefficient equation exactly.
    replayed_rows = 0
    for deficit in range(1, 16):
        cap = layer_cap(deficit)
        for index in range(cap + 1):
            residual = fixed_residual(deficit, index)
            p_coeff = p_layers[deficit][index]
            q_coeff = q_layers[deficit][index]
            row = residual + 5 * (4 * deficit - 5 * index) * (
                3 * q_coeff - 5 * p_coeff
            )
            assert sp.expand(row) == 0
            replayed_rows += 1

    print("total_degree_deficits=1..15")
    print("homogeneous_output_degrees=197..183")
    print("top_homogeneous_supports=P:{(15,60)},Q:{(25,100)}")
    print("standardized_common_top_form=x^5*y^20=(x*y^4)^5")
    print("support_layer_formula=0<=k<=floor(5d/6)")
    print("resonance_indices=5:4,10:8,15:12")
    print("resonance_compatibilities=all identically zero")
    print(f"replayed_supported_rows={replayed_rows}")
    print(f"checked_outside_support_rows={outside_support_rows}")
    print(f"free_parameter_count={expected_free_dimension}")
    print("PASS: exact positive-degree top rows are flexible through deficit 15")


if __name__ == "__main__":
    main()
