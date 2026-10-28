"""Continue the exact total-degree coefficient recursion through deficit 30.

The recursion is conditional on the theorem-facing F2 support envelope.  It
uses only positive-degree rows of [P,Q]=constant.  Nonresonant layers are
solved exactly for the new Q coefficients, with new P coefficients free; at
the first new resonance it reports the exact compatibility polynomial rather
than dividing by the vanishing principal coefficient.
"""

from __future__ import annotations

import sympy as sp


def layer_indices(deficit: int, *, p_side: bool) -> list[int]:
    """Indices k in a total-degree layer, using the F2 support envelope."""
    xmax, ymax, edge = (15, 60, 15) if p_side else (25, 100, 25)
    degree = 75 if p_side else 125
    return [
        i - (xmax - deficit)
        for i in range(xmax + 1)
        for j in range(ymax + 1)
        if i + j == degree - deficit and 5 * i - j <= edge
    ]


def main() -> None:
    p_layers: dict[int, dict[int, sp.Expr]] = {0: {0: sp.Integer(1)}}
    q_layers: dict[int, dict[int, sp.Expr]] = {0: {0: sp.Integer(1)}}
    free_symbols: list[sp.Symbol] = []
    resonances: list[tuple[int, int]] = []
    outside_rows = 0
    solved_rows = 0

    def fixed_residual(deficit: int, output_index: int) -> sp.Expr:
        value = sp.Integer(0)
        for p_deficit in range(1, deficit):
            q_deficit = deficit - p_deficit
            for p_index, p_coeff in p_layers.get(p_deficit, {}).items():
                for q_index, q_coeff in q_layers.get(q_deficit, {}).items():
                    if p_index + q_index != output_index:
                        continue
                    px = 15 - p_deficit + p_index
                    py = 60 - p_index
                    qx = 25 - q_deficit + q_index
                    qy = 100 - q_index
                    value += (px * qy - py * qx) * p_coeff * q_coeff
        return sp.expand(value)

    max_deficit = 30
    first_obstruction = None
    for deficit in range(1, max_deficit + 1):
        p_indices = layer_indices(deficit, p_side=True)
        q_indices = layer_indices(deficit, p_side=False)
        p_layers[deficit] = {}
        q_layers[deficit] = {}
        output_cap = max(p_indices + q_indices)

        # Coefficients beyond the support of either new layer have no new
        # variables. Their old-layer residual must already vanish.
        for output_index in range(output_cap + 1, deficit + 2):
            assert fixed_residual(deficit, output_index) == 0
            outside_rows += 1

        for output_index in range(output_cap + 1):
            residual = fixed_residual(deficit, output_index)
            principal = 4 * deficit - 5 * output_index
            has_p = output_index in p_indices
            has_q = output_index in q_indices

            if not has_p and not has_q:
                if residual != 0:
                    first_obstruction = (deficit, output_index, residual)
                    break
                solved_rows += 1
                continue

            if principal == 0:
                resonances.append((deficit, output_index))
                if residual != 0:
                    first_obstruction = (deficit, output_index, residual)
                    break
                p_coeff = sp.Symbol(f"p_{deficit}_{output_index}")
                q_coeff = sp.Symbol(f"q_{deficit}_{output_index}")
                p_layers[deficit][output_index] = p_coeff
                q_layers[deficit][output_index] = q_coeff
                free_symbols.extend((p_coeff, q_coeff))
            elif has_q:
                if has_p:
                    p_coeff = sp.Symbol(f"p_{deficit}_{output_index}")
                    p_layers[deficit][output_index] = p_coeff
                    free_symbols.append(p_coeff)
                    q_coeff = sp.expand(
                        (25 * principal * p_coeff - residual)
                        / (15 * principal)
                    )
                else:
                    q_coeff = sp.expand(-residual / (15 * principal))
                q_layers[deficit][output_index] = q_coeff
            elif has_p:
                p_layers[deficit][output_index] = sp.expand(
                    residual / (25 * principal)
                )
            else:
                assert residual == 0
            solved_rows += 1

        if first_obstruction is not None:
            break

    if first_obstruction is None:
        print(f"fully_solved_deficits=1..{max_deficit}")
    else:
        print(f"fully_solved_deficits=1..{deficit - 1}")
    print(f"supported_rows_solved={solved_rows}")
    print(f"outside_support_rows_checked={outside_rows}")
    print(f"resonance_indices={resonances}")
    print(f"free_parameters_before_first_obstruction={len(free_symbols)}")

    if first_obstruction is None:
        print(f"resonance_and_missing_slice_compatibilities=zero through {max_deficit}")
        print("scope=exact top-degree flexibility through deficit 30")
        return

    deficit, output_index, residual = first_obstruction
    variables = sorted(residual.free_symbols, key=str)
    polynomial = sp.Poly(residual, *variables, domain=sp.QQ)
    terms = polynomial.terms()
    print(f"first_nonzero_compatibility=deficit_{deficit}_output_{output_index}")
    print(f"compatibility_degree={polynomial.total_degree()}")
    print(f"compatibility_variables={len(variables)}")
    print(f"compatibility_monomials={len(terms)}")
    print(f"compatibility_factor={sp.factor(residual)}")
    print("scope=one necessary compatibility equation; not an exclusion")


if __name__ == "__main__":
    main()
