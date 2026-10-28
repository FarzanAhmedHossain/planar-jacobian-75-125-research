"""Exact rank audit for the reverse-y edge recurrence.

The top y-coefficient of a hypothetical (75,125) pair forces
p_60=a*R^3, q_100=b*R^5 with deg(R)=5. This script computes exact
rational ranks for a squarefree carrier and a pure fifth power, checking the
resonance formula derived in the accompanying note.
"""

from __future__ import annotations

import json

import sympy as sp


def profile(carrier: sp.Expr, x: sp.Symbol, exceptional_depths: set[int]) -> list[dict[str, int]]:
    p_top = sp.expand(carrier**3)
    q_top = sp.expand(carrier**5)
    rows: list[dict[str, int]] = []

    for depth in range(1, 61):
        p_degree = min(15, (60 - depth + 15) // 5)
        q_degree = min(25, (100 - depth + 25) // 5)
        output_degree = max(p_degree + 24, q_degree + 14)

        columns: list[sp.Expr] = []
        for power in range(p_degree + 1):
            u = x**power
            columns.append(
                sp.expand(
                    100 * sp.diff(u, x) * q_top
                    - (60 - depth) * u * sp.diff(q_top, x)
                )
            )
        for power in range(q_degree + 1):
            v = x**power
            columns.append(
                sp.expand(
                    (100 - depth) * sp.diff(p_top, x) * v
                    - 60 * p_top * sp.diff(v, x)
                )
            )

        matrix = sp.Matrix(
            [
                [sp.Poly(column, x).nth(power) for column in columns]
                for power in range(output_degree + 1)
            ]
        )
        rank = int(matrix.rank())
        domain_dimension = (p_degree + 1) + (q_degree + 1)
        target_dimension = output_degree + 1
        resonance = int(depth in exceptional_depths)
        expected_rank = q_degree + 1 - resonance
        assert rank == expected_rank, (depth, rank, expected_rank)
        rows.append(
            {
                "depth": depth,
                "deg_p_layer_cap": p_degree,
                "deg_q_layer_cap": q_degree,
                "domain_dimension": domain_dimension,
                "target_dimension": target_dimension,
                "rank": rank,
                "kernel_dimension": domain_dimension - rank,
                "cokernel_dimension": target_dimension - rank,
                "differential_kernel_dimension": resonance,
            }
        )
    return rows


def main() -> None:
    x = sp.symbols("x")
    squarefree = profile(x**5 + 1, x, {20, 40, 60})
    pure_fifth_power = profile(x**5, x, set(range(4, 61, 4)))
    print(
        json.dumps(
            {
                "squarefree_carrier": "R=x^5+1, a=b=1",
                "squarefree_resonance_depths": [20, 40, 60],
                "squarefree_rank_profile": squarefree,
                "pure_fifth_power_carrier": "R=x^5, a=b=1",
                "pure_fifth_power_resonance_depths": list(range(4, 61, 4)),
                "pure_fifth_power_rank_profile": pure_fifth_power,
                "depth_range": [1, 60],
                "scope": "exact rational confirmation of the derived operator-rank formula",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
