"""Audit the y-adic block structure of the exhaustive F2 support system.

Besides replaying the monomial Jacobian rows in (x,y), this records where a
new y-layer of P and Q enters the coefficient recurrence.  It is a structural
audit only; it does not solve the resulting polynomial system.
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict


def support(i_max: int, j_max: int, edge_max: int) -> list[tuple[int, int]]:
    return [
        (i, j)
        for i in range(i_max + 1)
        for j in range(j_max + 1)
        if 5 * i - j <= edge_max
    ]


def p_degree_bound(y_degree: int) -> int:
    if not 0 <= y_degree <= 60:
        return -1
    return min(15, (y_degree + 15) // 5)


def q_degree_bound(y_degree: int) -> int:
    if not 0 <= y_degree <= 100:
        return -1
    return min(25, (y_degree + 25) // 5)


def main() -> None:
    p_support = support(15, 60, 15)
    q_support = support(25, 100, 25)
    assert len(p_support) == 586
    assert len(q_support) == 1576

    # Store the number of nonzero monomial-pair contributions in each
    # coefficient [x^a y^n][P,Q].
    row_terms: dict[tuple[int, int], int] = defaultdict(int)
    row_witnesses: dict[tuple[int, int], list[tuple[int, int, int, int, int]]] = defaultdict(list)
    by_y_terms: Counter[int] = Counter()
    total_terms = 0
    for i, j in p_support:
        for u, v in q_support:
            coefficient = i * v - j * u
            if coefficient == 0:
                continue
            x_degree, y_degree = i + u - 1, j + v - 1
            assert x_degree >= 0 and y_degree >= 0
            row_terms[(y_degree, x_degree)] += 1
            row_witnesses[(y_degree, x_degree)].append((i, j, u, v, coefficient))
            by_y_terms[y_degree] += 1
            total_terms += 1

    assert total_terms == 912_713
    assert len(row_terms) == 3_790
    assert min(n for n, _ in row_terms) == 0
    assert max(n for n, _ in row_terms) == 159
    assert (159, 39) not in row_terms  # forced corner coefficient cancels

    # Exact y-adic layer sizes and the point at which new unknown layers stop.
    p_bounds = [p_degree_bound(j) for j in range(61)]
    q_bounds = [q_degree_bound(j) for j in range(101)]
    assert p_bounds[0] == 3 and p_bounds[-1] == 15
    assert q_bounds[0] == 5 and q_bounds[-1] == 25
    for j in range(1, 61):
        assert q_degree_bound(j) == p_degree_bound(j) + 2

    early_kernel_dimension = sum(p_degree_bound(j) - 1 for j in range(1, 61))
    early_edge_layers = [j for j in range(1, 61) if j % 5 == 0]
    # At each edge layer, the two leading coefficients of the new P and Q
    # polynomials are prescribed. Compatibility with the top edge equation
    # leaves one kernel direction to match them, reducing the affine freedom
    # by one in each of the 12 edge layers.
    early_after_edge_pinning = early_kernel_dimension - len(early_edge_layers)
    assert early_kernel_dimension == 462
    assert len(early_edge_layers) == 12
    assert early_after_edge_pinning == 450

    blocks = []
    for name, first, last, p_new, q_new in (
        ("both_new_layers", 0, 59, True, True),
        ("Q_only_new_layer", 60, 99, False, True),
        ("closed_tail", 100, 159, False, False),
    ):
        rows = [row for row in row_terms if first <= row[0] <= last]
        blocks.append(
            {
                "name": name,
                "y_output_range": [first, last],
                "slices": last - first + 1,
                "has_new_P_layer": p_new,
                "has_new_Q_layer": q_new,
                "nonzero_scalar_rows": len(rows),
                "nonzero_bilinear_terms": sum(row_terms[row] for row in rows),
            }
        )

    key_slices = []
    for n in (59, 60, 99, 100, 101, 159):
        xs = [x_degree for y_degree, x_degree in row_terms if y_degree == n]
        if not xs:
            continue
        x_degree = max(xs)
        key_slices.append(
            {
                "y_output": n,
                "max_x_output": x_degree,
                "max_x_row_terms": row_terms[(n, x_degree)],
                "max_x_row_monomial_pairs": row_witnesses[(n, x_degree)],
            }
        )

    for n in range(60, 100):
        assert p_degree_bound(n + 1) == -1
        assert q_degree_bound(n + 1) >= 0
    for n in range(100, 160):
        assert p_degree_bound(n + 1) == -1
        assert q_degree_bound(n + 1) == -1

    print(
        json.dumps(
            {
                "P_support_terms": len(p_support),
                "Q_support_terms": len(q_support),
                "nonzero_y_adic_rows": len(row_terms),
                "nonzero_bilinear_terms": total_terms,
                "y_output_range": [0, 159],
                "new_layer_transition": {
                    "both_P_and_Q": [0, 59],
                    "only_Q": [60, 99],
                    "neither_closed_tail": [100, 159],
                },
                "first_60_layers": {
                    "unconstrained_syzygy_dimension": early_kernel_dimension,
                    "edge_layers_with_prescribed_top_coefficients": len(early_edge_layers),
                    "remaining_freedom_after_edge_pinning": early_after_edge_pinning,
                    "edge_layer_indices": early_edge_layers,
                },
                "blocks": blocks,
                "key_slices": key_slices,
                "top_corner_row_159_x39_vanishes": True,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
