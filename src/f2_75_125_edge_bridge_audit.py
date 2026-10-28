"""Exact arithmetic audit for the (75,125) F2 support-edge bridge.

The source inputs are the standard-pair reduction and the F2, j=1 row
(A0,A0')=((5,20),(1,0)), (m,n)=(3,5).  The checker verifies the resulting
edge geometry and lattice arithmetic.  It also records that y -> y+1 does
not preserve the (5,-1) half-space, so that half-space must come from the
regular corner of the already standardized pair.

This is not a proof of the cited normal-form or complete-chain theorems and
does not eliminate the resulting coefficient system.
"""

from __future__ import annotations

import json
from math import gcd


def scaled(point: tuple[int, int], factor: int) -> tuple[int, int]:
    return factor * point[0], factor * point[1]


def support(i_max: int, j_max: int, edge_max: int) -> set[tuple[int, int]]:
    return {
        (i, j)
        for i in range(i_max + 1)
        for j in range(j_max + 1)
        if 5 * i - j <= edge_max
    }


def band_profile(points: set[tuple[int, int]]) -> tuple[int, int, int]:
    bands = {5 * i - j for i, j in points}
    return min(bands), max(bands), len(bands)


def edge_normal(
    start: tuple[int, int], end: tuple[int, int]
) -> tuple[int, int]:
    di = end[0] - start[0]
    dj = end[1] - start[1]
    divisor = gcd(abs(di), abs(dj))
    # For the upward-right F2 edge, (dj,-di) points outward.
    return dj // divisor, -di // divisor


def main() -> None:
    m, n = 3, 5
    a0 = (5, 20)
    a0_prime = (1, 0)

    p_edge = (scaled(a0_prime, m), scaled(a0, m))
    q_edge = (scaled(a0_prime, n), scaled(a0, n))
    assert p_edge == ((3, 0), (15, 60))
    assert q_edge == ((5, 0), (25, 100))
    assert edge_normal(*p_edge) == edge_normal(*q_edge) == (5, -1)

    p_edge_values = [5 * i - j for i, j in p_edge]
    q_edge_values = [5 * i - j for i, j in q_edge]
    assert p_edge_values == [15, 15]
    assert q_edge_values == [25, 25]

    p = support(15, 60, 15)
    q = support(25, 100, 25)
    assert len(p) == 586
    assert len(q) == 1576
    assert set(p_edge) <= p
    assert set(q_edge) <= q
    assert p <= q
    assert band_profile(p) == (-60, 15, 76)
    assert band_profile(q) == (-100, 25, 126)
    assert max(i + j for i, j in p) == 75
    assert max(i + j for i, j in q) == 125
    assert {(i, j) for i, j in p if i + j == 75} == {(15, 60)}
    assert {(i, j) for i, j in q if i + j == 125} == {(25, 100)}

    # A nonzero y-translation preserves the rectangular box but not this
    # half-space: x^15 y^60 -> x^15 (y+1)^60 contains x^15.
    translated_corner_term = (15, 0)
    assert 0 <= translated_corner_term[0] <= 15
    assert 0 <= translated_corner_term[1] <= 60
    assert 5 * translated_corner_term[0] - translated_corner_term[1] == 75
    assert 75 > 15

    print(
        json.dumps(
            {
                "source_inputs": {
                    "family": "F2",
                    "j": 1,
                    "A0": list(a0),
                    "A0_prime": list(a0_prime),
                    "m": m,
                    "n": n,
                },
                "P_edge_endpoints": [list(point) for point in p_edge],
                "Q_edge_endpoints": [list(point) for point in q_edge],
                "primitive_support_normal": [5, -1],
                "P_edge_support_value": 15,
                "Q_edge_support_value": 25,
                "P_support_points": len(p),
                "Q_support_points": len(q),
                "P_band_profile": list(band_profile(p)),
                "Q_band_profile": list(band_profile(q)),
                "P_unique_total_degree_vertex": [15, 60],
                "Q_unique_total_degree_vertex": [25, 100],
                "translation_preserves_rectangle": True,
                "translation_preserves_edge_halfspace": False,
                "translated_corner_witness": {
                    "monomial_exponent": list(translated_corner_term),
                    "5i_minus_j": 75,
                    "P_edge_bound": 15,
                },
                "scope": "exact arithmetic conditional on the cited F2 corner input",
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
