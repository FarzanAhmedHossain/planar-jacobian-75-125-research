"""Exact arithmetic audit for the standard F2 (75,125) support envelope.

This checks the lattice enumeration conditional on the standard-pair
reduction and the F2 corner of the already standardized pair. The edge
half-space is supplied by that regular corner; it is not transported through
the final y-translation. See f2_75_125_edge_bridge_audit.py for an exact
counterexample to such a transport argument.
"""

from __future__ import annotations

import json


def ceil_div(a: int, b: int) -> int:
    """Exact ceiling of a/b for positive b."""
    return -((-a) // b)


def make_support(i_max: int, j_max: int, edge_max: int) -> set[tuple[int, int]]:
    return {
        (i, j)
        for i in range(i_max + 1)
        for j in range(j_max + 1)
        if 5 * i - j <= edge_max
    }


def audit_bands(
    support: set[tuple[int, int]],
    i_max: int,
    j_max: int,
    band_min: int,
    band_max: int,
) -> list[dict[str, int]]:
    bands: dict[int, set[int]] = {k: set() for k in range(band_min, band_max + 1)}
    for i, j in support:
        bands[5 * i - j].add(i)

    profile = []
    for k in range(band_min, band_max + 1):
        lo = max(0, ceil_div(k, 5))
        hi = min(i_max, (k + j_max) // 5)
        expected = set(range(lo, hi + 1)) if lo <= hi else set()
        assert bands[k] == expected, (k, sorted(bands[k]), sorted(expected))
        assert expected, f"unexpected empty band {k}"
        profile.append({"band": k, "i_min": lo, "i_max": hi, "terms": len(expected)})
    return profile


def near_edge_pole_depths(
    top_band: int, top_t_power: int, last_distance: int, cap: int
) -> list[dict[str, int]]:
    """Audit the allowed t-pole depth after factoring x^top_t_power."""
    result = []
    for distance in range(last_distance + 1):
        band = top_band - distance
        i_min = max(0, ceil_div(band, 5))
        depth = top_t_power - i_min
        expected = min(distance // 5, cap)
        assert depth == expected, (distance, band, depth, expected)
        if distance % 5 == 0:
            result.append(
                {
                    "band_distance": distance,
                    "allowed_t_pole_depth": depth,
                }
            )
    return result


def main() -> None:
    p_support = make_support(i_max=15, j_max=60, edge_max=15)
    q_support = make_support(i_max=25, j_max=100, edge_max=25)

    p_count_formula = 4 * 61 + sum(76 - 5 * i for i in range(4, 16))
    q_count_formula = 6 * 101 + sum(126 - 5 * i for i in range(6, 26))
    assert len(p_support) == p_count_formula == 586
    assert len(q_support) == q_count_formula == 1576
    assert p_support <= q_support, "the polynomial target shear Q -> Q + lambda*P preserves S_Q"

    p_endpoints = {(3, 0), (15, 60)}
    q_endpoints = {(5, 0), (25, 100)}
    assert p_endpoints <= p_support
    assert q_endpoints <= q_support
    assert {5 * i - j for i, j in p_endpoints} == {15}
    assert {5 * i - j for i, j in q_endpoints} == {25}

    p_profile = audit_bands(p_support, 15, 60, -60, 15)
    q_profile = audit_bands(q_support, 25, 100, -100, 25)
    p_resonance_depths = near_edge_pole_depths(15, 3, 75, 3)
    q_resonance_depths = near_edge_pole_depths(25, 5, 125, 5)
    assert p_profile[-1] == {"band": 15, "i_min": 3, "i_max": 15, "terms": 13}
    assert q_profile[-1] == {"band": 25, "i_min": 5, "i_max": 25, "terms": 21}

    print(
        json.dumps(
            {
                "P_support_terms": len(p_support),
                "Q_support_terms": len(q_support),
                "P_support_contained_in_Q_support": p_support <= q_support,
                "P_band_range": [-60, 15],
                "Q_band_range": [-100, 25],
                "P_top_band_i_range": [3, 15],
                "Q_top_band_i_range": [5, 25],
                "P_resonance_t_pole_depths": p_resonance_depths,
                "Q_resonance_t_pole_depths": q_resonance_depths,
                "P_bands": p_profile,
                "Q_bands": q_profile,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
