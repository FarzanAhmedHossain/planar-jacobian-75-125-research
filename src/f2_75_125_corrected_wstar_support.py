"""Support-envelope audit for the corrected edge polynomial W*.

The high edge of W=Q^3-P^5 is cancelled by a finite list of lower power
products P^a Q^b.  This checker takes the union of support envelopes of the
two raw products and every listed correction, then computes safe t-degree
intervals on each (5,-1)-band.  It uses only exact integer support arithmetic.
"""

from __future__ import annotations

from collections import defaultdict


def source_support(i_max: int, j_max: int, edge_max: int) -> set[tuple[int, int]]:
    return {
        (i, j)
        for i in range(i_max + 1)
        for j in range(j_max + 1)
        if 5 * i - j <= edge_max
    }


def profile(support: set[tuple[int, int]]) -> dict[int, tuple[int, int]]:
    terms: dict[int, list[int]] = defaultdict(list)
    for i, j in support:
        terms[5 * i - j].append(i)
    return {w: (min(xs), max(xs)) for w, xs in terms.items()}


def convolve(
    left: dict[int, tuple[int, int]], right: dict[int, tuple[int, int]]
) -> dict[int, tuple[int, int]]:
    result: dict[int, tuple[int, int]] = {}
    for w, (lo, hi) in left.items():
        for v, (left_lo, left_hi) in right.items():
            weight = w + v
            interval = (lo + left_lo, hi + left_hi)
            if weight in result:
                old_lo, old_hi = result[weight]
                result[weight] = (min(old_lo, interval[0]), max(old_hi, interval[1]))
            else:
                result[weight] = interval
    return result


def power_product(
    p_profile: dict[int, tuple[int, int]],
    q_profile: dict[int, tuple[int, int]],
    p_power: int,
    q_power: int,
) -> dict[int, tuple[int, int]]:
    result = {0: (0, 0)}
    for _ in range(p_power):
        result = convolve(result, p_profile)
    for _ in range(q_power):
        result = convolve(result, q_profile)
    return result


def main() -> None:
    p_support = source_support(15, 60, 15)
    q_support = source_support(25, 100, 25)
    p_profile = profile(p_support)
    q_profile = profile(q_support)

    # Raw W terms followed by all possible corrections used to remove faces
    # 70, 65, 60, 55, 50, 45, and 40.  A zero correction only enlarges this
    # union, so it is safe to include it.
    products = {
        "Q^3": (0, 3),
        "P^5": (5, 0),
        "P^3Q": (3, 1),
        "PQ^2": (1, 2),
        "P^4": (4, 0),
        "P^2Q": (2, 1),
        "Q^2": (0, 2),
        "P^3": (3, 0),
        "PQ": (1, 1),
    }
    profiles = {
        name: power_product(p_profile, q_profile, *powers)
        for name, powers in products.items()
    }

    expected = {
        39: (8, 67),
        38: (8, 67),
        37: (8, 67),
        36: (8, 67),
        35: (7, 67),
        34: (7, 66),
        33: (7, 66),
        32: (7, 66),
        31: (7, 66),
        30: (6, 66),
    }
    union = {}
    for band in expected:
        intervals = [
            supports[band]
            for supports in profiles.values()
            if band in supports
        ]
        union[band] = (
            min(interval[0] for interval in intervals),
            max(interval[1] for interval in intervals),
        )
    assert union == expected, (union, expected)

    print("corrected_Wstar_support=PASS")
    print("products=" + ",".join(products))
    for band, (lo, hi) in union.items():
        print(f"band={band} t_degree={lo}..{hi}")
    print("scope=safe support envelope from the F2 outer support model")


if __name__ == "__main__":
    main()
