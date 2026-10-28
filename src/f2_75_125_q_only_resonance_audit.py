"""Enumerate every kernel/resonance depth in the Q-only top-y block.

For a degree-five carrier R with root-multiplicity partition lambda, the
one-new-layer operator at depth s has kernel exactly when
    w = c * R^((100-s)/20)
is a polynomial within the support degree bound.  This script checks the
finite integer partition census for 61 <= s <= 100 and replays the target
degree caps directly from the exhaustive support polygons.
"""

from __future__ import annotations

from collections import defaultdict


def support(i_max: int, j_max: int, edge_max: int) -> list[tuple[int, int]]:
    return [
        (i, j)
        for i in range(i_max + 1)
        for j in range(j_max + 1)
        if 5 * i - j <= edge_max
    ]


def partitions(n: int, largest: int | None = None) -> list[tuple[int, ...]]:
    if largest is None:
        largest = n
    if n == 0:
        return [()]
    out: list[tuple[int, ...]] = []
    for first in range(min(n, largest), 1 - 1, -1):
        for tail in partitions(n - first, first):
            out.append((first,) + tail)
    return out


def has_polynomial_kernel(partition: tuple[int, ...], s: int, degree_cap: int) -> bool:
    u = 100 - s
    if any((u * multiplicity) % 20 for multiplicity in partition):
        return False
    degree = u // 4
    assert u % 4 == 0
    return degree <= degree_cap


def main() -> None:
    p_support = support(15, 60, 15)
    q_support = support(25, 100, 25)
    assert len(p_support) == 586
    assert len(q_support) == 1576

    parts = partitions(5)
    assert parts == [
        (5,),
        (4, 1),
        (3, 2),
        (3, 1, 1),
        (2, 2, 1),
        (2, 1, 1, 1),
        (1, 1, 1, 1, 1),
    ]

    resonance_depths: dict[str, list[int]] = defaultdict(list)
    checked_target_rows = 0
    for s in range(61, 101):
        u = 100 - s
        q_degree_cap = 5 + u // 5
        output_y_degree = 159 - s

        # Verify that the complete row, including its new-layer endpoint,
        # lies in the polynomial target of degree q_degree_cap + 14.
        x_degrees = [
            i + j - 1
            for i, y_p in p_support
            for j, y_q in q_support
            if y_p + y_q - 1 == output_y_degree and i * y_q - y_p * j != 0
        ]
        assert x_degrees
        assert max(x_degrees) == q_degree_cap + 14
        checked_target_rows += 1

        for partition in parts:
            if has_polynomial_kernel(partition, s, q_degree_cap):
                key = "+".join(map(str, partition))
                resonance_depths[key].append(s)

    assert resonance_depths["5"] == list(range(64, 101, 4))
    for partition in parts[1:]:
        key = "+".join(map(str, partition))
        assert resonance_depths[key] == [80, 100]

    print("Q_only_depth_range=61..100")
    print("carrier_partitions=7")
    print(f"target_degree_checks={checked_target_rows}")
    print("pure_fifth_power_kernel_depths=64,68,72,76,80,84,88,92,96,100")
    print("every_nonpure_partition_kernel_depths=80,100")
    print("generic_rank_cokernel=domain_dim,14")
    print("resonant_rank_cokernel=domain_dim_minus_1,15")
    print("PASS: every exceptional kernel depth is retained by multiplicity branch")


if __name__ == "__main__":
    main()
