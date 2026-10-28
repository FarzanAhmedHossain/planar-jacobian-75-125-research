"""Exact modular rank certificate for both F2 2+2 terminal root values.

The weight-49 kernel has 64 coefficient columns: its first four are the
coefficients of S(t)=s_0+s_1*t+s_2*t^2+s_3*t^3.  The distinguished terminal
chart imposes S(1)=0.  The second eligible double-root chart imposes
S(alpha)=0, alpha=(3+sqrt(5))/2.  This script finds a split-prime
specialization whose augmented rank, together with the exact base rank 61,
determines the characteristic-zero rank of the combined system.

This is a finite-dimensional edge-kernel reduction only.  It is not an
exclusion of either full Newton-edge branch.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ

from f2_75_125_2plus2_terminal_s1_weight49_modrank_20261018 import (
    rank_mod_p,
    specialize_entry,
)
from f2_75_125_power_relation_second_band import build_matrix
from f2_75_125_power_relation_pq_nextbands import weight49_parameters


def main() -> None:
    t = sp.symbols("t")
    sqrt5 = sp.sqrt(5)
    alpha = (3 + sqrt5) / 2
    h = sp.expand((t - 1) ** 2 * (t - alpha) ** 2)
    slope = (4 * alpha - 11) / 3
    intercept = 3 - alpha
    v = sp.expand(t * (t - 1) * (t - alpha) * (slope * t + intercept))
    matrix = build_matrix(h, v)
    if matrix.shape != (72, 64):
        raise AssertionError(f"unexpected matrix shape: {matrix.shape}")

    parameters, _, _, exact_kernel = weight49_parameters(
        h, v, QQ.algebraic_field(sqrt5), "2+2"
    )
    exact_kernel_dimension = len(parameters)
    if exact_kernel.shape != (exact_kernel_dimension, matrix.cols):
        raise AssertionError(f"unexpected exact kernel shape: {exact_kernel.shape}")
    if exact_kernel_dimension != 3:
        raise AssertionError(f"expected exact base kernel dimension 3, got {exact_kernel_dimension}")
    exact_rank_base = matrix.cols - exact_kernel_dimension

    row_1 = [sp.Integer(1)] * 4 + [sp.Integer(0)] * 60
    best = None
    attempts = []
    for prime in (11, 19, 29, 31, 41, 59, 61, 71, 79, 89, 101):
        roots = [x for x in range(prime) if x * x % prime == 5 % prime]
        for root in roots:
            try:
                rows = [
                    [specialize_entry(entry, sqrt5, prime, root) for entry in row]
                    for row in matrix.tolist()
                ]
            except AssertionError:
                continue
            alpha_mod = ((3 + root) * pow(2, -1, prime)) % prime
            row_alpha = [pow(alpha_mod, j, prime) for j in range(4)] + [0] * 60
            ranks = (
                rank_mod_p(rows, prime),
                rank_mod_p(rows + [[int(x) for x in row_1]], prime),
                rank_mod_p(rows + [[int(x) for x in row_1], row_alpha], prime),
            )
            attempts.append((prime, root, alpha_mod, ranks))
            if best is None or ranks[2] > best[3][2]:
                best = (prime, root, alpha_mod, ranks)
            if ranks[2] == exact_rank_base + 2:
                break
        if best is not None and best[3][2] == exact_rank_base + 2:
            break

    if best is None:
        raise AssertionError("no valid split-prime specialization")
    prime, root, alpha_mod, ranks = best
    print(f"matrix_shape={matrix.shape}")
    print(f"exact_base_rank={exact_rank_base}")
    print(f"prime={prime}; sqrt5_image={root}; alpha_image={alpha_mod}")
    print(f"modular_ranks_base,S1,S1_and_Salpha={ranks}")
    print(f"exact_rank_upper_bound={exact_rank_base + 2}")
    print(f"exact_augmented_rank={ranks[2]}" if ranks[2] == exact_rank_base + 2
          else "exact_augmented_rank=not determined by this prime")
    print(f"exact_kernel_dimension={matrix.cols - ranks[2]}"
          if ranks[2] == exact_rank_base + 2
          else "exact_kernel_dimension=not determined")
    print("constraints=S(1)=0 and S(alpha)=0")
    print("scope=weight-49 kernel only; no full-branch exclusion")

    receipt = {
        "branch": "2+2",
        "matrix_shape": list(matrix.shape),
        "exact_base_kernel_dimension": exact_kernel_dimension,
        "exact_base_rank": exact_rank_base,
        "split_prime": prime,
        "sqrt5_image": root,
        "alpha_image": alpha_mod,
        "modular_ranks_base_S1_S1_and_Salpha": list(ranks),
        "exact_augmented_rank": ranks[2],
        "exact_augmented_kernel_dimension": matrix.cols - ranks[2],
        "constraints": ["S(1)=0", "S(alpha)=0"],
        "scope": "exact dimension reduction in weight-49 kernel; not a full-branch exclusion",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    Path("research/f2_75_125_dual_root_weight49_rank_receipt_20261018.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    print("receipt=research/f2_75_125_dual_root_weight49_rank_receipt_20261018.json")


if __name__ == "__main__":
    main()
