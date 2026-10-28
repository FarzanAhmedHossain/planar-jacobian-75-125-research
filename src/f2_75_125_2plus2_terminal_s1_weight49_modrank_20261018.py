"""Exact finite-field certificate for the F2 terminal contact S(1)=0.

The weight-49 matrix is the one in f2_75_125_power_relation_second_band.py.
For alpha=(3+sqrt(5))/2 we compare its rank with the matrix augmented by
S(1)=s0+s1+s2+s3=0.  The specialization sqrt(5)->4 in F_11 is a valid
field embedding because 4^2=5 mod 11.  Thus modular rank gives a rigorous
lower bound on the characteristic-zero algebraic rank.  The pre-existing
exact audit gives rank 61 for the unaugmented matrix; since adding one row
can raise rank by at most one, modular rank 62 for the augmented matrix
proves the exact augmented rank is 62 and the constrained kernel dimension
is 2.

This only imposes the necessary distinguished-root contact S(1)=0.  It is
not a branch or degree-pair exclusion.
"""

from __future__ import annotations

import sympy as sp

from f2_75_125_power_relation_second_band import build_matrix


def rank_mod_p(rows: list[list[int]], p: int) -> int:
    """RREF rank over the prime field, using only integer arithmetic."""
    if not rows:
        return 0
    matrix = [[entry % p for entry in row] for row in rows]
    height = len(matrix)
    width = len(matrix[0])
    pivot_row = 0
    for column in range(width):
        pivot = next((r for r in range(pivot_row, height)
                      if matrix[r][column] % p), None)
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], -1, p)
        matrix[pivot_row] = [(v * inverse) % p for v in matrix[pivot_row]]
        for r in range(height):
            if r == pivot_row:
                continue
            factor = matrix[r][column]
            if factor:
                matrix[r] = [(x - factor * y) % p
                             for x, y in zip(matrix[r], matrix[pivot_row])]
        pivot_row += 1
        if pivot_row == height:
            break
    return pivot_row


def specialize_entry(expr: sp.Expr, sqrt5: sp.Symbol, p: int, sqrt5_image: int) -> int:
    value = sp.cancel(expr.subs(sqrt5, sqrt5_image))
    if not value.is_Rational:
        raise AssertionError((expr, value))
    numerator = int(value.p) % p
    denominator = int(value.q) % p
    if denominator == 0:
        raise AssertionError(("bad denominator under specialization", expr, value))
    return numerator * pow(denominator, -1, p) % p


def main() -> None:
    t = sp.symbols("t")
    sqrt5 = sp.sqrt(5)
    alpha = (3 + sqrt5) / 2
    h = sp.expand((t - 1) ** 2 * (t - alpha) ** 2)
    slope = (4 * alpha - 11) / 3
    intercept = 3 - alpha
    v = sp.expand(t * (t - 1) * (t - alpha) * (slope * t + intercept))
    matrix = build_matrix(h, v)
    raw_rows = matrix.tolist()
    eval_row = [[sp.Integer(1)] * 4 + [sp.Integer(0)] * (matrix.cols - 4)]
    augmented_rows = raw_rows + eval_row
    assert matrix.cols == 64
    certificate = None
    # Choose a split prime where the chosen nonzero minors survive reduction.
    # Almost all split primes work; denominators are checked in specialize_entry.
    for p in (11, 19, 29, 31, 41, 59, 61, 71, 79, 89, 101):
        if any(pow(x, 2, p) == 5 % p for x in range(p)) is False:
            continue
        root = next(x for x in range(p) if pow(x, 2, p) == 5 % p)
        try:
            base_mod = [[specialize_entry(x, sqrt5, p, root) for x in row]
                        for row in raw_rows]
            aug_mod = [[specialize_entry(x, sqrt5, p, root) for x in row]
                       for row in augmented_rows]
        except AssertionError:
            continue
        base_rank_mod = rank_mod_p(base_mod, p)
        augmented_rank_mod = rank_mod_p(aug_mod, p)
        if base_rank_mod == 61 and augmented_rank_mod == 62:
            certificate = (p, root, base_rank_mod, augmented_rank_mod)
            break
    if certificate is None:
        raise AssertionError("no good prime among the checked split primes")
    p, root, base_rank_mod, augmented_rank_mod = certificate
    print(f"specialization=F_{p}, sqrt5->{root}")
    print(f"matrix_shape={matrix.rows}x{matrix.cols}")
    print(f"base_rank_mod={base_rank_mod}")
    print(f"augmented_rank_mod={augmented_rank_mod}")
    print("source_exact_base_rank=61 (f2_75_125_power_relation_second_band.py)")
    print("exact_augmented_rank=62 (modular lower bound and one-row upper bound)")
    print("terminal_contact=S(1)=0")
    print("exact_constrained_kernel_dimension=2")
    print("branch_excluded=False")


if __name__ == "__main__":
    main()
