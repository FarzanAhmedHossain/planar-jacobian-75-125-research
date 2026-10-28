#!/usr/bin/env python3
"""Find exact positive-grading candidates for one detached support slice.

For each polynomial, equal weighted degree requires w.(m-m0)=0 for every
monomial m relative to a reference monomial m0. This script computes that
integer linear system exactly. A positive solution would identify a torus
scaling of the slice; absence of a positive solution is only a statement
about this chosen support subsystem.
"""

from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix

from f2_75_125_211_degree_structure import streamed_array


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def exact_audit(source: Path, incidence_report: Path, nparameters: int,
                support_size: int) -> dict:
    incidence = json.loads(incidence_report.read_text(encoding="utf-8"))
    source_hash = sha256(source)
    if source_hash != incidence["source_sha256"]:
        raise ValueError("incidence report does not match the exact source hash")
    candidates = incidence["most_repeated_supports"][str(support_size)]
    if len(candidates) != 1:
        raise ValueError("support size is not represented by one exact variable set")
    support = tuple(candidates[0]["variables"])
    if len(support) != support_size:
        raise ValueError("support report has the wrong variable count")

    difference_rows: set[tuple[int, ...]] = set()
    selected = 0
    checked_monomials = 0
    for _tag, polynomial in streamed_array(source):
        variable_support = set()
        local_differences = set()
        base = None
        for sparse, _coefficient in polynomial:
            exponents = {i - 1: e for i, e in sparse
                         if 1 <= i <= nparameters and e}
            variable_support.update(exponents)
            vector = tuple(exponents.get(i, 0) for i in support)
            if base is None:
                base = vector
            else:
                delta = tuple(a - b for a, b in zip(vector, base))
                if any(delta):
                    local_differences.add(delta)
            checked_monomials += 1
        if tuple(sorted(variable_support)) == support:
            selected += 1
            difference_rows.update(local_differences)

    ordered_rows = sorted(difference_rows)
    matrix = DomainMatrix.from_Matrix(sp.Matrix(ordered_rows)).convert_to(QQ)
    reduced, pivots = matrix.rref()
    rank = len(pivots)
    ncols = len(support)
    reduced_matrix = reduced.to_Matrix()
    free_columns = [j for j in range(ncols) if j not in pivots]
    nullspace_basis = []
    for free in free_columns:
        vector = [sp.S.Zero] * ncols
        vector[free] = sp.S.One
        for row, pivot in enumerate(pivots):
            vector[pivot] = -reduced_matrix[row, free]
        nullspace_basis.append(vector)

    positive_vector = None
    search_status = "not_assessed"
    if not nullspace_basis:
        search_status = "no_nonzero_weight_vector"
    elif len(nullspace_basis) == 1:
        vector = nullspace_basis[0]
        if all(value > 0 for value in vector):
            positive_vector = vector
        elif all(value < 0 for value in vector):
            positive_vector = [-value for value in vector]
        search_status = "complete_one_dimensional_nullspace_check"
    elif len(nullspace_basis) <= 4:
        search_status = "finite_small_coefficient_search"
        for coefficients in product(range(-3, 4), repeat=len(nullspace_basis)):
            if not any(coefficients):
                continue
            vector = [sum(coefficients[k]*nullspace_basis[k][j]
                          for k in range(len(nullspace_basis)))
                      for j in range(ncols)]
            if all(value > 0 for value in vector):
                positive_vector = vector
                break

    return {
        "source": str(source),
        "source_sha256": source_hash,
        "scope": "exact equal-weight-degree test for the support-slice polynomials; not an exclusion",
        "support_size": support_size,
        "support_variables": [f"p{i}" for i in support],
        "selected_constraints": selected,
        "checked_monomials": checked_monomials,
        "distinct_exponent_difference_rows": len(ordered_rows),
        "weight_equation_rank_over_Q": rank,
        "weight_nullity": ncols-rank,
        "rational_nullspace_basis": [[str(value) for value in vector]
                                     for vector in nullspace_basis],
        "positive_weight_vector_found": positive_vector is not None,
        "positive_weight_vector": ([str(value) for value in positive_vector]
                                   if positive_vector is not None else None),
        "positive_search_status": search_status,
        "positive_grading_note": "A found vector is exact and makes every monomial in each selected equation have equal positive weight. A negative result from the finite search is not a proof that no positive grading exists.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("incidence_report", type=Path)
    parser.add_argument("nparameters", type=int)
    parser.add_argument("support_size", type=int)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    result = exact_audit(args.source, args.incidence_report,
                         args.nparameters, args.support_size)
    args.report.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
