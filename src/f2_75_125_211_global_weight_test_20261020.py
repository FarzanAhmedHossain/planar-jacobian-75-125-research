#!/usr/bin/env python3
"""Find exact common positive gradings for all detached 2+1+1 constraints.

This checks whether the full 754-equation depth-13 necessary subsystem is
weighted homogeneous under one positive torus action. It is structural
reconnaissance only: a grading is not a branch exclusion.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix

from f2_75_125_211_degree_structure import streamed_array


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_p43zero_depth13_constraints.jsonfrag"
INCIDENCE = HERE / "f2_75_125_211_p43zero_depth13_incidence_20261020.json"
REPORT = HERE / "f2_75_125_211_global_weight_test_20261020.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    incidence = json.loads(INCIDENCE.read_text(encoding="utf-8"))
    source_hash = sha256(SOURCE)
    if source_hash != incidence["source_sha256"]:
        raise SystemExit("constraint source hash differs from incidence receipt")
    active = tuple(
        int(name[1:])
        for name, count in incidence["variable_occurrence_counts"].items()
        if count
    )
    local = {global_i: local_i for local_i, global_i in enumerate(active)}
    integer_weights = {}
    for index in range(0, 5):
        integer_weights[index] = (1, 2, 3, 3, 3)[index]
    integer_weights.update({index: 4 for index in range(5, 13)})
    integer_weights.update({index: 5 for index in range(13, 20)})
    integer_weights.update({index: 6 for index in range(20, 25)})
    integer_weights.update({index: 7 for index in range(25, 30)})
    integer_weights.update({index: 8 for index in range(30, 36)})
    integer_weights.update({index: 9 for index in range(37, 43)})
    differences: set[tuple[tuple[int, int], ...]] = set()
    weighted_degree_histogram: Counter[int] = Counter()
    constraint_count = term_count = 0
    for tag, polynomial in streamed_array(SOURCE):
        constraint_count += 1
        base = None
        constraint_degree = None
        for sparse, _values in polynomial:
            vector = {}
            for generator, exponent in sparse:
                parameter = generator - 1
                if 0 <= parameter < 75 and exponent:
                    if parameter not in local:
                        raise ValueError(f"unexpected active parameter p{parameter} in {tag}")
                    vector[local[parameter]] = exponent
            degree = sum(integer_weights[active[index]] * exponent
                         for index, exponent in vector.items())
            if constraint_degree is None:
                constraint_degree = degree
            elif degree != constraint_degree:
                raise ValueError(f"global weights do not grade {tag}")
            if base is None:
                base = vector
            else:
                delta = tuple(
                    (i, value)
                    for i in sorted(set(base) | set(vector))
                    if (value := vector.get(i, 0) - base.get(i, 0))
                )
                if delta:
                    differences.add(delta)
            term_count += 1
        if constraint_degree is not None:
            weighted_degree_histogram[constraint_degree] += 1
    rows = []
    for sparse in sorted(differences):
        row = [0] * len(active)
        for index, value in sparse:
            row[index] = value
        rows.append(row)
    matrix = DomainMatrix(rows, (len(rows), len(active)), QQ)
    reduced, pivots = matrix.rref()
    free = [j for j in range(len(active)) if j not in pivots]
    reduced_rows = reduced.to_list()
    basis = []
    for column in free:
        vector = [sp.S.Zero] * len(active)
        vector[column] = sp.S.One
        for row, pivot in enumerate(pivots):
            vector[pivot] = -reduced_rows[row][column]
        basis.append(vector)

    positive = None
    if len(basis) == 1:
        candidate = basis[0]
        if all(value > 0 for value in candidate):
            positive = candidate
        elif all(value < 0 for value in candidate):
            positive = [-value for value in candidate]
    result = {
        "source": SOURCE.name,
        "source_sha256": source_hash,
        "branch": "2+1+1",
        "scope": "common torus grading test for the entire detached depth-13 necessary subsystem; not exclusion",
        "constraints": constraint_count,
        "terms_scanned": term_count,
        "active_variables": [f"p{i}" for i in active],
        "positive_integer_weights_by_variable": {
            f"p{i}": integer_weights[i] for i in active
        },
        "weighted_degree_histogram": dict(sorted(weighted_degree_histogram.items())),
        "distinct_exponent_difference_rows": len(differences),
        "weight_equation_rank_over_Q": len(pivots),
        "weight_nullity": len(free),
        "nullspace_basis": [[str(value) for value in vector] for vector in basis],
        "common_positive_weight_found": positive is not None,
        "positive_weights": ([str(value) for value in positive] if positive else None),
        "warning": "A common positive grading can enable projective methods but does not prove the projective zero set empty or exclude the F2 branch.",
    }
    REPORT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
