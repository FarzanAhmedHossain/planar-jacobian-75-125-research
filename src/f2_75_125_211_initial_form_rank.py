#!/usr/bin/env python3
"""Exact rank of the lowest ordinary-degree forms in a support slice.

The output is a tangent-cone diagnostic for a detached necessary ideal, not
an exclusion of the F2 branch. In particular, an isolated origin in a
nonhomogeneous affine system does not rule out other components.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
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


def coefficient(values):
    root = sp.sqrt(5)
    return sp.Add(*(sp.Rational(value) * root ** (len(values) - 1 - i)
                    for i, value in enumerate(values)))


def audit(source: Path, nparameters: int, target_support_size: int) -> dict:
    target_support = None
    rows = []
    tags = []
    source_is_quadratic = False
    for tag, polynomial in streamed_array(source):
        support = set()
        quadratic = defaultdict(lambda: sp.S.Zero)
        for sparse, values in polynomial:
            exponents = [(i - 1, e) for i, e in sparse
                         if 1 <= i <= nparameters and e]
            support.update(i for i, _ in exponents)
            if sum(e for _, e in exponents) != 2:
                continue
            if len(values) > 1:
                source_is_quadratic = True
            if len(exponents) == 1 and exponents[0][1] == 2:
                monomial = (exponents[0][0], exponents[0][0])
            elif len(exponents) == 2 and all(e == 1 for _, e in exponents):
                monomial = tuple(sorted(i for i, _ in exponents))
            else:
                raise ValueError(f"unexpected degree-two monomial in {tag}: {exponents}")
            quadratic[monomial] += coefficient(values)
        if len(support) != target_support_size:
            continue
        support_tuple = tuple(sorted(support))
        if target_support is None:
            target_support = support_tuple
        elif support_tuple != target_support:
            raise ValueError("more than one support set has the requested size")
        rows.append(quadratic)
        tags.append(tag)

    if target_support is None:
        raise ValueError("no constraints matched target support size")
    local = {global_i: local_i for local_i, global_i in enumerate(target_support)}
    monomials = [(i, j) for i in range(len(target_support))
                 for j in range(i, len(target_support))]
    columns = {monomial: j for j, monomial in enumerate(monomials)}
    matrix = []
    for row in rows:
        values = [sp.S.Zero] * len(monomials)
        for (global_i, global_j), value in row.items():
            if value:
                pair = tuple(sorted((local[global_i], local[global_j])))
                values[columns[pair]] = value
        matrix.append(values)

    field = QQ.algebraic_field(sp.sqrt(5)) if source_is_quadratic else QQ
    exact = DomainMatrix.from_Matrix(sp.Matrix(matrix)).convert_to(field)
    rref, pivots = exact.rref()
    selected_rows = list(exact.transpose().rref()[1])
    minor_det = None
    if len(selected_rows) == len(monomials):
        minor_det = str(exact.extract(selected_rows, list(range(len(monomials)))).det())

    return {
        "source": str(source),
        "source_sha256": sha256(source),
        "branch": "2+1+1",
        "scope": "rank of quadratic initial forms in one exact-support slice; not a branch exclusion",
        "nparameters": nparameters,
        "support_size": target_support_size,
        "support_variables": [f"p{i}" for i in target_support],
        "constraints_in_slice": len(rows),
        "quadratic_monomial_count": len(monomials),
        "quadratic_initial_form_rank": len(pivots),
        "tangent_space_dimension_in_slice": len(monomials) - len(pivots),
        "exact_coefficient_field": str(field),
        "independent_row_indices": selected_rows,
        "full_rank_minor_determinant": minor_det,
        "constraint_tag_samples": tags[:12],
        "warning": "These constraints have mixed ordinary degrees; tangent-cone data alone does not exclude nonzero affine solutions.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("nparameters", type=int)
    parser.add_argument("support_size", type=int)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    result = audit(args.source, args.nparameters, args.support_size)
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
