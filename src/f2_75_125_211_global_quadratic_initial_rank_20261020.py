#!/usr/bin/env python3
"""Exact quadratic initial-form rank for the full detached 2+1+1 slice.

This extracts the lowest ordinary-degree forms from all 754 necessary
constraints. Rank is only a tangent-cone diagnostic; it is not a branch
exclusion.
"""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix

from f2_75_125_211_degree_structure import streamed_array


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_p43zero_depth13_constraints.jsonfrag"
INCIDENCE = HERE / "f2_75_125_211_p43zero_depth13_incidence_20261020.json"
REPORT = HERE / "f2_75_125_211_global_quadratic_initial_rank_20261020.json"


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
    monomials = [(i, j) for i in range(len(active)) for j in range(i, len(active))]
    columns = {pair: index for index, pair in enumerate(monomials)}
    rows = []
    tags = []
    terms_scanned = 0
    for tag, polynomial in streamed_array(SOURCE):
        row = [QQ.zero] * len(monomials)
        for sparse, values in polynomial:
            powers = {
                local[generator - 1]: exponent
                for generator, exponent in sparse
                if 1 <= generator <= 75 and exponent
            }
            terms_scanned += 1
            if sum(powers.values()) != 2:
                continue
            if len(values) != 1:
                raise ValueError(f"expected rational coefficient in {tag}")
            if len(powers) == 1:
                (i, exponent), = powers.items()
                if exponent != 2:
                    raise ValueError(f"unexpected quadratic monomial in {tag}")
                pair = (i, i)
            elif len(powers) == 2 and all(exponent == 1 for exponent in powers.values()):
                pair = tuple(sorted(powers))
            else:
                raise ValueError(f"unexpected quadratic monomial in {tag}")
            row[columns[pair]] += QQ.convert(sp.Rational(values[0]))
        if any(row):
            rows.append(row)
            tags.append(tag)
    matrix = DomainMatrix(rows, (len(rows), len(monomials)), QQ)
    _, pivots = matrix.rref()
    selected_rows = list(matrix.transpose().rref()[1])
    result = {
        "source": SOURCE.name,
        "source_sha256": source_hash,
        "branch": "2+1+1",
        "scope": "ordinary quadratic initial forms from the full detached depth-13 necessary subsystem; tangent-cone diagnostic only",
        "constraints": incidence["constraints"],
        "constraints_with_nonzero_quadratic_initial_form": len(rows),
        "terms_scanned": terms_scanned,
        "active_variables": [f"p{i}" for i in active],
        "quadratic_monomials": len(monomials),
        "quadratic_initial_form_rank_over_Q": len(pivots),
        "tangent_cone_quadratic_span_codimension": len(monomials) - len(pivots),
        "independent_input_row_indices": selected_rows,
        "independent_input_row_tags": [tags[i] for i in selected_rows],
        "warning": "A tangent-cone rank result alone does not exclude affine solutions or establish emptiness of the full 2+1+1 Keller branch.",
    }
    REPORT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
