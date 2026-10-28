#!/usr/bin/env python3
"""Exact rank/factor audit of the independent quadratic tangent forms."""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix

from f2_75_125_211_degree_structure import streamed_array


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_p43zero_depth13_constraints.jsonfrag"
INPUT = HERE / "f2_75_125_211_global_quadratic_initial_rank_20261020.json"
OUTPUT = HERE / "f2_75_125_211_global_quad_form_rank_20261021.json"
LOW_ZERO = {0, 1, 2, 3, 4}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def quadratic_part(saved, *, boundary: bool):
    result = defaultdict(Fraction)
    for sparse, values in saved:
        exponents = {}
        vanishes = False
        for raw_index, exponent in sparse:
            if not exponent or not 1 <= raw_index <= 75:
                continue
            variable = raw_index - 1
            if boundary and variable in LOW_ZERO:
                vanishes = True
                break
            exponents[variable] = exponent
        if vanishes or sum(exponents.values()) != 2:
            continue
        if len(values) != 1:
            raise ValueError("expected rational quadratic coefficients")
        monomial = tuple(sorted((i, e) for i, e in exponents.items() if e))
        result[monomial] += Fraction(values[0])
    return {monomial: coefficient for monomial, coefficient in result.items() if coefficient}


def analyze_form(tag, form, variables, *, boundary):
    local = {variable: i for i, variable in enumerate(variables)}
    entries = {}
    symbols = sp.symbols(" ".join(f"p{i}" for i in variables))
    expression = sp.S.Zero
    for monomial, coefficient in form.items():
        expanded = []
        term = sp.Rational(coefficient.numerator, coefficient.denominator)
        for variable, exponent in monomial:
            if variable not in local:
                raise AssertionError(f"unexpected variable p{variable} in restricted form")
            expanded.extend([local[variable]] * exponent)
            term *= symbols[local[variable]] ** exponent
        expression += term
        if len(expanded) != 2:
            raise AssertionError("unexpected quadratic degree")
        i, j = sorted(expanded)
        if i == j:
            entries[(i, i)] = entries.get((i, i), 0) + sp.Rational(
                coefficient.numerator, coefficient.denominator
            )
        else:
            half = sp.Rational(coefficient.numerator, 2 * coefficient.denominator)
            entries[(i, j)] = entries.get((i, j), 0) + half
            entries[(j, i)] = entries.get((j, i), 0) + half
    matrix = sp.SparseMatrix(len(variables), len(variables), entries)
    rank = len(DomainMatrix.from_Matrix(matrix).convert_to(QQ).rref()[1])
    factorization = str(sp.factor(expression)) if rank <= 2 and expression else None
    return {
        "tag": tag,
        "quadratic_term_count": len(form),
        "support": [f"p{i}" for i in variables
                    if any(any(j == i for j, _ in monomial) for monomial in form)],
        "symmetric_matrix_rank_over_Q": rank,
        "factors_if_rank_at_most_2": factorization,
        "restricted_boundary": boundary,
    }


def main() -> None:
    if digest(SOURCE) != json.loads(INPUT.read_text(encoding="utf-8"))["source_sha256"]:
        raise SystemExit("quadratic input does not match its source hash")
    tags = json.loads(INPUT.read_text(encoding="utf-8"))["independent_input_row_tags"]
    wanted = set(tags)
    selected = {}
    for tag, saved in streamed_array(SOURCE):
        if tag in wanted:
            if tag in selected:
                raise ValueError(f"duplicate input tag {tag}")
            selected[tag] = saved
    if set(selected) != wanted:
        raise ValueError("one or more independent quadratic rows were not found")

    variables_full = [i for i in range(43) if i != 36]
    variables_boundary = [i for i in variables_full if i not in LOW_ZERO]
    full_forms = [analyze_form(tag, quadratic_part(selected[tag], boundary=False),
                               variables_full, boundary=False) for tag in tags]
    boundary_forms = [analyze_form(tag, quadratic_part(selected[tag], boundary=True),
                                   variables_boundary, boundary=True) for tag in tags]
    receipt = {
        "source": SOURCE.name,
        "source_sha256": digest(SOURCE),
        "independent_row_tags": tags,
        "full_variable_count": len(variables_full),
        "boundary": "p0=p1=p2=p3=p4=0",
        "boundary_variable_count": len(variables_boundary),
        "full_forms": full_forms,
        "boundary_forms": boundary_forms,
        "scope": "individual quadratic tangent-form ranks only; does not exclude any affine or projective branch",
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "full_form_ranks": [item["symmetric_matrix_rank_over_Q"] for item in full_forms],
        "boundary_form_ranks": [item["symmetric_matrix_rank_over_Q"] for item in boundary_forms],
        "factorizations": [item["factors_if_rank_at_most_2"] for item in full_forms
                           if item["factors_if_rank_at_most_2"] is not None],
        "receipt": str(OUTPUT),
        "scope": receipt["scope"],
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
