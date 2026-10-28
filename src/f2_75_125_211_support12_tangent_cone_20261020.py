#!/usr/bin/env python3
"""Exact FLINT-backed tangent-cone split of the 2+1+1 support-12 slice.

This replays the detached rational equations against their SHA-256 receipt,
extracts the ordinary quadratic initial forms on the unique support-12 slice,
and computes the two linear-factor charts of the first quadratic form.  These
are tangent-cone data for a necessary subideal, not an exclusion certificate.
"""

from __future__ import annotations

import hashlib
import json
import os
from collections import defaultdict
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix

from f2_75_125_211_degree_structure import streamed_array


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_p43zero_depth13_constraints.jsonfrag"
PROFILE = HERE / "f2_75_125_211_p43zero_depth13_support12_quadratic_rank_20261020.json"
RECEIPT = HERE / "f2_75_125_211_support12_tangent_cone_20261020.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def quadratic_matrix(expr: sp.Expr, variables: tuple[sp.Symbol, ...]) -> sp.Matrix:
    poly = sp.Poly(expr, *variables, domain=QQ)
    matrix = sp.zeros(len(variables))
    for monomial, coefficient in poly.terms():
        if sum(monomial) != 2:
            raise ValueError("quadratic matrix received a non-quadratic expression")
        indices = [i for i, exponent in enumerate(monomial) for _ in range(exponent)]
        i, j = indices
        if i == j:
            matrix[i, i] = coefficient
        else:
            matrix[i, j] = matrix[j, i] = coefficient / 2
    return matrix


def main() -> None:
    profile = json.loads(PROFILE.read_text(encoding="utf-8"))
    source_hash = sha256(SOURCE)
    if source_hash != profile["source_sha256"]:
        raise SystemExit("detached source hash does not match the quadratic-rank report")

    support = tuple(int(name[1:]) for name in profile["support_variables"])
    variables = tuple(sp.Symbol(f"p{i}") for i in support)
    local_index = {global_index: local for local, global_index in enumerate(support)}
    monomials = tuple((i, j) for i in range(len(support)) for j in range(i, len(support)))
    column_index = {monomial: column for column, monomial in enumerate(monomials)}

    rows: list[list[sp.Expr]] = []
    tags: list[str] = []
    expressions: list[sp.Expr] = []
    for tag, saved in streamed_array(SOURCE):
        active: set[int] = set()
        quadratic: defaultdict[tuple[int, int], sp.Expr] = defaultdict(lambda: sp.S.Zero)
        for sparse, values in saved:
            exponents = {
                index - 1: exponent
                for index, exponent in sparse
                if 1 <= index <= 75 and exponent
            }
            active.update(exponents)
            if sum(exponents.values()) != 2:
                continue
            if len(values) != 1:
                raise ValueError("support-12 source is expected to be rational")
            active_items = list(exponents.items())
            if len(active_items) == 1 and active_items[0][1] == 2:
                pair = (active_items[0][0], active_items[0][0])
            elif len(active_items) == 2 and all(exponent == 1 for _, exponent in active_items):
                pair = tuple(sorted((active_items[0][0], active_items[1][0])))
            else:
                raise ValueError(f"unexpected quadratic monomial in {tag}")
            quadratic[pair] += sp.Rational(values[0])

        if tuple(sorted(active)) != support:
            continue
        vector = [sp.S.Zero] * len(monomials)
        expression = sp.S.Zero
        for global_pair, coefficient in quadratic.items():
            if coefficient:
                pair = tuple(sorted(local_index[i] for i in global_pair))
                vector[column_index[pair]] = coefficient
                i, j = pair
                expression += coefficient * variables[i] * variables[j]
        rows.append(vector)
        tags.append(tag)
        expressions.append(sp.expand(expression))

    if len(rows) != 147:
        raise SystemExit(f"expected 147 support-12 rows, got {len(rows)}")

    exact = DomainMatrix.from_Matrix(sp.Matrix(rows)).convert_to(QQ)
    _, pivots = exact.rref()
    independent_rows = list(exact.transpose().rref()[1])
    if len(pivots) != 2 or independent_rows != [0, 1]:
        raise SystemExit(f"unexpected tangent-cone rank/chart: {len(pivots)}, {independent_rows}")

    q1 = expressions[independent_rows[0]]
    q2 = expressions[independent_rows[1]]
    factor1 = sp.factor_list(q1)
    if len(factor1[1]) != 2 or any(power != 1 for _, power in factor1[1]):
        raise SystemExit("first quadratic does not have the expected two distinct linear factors")
    linear_factors = [sp.expand(factor) for factor, _ in factor1[1]]

    charts = []
    for linear in linear_factors:
        poly = sp.Poly(linear, *variables, domain=QQ)
        pivot_index = next(i for i, variable in enumerate(variables) if poly.coeff_monomial(variable))
        pivot = variables[pivot_index]
        pivot_coefficient = poly.coeff_monomial(pivot)
        substitution = sp.solve(linear, pivot)[0]
        remaining = tuple(variable for variable in variables if variable != pivot)
        restricted = sp.factor(q2.subs(pivot, substitution))
        restricted_matrix = quadratic_matrix(restricted, remaining)
        rank = restricted_matrix.rank()
        charts.append(
            {
                "linear_factor": str(linear),
                "pivot_variable": str(pivot),
                "pivot_coefficient": str(pivot_coefficient),
                "pivot_substitution": str(substitution),
                "remaining_variables": [str(variable) for variable in remaining],
                "restricted_second_quadratic": str(restricted),
                "restricted_quadratic_matrix_rank": int(rank),
                "restricted_quadratic_matrix_size": len(remaining),
                "restricted_factorization": [
                    [str(factor), power] for factor, power in sp.factor_list(restricted)[1]
                ],
            }
        )

    result = {
        "source": SOURCE.name,
        "source_sha256": source_hash,
        "quadratic_rank_report_sha256": sha256(PROFILE),
        "scope": "ordinary tangent cone of the support-12 subideal only; no affine or branch exclusion",
        "branch": "2+1+1",
        "coefficient_field": "QQ",
        "support_variables": [str(variable) for variable in variables],
        "constraints_on_support": len(rows),
        "quadratic_matrix_shape": [len(rows), len(monomials)],
        "quadratic_monomials": len(monomials),
        "quadratic_initial_form_rank": len(pivots),
        "independent_row_indices": independent_rows,
        "independent_row_tags": [tags[i] for i in independent_rows],
        "first_quadratic": str(q1),
        "first_quadratic_factorization": [
            [str(factor), power] for factor, power in factor1[1]
        ],
        "second_quadratic": str(q2),
        "linear_factor_charts": charts,
        "arithmetic_backend": os.environ.get("SYMPY_GROUND_TYPES", "default"),
        "warning": "The higher ordinary-degree terms may bend or remove these tangent branches; each chart is not an exact branch of the affine support-12 ideal.",
    }
    encode = lambda value: int(value) if isinstance(value, sp.Integer) else str(value)
    rendered = json.dumps(result, indent=2, default=encode)
    RECEIPT.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
