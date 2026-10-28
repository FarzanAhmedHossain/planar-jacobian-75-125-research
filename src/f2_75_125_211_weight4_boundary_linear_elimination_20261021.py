#!/usr/bin/env python3
"""Exact linear elimination on one weighted boundary of the 2+1+1 lift.

We set p0,...,p4=0 (the variables of weights 1,2,3) in the saved depth-13
necessary equations.  Because every equation has weighted degree 8 through
13, p30,...,p35 can occur linearly only in weight 8, and p37,...,p42 only
linearly in weight 9.  This script solves those two affine-linear blocks over
Q, substitutes the pivot variables into the complete saved equation array,
and profiles the exact residual system.  It is a boundary-pruning calculation,
not a whole-branch exclusion.
"""

from __future__ import annotations

from collections import Counter
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
OUTPUT = HERE / "f2_75_125_211_weight4_boundary_linear_elimination_20261021.json"
WEIGHT = {0: 1, 1: 2, 2: 3, 3: 3, 4: 3}
WEIGHT.update({i: 4 for i in range(5, 13)})
WEIGHT.update({i: 5 for i in range(13, 20)})
WEIGHT.update({i: 6 for i in range(20, 25)})
WEIGHT.update({i: 7 for i in range(25, 30)})
WEIGHT.update({i: 8 for i in range(30, 36)})
WEIGHT.update({i: 9 for i in range(37, 43)})
ZERO = {0, 1, 2, 3, 4}
HIGH_BY_DEGREE = {8: tuple(range(30, 36)), 9: tuple(range(37, 43))}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def add_term(poly: dict[tuple[tuple[int, int], ...], Fraction],
             mon: tuple[tuple[int, int], ...], value: Fraction) -> None:
    if not value:
        return
    updated = poly.get(mon, Fraction(0)) + value
    if updated:
        poly[mon] = updated
    else:
        poly.pop(mon, None)


def merge_monomials(left, right):
    result = dict(left)
    for variable, exponent in right:
        result[variable] = result.get(variable, 0) + exponent
    return tuple(sorted((i, e) for i, e in result.items() if e))


def parse_boundary(saved):
    poly = {}
    for sparse, values in saved:
        exponents = {}
        vanishes = False
        for raw_index, exponent in sparse:
            if not exponent:
                continue
            if 1 <= raw_index <= 5:
                vanishes = True
                break
            if 6 <= raw_index <= 43:
                variable = raw_index - 1
                exponents[variable] = exponent
            elif raw_index > 43:
                raise ValueError(f"unexpected active parameter index {raw_index}")
        if vanishes:
            continue
        if len(values) != 1:
            raise ValueError("2+1+1 input is expected to have rational coefficients")
        mon = tuple(sorted(exponents.items()))
        add_term(poly, mon, Fraction(values[0]))
    return poly


def weighted_degree(poly):
    degrees = {sum(WEIGHT[i] * e for i, e in mon) for mon in poly}
    if len(degrees) > 1:
        raise ValueError(f"boundary restriction lost homogeneity: {degrees}")
    return next(iter(degrees)) if degrees else None


def add_scaled(left, right, scalar):
    result = dict(left)
    for mon, coefficient in right.items():
        add_term(result, mon, scalar * coefficient)
    return result


def eliminate_block(polynomials, variables, degree):
    """Gauss-Jordan eliminate A*y+b=0, with polynomial b over Q."""
    columns = {variable: j for j, variable in enumerate(variables)}
    matrix = []
    forcing = []
    for poly in polynomials:
        row = [Fraction(0)] * len(variables)
        other = {}
        for mon, coefficient in poly.items():
            if len(mon) == 1 and mon[0][1] == 1 and mon[0][0] in columns:
                row[columns[mon[0][0]]] += coefficient
            else:
                other[mon] = coefficient
        matrix.append(row)
        forcing.append(other)

    pivot_row = 0
    pivot_columns = []
    for column in range(len(variables)):
        selected = next((r for r in range(pivot_row, len(matrix))
                         if matrix[r][column]), None)
        if selected is None:
            continue
        matrix[pivot_row], matrix[selected] = matrix[selected], matrix[pivot_row]
        forcing[pivot_row], forcing[selected] = forcing[selected], forcing[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [value / scale for value in matrix[pivot_row]]
        forcing[pivot_row] = {m: c / scale for m, c in forcing[pivot_row].items()}
        for row in range(len(matrix)):
            if row == pivot_row:
                continue
            scale = matrix[row][column]
            if not scale:
                continue
            matrix[row] = [a - scale * b for a, b in zip(matrix[row], matrix[pivot_row])]
            forcing[row] = add_scaled(forcing[row], forcing[pivot_row], -scale)
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break

    pivot_set = set(pivot_columns)
    free_columns = [j for j in range(len(variables)) if j not in pivot_set]
    substitution = {}
    for r, column in enumerate(pivot_columns):
        expression = {m: -c for m, c in forcing[r].items()}
        for free in free_columns:
            coefficient = -matrix[r][free]
            if coefficient:
                add_term(expression, ((variables[free], 1),), coefficient)
        substitution[variables[column]] = expression

    compatibility = []
    for row in range(pivot_row, len(matrix)):
        if any(matrix[row]):
            raise AssertionError("nonzero row remained below RREF pivots")
        if forcing[row]:
            compatibility.append(forcing[row])

    return {
        "degree": degree,
        "rows": len(matrix),
        "variables": list(variables),
        "rank": len(pivot_columns),
        "pivot_variables": [variables[j] for j in pivot_columns],
        "free_variables": [variables[j] for j in free_columns],
        "substitution": substitution,
        "compatibility": compatibility,
    }


def substitute(poly, substitutions):
    result = {}
    for mon, coefficient in poly.items():
        pivot = [(variable, exponent) for variable, exponent in mon
                 if variable in substitutions]
        if len(pivot) > 1 or (pivot and pivot[0][1] != 1):
            raise ValueError("weighted-degree bound permits at most one linear pivot variable")
        if not pivot:
            add_term(result, mon, coefficient)
            continue
        variable, _ = pivot[0]
        base = tuple((i, e) for i, e in mon if i != variable)
        for replacement_mon, replacement_coefficient in substitutions[variable].items():
            add_term(result, merge_monomials(base, replacement_mon),
                     coefficient * replacement_coefficient)
    return result


def encode_poly(poly):
    return [
        {"monomial": [[i, e] for i, e in mon], "coefficient": str(coefficient)}
        for mon, coefficient in sorted(poly.items())
    ]


def main() -> None:
    source_hash = digest(SOURCE)
    blocks = {d: [] for d in HIGH_BY_DEGREE}
    for tag, saved in streamed_array(SOURCE):
        poly = parse_boundary(saved)
        degree = weighted_degree(poly)
        if degree in blocks:
            blocks[degree].append((tag, poly))

    eliminated = {}
    substitutions = {}
    for degree, variables in HIGH_BY_DEGREE.items():
        rows = blocks[degree]
        result = eliminate_block([poly for _, poly in rows], variables, degree)
        result["tags"] = [tag for tag, _ in rows]
        substitutions.update(result["substitution"])
        eliminated[str(degree)] = result

    counts = Counter()
    support_counts = Counter()
    residual_terms = 0
    zero_count = 0
    quadratic_rows = []
    retained_tags = []
    for tag, saved in streamed_array(SOURCE):
        original = parse_boundary(saved)
        reduced = substitute(original, substitutions)
        if not reduced:
            zero_count += 1
            continue
        retained_tags.append(tag)
        residual_terms += len(reduced)
        support = {i for mon in reduced for i, _ in mon}
        support_counts[len(support)] += 1
        ordinary_degrees = {sum(e for _, e in mon) for mon in reduced}
        counts[f"ordinary_degree_span_{min(ordinary_degrees)}_{max(ordinary_degrees)}"] += 1
        quadratic = {mon: c for mon, c in reduced.items()
                     if sum(e for _, e in mon) == 2}
        if quadratic:
            quadratic_rows.append((tag, quadratic))

    remaining_variables = sorted(set(WEIGHT) - ZERO)
    quadratic_monomials = [(i, j) for a, i in enumerate(remaining_variables)
                           for j in remaining_variables[a:]]
    quadratic_column = {pair: column for column, pair in enumerate(quadratic_monomials)}
    sparse_entries = {}
    for row, (_, poly) in enumerate(quadratic_rows):
        for mon, coefficient in poly.items():
            expanded = []
            for variable, exponent in mon:
                expanded.extend([variable] * exponent)
            if len(expanded) != 2:
                raise AssertionError("quadratic initial form has wrong ordinary degree")
            pair = tuple(sorted(expanded))
            sparse_entries[(row, quadratic_column[pair])] = sp.Rational(
                coefficient.numerator, coefficient.denominator
            )
    quadratic_matrix = sp.SparseMatrix(
        len(quadratic_rows), len(quadratic_monomials), sparse_entries
    )
    if quadratic_rows:
        _, quadratic_pivots = DomainMatrix.from_Matrix(quadratic_matrix).convert_to(QQ).rref()
        quadratic_rank = len(quadratic_pivots)
    else:
        quadratic_rank = 0

    receipt = {
        "source": SOURCE.name,
        "source_sha256": source_hash,
        "boundary": "p0=p1=p2=p3=p4=0 (all weights 1,2,3 vanish)",
        "coefficient_field": "Q",
        "weighted_variable_count_before_elimination": len(WEIGHT) - len(ZERO),
        "weighted_variable_weights": {f"p{i}": WEIGHT[i] for i in sorted(WEIGHT)
                                       if i not in ZERO},
        "linear_blocks": {
            degree: {
                "rows": info["rows"],
                "rank": info["rank"],
                "pivot_variables": [f"p{i}" for i in info["pivot_variables"]],
                "free_variables": [f"p{i}" for i in info["free_variables"]],
                "compatibility_polynomial_count": len(info["compatibility"]),
                "pivot_substitutions": {
                    f"p{variable}": encode_poly(expression)
                    for variable, expression in info["substitution"].items()
                },
                "input_tags": info["tags"],
            }
            for degree, info in eliminated.items()
        },
        "equations_before": 754,
        "equations_after_linear_substitution": 754 - zero_count,
        "equations_becoming_zero": zero_count,
        "residual_term_count": residual_terms,
        "residual_support_size_histogram": dict(sorted(support_counts.items())),
        "residual_ordinary_degree_span_histogram": dict(sorted(counts.items())),
        "quadratic_initial_form_count": len(quadratic_rows),
        "quadratic_initial_form_tags": [tag for tag, _ in quadratic_rows],
        "quadratic_monomial_count": len(quadratic_monomials),
        "quadratic_initial_form_rank_over_Q": quadratic_rank,
        "quadratic_span_codimension": len(quadratic_monomials) - quadratic_rank,
        "scope": "exact linear elimination on one weighted boundary of the detached depth-13 necessary subsystem; not a complete branch exclusion",
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "source_sha256": source_hash,
        "linear_blocks": {d: {"rows": x["rows"], "rank": x["rank"],
                              "pivots": x["pivot_variables"],
                              "free": x["free_variables"],
                              "compatibility": len(x["compatibility"])}
                          for d, x in eliminated.items()},
        "equations_after": receipt["equations_after_linear_substitution"],
        "equations_zero": zero_count,
        "residual_terms": residual_terms,
        "support_histogram": receipt["residual_support_size_histogram"],
        "quadratic_initial_forms": len(quadratic_rows),
        "quadratic_rank_over_Q": quadratic_rank,
        "quadratic_monomials": len(quadratic_monomials),
        "quadratic_span_codimension": len(quadratic_monomials) - quadratic_rank,
        "receipt": str(OUTPUT),
        "scope": receipt["scope"],
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
