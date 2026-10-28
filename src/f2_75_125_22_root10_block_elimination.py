#!/usr/bin/env python3
"""Exact elimination of the first small 2+2 F3 root-10 coefficient block.

This uses six source equations over Q(sqrt(5)) and the source equation that
forces p57=0. It certifies p0=p5=p6=p7=p8=0 on this necessary root-10 locus.
The calculation does not exclude the full 2+2 F2 branch.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix

from f2_75_125_22_root10_support_profile import constraint_array


TARGET_TAGS = (
    "terminal-alpha-P-distance3:u^2",
    "terminal-alpha-P-distance3:u^3",
    "terminal-alpha-P-distance3:u^4",
    "terminal-alpha-Q-distance3:u^6",
    "terminal-alpha-Q-distance3:u^7",
    "terminal-alpha-Q-distance3:u^8",
)
SECOND_BLOCK_TAGS = (
    "terminal-alpha-P-distance5:u^0",
    "terminal-alpha-P-distance5:u^1",
    "terminal-alpha-P-distance5:u^2",
    "terminal-alpha-P-distance5:u^3",
    "terminal-alpha-Q-distance5:u^4",
    "terminal-alpha-Q-distance5:u^5",
    "terminal-alpha-Q-distance5:u^6",
    "terminal-alpha-Q-distance5:u^7",
)
THIRD_BLOCK_TAGS = (
    "terminal-alpha-P-distance7:u^0",
    "terminal-alpha-P-distance7:u^1",
    "terminal-alpha-P-distance7:u^2",
    "terminal-alpha-P-distance7:u^3",
    "terminal-alpha-Q-distance7:u^4",
    "terminal-alpha-Q-distance7:u^5",
    "terminal-alpha-Q-distance7:u^6",
    "terminal-alpha-Q-distance7:u^7",
)
P57_WITNESS_TAG = "10:1:t^38"
N_PARAMETERS = 69


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def add_pair(a: tuple[Fraction, Fraction], b: tuple[Fraction, Fraction]):
    return a[0] + b[0], a[1] + b[1]


def parse_source_polynomial(encoded):
    poly: dict[tuple[int, ...], tuple[Fraction, Fraction]] = {}
    for sparse, values in encoded:
        exp_map = {int(i): int(e) for i, e in sparse if e}
        if 0 in exp_map:
            raise AssertionError("root-10 constraints must not retain t")
        if any(i > N_PARAMETERS for i in exp_map):
            raise AssertionError("root-10 constraint contains a non-parameter generator")
        exponent = tuple(exp_map.get(i + 1, 0) for i in range(N_PARAMETERS))
        if len(values) != 2:
            raise AssertionError("expected Q(sqrt(5)) coefficient pair")
        coefficient = (Fraction(values[0]), Fraction(values[1]))
        poly[exponent] = add_pair(
            poly.get(exponent, (Fraction(0), Fraction(0))), coefficient
        )
    return {e: c for e, c in poly.items() if c != (0, 0)}


def field_expr(coefficient, sqrt5):
    rational, radical = coefficient
    return (sp.Rational(rational.numerator, rational.denominator)
            + sp.Rational(radical.numerator, radical.denominator) * sqrt5)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("state", type=Path)
    parser.add_argument("receipt", type=Path)
    args = parser.parse_args()
    source_hash = sha256(args.state)

    selected = {}
    second_selected = {}
    third_selected = {}
    p57_poly = None
    for tag, encoded in constraint_array(args.state):
        if (tag not in TARGET_TAGS and tag not in SECOND_BLOCK_TAGS
                and tag not in THIRD_BLOCK_TAGS and tag != P57_WITNESS_TAG):
            continue
        poly = parse_source_polynomial(encoded)
        if tag in TARGET_TAGS:
            if tag in selected:
                raise AssertionError(f"duplicate source tag {tag}")
            selected[tag] = poly
        elif tag in SECOND_BLOCK_TAGS:
            if tag in second_selected:
                raise AssertionError(f"duplicate source tag {tag}")
            second_selected[tag] = poly
        elif tag in THIRD_BLOCK_TAGS:
            if tag in third_selected:
                raise AssertionError(f"duplicate source tag {tag}")
            third_selected[tag] = poly
        else:
            p57_poly = poly
    if set(selected) != set(TARGET_TAGS):
        raise ValueError("not all six exact root-distance-3 source rows were found")
    if p57_poly is None or len(p57_poly) != 1:
        raise ValueError("p57 witness row is not a single nonzero monomial")
    p57_exponent, p57_coefficient = next(iter(p57_poly.items()))
    if p57_exponent[57] != 1 or sum(p57_exponent) != 1:
        raise ValueError("p57 source witness does not have the claimed form")

    sqrt5 = sp.sqrt(5)
    p0, p3, p5, p6, p7, p8 = sp.symbols("p0 p3 p5 p6 p7 p8")
    variable = {0: p0, 3: p3, 5: p5, 6: p6, 7: p7, 8: p8}
    equations = []
    for tag in TARGET_TAGS:
        expression = sp.S.Zero
        for exponent, coefficient in selected[tag].items():
            monomial = sp.prod(variable[i] ** exponent[i]
                               for i in variable if exponent[i])
            if any(exponent[i] for i in range(N_PARAMETERS) if i not in variable):
                raise AssertionError(f"unexpected variable in selected row {tag}")
            expression += field_expr(coefficient, sqrt5) * monomial
        equations.append(sp.expand(expression))

    unknowns = (p5, p6, p7, p8)
    matrix_expr, forcing_expr = sp.linear_eq_to_matrix(equations, unknowns)
    field = QQ.algebraic_field(sqrt5)
    matrix_rows = [
        [field.convert(entry) for entry in matrix_expr.row(i)]
        for i in range(matrix_expr.rows)
    ]
    matrix = DomainMatrix(matrix_rows, matrix_expr.shape, field)
    _, pivots = matrix.rref()
    rank = len(pivots)
    if rank != 4:
        raise AssertionError(f"the coefficient block has rank {rank}, expected 4")

    transpose_rows = [
        [matrix_rows[i][j] for i in range(len(matrix_rows))]
        for j in range(4)
    ]
    _, independent_row_indices = DomainMatrix(
        transpose_rows, (4, len(matrix_rows)), field
    ).rref()
    if len(independent_row_indices) != 4:
        raise AssertionError("failed to extract four independent source rows")
    minor_expr = matrix_expr.extract(independent_row_indices, [0, 1, 2, 3])
    minor_det = sp.simplify(sp.expand(minor_expr.det()))
    if field.convert(minor_det) == field.zero:
        raise AssertionError("reported 4x4 minor determinant is zero")

    left_kernel = sp.Matrix(matrix_expr.T).nullspace()
    if len(left_kernel) != 2:
        raise AssertionError(f"expected 2 compatibility rows, got {len(left_kernel)}")
    compatibility = []
    compatibility_coefficients = []
    for vector in left_kernel:
        residual = sp.factor(sum(vector[i] * forcing_expr[i]
                                 for i in range(len(equations))))
        polynomial = sp.Poly(residual, p0, p3, extension=sqrt5)
        allowed = {(3, 0), (1, 1)}
        if set(polynomial.monoms()) - allowed:
            raise AssertionError("compatibility row is not in span(p0^3,p0*p3)")
        coefficient_cubic = sp.simplify(polynomial.coeff_monomial(p0 ** 3))
        coefficient_mixed = sp.simplify(polynomial.coeff_monomial(p0 * p3))
        if coefficient_cubic == 0 or coefficient_mixed == 0:
            raise AssertionError("compatibility row lost one of the two base monomials")
        compatibility.append(residual)
        compatibility_coefficients.append((coefficient_cubic, coefficient_mixed))

    compat_det = sp.simplify(
        compatibility_coefficients[0][0] * compatibility_coefficients[1][1]
        - compatibility_coefficients[1][0] * compatibility_coefficients[0][1]
    )
    if field.convert(compat_det) == field.zero:
        raise AssertionError("compatibility equations do not separate p0^3 and p0*p3")

    zero_forcing = [sp.simplify(value.subs(p0, 0)) for value in forcing_expr]
    if any(value != 0 for value in zero_forcing):
        raise AssertionError("setting p0=0 does not annihilate the coefficient block forcing")

    if set(second_selected) != set(SECOND_BLOCK_TAGS):
        raise ValueError("not all eight distance-5 source rows were found")
    second_unknown_indices = (18, 19, 20, 21, 22, 23)
    second_rows_expr = []
    for tag in SECOND_BLOCK_TAGS:
        poly = second_selected[tag]
        row = [sp.S.Zero] * len(second_unknown_indices)
        for exponent, coefficient in poly.items():
            if any(exponent[i] for i in (0, 5, 6, 7, 8, 57)):
                continue
            degree = sum(exponent)
            if degree == 0:
                raise AssertionError(f"distance-5 row {tag} has a nonzero constant")
            if degree != 1:
                raise AssertionError(f"distance-5 row {tag} is not linear after first-block reduction")
            support = [i for i, power in enumerate(exponent) if power]
            if len(support) != 1 or support[0] not in second_unknown_indices:
                raise AssertionError(f"unexpected distance-5 variable in {tag}: {support}")
            column = second_unknown_indices.index(support[0])
            row[column] += field.convert(field_expr(coefficient, sqrt5))
        second_rows_expr.append(row)
    second_matrix = DomainMatrix(
        second_rows_expr,
        (len(second_rows_expr), len(second_unknown_indices)),
        field,
    )
    _, second_pivots = second_matrix.rref()
    second_rank = len(second_pivots)
    if second_rank != len(second_unknown_indices):
        raise AssertionError(f"distance-5 block has rank {second_rank}, expected 6")
    second_transpose_rows = [
        [second_rows_expr[i][j] for i in range(len(second_rows_expr))]
        for j in range(len(second_unknown_indices))
    ]
    _, second_independent_rows = DomainMatrix(
        second_transpose_rows,
        (len(second_unknown_indices), len(second_rows_expr)),
        field,
    ).rref()
    second_matrix_expr = sp.Matrix([
        [sp.simplify(field.to_sympy(entry)) for entry in row]
        for row in second_rows_expr
    ])
    second_minor = second_matrix_expr.extract(
        second_independent_rows, list(range(len(second_unknown_indices)))
    )
    second_minor_det = sp.simplify(sp.expand(second_minor.det()))
    if field.convert(second_minor_det) == field.zero:
        raise AssertionError("reported 6x6 distance-5 minor determinant is zero")

    if set(third_selected) != set(THIRD_BLOCK_TAGS):
        raise ValueError("not all eight distance-7 source rows were found")
    prior_zero_indices = {0, 5, 6, 7, 8, 18, 19, 20, 21, 22, 23, 57}
    third_unknown_indices = (30, 31, 32, 33, 34, 35)
    third_rows = []
    for tag in THIRD_BLOCK_TAGS:
        poly = third_selected[tag]
        row = [sp.S.Zero] * len(third_unknown_indices)
        for exponent, coefficient in poly.items():
            if any(exponent[i] for i in prior_zero_indices):
                continue
            degree = sum(exponent)
            if degree == 0:
                raise AssertionError(f"distance-7 row {tag} has a nonzero constant")
            if degree != 1:
                raise AssertionError(f"distance-7 row {tag} is not linear after prior reductions")
            support = [i for i, power in enumerate(exponent) if power]
            if len(support) != 1 or support[0] not in third_unknown_indices:
                raise AssertionError(f"unexpected distance-7 variable in {tag}: {support}")
            column = third_unknown_indices.index(support[0])
            row[column] += field.convert(field_expr(coefficient, sqrt5))
        third_rows.append(row)
    third_matrix = DomainMatrix(
        third_rows, (len(third_rows), len(third_unknown_indices)), field
    )
    _, third_pivots = third_matrix.rref()
    third_rank = len(third_pivots)
    if third_rank != len(third_unknown_indices):
        raise AssertionError(f"distance-7 block has rank {third_rank}, expected 6")
    third_transpose_rows = [
        [third_rows[i][j] for i in range(len(third_rows))]
        for j in range(len(third_unknown_indices))
    ]
    _, third_independent_rows = DomainMatrix(
        third_transpose_rows,
        (len(third_unknown_indices), len(third_rows)),
        field,
    ).rref()
    third_matrix_expr = sp.Matrix([
        [sp.simplify(field.to_sympy(entry)) for entry in row]
        for row in third_rows
    ])
    third_minor = third_matrix_expr.extract(
        third_independent_rows, list(range(len(third_unknown_indices)))
    )
    third_minor_det = sp.simplify(sp.expand(third_minor.det()))
    if field.convert(third_minor_det) == field.zero:
        raise AssertionError("reported 6x6 distance-7 minor determinant is zero")

    result = {
        "source": str(args.state),
        "source_sha256": source_hash,
        "scope": "exact elimination inside one 2+2 F3 dual-root distance-10 necessary system; not a whole-branch exclusion",
        "coefficient_field": "Q(sqrt(5))",
        "p57_witness": {
            "tag": P57_WITNESS_TAG,
            "equation": "c*p57=0",
            "coefficient": [str(p57_coefficient[0]), str(p57_coefficient[1])],
            "interpretation": "c=a+b*sqrt(5) is nonzero, hence p57=0",
        },
        "block_tags": list(TARGET_TAGS),
        "block_shape": list(matrix_expr.shape),
        "block_coefficient_rank": rank,
        "independent_row_indices_zero_based": list(independent_row_indices),
        "independent_4x4_minor_determinant": str(minor_det),
        "left_nullspace_basis": [[str(sp.simplify(entry)) for entry in vector]
                                  for vector in left_kernel],
        "compatibility_equations": [str(expression) for expression in compatibility],
        "compatibility_coefficient_matrix_rows_p0cubed_p0p3": [
            [str(cubic), str(mixed)] for cubic, mixed in compatibility_coefficients
        ],
        "compatibility_determinant": str(compat_det),
        "forced_zero_parameters": ["p57", "p0", "p5", "p6", "p7", "p8"],
        "second_block": {
            "tags": list(SECOND_BLOCK_TAGS),
            "matrix_shape": [len(second_rows_expr), len(second_unknown_indices)],
            "variables": [f"p{i}" for i in second_unknown_indices],
            "rank": second_rank,
            "independent_row_indices_zero_based": list(second_independent_rows),
            "independent_6x6_minor_determinant": str(second_minor_det),
            "consequence": "The block is homogeneous, square full rank after choosing six rows, so p18 through p23 all vanish.",
        },
        "forced_zero_parameters_after_second_block": [
            "p57", "p0", "p5", "p6", "p7", "p8",
            "p18", "p19", "p20", "p21", "p22", "p23",
        ],
        "third_block": {
            "tags": list(THIRD_BLOCK_TAGS),
            "matrix_shape": [len(third_rows), len(third_unknown_indices)],
            "variables": [f"p{i}" for i in third_unknown_indices],
            "rank": third_rank,
            "independent_row_indices_zero_based": list(third_independent_rows),
            "independent_6x6_minor_determinant": str(third_minor_det),
            "consequence": "After the first two block substitutions, these homogeneous rows force p30 through p35 to vanish.",
        },
        "forced_zero_parameters_after_third_block": [
            "p57", "p0", "p5", "p6", "p7", "p8",
            "p18", "p19", "p20", "p21", "p22", "p23",
            "p30", "p31", "p32", "p33", "p34", "p35",
        ],
        "derivation": [
            "The source row c*p57=0 has c nonzero in Q(sqrt(5)), so p57=0.",
            "The six tagged equations are linear in p5,p6,p7,p8 with coefficient rank 4.",
            "The two left-kernel compatibility equations are linear combinations of p0^3 and p0*p3 with nonzero 2x2 determinant.",
            "Thus p0^3=0, hence p0=0; the block forcing then vanishes and rank 4 forces p5=p6=p7=p8=0.",
            "After that substitution, the eight distance-5 source rows are homogeneous linear equations in p18,...,p23; a nonzero 6x6 minor forces all six to zero.",
            "After the distance-5 zeros, eight homogeneous distance-7 rows have a nonzero 6x6 minor in p30,...,p35, forcing those six to zero as well.",
        ],
        "warning": "The remaining root-10 constraints and other 2+2 sheets remain. This exact local reduction is not an exclusion of the 2+2 branch.",
    }
    args.receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
