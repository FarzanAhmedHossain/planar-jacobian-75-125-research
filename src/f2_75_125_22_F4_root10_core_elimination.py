#!/usr/bin/env python3
"""Exact small-support elimination in the F4 dual-root distance-10 sheet."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix

from f2_75_125_22_root10_support_profile import constraint_array


TAGS = (
    "terminal-alpha-P-distance3:u^1",
    "terminal-alpha-P-distance3:u^2",
    "terminal-alpha-P-distance3:u^3",
    "terminal-alpha-P-distance3:u^4",
    "terminal-alpha-Q-distance3:u^5",
    "terminal-alpha-Q-distance3:u^6",
    "terminal-alpha-Q-distance3:u^7",
    "terminal-alpha-Q-distance3:u^8",
)
WITNESS_TAG = "10:1:t^38"
UNKNOWN = (4, 5, 6, 7, 8)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def coeff_pair(values):
    if len(values) == 1:
        return Fraction(values[0]), Fraction(0)
    if len(values) == 2:
        return Fraction(values[0]), Fraction(values[1])
    raise ValueError(f"unexpected Q(sqrt(5)) coefficient encoding: {values!r}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("state", type=Path)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    sqrt5 = sp.sqrt(5)
    field = QQ.algebraic_field(sqrt5)
    selected = {}
    witness = None
    for tag, encoded in constraint_array(args.state):
        if tag not in TAGS and tag != WITNESS_TAG:
            continue
        polynomial = {}
        for sparse, values in encoded:
            exp_map = {int(i): int(e) for i, e in sparse if e}
            if 0 in exp_map:
                raise AssertionError("root-10 F4 constraint unexpectedly retains t")
            exponent = tuple(exp_map.get(i + 1, 0) for i in range(69))
            rational, radical = coeff_pair(values)
            old = polynomial.get(exponent, (Fraction(0), Fraction(0)))
            polynomial[exponent] = (old[0] + rational, old[1] + radical)
        polynomial = {e: c for e, c in polynomial.items() if c != (0, 0)}
        if tag == WITNESS_TAG:
            witness = polynomial
        else:
            selected[tag] = polynomial
    if witness is None or len(witness) != 1:
        raise ValueError("p57 one-term source witness is missing")
    (witness_exp, witness_coeff), = witness.items()
    if witness_exp[57] != 1 or sum(witness_exp) != 1:
        raise AssertionError("the one-term witness is not c*p57=0")
    if set(selected) != set(TAGS):
        raise ValueError("not all eight distance-3 source rows were found")

    coefficient_rows = []
    forcing_rows = []
    for tag in TAGS:
        coefficients = [field.zero] * len(UNKNOWN)
        forcing = [field.zero, field.zero]  # p0*p3, p0^3
        for exponent, (rational, radical) in selected[tag].items():
            if exponent[57]:
                continue  # justified by the exact witness c*p57=0
            support = [i for i, power in enumerate(exponent) if power]
            coefficient_expr = sp.Rational(rational.numerator, rational.denominator)
            coefficient_expr += sp.Rational(radical.numerator, radical.denominator) * sqrt5
            value = field.convert(coefficient_expr)
            if support == [0, 3] and exponent[0] == exponent[3] == 1:
                forcing[0] += value
            elif support == [0] and exponent[0] == 3:
                forcing[1] += value
            elif len(support) == 1 and support[0] in UNKNOWN and exponent[support[0]] == 1:
                coefficients[UNKNOWN.index(support[0])] += value
            else:
                raise AssertionError(f"unexpected monomial in {tag}: {support}, {exponent}")
        coefficient_rows.append(coefficients)
        forcing_rows.append(forcing)

    matrix = DomainMatrix(coefficient_rows, (len(TAGS), len(UNKNOWN)), field)
    _, pivots = matrix.rref()
    rank = len(pivots)
    if rank != len(UNKNOWN):
        raise AssertionError(f"degree-3 coefficient block has rank {rank}, expected 5")
    transpose_rows = [[coefficient_rows[i][j] for i in range(len(TAGS))]
                      for j in range(len(UNKNOWN))]
    _, independent_rows = DomainMatrix(
        transpose_rows, (len(UNKNOWN), len(TAGS)), field
    ).rref()
    minor = DomainMatrix(
        [coefficient_rows[i] for i in independent_rows],
        (len(UNKNOWN), len(UNKNOWN)), field,
    )
    minor_det = minor.det()
    if minor_det == field.zero:
        raise AssertionError("selected 5x5 minor is singular")

    left_nullspace = matrix.transpose().nullspace()
    compatibility_rows = []
    for vector in left_nullspace.to_list():
        compatibility_rows.append([
            sum((vector[i] * forcing_rows[i][j] for i in range(len(TAGS))), field.zero)
            for j in range(2)
        ])
    compatibility = DomainMatrix(
        compatibility_rows, (len(compatibility_rows), 2), field
    )
    _, compatibility_pivots = compatibility.rref()
    compatibility_rank = len(compatibility_pivots)
    if compatibility_rank == 2:
        compat_transpose = [
            [compatibility_rows[i][j] for i in range(len(compatibility_rows))]
            for j in range(2)
        ]
        _, compat_independent = DomainMatrix(
            compat_transpose, (2, len(compatibility_rows)), field
        ).rref()
        compat_minor = DomainMatrix(
            [compatibility_rows[i] for i in compat_independent], (2, 2), field
        )
        compat_det = compat_minor.det()
        if compat_det == field.zero:
            raise AssertionError("selected compatibility determinant vanished")
        forced = ["p57", "p0", *[f"p{i}" for i in UNKNOWN]]
        consequence = "The two independent compatibility equations force p0*p3=p0^3=0, hence p0=0; rank 5 then forces p4=...=p8=0."
    else:
        compat_independent = []
        compat_det = field.zero
        forced = ["p57"]
        consequence = "The compatibility rank is below 2; retain the exact compatibility equations and analyze them with later rows."

    result = {
        "source": str(args.state),
        "source_sha256": sha256(args.state),
        "scope": "exact small-support elimination inside one 2+2 F4 dual-root distance-10 necessary system",
        "coefficient_field": "Q(sqrt(5))",
        "p57_witness": {
            "tag": WITNESS_TAG,
            "equation": "c*p57=0",
            "coefficient_pair": [str(witness_coeff[0]), str(witness_coeff[1])],
        },
        "source_tags": list(TAGS),
        "variables": [f"p{i}" for i in UNKNOWN],
        "forcing_monomials": ["p0*p3", "p0^3"],
        "coefficient_matrix_shape": [len(TAGS), len(UNKNOWN)],
        "coefficient_rank": rank,
        "independent_source_row_indices_zero_based": list(independent_rows),
        "independent_5x5_minor_determinant": str(field.to_sympy(minor_det)),
        "left_nullity": len(compatibility_rows),
        "compatibility_matrix_rows_p0p3_p0cubed": [
            [str(field.to_sympy(entry)) for entry in row]
            for row in compatibility_rows
        ],
        "compatibility_rank": compatibility_rank,
        "compatibility_independent_row_indices_zero_based": list(compat_independent),
        "compatibility_minor_determinant": str(field.to_sympy(compat_det)),
        "consequence": consequence,
        "forced_zero_parameters_after_core_block": forced,
        "warning": "This local reduction concerns one F4 dual-root root-10 necessary sheet and is not a 2+2 branch exclusion.",
    }
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
