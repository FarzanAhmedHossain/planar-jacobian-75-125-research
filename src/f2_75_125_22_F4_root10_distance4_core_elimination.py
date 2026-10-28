#!/usr/bin/env python3
"""Certify the F4 distance-4 core block after preceding zero reductions."""

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
    "terminal-alpha-P-distance4:u^1",
    "terminal-alpha-P-distance4:u^2",
    "terminal-alpha-P-distance4:u^3",
    "terminal-alpha-P-distance4:u^4",
    "terminal-alpha-Q-distance4:u^5",
    "terminal-alpha-Q-distance4:u^6",
    "terminal-alpha-Q-distance4:u^7",
    "terminal-alpha-Q-distance4:u^8",
)
UNKNOWN = (11, 12, 13, 14, 15, 16)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("state", type=Path)
    parser.add_argument("prior_zero_certificate", type=Path)
    parser.add_argument("receipt", type=Path)
    args = parser.parse_args()
    prior = json.loads(args.prior_zero_certificate.read_text(encoding="utf-8"))
    source_hash = sha256(args.state)
    if prior["source_sha256"] != source_hash:
        raise ValueError("prior zero certificate refers to another source state")
    prior_zeros = prior["forced_zero_parameters_after_outer_linear_block"]
    prior_zero_indices = {int(name[1:]) for name in prior_zeros}

    sqrt5 = sp.sqrt(5)
    field = QQ.algebraic_field(sqrt5)
    selected = {}
    for tag, encoded in constraint_array(args.state):
        if tag not in TAGS:
            continue
        row = [field.zero] * len(UNKNOWN)
        forcing = field.zero  # coefficient of p3^2
        for sparse, values in encoded:
            exp_map = {int(i): int(e) for i, e in sparse if e}
            if 0 in exp_map:
                raise AssertionError("root-10 row unexpectedly retains t")
            exponent = tuple(exp_map.get(i + 1, 0) for i in range(69))
            if any(exponent[i] for i in prior_zero_indices):
                continue
            if len(values) == 1:
                rational, radical = Fraction(values[0]), Fraction(0)
            elif len(values) == 2:
                rational, radical = Fraction(values[0]), Fraction(values[1])
            else:
                raise AssertionError("unexpected Q(sqrt(5)) coefficient encoding")
            coefficient = sp.Rational(rational.numerator, rational.denominator)
            coefficient += sp.Rational(radical.numerator, radical.denominator) * sqrt5
            value = field.convert(coefficient)
            support = [i for i, power in enumerate(exponent) if power]
            if support == [3] and exponent[3] == 2:
                forcing += value
            elif len(support) == 1 and support[0] in UNKNOWN and exponent[support[0]] == 1:
                row[UNKNOWN.index(support[0])] += value
            else:
                raise AssertionError(f"unexpected monomial in {tag}: {support}, {exponent}")
        selected[tag] = (row, forcing)

    if set(selected) != set(TAGS):
        raise ValueError("not all eight distance-4 source rows were found")
    rows = [selected[tag][0] for tag in TAGS]
    forcing = [selected[tag][1] for tag in TAGS]
    matrix = DomainMatrix(rows, (len(rows), len(UNKNOWN)), field)
    _, pivots = matrix.rref()
    rank = len(pivots)
    if rank != len(UNKNOWN):
        raise AssertionError(f"distance-4 core coefficient rank is {rank}, expected 6")
    transpose_rows = [[rows[i][j] for i in range(len(rows))]
                      for j in range(len(UNKNOWN))]
    _, independent_rows = DomainMatrix(
        transpose_rows, (len(UNKNOWN), len(rows)), field
    ).rref()
    minor = DomainMatrix([rows[i] for i in independent_rows],
                         (len(UNKNOWN), len(UNKNOWN)), field)
    determinant = minor.det()
    if determinant == field.zero:
        raise AssertionError("selected exact 6x6 determinant is zero")

    left_nullspace = matrix.transpose().nullspace().to_list()
    compatibility = [
        sum((vector[i] * forcing[i] for i in range(len(TAGS))), field.zero)
        for vector in left_nullspace
    ]
    nonzero_compatibility = [value for value in compatibility if value != field.zero]
    if not nonzero_compatibility:
        raise AssertionError("left-kernel compatibility does not force p3^2=0")
    forced = list(prior_zeros) + ["p3", *[f"p{i}" for i in UNKNOWN]]

    result = {
        "source": str(args.state),
        "source_sha256": source_hash,
        "prior_zero_certificate": str(args.prior_zero_certificate),
        "prior_zero_certificate_sha256": sha256(args.prior_zero_certificate),
        "scope": "exact distance-4 core elimination inside one 2+2 F4 dual-root distance-10 necessary system",
        "coefficient_field": "Q(sqrt(5))",
        "source_tags": list(TAGS),
        "variables": [f"p{i}" for i in UNKNOWN],
        "forcing_monomial": "p3^2",
        "matrix_shape": [len(rows), len(UNKNOWN)],
        "coefficient_rank": rank,
        "independent_source_row_indices_zero_based": list(independent_rows),
        "independent_6x6_minor_determinant": str(field.to_sympy(determinant)),
        "left_nullity": len(left_nullspace),
        "compatibility_coefficients_of_p3_squared": [str(field.to_sympy(value)) for value in compatibility],
        "nonzero_compatibility_witness": str(field.to_sympy(nonzero_compatibility[0])),
        "consequence": "A left-kernel combination gives c*p3^2=0 with c nonzero, so p3=0; the full-rank homogeneous block then forces p11=...=p16=0.",
        "forced_zero_parameters_after_core_block": forced,
        "warning": "This remains a local reduction of one F4 root-10 necessary sheet, not an exclusion of the F4 component or 2+2 branch.",
    }
    args.receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
