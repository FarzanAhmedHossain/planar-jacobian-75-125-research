#!/usr/bin/env python3
"""Certify the final six-variable linear block after core elimination."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix

from f2_75_125_22_root10_block_elimination import N_PARAMETERS
from f2_75_125_22_root10_support_profile import constraint_array


TAGS = (
    "terminal-alpha-P-distance6:u^0",
    "terminal-alpha-P-distance6:u^1",
    "terminal-alpha-P-distance6:u^2",
    "terminal-alpha-P-distance6:u^3",
    "terminal-alpha-Q-distance6:u^4",
    "terminal-alpha-Q-distance6:u^5",
    "terminal-alpha-Q-distance6:u^6",
    "terminal-alpha-Q-distance6:u^7",
)
UNKNOWN_INDICES = (24, 25, 26, 27, 28, 29)


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
    prior_zeros = prior["forced_zero_parameters_after_core_degree2"]
    prior_zero_indices = {int(name[1:]) for name in prior_zeros}

    sqrt5 = sp.sqrt(5)
    field = QQ.algebraic_field(sqrt5)
    selected = {}
    for tag, encoded in constraint_array(args.state):
        if tag not in TAGS:
            continue
        row = [field.zero] * len(UNKNOWN_INDICES)
        for sparse, values in encoded:
            exp_map = {int(i): int(e) for i, e in sparse if e}
            if 0 in exp_map:
                raise AssertionError("root-10 row unexpectedly contains t")
            exponent = tuple(exp_map.get(i + 1, 0) for i in range(N_PARAMETERS))
            if any(exponent[i] for i in prior_zero_indices):
                continue
            if len(values) == 1:
                rational, radical = Fraction(values[0]), Fraction(0)
            elif len(values) == 2:
                rational, radical = Fraction(values[0]), Fraction(values[1])
            else:
                raise AssertionError("unexpected Q(sqrt(5)) coefficient encoding")
            degree = sum(exponent)
            if degree != 1:
                raise AssertionError(f"distance-6 row {tag} is not homogeneous linear")
            support = [i for i, power in enumerate(exponent) if power]
            if len(support) != 1 or support[0] not in UNKNOWN_INDICES:
                raise AssertionError(f"unexpected variable in {tag}: {support}")
            column = UNKNOWN_INDICES.index(support[0])
            coefficient = sp.Rational(rational.numerator, rational.denominator)
            coefficient += sp.Rational(radical.numerator, radical.denominator) * sqrt5
            row[column] += field.convert(coefficient)
        selected[tag] = row

    if set(selected) != set(TAGS):
        raise ValueError("not all eight distance-6 source rows were found")
    rows = [selected[tag] for tag in TAGS]
    matrix = DomainMatrix(rows, (len(rows), len(UNKNOWN_INDICES)), field)
    _, pivots = matrix.rref()
    rank = len(pivots)
    transpose_rows = [[rows[i][j] for i in range(len(rows))]
                      for j in range(len(UNKNOWN_INDICES))]
    _, independent_rows = DomainMatrix(
        transpose_rows, (len(UNKNOWN_INDICES), len(rows)), field
    ).rref()
    if rank != len(UNKNOWN_INDICES) or len(independent_rows) != len(UNKNOWN_INDICES):
        raise AssertionError(f"distance-6 block rank is {rank}, expected 6")
    minor = DomainMatrix(
        [rows[i] for i in independent_rows],
        (len(UNKNOWN_INDICES), len(UNKNOWN_INDICES)),
        field,
    )
    determinant = minor.det()
    if determinant == field.zero:
        raise AssertionError("selected exact 6x6 minor determinant is zero")

    result = {
        "source": str(args.state),
        "source_sha256": source_hash,
        "prior_zero_certificate": str(args.prior_zero_certificate),
        "prior_zero_certificate_sha256": sha256(args.prior_zero_certificate),
        "scope": "exact homogeneous outer linear elimination inside one 2+2 F3 root-10 necessary system",
        "coefficient_field": "Q(sqrt(5))",
        "source_tags": list(TAGS),
        "matrix_shape": [len(rows), len(UNKNOWN_INDICES)],
        "variables": [f"p{i}" for i in UNKNOWN_INDICES],
        "rank": rank,
        "independent_source_row_indices_zero_based": list(independent_rows),
        "independent_6x6_minor_determinant": str(field.to_sympy(determinant)),
        "consequence": "The eight source equations are homogeneous linear in p24,...,p29 and have coefficient rank 6, so p24=...=p29=0.",
        "forced_zero_parameters_after_outer_linear_block": list(prior_zeros) + [f"p{i}" for i in UNKNOWN_INDICES],
        "warning": "This solves one necessary root-10 sheet locally; it does not exclude the full 2+2 branch or establish a contradiction with a required nonzero parameter.",
    }
    args.receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
