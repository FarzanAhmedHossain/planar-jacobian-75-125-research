#!/usr/bin/env python3
"""Replay the full-rank degree-2 core block in the F3 root-10 system."""

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
    "terminal-alpha-P-distance4:u^0",
    "terminal-alpha-P-distance4:u^1",
    "terminal-alpha-P-distance4:u^2",
    "terminal-alpha-P-distance4:u^3",
    "terminal-alpha-P-distance4:u^4",
    "terminal-alpha-Q-distance4:u^4",
    "terminal-alpha-Q-distance4:u^5",
    "terminal-alpha-Q-distance4:u^6",
)
MONOMIALS = ((3, 2), (10, 1), (11, 1), (12, 1), (13, 1), (14, 1), (15, 1), (16, 1))
ZERO_INDICES = (0, 5, 6, 7, 8, 18, 19, 20, 21, 22, 23, 30, 31, 32, 33, 34, 35, 57)


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
    if prior["source_sha256"] != sha256(args.state):
        raise ValueError("prior zero certificate refers to another source state")
    expected_prior = {f"p{i}" for i in ZERO_INDICES}
    if not expected_prior <= set(prior["forced_zero_parameters_after_third_block"]):
        raise ValueError("prior certificate does not prove every required zero substitution")

    sqrt5 = sp.sqrt(5)
    field = QQ.algebraic_field(sqrt5)
    selected = {}
    for tag, encoded in constraint_array(args.state):
        if tag not in TAGS:
            continue
        row = [field.zero] * len(MONOMIALS)
        for sparse, values in encoded:
            exp_map = {int(i): int(e) for i, e in sparse if e}
            if 0 in exp_map:
                raise AssertionError("root-10 row unexpectedly contains t")
            exponent = tuple(exp_map.get(i + 1, 0) for i in range(N_PARAMETERS))
            if any(exponent[i] for i in ZERO_INDICES):
                continue
            if len(values) == 1:
                rational, radical = Fraction(values[0]), Fraction(0)
            elif len(values) == 2:
                rational, radical = Fraction(values[0]), Fraction(values[1])
            else:
                raise AssertionError("unexpected Q(sqrt(5)) coefficient encoding")
            support = [i for i, power in enumerate(exponent) if power]
            if len(support) != 1:
                raise AssertionError(f"unexpected monomial support in {tag}: {support}")
            parameter = support[0]
            power = exponent[parameter]
            if (parameter, power) not in MONOMIALS:
                raise AssertionError(f"unexpected core monomial p{parameter}^{power} in {tag}")
            column = MONOMIALS.index((parameter, power))
            coefficient = sp.Rational(rational.numerator, rational.denominator)
            coefficient += sp.Rational(radical.numerator, radical.denominator) * sqrt5
            row[column] += field.convert(coefficient)
        if tag in selected:
            raise AssertionError(f"duplicate source row {tag}")
        selected[tag] = row

    if set(selected) != set(TAGS):
        raise ValueError("not all eight selected source rows were found")
    matrix = DomainMatrix([selected[tag] for tag in TAGS], (8, 8), field)
    _, pivots = matrix.rref()
    rank = len(pivots)
    determinant = matrix.det()
    if rank != 8 or determinant == field.zero:
        raise AssertionError(f"the core degree-2 block is not full rank: rank={rank}")
    result = {
        "source": str(args.state),
        "source_sha256": sha256(args.state),
        "prior_zero_certificate": str(args.prior_zero_certificate),
        "prior_zero_certificate_sha256": sha256(args.prior_zero_certificate),
        "scope": "exact degree-2 core elimination inside one 2+2 F3 root-10 necessary system",
        "coefficient_field": "Q(sqrt(5))",
        "source_tags": list(TAGS),
        "unknown_monomials": ["p3^2", "p10", "p11", "p12", "p13", "p14", "p15", "p16"],
        "matrix_shape": [8, 8],
        "rank": rank,
        "determinant": str(field.to_sympy(determinant)),
        "consequence": "The only solution for the eight monomial values is zero. Hence p3^2=0 and p10=...=p16=0; over the coefficient field's point values this forces p3=p10=...=p16=0.",
        "forced_zero_parameters_this_step": ["p3", "p10", "p11", "p12", "p13", "p14", "p15", "p16"],
        "forced_zero_parameters_after_core_degree2": list(prior["forced_zero_parameters_after_third_block"]) + ["p3", "p10", "p11", "p12", "p13", "p14", "p15", "p16"],
        "warning": "This eliminates the core variables only on the selected F3 root-10 necessary sheet; the remaining p17,p24,...,p29 equations and other sheets remain.",
    }
    args.receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
