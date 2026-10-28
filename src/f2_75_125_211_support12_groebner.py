#!/usr/bin/env python3
"""Exact Gröbner study of the positive-graded support-12 subideal.

The input is the detached depth-13 necessary-equation array. The script
selects every polynomial on the unique 12-parameter support, verifies the
positive grading, computes an exact characteristic-zero grevlex F5B basis,
and checks whether each support variable is nilpotent modulo the ideal.
Even a result forcing these 12 variables to zero would be only a necessary
branch reduction; the other parameters and full Keller lift remain.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from collections import Counter, defaultdict
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ

from f2_75_125_211_degree_structure import streamed_array


WEIGHTS = (1, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("incidence_report", type=Path)
    parser.add_argument("nparameters", type=int)
    parser.add_argument("support_size", type=int)
    parser.add_argument("receipt", type=Path)
    args = parser.parse_args()

    source_hash = sha256(args.source)
    incidence = json.loads(args.incidence_report.read_text(encoding="utf-8"))
    if source_hash != incidence["source_sha256"]:
        raise ValueError("incidence report does not match the exact source hash")
    supports = incidence["most_repeated_supports"][str(args.support_size)]
    if len(supports) != 1:
        raise ValueError("requested support size is not a unique support set")
    support = tuple(supports[0]["variables"])
    if len(support) != len(WEIGHTS):
        raise ValueError("the verified weight vector is for exactly 12 variables")
    symbols = sp.symbols(" ".join(f"p{i}" for i in support))

    polynomials = []
    tags = []
    weighted_degrees = Counter()
    total_terms = 0
    start = time.monotonic()
    print(f"source_sha256={source_hash}", flush=True)
    print(f"support={[f'p{i}' for i in support]} weights={WEIGHTS}", flush=True)

    for tag, saved in streamed_array(args.source):
        term_map = defaultdict(lambda: sp.S.Zero)
        parameter_support = set()
        for sparse, values in saved:
            exponent_map = {i - 1: e for i, e in sparse
                            if 1 <= i <= args.nparameters and e}
            parameter_support.update(exponent_map)
            exponent = tuple(exponent_map.get(i, 0) for i in support)
            coefficient = sp.Rational(values[0]) if len(values) == 1 else None
            if coefficient is None:
                raise ValueError("support-12 input is not rational over Q")
            term_map[exponent] += coefficient
            total_terms += 1
        if tuple(sorted(parameter_support)) != support:
            continue
        term_map = {exponent: coefficient for exponent, coefficient in term_map.items()
                    if coefficient}
        degrees = {sum(weight*power for weight, power in zip(WEIGHTS, exponent))
                   for exponent in term_map}
        if len(degrees) != 1:
            raise AssertionError(f"weight vector fails on {tag}: {sorted(degrees)}")
        weighted_degrees[next(iter(degrees))] += 1
        polynomial = sp.Poly.from_dict(term_map, *symbols, domain=QQ)
        if polynomial:
            polynomials.append(polynomial.as_expr())
            tags.append(tag)

    print(f"selected={len(polynomials)} terms_scanned={total_terms} "
          f"weighted_degrees={dict(sorted(weighted_degrees.items()))} "
          f"parse_seconds={time.monotonic()-start:.2f}", flush=True)
    if not polynomials:
        raise ValueError("selected support subideal is empty")

    gb_start = time.monotonic()
    print("starting exact QQ grevlex F5B", flush=True)
    basis = sp.groebner(polynomials, *symbols, domain=QQ,
                        order="grevlex", method="f5b")
    gb_seconds = time.monotonic() - gb_start
    print(f"groebner_done basis_size={len(basis.polys)} "
          f"zero_dimensional={basis.is_zero_dimensional} "
          f"seconds={gb_seconds:.2f}", flush=True)

    # A homogeneous ideal whose radical is the origin contains a power of
    # every variable. Check those ideal-membership claims by exact reduction.
    nilpotence_powers = {}
    for symbol in symbols:
        found = None
        for power in range(1, 65):
            remainder = basis.reduce(symbol**power)[1]
            if remainder == 0:
                found = power
                break
        nilpotence_powers[str(symbol)] = found
        print(f"variable_power={symbol}^{found if found else 'not_found_through_64'}",
              flush=True)

    input_remainders_zero = all(basis.reduce(polynomial)[1] == 0
                                for polynomial in polynomials)
    result = {
        "source": str(args.source),
        "source_sha256": source_hash,
        "branch": "2+1+1",
        "scope": "exact support-12 subideal from detached depth-13 necessary equations; not a whole-branch exclusion",
        "coefficient_field": "Q",
        "monomial_order": "grevlex",
        "algorithm": "SymPy F5B exact Buchberger-family computation",
        "support_variables": [str(symbol) for symbol in symbols],
        "positive_integer_weights": list(WEIGHTS),
        "selected_constraints": len(polynomials),
        "selected_tags": tags,
        "weighted_degree_histogram": dict(sorted(weighted_degrees.items())),
        "basis_size": len(basis.polys),
        "zero_dimensional": basis.is_zero_dimensional,
        "basis": [str(poly.as_expr()) for poly in basis.polys],
        "variable_nilpotence_powers": nilpotence_powers,
        "all_input_generators_reduce_to_zero": input_remainders_zero,
        "groebner_seconds": gb_seconds,
        "warning": "The 12 support variables are only part of the 75 carried parameters. The remaining necessary equations and nonzero/open guards still have to be combined before any branch percentage can change.",
    }
    args.receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"receipt={args.receipt}", flush=True)


if __name__ == "__main__":
    main()
