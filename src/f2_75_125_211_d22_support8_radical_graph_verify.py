#!/usr/bin/env python3
"""Verify the exact radical graph of the saved d22 support-eight subsystem.

This checks the three stored Groebner generators over QQ, the membership
relations L**3 in I and M**2 in I+(L), and the graph parametrization
L=M=0.  Its scope is only the support-eight subsystem of the saved
2+1+1, U1=0 checkpoint; it is not a branch-exclusion certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ


NAMES = ("p5", "p15", "p34", "p35", "p36", "p37", "p38", "p39")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("support8_report", type=Path)
    parser.add_argument("radical_report", type=Path)
    args = parser.parse_args()

    support = json.loads(args.support8_report.read_text(encoding="utf-8"))
    radical = json.loads(args.radical_report.read_text(encoding="utf-8"))
    if support.get("branch") != "2+1+1":
        raise ValueError("support-eight report is not for the 2+1+1 branch")
    if tuple(support.get("selected_variables", ())) != NAMES:
        raise ValueError("unexpected support-eight variable order")
    if support.get("monomial_order") != "grevlex":
        raise ValueError("unexpected monomial order")
    if radical.get("source_sha256") != support.get("source_sha256"):
        raise ValueError("radical scan and support-eight report use different checkpoints")
    if radical.get("support8_report_sha256") != sha256(args.support8_report):
        raise ValueError("support-eight report hash does not match the radical scan")

    variables = sp.symbols(" ".join(NAMES))
    symbols = dict(zip(NAMES, variables))
    basis = [
        sp.sympify(text, locals=symbols, rational=True)
        for text in support["groebner_basis"]
    ]
    if len(basis) != support.get("groebner_basis_length"):
        raise ValueError("stored Groebner basis length is inconsistent")

    p5, p15, p34, p35, p36, p37, p38, p39 = variables
    linear = (
        224725 * p34 + 1639725 * p35 + 6284625 * p36
        + 15845845 * p37 + 24712395 * p38 - 163942748 * p39
    )
    second = (
        375 * p35 + 3750 * p36 + 20625 * p37 + 82650 * p38
        - 3016650 * p39 + 251692 * p5 * p15
    )
    p35_value = (
        -sp.Rational(251692, 375) * p15 * p5 - 10 * p36 - 55 * p37
        - sp.Rational(1102, 5) * p38 + sp.Rational(40222, 5) * p39
    )
    p34_value = (
        sp.Rational(612164, 125) * p15 * p5 + 45 * p36
        + sp.Rational(1654, 5) * p37 + sp.Rational(7491, 5) * p38
        - sp.Rational(1449178, 25) * p39
    )
    graph = {p34: p34_value, p35: p35_value}

    ideal_gb = sp.groebner(basis, *variables, order="grevlex", domain=QQ)
    ideal_plus_linear_gb = sp.groebner(
        [*basis, linear], *variables, order="grevlex", domain=QQ
    )
    checks = {
        "L_cubed_in_I": ideal_gb.reduce(linear**3)[1] == 0,
        "M_squared_in_I_plus_L": ideal_plus_linear_gb.reduce(second**2)[1] == 0,
        "L_vanishes_on_graph": sp.expand(linear.subs(graph, simultaneous=True)) == 0,
        "M_vanishes_on_graph": sp.expand(second.subs(graph, simultaneous=True)) == 0,
        "all_basis_generators_vanish_on_graph": all(
            sp.expand(poly.subs(graph, simultaneous=True)) == 0 for poly in basis
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"exact check failed: {checks}")

    result = {
        "branch": support["branch"],
        "source_sha256": support["source_sha256"],
        "support8_report_sha256": sha256(args.support8_report),
        "basis_length": len(basis),
        "recomputed_groebner_basis_length": len(ideal_gb.polys),
        "groebner_basis_with_L_length": len(ideal_plus_linear_gb.polys),
        "checks": checks,
        "scope": (
            "The listed support-eight ideal has zero set exactly the graph "
            "p34=p34_value, p35=p35_value over an algebraically closed field "
            "of characteristic zero. This is a nonempty necessary subsystem, "
            "not a branch exclusion."
        ),
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
