#!/usr/bin/env python3
"""Intersect the exact d19 support-10 ideal with simple-root equations.

This is a scoped coefficient-compatibility probe for the 2+1+1, U1=0
branch.  It pulls back the complete-support support-10 ideal along the exact
support-3 line and the two added simple-root equations, solving only for
p34,p35.  The remaining line parameter r is unrestricted.  The result is a
necessary subsystem, not a whole-branch exclusion.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ


HERE = Path(__file__).resolve().parent
SUPPORT10 = HERE / "f2_75_125_sparse_ring_lift_211_roots19_u1zero_support10_gb.json"
MEMBERSHIP = HERE / "f2_75_125_211_simple_root_reduction_membership_20261025.json"
EXPECTED_D19_SHA256 = "0f293c437b1cc18f379e3a8d0bb2e98647616e062fa321e91cfdc45ebe67926d"
OUTPUT = HERE / "f2_75_125_211_simple_root_support10_lift_20261028.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    support_hash = digest(SUPPORT10)
    membership_hash = digest(MEMBERSHIP)
    membership = json.loads(MEMBERSHIP.read_text(encoding="utf-8"))
    if membership["d19_source_sha256"] != EXPECTED_D19_SHA256:
        raise ValueError("simple-root report is not bound to the expected d19 source")
    if membership["support10_report_sha256"] != support_hash:
        raise ValueError("support-10 report changed since the membership analysis")
    if membership["scope"] != "exact ideal-membership check for added simple-root equations only":
        raise ValueError("unexpected membership-report scope")

    support = json.loads(SUPPORT10.read_text(encoding="utf-8"))
    if support["source_sha256"] != EXPECTED_D19_SHA256:
        raise ValueError("support-10 report is not bound to the expected d19 source")
    variables = [sp.Symbol(name) for name in
                 ("p5", "r", "p36", "p37", "p38", "p39")]
    p5, r, p36, p37, p38, p39 = variables
    names = {name: sp.Symbol(name) for name in
             ("p5", "p13", "p14", "p15", "p34", "p35", "p36", "p37", "p38", "p39", "r")}
    line = {
        names["p13"]: sp.Rational(452, 5) * names["r"],
        names["p14"]: -14 * names["r"],
        names["p15"]: names["r"],
    }
    p34_solution = sp.sympify(
        membership["restriction_to_known_support3_line"]["solution"]["p34"],
        locals=names,
    )
    p35_solution = sp.sympify(
        membership["restriction_to_known_support3_line"]["solution"]["p35"],
        locals=names,
    )
    substitutions = {
        **line,
        names["p34"]: p34_solution,
        names["p35"]: p35_solution,
    }
    pulled_back = [
        sp.Poly(sp.sympify(text, locals=names).subs(substitutions), *variables, domain=QQ)
        for text in support["groebner_basis"]
    ]
    nonzero = [poly.as_expr() for poly in pulled_back if not poly.is_zero]
    gb = sp.groebner(nonzero, *variables, order="grevlex", domain=QQ, method="f5b")
    basis = [sp.Poly(poly.as_expr(), *variables, domain=QQ) for poly in gb.polys]
    unit = any(poly.as_expr() == 1 for poly in gb.polys)
    result = {
        "scope": "pullback of the exact d19 complete-support support-10 ideal along the support-3 line and added simple-root equations; necessary subsystem only",
        "branch": "2+1+1, U1=0",
        "d19_source_sha256": EXPECTED_D19_SHA256,
        "support10_sha256": support_hash,
        "simple_root_membership_report_sha256": membership_hash,
        "support10_basis_length": len(pulled_back),
        "nonzero_pulled_back_generators": len(nonzero),
        "variables": [str(variable) for variable in variables],
        "substitutions": {
            "p13": "452*r/5",
            "p14": "-14*r",
            "p15": "r",
            "p34": str(p34_solution),
            "p35": str(p35_solution),
        },
        "groebner_order": "grevlex over QQ via SymPy F5B",
        "groebner_basis_length": len(basis),
        "unit_ideal": unit,
        "zero_dimensional": gb.is_zero_dimensional,
        "groebner_basis": [str(poly.as_expr()) for poly in gb.polys],
        "whole_2+1+1_u1zero_branch_excluded_percent": 0.0,
        "conclusion": (
            "the pulled-back necessary ideal is the unit ideal"
            if unit else
            "the pulled-back necessary ideal remains nonunit; no exclusion follows"
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "scope", "support10_basis_length", "nonzero_pulled_back_generators",
        "variables", "groebner_basis_length", "unit_ideal", "zero_dimensional",
        "groebner_basis", "conclusion", "whole_2+1+1_u1zero_branch_excluded_percent"
    )}, indent=2), flush=True)
    print(f"report={OUTPUT}", flush=True)


if __name__ == "__main__":
    main()
