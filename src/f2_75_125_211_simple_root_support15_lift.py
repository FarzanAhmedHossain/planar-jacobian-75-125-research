#!/usr/bin/env python3
"""Pull the exact d19 support-15 ideal back to the simple-root line.

This intersects a previously certified support-15 necessary ideal with the
support-3 radical line and the two exact simple-root equations solving for
p34,p35.  It is a scoped compatibility calculation, not a branch exclusion.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ


HERE = Path(__file__).resolve().parent
SUPPORT15 = HERE / "f2_75_125_sparse_ring_lift_211_roots19_u1zero_support15_gb.json"
MEMBERSHIP = HERE / "f2_75_125_211_simple_root_reduction_membership_20261025.json"
SUPPORT3_LINE = HERE / "f2_75_125_211_support3_radical_line_20261025.json"
OUTPUT = HERE / "f2_75_125_211_simple_root_support15_lift_20261028.json"
EXPECTED_D19_SHA256 = "0f293c437b1cc18f379e3a8d0bb2e98647616e062fa321e91cfdc45ebe67926d"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    support_hash = sha256(SUPPORT15)
    membership_hash = sha256(MEMBERSHIP)
    line_hash = sha256(SUPPORT3_LINE)
    support = json.loads(SUPPORT15.read_text(encoding="utf-8"))
    membership = json.loads(MEMBERSHIP.read_text(encoding="utf-8"))
    line_report = json.loads(SUPPORT3_LINE.read_text(encoding="utf-8"))

    if support.get("source_sha256") != EXPECTED_D19_SHA256:
        raise ValueError("support-15 report is not bound to the expected d19 source")
    if membership.get("d19_source_sha256") != EXPECTED_D19_SHA256:
        raise ValueError("simple-root report is not bound to the expected d19 source")
    if membership.get("support15_report_sha256") not in (None, support_hash):
        raise ValueError("simple-root report is bound to a different support-15 report")
    if line_report.get("source_sha256") != EXPECTED_D19_SHA256:
        raise ValueError("support-3 line report has a different d19 source")
    if not line_report.get("line_containment_verified") or not line_report.get("basis_restrictions_verified"):
        raise ValueError("support-3 line report lacks exact verification flags")
    if support.get("branch") != "2+1+1" or support.get("unit_ideal"):
        raise ValueError("expected the exact nonunit 2+1+1 support-15 ideal")
    if support.get("groebner_basis_length") != 4 or len(support.get("groebner_basis", [])) != 4:
        raise ValueError("support-15 Groebner basis report is incomplete")

    names = tuple(support["selected_variables"])
    symbols = {name: sp.Symbol(name) for name in names}
    p5, p13, p14, p15 = (symbols[name] for name in ("p5", "p13", "p14", "p15"))
    r = sp.Symbol("r")
    target_names = ("p5", "r", "p24", "p25", "p26", "p27", "p28",
                    "p36", "p37", "p38", "p39")
    target_symbols = tuple(sp.Symbol(name) for name in target_names)
    target = dict(zip(target_names, target_symbols))

    solution = membership["restriction_to_known_support3_line"]["solution"]
    p34_solution = sp.sympify(solution["p34"], locals={**symbols, "r": r})
    p35_solution = sp.sympify(solution["p35"], locals={**symbols, "r": r})
    substitutions = {
        p13: sp.Rational(452, 5) * r,
        p14: -14 * r,
        p15: r,
        symbols["p34"]: p34_solution,
        symbols["p35"]: p35_solution,
    }
    substitutions.update({symbols[name]: target[name] for name in target_names
                          if name in symbols and name != "r"})
    substitutions[r] = target["r"]

    pulled_back = [
        sp.Poly(sp.sympify(expr, locals=symbols).subs(substitutions),
                *target_symbols, domain=QQ)
        for expr in support["groebner_basis"]
    ]
    nonzero = [poly.as_expr() for poly in pulled_back if not poly.is_zero]
    gb = sp.groebner(nonzero, *target_symbols, order="grevlex", domain=QQ,
                     method="f5b")
    basis = [sp.Poly(poly.as_expr(), *target_symbols, domain=QQ) for poly in gb.polys]
    unit = any(poly.as_expr() == 1 for poly in gb.polys)

    result = {
        "scope": (
            "pullback of the exact d19 support-15 ideal along the support-3 radical line "
            "and the exact added simple-root equations; necessary subsystem only"
        ),
        "branch": "2+1+1, U1=0",
        "d19_source_sha256": EXPECTED_D19_SHA256,
        "support15_report_sha256": support_hash,
        "simple_root_membership_report_sha256": membership_hash,
        "support3_line_report_sha256": line_hash,
        "support15_basis_length": len(pulled_back),
        "nonzero_pulled_back_generators": len(nonzero),
        "variables": list(target_names),
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
        "scope", "support15_basis_length", "nonzero_pulled_back_generators",
        "variables", "groebner_basis_length", "unit_ideal", "zero_dimensional",
        "groebner_basis", "conclusion", "whole_2+1+1_u1zero_branch_excluded_percent"
    )}, indent=2), flush=True)
    print(f"report={OUTPUT}", flush=True)


if __name__ == "__main__":
    main()
