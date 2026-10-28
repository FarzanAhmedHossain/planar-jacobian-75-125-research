#!/usr/bin/env python3
"""Restrict the exact support-eight basis to its forced linear section."""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

import sympy as sp


SOURCE_SHA256 = "e2e35061e883437806d6486d8f61fdc4cc6b963c49655366a6795db210d874e7"
REPORT = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support8_gb.json")
NAMES = ("p5", "p34", "p35", "p36", "p37", "p38", "p39", "r")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    data = json.loads(REPORT.read_text(encoding="utf-8"))
    if data.get("source_sha256") != SOURCE_SHA256:
        raise ValueError("support-eight report has the wrong checkpoint hash")
    if data.get("selected_variables") != list(NAMES) or len(data.get("groebner_basis", [])) != 3:
        raise ValueError("unexpected support-eight report")
    p5, p34, p35, p36, p37, p38, p39, r = sp.symbols(" ".join(NAMES))
    local = dict(zip(NAMES, (p5, p34, p35, p36, p37, p38, p39, r)))
    basis = [sp.sympify(expression, locals=local) for expression in data["groebner_basis"]]
    p34_section = (
        -1639725*p35 - 6284625*p36 - 15845845*p37
        - 24712395*p38 + 163942748*p39
    ) / sp.Integer(224725)
    substitutions = {p34: p34_section}
    restricted = [sp.cancel(polynomial.subs(substitutions)) for polynomial in basis]
    linear = 224725*p34 + 1639725*p35 + 6284625*p36 + 15845845*p37 + 24712395*p38 - 163942748*p39
    if sp.expand(basis[1] - linear**3/sp.Integer(11348910401078125)) != 0:
        raise ValueError("basis cube identity failed exact replay")
    if sp.expand(restricted[1]) != 0 or sp.expand(restricted[2]) != 0:
        raise ValueError("the remaining basis elements did not vanish on L=0")
    variables = (p5, p35, p36, p37, p38, p39, r)
    q = sp.Poly(restricted[0], *variables, domain=sp.QQ)
    if q.is_zero:
        raise ValueError("unexpectedly lost the quadratic support-eight relation")
    report = {
        "source_sha256": SOURCE_SHA256,
        "support8_report_sha256": sha256(REPORT),
        "support8_groebner_basis_length": len(basis),
        "linear_radical_relation": str(linear),
        "line_section": {
            "p34": "(-1639725*p35-6284625*p36-15845845*p37-24712395*p38+163942748*p39)/224725"
        },
        "basis_restrictions": [str(value) for value in restricted],
        "remaining_equation": str(q.as_expr()),
        "remaining_equation_terms": len(q.terms()),
        "remaining_equation_total_degree": q.total_degree(),
        "remaining_equation_grevlex_leading_monomial": str(q.terms(order="grevlex")[0][0]),
        "factorization": str(sp.factor(q.as_expr())),
        "interpretation": "On the forced linear section L=0, the complete support-eight ideal restricts to the single displayed nonzero equation Q=0; the other two Groebner basis elements vanish there.",
    }
    output = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support8_line_section.json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "source_sha256": SOURCE_SHA256,
        "support8_report_sha256": report["support8_report_sha256"],
        "remaining_equation_terms": report["remaining_equation_terms"],
        "remaining_equation_total_degree": report["remaining_equation_total_degree"],
        "remaining_equation_grevlex_leading_monomial": report["remaining_equation_grevlex_leading_monomial"],
        "factorization": report["factorization"],
        "report": str(output),
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
