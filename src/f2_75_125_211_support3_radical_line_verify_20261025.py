"""Replay the radical line forced by the exact 2+1+1 support-three basis."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "f2_75_125_sparse_ring_lift_211_roots19_u1zero_streamed_state.json"
BASIS = ROOT / "f2_75_125_sparse_ring_lift_211_roots19_u1zero_support3_gb.json"
OUTPUT = ROOT / "f2_75_125_211_support3_radical_line_20261025.json"
EXPECTED_SOURCE_SHA256 = "0f293c437b1cc18f379e3a8d0bb2e98647616e062fa321e91cfdc45ebe67926d"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    report = json.loads(BASIS.read_text(encoding="utf-8"))
    if report.get("source_sha256") != EXPECTED_SOURCE_SHA256:
        raise ValueError("support-three basis is not bound to the d19 source")
    if sha256(SOURCE) != EXPECTED_SOURCE_SHA256:
        raise ValueError("d19 source hash mismatch")
    if report.get("branch") != "2+1+1":
        raise ValueError("wrong factor branch")
    if report.get("selected_variables") != ["p13", "p14", "p15"]:
        raise ValueError("wrong support-three variable set")
    if report.get("groebner_method") != "f5b" or report.get("monomial_order") != "grevlex":
        raise ValueError("unexpected recorded Groebner computation")

    p13, p14, p15 = sp.symbols("p13 p14 p15")
    basis = [sp.sympify(value, locals={"p13": p13, "p14": p14, "p15": p15})
             for value in report["groebner_basis"]]
    if len(basis) != 4:
        raise ValueError("expected the recorded four-element basis")
    if sp.expand(basis[0] - (p14 + 14*p15)**5) != 0:
        raise AssertionError("first basis element is not the exact fifth power")

    line_substitution = {p14: -14*p15}
    restricted = [sp.factor(poly.subs(line_substitution)) for poly in basis]
    expected = [
        sp.S.Zero,
        sp.S.Zero,
        (5*p13 - 452*p15)**3 / 125,
        sp.S.Zero,
    ]
    if any(sp.expand(a-b) != 0 for a, b in zip(restricted, expected)):
        raise AssertionError("basis restriction does not give the claimed cube")

    line_substitution[p13] = sp.Rational(452, 5)*p15
    if any(sp.expand(poly.subs(line_substitution)) != 0 for poly in basis):
        raise AssertionError("the claimed line is not contained in the basis zero set")

    result = {
        "scope": (
            "Exact set-theoretic consequence of the recorded complete support-three "
            "Groebner basis only; it is not a 2+1+1 or degree-pair exclusion."
        ),
        "branch": "2+1+1, U1=0, d19 necessary system",
        "source_sha256": EXPECTED_SOURCE_SHA256,
        "basis_report_sha256": sha256(BASIS),
        "variables": ["p13", "p14", "p15"],
        "support_three_input_equations": report["input_equations"],
        "unique_equations_up_to_rational_scale": report["unique_equations_up_to_nonzero_rational_scale"],
        "exact_relations": [
            "(p14 + 14*p15)^5 = 0",
            "(5*p13 - 452*p15)^3 = 0 after p14 = -14*p15",
        ],
        "geometric_zero_set": "(p13,p14,p15) = (452*s/5,-14*s,s), s in algebraic closure",
        "line_containment_verified": True,
        "basis_restrictions_verified": True,
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
