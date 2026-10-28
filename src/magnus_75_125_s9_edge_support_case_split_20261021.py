"""Replay exact case splits from the support-tight Magnus s=9 log.

This extracts only necessary local coefficient rows. It does not claim a
Keller-pair or F2-branch exclusion.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parent
LOG = ROOT / "magnus_75_125_tail_remainder_depth8_s9_edge_support_refinement_20261021.log"
SOURCE = ROOT / "magnus_75_125_tail_remainder_depth8_s9_edge_support_refinement_20261021.py"
RECEIPT = ROOT / "magnus_75_125_s9_edge_support_case_split_20261021.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_rhs(line: str) -> sp.Expr:
    rhs = line.split("=", 1)[1]
    names = set(re.findall(r"\b[A-Za-z_]\w*\b", rhs))
    local = {name: sp.Symbol(name) for name in names}
    return sp.sympify(rhs, locals=local)


def one_line(lines: list[str], prefix: str) -> str:
    matches = [line for line in lines if line.startswith(prefix)]
    if len(matches) != 1:
        raise AssertionError(f"expected one {prefix}, found {len(matches)}")
    return matches[0]


def main() -> None:
    lines = LOG.read_text(encoding="utf-8").splitlines()
    row7_lines = sorted(
        (line for line in lines if re.match(r"^depth7_fixed_x\d+=", line)),
        key=lambda line: int(re.search(r"x(\d+)=", line).group(1)),
    )
    row8_lines = sorted(
        (line for line in lines if re.match(r"^depth8_fixed_x\d+=", line)),
        key=lambda line: int(re.search(r"x(\d+)=", line).group(1)),
    )
    if len(row7_lines) != 12 or len(row8_lines) != 8:
        raise AssertionError("exact-support fixed-row census changed")
    row7_by_degree: dict[int, sp.Expr] = {}
    for line in row7_lines:
        degree = int(re.search(r"x(\d+)=", line).group(1))
        value = parse_rhs(line)
        if degree in row7_by_degree and sp.expand(row7_by_degree[degree]-value) != 0:
            raise AssertionError("the repeated depth-7 row blocks differ")
        row7_by_degree[degree] = value
    if len(row7_by_degree) != 6:
        raise AssertionError("expected six distinct depth-7 fixed rows")

    mu, sigma0, sigma1, rho2, a4, c2, c3, u2, d2, u3 = sp.symbols(
        "mu sigma0 sigma1 rho2 a4 c2 c3 u2 d2 u3"
    )
    rows7 = row7_by_degree
    rows8 = {int(re.search(r"x(\d+)=", line).group(1)): parse_rhs(line)
             for line in row8_lines}

    z6_line = one_line(lines, "z6_remaining_commutator_rows=")
    if z6_line != "z6_remaining_commutator_rows=[(3, 612*mu*sigma0**2/5)]":
        raise AssertionError("z6 residual changed")
    z7_line = one_line(lines, "z7_remaining_commutator_rows=")
    z7_expr = sp.sympify(
        z7_line.split("=", 1)[1],
        locals={str(symbol): symbol for symbol in (
            mu, sigma0, sigma1, rho2, a4, sp.Symbol("rho3"), sp.Symbol("tau0")
        )},
    )
    z7_rows = dict(z7_expr)
    if set(z7_rows) != {2, 3}:
        raise AssertionError("support-tight z7 residual degrees changed")

    # Since mu is the nonzero chart scale, the z6 residual forces sigma0=0.
    base = {sigma0: 0}
    z7_x3_zero_sigma = sp.factor(z7_rows[3].subs(base))
    expected_z7 = sp.Rational(594, 125) * mu * rho2**2 * (2*a4*rho2 - 5*sigma1)
    if sp.expand(z7_x3_zero_sigma - expected_z7) != 0:
        raise AssertionError("zero-sigma z7 factorization changed")

    # In the rho2 != 0 stratum, the following deductions use only exact rows:
    # z7 fixes sigma1; depth8 x1 fixes c2; depth8 x2 fixes c3; depth7 x3 fixes u2.
    nonzero_rho = {
        **base,
        sigma1: sp.Rational(2, 5) * a4 * rho2,
        c2: 0,
        c3: 0,
        u2: 0,
    }
    chain = [
        sp.factor(z7_x3_zero_sigma.subs(sigma1, sp.Rational(2, 5)*a4*rho2)),
        sp.factor(rows8[1].subs(base)),
        sp.factor(rows7[3].subs({**base, c2: 0})),
        sp.factor(rows8[2].subs({
            **base, sigma1: sp.Rational(2, 5)*a4*rho2, c2: 0, u2: 0
        })),
    ]
    expected_chain = [
        0,
        sp.Rational(72, 25) * mu * c2 * rho2**3,
        sp.Rational(10944, 125) * mu * rho2**3 * u2,
        sp.Rational(12312, 125) * mu * c3 * rho2**3,
    ]
    if any(sp.expand(a-b) != 0 for a, b in zip(chain, expected_chain)):
        raise AssertionError(
            f"rho2-nonzero deduction chain changed: actual={chain}; expected={expected_chain}"
        )

    # The next fixed row gives one more exact graph relation for d2.
    x4_after = sp.factor(rows7[4].subs(nonzero_rho))
    d2_relation = sp.solve(sp.together(x4_after), d2)
    if len(d2_relation) != 1 or sp.simplify(d2_relation[0] - sp.Rational(107, 103)*rho2*u3) != 0:
        raise AssertionError("depth7 x4 relation changed")
    graph = {**nonzero_rho, d2: d2_relation[0]}

    remaining = []
    for band, row_map in (("depth7", rows7), ("depth8", rows8)):
        for degree, equation in row_map.items():
            reduced = sp.factor(equation.subs(graph))
            if reduced != 0:
                remaining.append({
                    "band": band,
                    "x_degree": degree,
                    "factorized_equation": str(reduced),
                })

    result = {
        "source_sha256": sha256(SOURCE),
        "log_sha256": sha256(LOG),
        "checked_depth7_fixed_rows": len(rows7),
        "repeated_depth7_output_blocks_identical": True,
        "checked_depth8_fixed_rows": len(rows8),
        "exact_branch": {
            "guards": ["mu != 0", "rho2 != 0"],
            "forced": [
                "sigma0=0",
                "sigma1=2*a4*rho2/5",
                "c2=0",
                "c3=0",
                "u2=0",
                "d2=107*rho2*u3/103",
            ],
            "remaining_nonzero_rows_after_graph_substitution": remaining,
        },
        "scope": "exact necessary local reduction only; no complete s=9 leaf or F2 branch excluded",
    }
    RECEIPT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
