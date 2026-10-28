"""Replay the saved depth-7/8 rows on the s=9, sigma0=rho2=0 leaf.

This is a necessary-equation reduction only. It does not claim that the leaf,
the two-factor carrier, or either F2 quartic branch is empty.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "magnus_75_125_tail_remainder_depth8_s9_a3rho1zero_20261020.py"
LOG = ROOT / "magnus_75_125_tail_remainder_depth8_s9_a3rho1zero_20261020.log"
RECEIPT = ROOT / "magnus_75_125_depth8_s9_rho2zero_reduction_20261021.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rows(prefix: str, degrees: range) -> dict[int, sp.Expr]:
    lines = LOG.read_text(encoding="utf-8").splitlines()
    rows: dict[int, sp.Expr] = {}
    for degree in degrees:
        marker = f"{prefix}{degree}="
        hits = [line for line in lines if line.startswith(marker)]
        unique_hits = sorted(set(hits))
        if len(unique_hits) != 1:
            raise AssertionError(
                f"expected one unique row {marker}, found {len(unique_hits)}"
            )
        source = unique_hits[0].split("=", 1)[1]
        names = set(re.findall(r"\b[A-Za-z_]\w*\b", source))
        local = {name: sp.Symbol(name, nonzero=True) if name == "mu"
                 else sp.Symbol(name) for name in names}
        rows[degree] = sp.sympify(source, locals=local)
    return rows


def main() -> None:
    if "mu = sp.symbols(\"mu\", nonzero=True)" not in SOURCE.read_text(encoding="utf-8"):
        raise AssertionError("source does not declare the chart scale mu nonzero")

    substitutions = {sp.Symbol("sigma0"): 0, sp.Symbol("rho2"): 0}
    specialized: dict[str, list[dict[str, object]]] = {}
    exact_rows: dict[tuple[str, int], sp.Expr] = {}
    for prefix, degrees in (("depth7_fixed_x", range(2, 8)),
                            ("depth8_fixed_x", range(0, 8))):
        items = []
        for degree, source_row in read_rows(prefix, degrees).items():
            reduced = sp.factor_terms(sp.expand(source_row.subs(substitutions)))
            exact_rows[(prefix, degree)] = reduced
            if reduced == 0:
                continue
            variables = sorted(str(symbol) for symbol in reduced.free_symbols)
            items.append({
                "x_degree": degree,
                "equation": str(reduced),
                "term_count": len(sp.Add.make_args(sp.expand(reduced))),
                "variables": variables,
            })
        specialized[prefix] = items

    mu = sp.Symbol("mu", nonzero=True)
    c2, c3, sigma1, u2, a4, tau0 = sp.symbols(
        "c2 c3 sigma1 u2 a4 tau0"
    )
    expected_x2 = -sp.Rational(11493, 25) * mu * c2 * sigma1**2
    expected_x4 = sp.Rational(9, 25) * mu * (
        a4**2 * c2 * sigma1 - 1690 * a4 * c2 * tau0
        - 1265 * c2 * sp.Symbol("rho3") * sigma1 + 140 * sigma1**2 * u2
    )
    if exact_rows[("depth8_fixed_x", 2)] != expected_x2:
        raise AssertionError("the rho2=0 depth-eight x^2 factor changed")
    if exact_rows[("depth7_fixed_x", 4)] != expected_x4:
        raise AssertionError("the rho2=0 depth-seven x^4 factor changed")
    x4_c2_zero = sp.factor_terms(expected_x4.subs(c2, 0))
    x4_sigma1_zero = sp.factor_terms(expected_x4.subs(sigma1, 0))
    if x4_c2_zero != sp.Rational(252, 5) * mu * sigma1**2 * u2:
        raise AssertionError("the c2=0 x^4 split changed")
    if x4_sigma1_zero != -sp.Rational(3042, 5) * mu * a4 * c2 * tau0:
        raise AssertionError("the sigma1=0 x^4 split changed")

    row8_x3 = exact_rows[("depth8_fixed_x", 3)]
    x3_c2_u2_zero = sp.factor_terms(row8_x3.subs({c2: 0, u2: 0}))
    expected_x3_c2_u2_zero = -sp.Rational(11133, 25) * mu * c3 * sigma1**2
    if x3_c2_u2_zero != expected_x3_c2_u2_zero:
        raise AssertionError("the c2=u2=0 x^3 split changed")

    result = {
        "source_sha256": digest(SOURCE),
        "log_sha256": digest(LOG),
        "specialization": ["sigma0=0", "rho2=0"],
        "exact_splits": [
            "depth8 x^2: c2*sigma1^2=0, hence c2=0 or sigma1=0",
            "depth7 x^4 with c2=0: sigma1^2*u2=0",
            "depth7 x^4 with sigma1=0: a4*c2*tau0=0",
            "depth8 x^3 with c2=u2=0: c3*sigma1^2=0",
        ],
        "nonzero_depth7_rows": specialized["depth7_fixed_x"],
        "nonzero_depth8_rows": specialized["depth8_fixed_x"],
        "scope": "necessary local equations only; no leaf or F2 branch excluded",
    }
    RECEIPT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    summary = {
        "specialization": result["specialization"],
        "nonzero_depth7_degrees": [r["x_degree"] for r in specialized["depth7_fixed_x"]],
        "nonzero_depth8_degrees": [r["x_degree"] for r in specialized["depth8_fixed_x"]],
        "depth7_term_counts": [r["term_count"] for r in specialized["depth7_fixed_x"]],
        "depth8_term_counts": [r["term_count"] for r in specialized["depth8_fixed_x"]],
        "receipt": RECEIPT.name,
        "scope": result["scope"],
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
