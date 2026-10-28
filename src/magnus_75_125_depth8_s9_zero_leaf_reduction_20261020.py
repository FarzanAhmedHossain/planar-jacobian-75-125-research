"""Replay one exact branch split from the Magnus depth-8 log.

This checks only necessary equations on the s=9 zero-vector leaf. It is not
an exclusion certificate for the full leaf or either F2 branch.
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
RECEIPT = ROOT / "magnus_75_125_depth8_s9_zero_leaf_reduction_20261020.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_equation(prefix: str) -> sp.Expr:
    hits = [line for line in LOG.read_text(encoding="utf-8").splitlines()
            if line.startswith(prefix)]
    if not hits:
        raise AssertionError(f"missing output row {prefix}")
    expression = hits[0].split("=", 1)[1]
    identifiers = set(re.findall(r"\b[A-Za-z_]\w*\b", expression))
    local_symbols = {name: sp.Symbol(name) for name in identifiers}
    return sp.sympify(expression, locals=local_symbols)


def main() -> None:
    if "mu = sp.symbols(\"mu\", nonzero=True)" not in SOURCE.read_text(encoding="utf-8"):
        raise AssertionError("source does not declare the chart scale mu nonzero")

    mu, c2, c3, rho2, sigma0, sigma1, u2, a4 = sp.symbols(
        "mu c2 c3 rho2 sigma0 sigma1 u2 a4"
    )
    row7 = read_equation("depth7_fixed_x2=")
    row7_3 = read_equation("depth7_fixed_x3=")
    row8_0 = read_equation("depth8_fixed_x0=")
    row8_1 = read_equation("depth8_fixed_x1=")
    z6_rows = [line for line in LOG.read_text(encoding="utf-8").splitlines()
               if line.startswith("z6_remaining_commutator_rows=")]
    expected_z6 = "z6_remaining_commutator_rows=[(3, 612*mu*sigma0**2/5)]"
    if z6_rows != [expected_z6]:
        raise AssertionError("the exact z6 residual row changed")
    z7_rows = read_equation("z7_remaining_commutator_rows=")

    expected7 = -18 * mu * sigma0 * (sp.Rational(1281, 25) * c2 * rho2
                                      + 22 * sigma0 * u2)
    expected8_0 = -459 * c2 * mu * sigma0**2
    expected8_1 = 9 * mu * (
        sp.Rational(456, 125) * a4 * c2 * rho2 * sigma0
        + 20 * a4 * sigma0**2 * u2
        + sp.Rational(8, 25) * c2 * rho2**3
        - sp.Rational(558, 5) * c2 * sigma0 * sigma1
        - 47 * c3 * sigma0**2
        + sp.Rational(432, 125) * rho2**2 * sigma0 * u2
    )
    symbols = {str(s): s for s in (mu, c2, c3, rho2, sigma0,
                                    sigma1, u2, a4)}
    for actual, expected in ((row7, expected7),
                             (row8_0, expected8_0),
                             (row8_1, expected8_1)):
        if sp.expand(actual - expected) != 0:
            raise AssertionError("logged exact row does not match replay formula")

    # Branch sigma0 != 0: row 8, x^0 gives c2=0; row 7, x^2 then gives
    # u2=0; row 8, x^1 then gives c3=0.
    nonzero_sigma = [sp.factor(row8_0.subs(mu, 1)),
                     sp.factor(row7.subs({mu: 1, c2: 0})),
                     sp.factor(row8_1.subs({mu: 1, c2: 0, u2: 0}))]
    if nonzero_sigma != [-459 * c2 * sigma0**2,
                          -396 * sigma0**2 * u2,
                          -423 * c3 * sigma0**2]:
        raise AssertionError("nonzero-sigma0 branch reduction changed")

    # Branch sigma0 = 0: row 8, x^1 reduces to (72/25) mu*c2*rho2^3=0.
    zero_sigma = sp.factor(row8_1.subs({mu: 1, sigma0: 0}))
    if zero_sigma != sp.Rational(72, 25) * c2 * rho2**3:
        raise AssertionError("zero-sigma0 branch reduction changed")

    r13_8, r13_9, r13_10, rho3, tau0 = sp.symbols(
        "r13_8 r13_9 r13_10 rho3 tau0"
    )
    expected_z7 = [
        (18, -243 * mu * r13_10),
        (17, -207 * mu * r13_9),
        (16, -171 * mu * r13_8),
        (3, sp.Rational(594, 125) * mu * (
            6 * a4**2 * rho2 * sigma0 + 2 * a4 * rho2**3
            - 10 * a4 * sigma0 * sigma1 - 5 * rho2**2 * sigma1
            - 10 * rho2 * rho3 * sigma0 + 50 * sigma0 * tau0
        )),
        (2, sp.Rational(18, 25) * mu * sigma0
         * (102 * a4 * sigma0 - 37 * rho2**2)),
    ]
    if z7_rows != expected_z7:
        raise AssertionError("the exact z7 residual rows changed")
    z7_x3_zero = sp.factor(dict(z7_rows)[3].subs(sigma0, 0))
    expected_z7_x3_zero = (sp.Rational(594, 125) * mu * rho2**2
                             * (2 * a4 * rho2 - 5 * sigma1))
    if z7_x3_zero != expected_z7_x3_zero:
        raise AssertionError("zero-sigma0 z7 x^3 factorization changed")

    # The z6 commutator residual is unconditional and mu is nonzero, so only
    # sigma0=0 remains. On its c2=0, rho2!=0 subcase, depth7 x^3 forces u2=0.
    row7_3_zero = sp.factor(row7_3.subs({sp.Symbol("sigma0"): 0,
                                         sp.Symbol("c2"): 0}))
    expected_row7_3_zero = (sp.Rational(10944, 125) * sp.Symbol("mu")
                            * sp.Symbol("rho2")**3 * sp.Symbol("u2"))
    if row7_3_zero != expected_row7_3_zero:
        raise AssertionError("the c2=0 depth7 x^3 leaf changed")

    result = {
        "source_sha256": digest(SOURCE),
        "log_sha256": digest(LOG),
        "rows_checked": ["depth7_fixed_x2", "depth8_fixed_x0",
                         "depth8_fixed_x1", "depth7_fixed_x3",
                         "z6_remaining_commutator_rows",
                         "z7_remaining_commutator_rows"],
        "exact_implications": [
            "mu != 0 and z6 residual imply sigma0=0",
            "sigma0=0 implies c2*rho2^3=0",
            "sigma0=0, c2=0, rho2!=0 imply u2=0",
            "z7 residual rows at x^16,x^17,x^18 force r13_8=r13_9=r13_10=0",
            "sigma0=0 and z7 x^3 row imply rho2=0 or 5*sigma1=2*a4*rho2",
        ],
        "scope": "necessary local reduction only; no complete leaf or F2 branch excluded",
    }
    RECEIPT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
