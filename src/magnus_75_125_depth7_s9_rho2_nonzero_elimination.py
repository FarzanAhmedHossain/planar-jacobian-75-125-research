"""Replay a small exact elimination on the s=9, rho2!=0 Magnus leaf.

The input equations are extracted from the saved depth-7 local-split report.
The calculation uses only the explicit branch guard rho2!=0 and the already
stated mu!=0; it is a further necessary reduction, not a branch exclusion.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
REPORT = HERE / "magnus_75_125_tail_remainder_depth7_s9_local_split_20261020.json"
RECEIPT = HERE / "magnus_75_125_depth7_s9_rho2_nonzero_elimination_20261020.json"


def main() -> None:
    report_bytes = REPORT.read_bytes()
    report = json.loads(report_bytes)
    source_log = HERE / report["source_log"]
    source_sha256 = hashlib.sha256(source_log.read_bytes()).hexdigest()
    if source_sha256 != report["source_sha256"]:
        raise SystemExit("depth-7 source log digest does not match its local-split report")

    mu, rho2, d2, u3, c4, d3, u4, sigma1 = sp.symbols(
        "mu rho2 d2 u3 c4 d3 u4 sigma1"
    )
    local = {str(s): s for s in (mu, rho2, d2, u3, c4, d3, u4, sigma1)}
    face = report["remaining_rows_after_split"]["rho2_nonzero_face"]

    row_d2_u3_a = sp.sympify(face["depth5_x7"], locals=local)
    row_d2_u3_b = sp.sympify(face["depth7_x4"], locals=local)
    row_c4 = sp.sympify(face["depth6_x6"], locals=local)

    expected_a = 36 * mu * rho2 * (-11 * d2 + sp.Rational(57, 5) * rho2 * u3)
    expected_b = sp.Rational(18, 5) * mu * rho2**2 * (-103 * d2 + 107 * rho2 * u3)
    if sp.expand(row_d2_u3_a - expected_a) != 0:
        raise SystemExit("saved depth5 row does not match the audited linear system")
    if sp.expand(row_d2_u3_b - expected_b) != 0:
        raise SystemExit("saved depth7 row does not match the audited linear system")

    matrix = sp.Matrix([
        [-11, sp.Rational(57, 5)],
        [-103, 107],
    ])
    determinant = sp.factor(matrix.det())
    if determinant != -sp.Rational(14, 5):
        raise SystemExit(f"unexpected determinant: {determinant}")

    # The equations are in d2 and v=rho2*u3.  A nonzero determinant gives
    # d2=v=0; the branch guard rho2!=0 then gives u3=0.
    reduced_c4_row = sp.factor(row_c4.subs({d2: 0, u3: 0}))
    expected_c4_row = sp.Rational(459, 5) * mu * c4 * rho2**2
    if sp.expand(reduced_c4_row - expected_c4_row) != 0:
        raise SystemExit("saved depth6 row does not reduce to the claimed c4 equation")

    next_row = sp.sympify(
        report["rho2_nonzero_refinement"]["next_exact_row"], locals=local
    )
    expected_next_row = sp.Rational(18, 125) * mu * rho2**2 * (
        725 * d3 - 738 * rho2 * u4
    )
    if sp.expand(next_row - expected_next_row) != 0:
        raise SystemExit("saved next row does not match the claimed d3/u4 relation")
    next_relation = sp.solve(next_row, d3)
    if next_relation != [sp.Rational(738, 725) * rho2 * u4]:
        raise SystemExit(f"unexpected next-row solution: {next_relation}")

    receipt = {
        "source_report": REPORT.name,
        "source_report_sha256": hashlib.sha256(report_bytes).hexdigest(),
        "source_log": source_log.name,
        "source_log_sha256": source_sha256,
        "leaf": report["leaf"],
        "subleaf": "rho2!=0, c2=u2=c3=0",
        "coefficient_field": "Q",
        "assumptions": ["mu!=0", "rho2!=0", "prior exact branch equations in the source report"],
        "rows": [
            "depth5_x7 = 36*mu*rho2*(-11*d2 + 57*rho2*u3/5)",
            "depth7_x4 = 18*mu*rho2^2*(-103*d2 + 107*rho2*u3)/5",
            "depth6_x6 after d2=u3=0 = 459*mu*c4*rho2^2/5",
        ],
        "linear_system_matrix_in_(d2,rho2_u3)": [[-11, "57/5"], [-103, 107]],
        "determinant": str(determinant),
        "forced_zero_parameters": ["d2", "u3", "c4"],
        "next_row": "18*mu*rho2^2*(725*d3 - 738*rho2*u4)/125 = 0",
        "next_forced_relation": "d3 = (738/725)*rho2*u4",
        "next_relation_assumptions": ["mu!=0", "rho2!=0"],
        "scope": "Exact necessary-equation reduction on one s=9 root-10 leaf; not an exclusion of F3 or either full F2 branch.",
    }
    RECEIPT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
