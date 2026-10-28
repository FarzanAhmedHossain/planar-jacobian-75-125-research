"""Replay the exact local split from the Magnus depth-seven log.

This verifies only the leaf
  a3=rho1=u0=u1=c0=c1=0, with u2,c2,u3,c3 initially free.
It is not a certificate for either full F2 root-factor branch.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
LOG = HERE / "magnus_75_125_tail_remainder_depth7_s9_a3rho1zero_20261020_replay.log"
OUT = HERE / "magnus_75_125_tail_remainder_depth7_s9_local_split_20261020.json"


def get_rhs(lines: list[str], prefix: str) -> str:
    matches = [line[len(prefix):] for line in lines if line.startswith(prefix)]
    if len(matches) != 1:
        raise AssertionError(f"expected one {prefix!r} row, got {len(matches)}")
    return matches[0]


def main() -> None:
    raw = LOG.read_bytes()
    lines = raw.decode("utf-8").splitlines()
    branch = get_rhs(lines, "branch=")
    if branch != "s9,a3=rho1=u0=u1=c0=c1=0;u2,c2,u3,c3 free":
        raise AssertionError(f"unexpected source leaf: {branch}")

    z6_text = get_rhs(lines, "z6_remaining_commutator_rows=")
    expected_z6_text = "[(3, 612*mu*sigma0**2/5)]"
    if z6_text != expected_z6_text:
        raise AssertionError(f"unexpected z6 commutator obstruction: {z6_text}")

    k5_text = get_rhs(lines, "prior_k5_fixed_x7=")
    x3_text = get_rhs(lines, "depth7_fixed_x3=")
    x5_text = get_rhs(lines, "fixed_x5=")

    row_prefixes = {
        "depth5_x5": "depth5_fixed_x5=",
        "depth5_x6": "depth5_fixed_x6=",
        "depth5_x7": "depth5_fixed_x7=",
        "depth6_x4": "fixed_x4=",
        "depth6_x5": "fixed_x5=",
        "depth6_x6": "fixed_x6=",
        "depth6_x7": "fixed_x7=",
        "depth7_x2": "depth7_fixed_x2=",
        "depth7_x3": "depth7_fixed_x3=",
        "depth7_x4": "depth7_fixed_x4=",
        "depth7_x5": "depth7_fixed_x5=",
        "depth7_x6": "depth7_fixed_x6=",
        "depth7_x7": "depth7_fixed_x7=",
    }
    row_texts = {name: get_rhs(lines, prefix) for name, prefix in row_prefixes.items()}
    token_names = sorted(set(re.findall(r"\b[a-z][a-zA-Z0-9_]*\b", " ".join(
        (k5_text, x3_text, x5_text, "mu sigma0 c2 rho2 u2 c3", *row_texts.values())
    ))))
    syms = {name: sp.Symbol(name) for name in token_names}
    parse = lambda value: sp.sympify(value, locals=syms)
    mu, sigma0, c2 = syms["mu"], syms["sigma0"], syms["c2"]
    rho2, u2, c3 = syms["rho2"], syms["u2"], syms["c3"]

    k5 = parse(k5_text)
    x3 = parse(x3_text)
    x5 = parse(x5_text)
    row_exprs = {name: parse(value) for name, value in row_texts.items()}

    # Since mu is nonzero in the coefficient field, the z6 obstruction
    # 612*mu*sigma0^2/5 = 0 forces sigma0=0.
    sigma0_row = parse("612*mu*sigma0**2/5")
    if sp.expand(sigma0_row - sp.Rational(612, 5) * mu * sigma0**2) != 0:
        raise AssertionError("z6 row failed exact reconstruction")

    k5_after_sigma0 = sp.factor(k5.subs(sigma0, 0))
    expected_k5 = -sp.Rational(2016, 5) * mu * c2 * rho2
    if sp.expand(k5_after_sigma0 - expected_k5) != 0:
        raise AssertionError("K5 row did not reduce to c2*rho2=0")

    x3_c2_zero = sp.factor(x3.subs({sigma0: 0, c2: 0}))
    expected_x3 = sp.Rational(18 * 608, 125) * mu * rho2**3 * u2
    if sp.expand(x3_c2_zero - expected_x3) != 0:
        raise AssertionError("depth-seven x3 row did not reduce as expected")

    x5_c2_u2_zero = sp.factor(x5.subs({sigma0: 0, c2: 0, u2: 0}))
    expected_x5 = sp.Rational(9 * 1475, 125) * mu * c3 * rho2**2
    if sp.expand(x5_c2_u2_zero - expected_x5) != 0:
        raise AssertionError("depth-six x5 row did not reduce as expected")

    # Profile every logged residual on the two exact subbranches just found.
    # This is a specialization audit only; nonzero polynomial rows remain
    # constraints and are never divided by their parameter factors.
    branches = {
        "rho2_zero": {sigma0: 0, rho2: 0},
        "rho2_nonzero_face": {sigma0: 0, c2: 0, u2: 0, c3: 0},
    }
    branch_rows = {}
    for branch_name, substitutions in branches.items():
        residuals = {}
        for row_name, expression in row_exprs.items():
            value = sp.factor_terms(expression.subs(substitutions))
            if value != 0:
                residuals[row_name] = str(value)
        branch_rows[branch_name] = residuals

    def residual_rows(substitutions: dict[sp.Symbol, sp.Expr]) -> dict[str, str]:
        result = {}
        for row_name, expression in row_exprs.items():
            value = sp.factor_terms(expression.subs(substitutions))
            if value != 0:
                result[row_name] = str(value)
        return result

    # A first exhaustive split on the rho2=0 branch. The implications below
    # are kept as locally closed subcases; no factor is divided away silently.
    rho2_zero_base = {sigma0: 0, rho2: 0}
    a1_base = {**rho2_zero_base, sigma0: 0, syms["sigma1"]: 0,
               syms["tau0"]: 0}
    sigma2 = syms["sigma2"]
    a1_sigma2_equation = sp.factor_terms(
        row_exprs["depth5_x7"].subs(a1_base)
    )
    a1_sigma2_solutions = sp.solve(a1_sigma2_equation, sigma2, dict=True)
    if len(a1_sigma2_solutions) != 1:
        raise AssertionError(f"expected unique sigma2 relation, got {a1_sigma2_solutions}")
    a1_subs = {**a1_base, sigma2: a1_sigma2_solutions[0][sigma2]}

    a2_base = {**rho2_zero_base, c2: 0, syms["sigma1"]: 0}
    a2_tau0_zero = {**a2_base, syms["tau0"]: 0}
    a2_tau0_nonzero = {
        **a2_base, u2: 0, syms["tau0"]: syms["tau0"], c3: 0,
    }
    # On this last leaf tau0!=0 is assumed, and the depth-six x6 row forces
    # c3=0. The identity is verified without inverting tau0 in the checker.
    a2_pre = residual_rows({**a2_base, u2: 0})
    a2_x6 = parse(a2_pre["depth6_x6"])
    if sp.expand(a2_x6 + 594 * mu * c3 * syms["tau0"]) != 0:
        raise AssertionError("the 2+1+1 local row c3*tau0=0 changed")

    sigma1, u3 = syms["sigma1"], syms["u3"]
    a3_subs = {
        **rho2_zero_base, c2: 0, u2: 0, c3: 0, syms["d2"]: 0,
    }
    a3_rows = residual_rows(a3_subs)
    a3_x7 = parse(a3_rows["depth6_x7"])
    a3_d3_solutions = sp.solve(a3_x7, syms["d3"], dict=True)
    if len(a3_d3_solutions) != 1:
        raise AssertionError(f"expected one d3 relation, got {a3_d3_solutions}")
    a3_subs[syms["d3"]] = a3_d3_solutions[0][syms["d3"]]
    a3_later_rows = residual_rows(a3_subs)
    a3_u3_row = parse(a3_later_rows["depth7_x5"])
    expected_a3_u3_row = sp.Rational(144, 5) * mu * sigma1**2 * u3
    if sp.expand(a3_u3_row - expected_a3_u3_row) != 0:
        raise AssertionError("the rho2=0, sigma1!=0 branch did not force u3=0")
    a3_subs[u3] = 0
    a3_final_rows = residual_rows(a3_subs)

    # Refine the nonzero-rho2 face without dividing by any parameter that is
    # allowed to vanish. Here rho2!=0 is the explicit branch condition.
    b_subs = {sigma0: 0, c2: 0, u2: 0, c3: 0}
    b_first = sp.factor_terms(row_exprs["depth5_x7"].subs(b_subs))
    b_second = sp.factor_terms(row_exprs["depth7_x4"].subs(b_subs))
    d2 = syms["d2"]
    c4, d3, u4 = syms["c4"], syms["d3"], syms["u4"]
    # Divide only by the nonzero coefficients mu and rho2 on this face.
    b1 = sp.cancel(b_first / (36 * mu * rho2))
    b2 = sp.cancel(b_second / (sp.Rational(18, 5) * mu * rho2**2))
    b_matrix = sp.Matrix([
        [sp.Poly(row, d2, u3).coeff_monomial(d2),
         sp.Poly(row, d2, u3).coeff_monomial(u3)]
        for row in (b1, b2)
    ])
    b_det = sp.factor(b_matrix.det())
    if sp.expand(b_det + sp.Rational(14, 5) * rho2) != 0:
        raise AssertionError(f"unexpected branch-B determinant: {b_det}")
    b_solution = sp.solve([b1, b2], [d2, u3], dict=True)
    if b_solution != [{d2: 0, u3: 0}]:
        raise AssertionError(f"unexpected branch-B linear solution: {b_solution}")

    b_after_linear = {**b_subs, d2: 0, u3: 0}
    b_c4_row = sp.factor_terms(row_exprs["depth6_x6"].subs(b_after_linear))
    expected_c4_row = sp.Rational(459, 5) * mu * c4 * rho2**2
    if sp.expand(b_c4_row - expected_c4_row) != 0:
        raise AssertionError("branch-B depth-six x6 row did not force c4=0")
    b_after_c4 = {**b_after_linear, c4: 0}
    b_d3_u4_row = sp.factor_terms(row_exprs["depth7_x5"].subs(b_after_c4))
    expected_d3_u4_row = sp.Rational(18, 125) * mu * rho2**2 * (
        725 * d3 - 738 * rho2 * u4
    )
    if sp.expand(b_d3_u4_row - expected_d3_u4_row) != 0:
        raise AssertionError("branch-B depth-seven x5 row did not reduce as expected")
    b_after_d3 = {
        **b_after_c4,
        d3: sp.Rational(738, 725) * rho2 * u4,
    }
    b_later_rows = {}
    for row_name in ("depth6_x7", "depth7_x6", "depth7_x7"):
        value = sp.factor_terms(row_exprs[row_name].subs(b_after_d3))
        if value != 0:
            b_later_rows[row_name] = str(value)

    receipt = {
        "source_log": LOG.name,
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "leaf": branch,
        "exact_rows": {
            "z6_commutator": z6_text,
            "K5_x7": k5_text,
            "depth7_x3": x3_text,
            "depth6_x5": x5_text,
        },
        "reductions": {
            "z6_row_forces": "sigma0=0 (mu is a nonzero field element)",
            "K5_after_sigma0": str(k5_after_sigma0),
            "depth7_x3_after_sigma0_c2zero": str(x3_c2_zero),
            "depth6_x5_after_sigma0_c2_u2zero": str(x5_c2_u2_zero),
        },
        "covering_split": [
            "rho2=0",
            "rho2!=0, c2=0, u2=0, c3=0",
        ],
        "remaining_rows_after_split": branch_rows,
        "rho2_zero_refinement": {
            "first_row": str(sp.factor_terms(row_exprs["depth5_x6"].subs(rho2_zero_base))),
            "covering_cases": [
                "c2!=0, sigma1=0, tau0=0, sigma2=(-444*a4^3+1665*a4*rho3)/2775",
                "c2=0, sigma1=0, tau0=0",
                "c2=0, sigma1=0, tau0!=0, u2=0, c3=0",
                "c2=0, sigma1!=0, u2=0, c3=0, d2=0, u3=0, d3 given by the exact linear row",
            ],
            "case_residual_rows": {
                "c2_nonzero_case": residual_rows(a1_subs),
                "c2_zero_tau0_zero_case": residual_rows(a2_tau0_zero),
                "c2_zero_tau0_nonzero_case": residual_rows(a2_tau0_nonzero),
                "sigma1_nonzero_case": a3_final_rows,
            },
            "sigma1_nonzero_d3_relation": str(a3_d3_solutions[0][syms["d3"]]),
            "sigma1_nonzero_u3_row": str(a3_u3_row),
            "scope": "exact necessary-equation refinement on rho2=0; cases remain and this does not exclude the full Magnus or F2 branches",
        },
        "rho2_nonzero_refinement": {
            "conditions": "sigma0=c2=u2=c3=0, rho2!=0",
            "normalized_linear_rows": [str(b1), str(b2)],
            "linear_coefficient_matrix": [[str(x) for x in row] for row in b_matrix.tolist()],
            "determinant": str(b_det),
            "forced": ["d2=0", "u3=0", "c4=0"],
            "next_exact_row": str(b_d3_u4_row),
            "next_relation": "d3=(738/725)*rho2*u4",
            "later_rows_after_relation": b_later_rows,
            "scope": "necessary constraints on the explicit rho2!=0 subbranch; not an exclusion",
        },
        "derivation": [
            "z6 row and mu!=0 imply sigma0=0.",
            "K5 x7 then gives c2*rho2=0.",
            "If rho2=0, retain that branch.",
            "If rho2!=0, K5 gives c2=0; depth7 x3 gives u2=0; depth6 x5 gives c3=0.",
        ],
        "scope": "exact necessary-equation split on this selected Magnus leaf; not an exclusion of all F2 cases or of either full 2+2/2+1+1 branch",
        "checks": "exact SymPy rational identities; no parameter division except the stated nonzero field coefficient mu and the explicit rho2!=0 case",
    }
    OUT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "receipt": OUT.name,
        "source_sha256": receipt["source_sha256"],
        "covering_split": receipt["covering_split"],
        "scope": receipt["scope"],
    }, indent=2))


if __name__ == "__main__":
    main()
