"""Test the alpha-local witness against the licensed t=1 terminal-P chart.

The saved depth-six state is decoded directly after exact specialization, so
all terminal chart computations are univariate over Q(sqrt(5)).
"""

from __future__ import annotations

import hashlib
import json
from math import comb, factorial
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parent
STATE = ROOT / "f2_75_125_sparse_ring_lift_22_terminal_roots6_q_replayed_state.json"
REPLAY = ROOT / "f2_75_125_sparse_ring_lift_22_terminal_roots6_q_replayed.json"
OUTPUT = ROOT / "f2_75_125_alpha_witness_terminal_P_receipt_20261018.json"
ZERO = sp.Integer(0)


def witness_values():
    return {
        "p1": sp.Integer(1), "p2": ZERO,
        "p7": sp.Rational(11108345, 43666864)
        - sp.Rational(9423059, 131000592) * sp.sqrt(5),
        "p8": ZERO, "p9": ZERO,
        "p16": -sp.Rational(1797920925871, 24427418389056)
        + sp.Rational(3400480314037, 122137091945280) * sp.sqrt(5),
        "p17": ZERO, "p18": ZERO, "p19": ZERO, "p20": ZERO,
        "p27": sp.Rational(3748511645183183, 469071572852246016)
        - sp.Rational(78142381152123083, 21108220778351070720) * sp.sqrt(5),
        "p28": ZERO, "p29": ZERO, "p30": ZERO, "p31": ZERO,
        "p32": ZERO, "p33": ZERO,
    }


def decode_specialized(terms, generator_names, values):
    """Return a univariate t-polynomial as {degree: exact coefficient}."""
    out = {}
    for sparse, encoded_coeff in terms:
        if len(encoded_coeff) == 2:
            coeff = sp.Rational(encoded_coeff[0]) * sp.sqrt(5) + sp.Rational(encoded_coeff[1])
        elif len(encoded_coeff) == 1:
            coeff = sp.Rational(encoded_coeff[0])
        else:
            raise ValueError("unexpected coefficient encoding")
        tdegree = 0
        for gen_index, exponent in sparse:
            name = generator_names[gen_index]
            if name == "t":
                tdegree += exponent
            elif name.startswith("p"):
                coeff *= values.get(name, ZERO) ** exponent
            elif name.startswith("x"):
                # No unsolved x-coordinate is assigned a nonzero value.
                coeff = ZERO
            else:
                raise ValueError(f"unexpected checkpoint generator {name}")
        if coeff:
            out[tdegree] = sp.simplify(out.get(tdegree, ZERO) + coeff)
    return {degree: coeff for degree, coeff in out.items() if coeff != 0}


def add_scaled(out, source, scale):
    for degree, coeff in source.items():
        out[degree] = sp.simplify(out.get(degree, ZERO) + scale * coeff)


def local_slice(p_poly, distance, bound):
    """The checkpoint chart X^-15 P(X^5,X^-1(1+u)), truncated in u."""
    out = {}
    for i, coeff in p_poly.items():
        exponent = 5 * i - 15 + distance
        if exponent < 0:
            raise AssertionError("source monomial fell outside chart support")
        for degree in range(min(exponent, bound) + 1):
            out[degree] = sp.simplify(
                out.get(degree, ZERO) + coeff * comb(exponent, degree)
            )
    return {degree: coeff for degree, coeff in out.items() if coeff != 0}


def derivative(poly, order):
    out = dict(poly)
    for _ in range(order):
        out = {degree - 1: degree * coeff
               for degree, coeff in out.items() if degree > 0}
    return {degree: sp.simplify(coeff) for degree, coeff in out.items() if coeff != 0}


def centered_slice(p_list, distance, bound, center_a, center_b):
    out = {}
    for i in range(distance + 1):
        deficit = distance - i
        for n in range(deficit // 2 + 1):
            m = deficit - 2 * n
            chart = local_slice(p_list[i], i, bound + m + n)
            chart = derivative(chart, m + n)
            scale = center_a ** m * center_b ** n / (factorial(m) * factorial(n))
            add_scaled(out, chart, scale)
    return {degree: sp.simplify(coeff) for degree, coeff in out.items()
            if degree <= bound and sp.simplify(coeff) != 0}


def poly_at_one(poly):
    return sp.simplify(sum(poly.values(), ZERO))


def main():
    state = json.loads(STATE.read_text(encoding="utf-8"))
    replay = json.loads(REPLAY.read_text(encoding="utf-8"))
    values = witness_values()
    p_list = [decode_specialized(terms, state["generators"], values)
              for terms in state["p"]]
    S = decode_specialized(state["S"], state["generators"], values)

    # Terminal chart initialization and recentering at t=1, exactly as in the
    # saved lift: c0=[u^6]P0, a=-[u^5]P1/(6c0), then determine b from P2.
    local_c = local_slice(p_list[0], 0, 6)
    c0 = sp.simplify(local_c.get(6, ZERO))
    a = sp.simplify(-local_slice(p_list[1], 1, 5).get(5, ZERO) / (6 * c0))
    b0_slice = centered_slice(p_list, 2, 5, a, ZERO)
    b = sp.simplify(-b0_slice.get(5, ZERO) / (6 * c0))

    init0 = centered_slice(p_list, 0, 5, a, b)
    init1 = centered_slice(p_list, 1, 5, a, b)
    terminal_rows = []
    for distance in range(len(p_list)):
        bound = (71 - 5 * distance) // 12
        if bound < 0:
            continue
        residual = centered_slice(p_list, distance, bound, a, b)
        terminal_rows.append({
            "distance": distance,
            "bound": bound,
            "positions_checked": bound + 1,
            "nonzero_residuals": [
                {"u_degree": degree, "value": str(sp.simplify(coeff))}
                for degree, coeff in sorted(residual.items())
            ],
            "zero": not bool(residual),
        })

    S_one = poly_at_one(S)
    assert c0 != 0
    assert S_one == 0
    assert not init0 and not init1
    assert all(row["zero"] for row in terminal_rows)
    assert sum(row["positions_checked"] for row in terminal_rows if row["distance"] <= 5) == 32
    assert sum(row["positions_checked"] for row in terminal_rows) == 36

    receipt = {
        "source_state": STATE.name,
        "source_state_sha256": hashlib.sha256(STATE.read_bytes()).hexdigest(),
        "checkpoint_depth": len(p_list) - 1,
        "root_condition": "S(1)=0",
        "root_condition_residual": str(S_one),
        "terminal_chart_initialization": {
            "P0_centered_mod_u6_zero": not bool(init0),
            "P1_centered_mod_u6_zero": not bool(init1),
            "local_leading_coefficient_c0": str(c0),
            "center_a": str(a),
            "center_b": str(b),
        },
        "terminal_P_rows": terminal_rows,
        "original_32_rows_zero": all(row["zero"] for row in terminal_rows
                                      if row["distance"] <= 5),
        "additional_depth6_rows_zero": all(row["zero"] for row in terminal_rows
                                            if row["distance"] == 6),
        "Q_replay_coverage": {
            "saved_corrected_power_relations_recovered": all(
                row["saved_Q_slice_recovered"]
                for row in replay["corrected_power_replay"]
            ),
            "terminal_Q_positions_checked": replay["terminal_Q_positions_checked"],
            "retained_constraint_count": replay["retained_constraint_count"],
        },
        "result": "the alpha-local witness survives every checked licensed t=1 terminal-P chart/root condition",
        "scope": "exact specialization of the depth-six checkpoint; no alpha terminal chart and no global K6-K8 ideal claim",
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
