"""Replay an exact depth-six witness in the rootwise F2 terminal chart at t=alpha.

The rootwise licensing argument is documented in F2-75-125-band-model.md.
This probe does not exclude a branch; it checks whether one saved local witness
survives the additional alpha-centred terminal support conditions.
"""

from __future__ import annotations

import hashlib
import json
from math import comb, factorial
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ


ROOT = Path(__file__).resolve().parent
STATE = ROOT / "f2_75_125_sparse_ring_lift_22_terminal_roots6_q_replayed_state.json"
ROOTWISE = ROOT / "f2_75_125_alpha_tuned_rootwise_receipt_20261018.json"
OUTPUT = ROOT / "f2_75_125_alpha_terminal_chart_probe_20261018.json"
ZERO = QQ.zero


def witness_values(theta_field):
    sqrt5 = 2 * theta_field**5 - 3
    return {
        1: QQ.one,
        7: (theta_field**0) * QQ.from_sympy(sp.Rational(11108345, 43666864))
           - sqrt5 * QQ.from_sympy(sp.Rational(9423059, 131000592)),
        16: -QQ.from_sympy(sp.Rational(1797920925871, 24427418389056))
            + sqrt5 * QQ.from_sympy(sp.Rational(3400480314037, 122137091945280)),
        27: QQ.from_sympy(sp.Rational(3748511645183183, 469071572852246016))
            - sqrt5 * QQ.from_sympy(sp.Rational(78142381152123083, 21108220778351070720)),
    }


def decode_coefficient(encoded, field, sqrt5):
    if len(encoded) == 1:
        return field.from_sympy(sp.Rational(encoded[0]))
    if len(encoded) == 2:
        return (field.from_sympy(sp.Rational(encoded[0])) * sqrt5
                + field.from_sympy(sp.Rational(encoded[1])))
    raise ValueError(f"unexpected coefficient encoding: {encoded}")


def specialized_polynomial(terms, generators, values, field, sqrt5):
    """Decode a saved sparse polynomial and evaluate its parameter variables."""
    result = {}
    for sparse, encoded in terms:
        coefficient = decode_coefficient(encoded, field, sqrt5)
        tdegree = 0
        for index, exponent in sparse:
            name = generators[index]
            if name == "t":
                tdegree += exponent
            elif name.startswith("p"):
                pindex = int(name[1:])
                coefficient *= values.get(pindex, ZERO) ** exponent
            elif name.startswith("x"):
                # This witness assigns every uneliminated affine coordinate 0.
                coefficient = ZERO
            else:
                raise ValueError(f"unexpected generator {name}")
        if coefficient:
            result[tdegree] = result.get(tdegree, ZERO) + coefficient
    return {degree: value for degree, value in result.items() if value}


def local_slice(poly, distance, bound, root, base_degree, field, q_support=False):
    """Expand a band in y=X^-1(root+u), keeping u^0..u^bound."""
    out = [ZERO for _ in range(bound + 1)]
    if q_support:
        low = -(-(25 - distance) // 5)
        high = min(25, (125 - distance) // 5)
    for tdegree, coefficient in poly.items():
        if q_support and not (low <= tdegree <= high):
            continue
        exponent = 5 * tdegree - base_degree + distance
        if exponent < 0:
            raise AssertionError("source monomial fell outside the local chart")
        for degree in range(min(exponent, bound) + 1):
            out[degree] += coefficient * field.from_sympy(sp.Integer(comb(exponent, degree))) \
                           * root ** (exponent - degree)
    return out


def derivative_coefficients(poly, order, bound, field):
    out = [ZERO for _ in range(bound + 1)]
    for degree in range(bound + 1):
        source_degree = degree + order
        if source_degree < len(poly):
            falling = factorial(source_degree) // factorial(degree)
            out[degree] = poly[source_degree] * field.from_sympy(sp.Integer(falling))
    return out


def centered_slice(poly_list, distance, bound, center_a, center_b, root,
                   base_degree, field, q_support=False):
    out = [ZERO for _ in range(bound + 1)]
    for i in range(distance + 1):
        deficit = distance - i
        for n in range(deficit // 2 + 1):
            m = deficit - 2 * n
            local = local_slice(poly_list[i], i, bound + m + n, root,
                                base_degree, field, q_support=q_support)
            derivative = derivative_coefficients(local, m + n, bound, field)
            scale = (center_a**m * center_b**n
                     / field.from_sympy(sp.Integer(factorial(m) * factorial(n))))
            for degree in range(bound + 1):
                out[degree] += derivative[degree] * scale
    return [value for value in out]


def nonzero_terms(poly):
    return [(degree, value) for degree, value in enumerate(poly) if value]


def main():
    state = json.loads(STATE.read_text(encoding="utf-8"))
    rootwise = json.loads(ROOTWISE.read_text(encoding="utf-8"))
    state_hash = hashlib.sha256(STATE.read_bytes()).hexdigest()
    if rootwise["source_state_sha256"] != state_hash:
        raise AssertionError("alpha-local receipt does not match the selected checkpoint")

    x = sp.Symbol("x")
    minpoly = x**10 - 3*x**5 + 1
    theta_expr = sp.CRootOf(minpoly, 1)  # positive real fifth root of alpha
    field = QQ.algebraic_field(theta_expr)
    theta = field.from_sympy(theta_expr)
    alpha = theta**5
    sqrt5 = 2 * alpha - 3
    values = witness_values(theta)

    p_list = [specialized_polynomial(poly, state["generators"], values,
                                     field, sqrt5) for poly in state["p"]]
    q_list = [specialized_polynomial(poly, state["generators"], values,
                                     field, sqrt5) for poly in state["q"]]
    S = specialized_polynomial(state["S"], state["generators"], values,
                               field, sqrt5)
    S_alpha = sum((value * alpha**degree for degree, value in S.items()), ZERO)

    c0 = local_slice(p_list[0], 0, 6, theta, 15, field)[6]
    first_band = local_slice(p_list[1], 1, 5, theta, 15, field)
    center_a = -first_band[5] / (field.from_sympy(sp.Integer(6)) * c0)
    distance2_at_a = centered_slice(p_list, 2, 5, center_a, ZERO, theta,
                                    15, field)
    center_b = -distance2_at_a[5] / (field.from_sympy(sp.Integer(6)) * c0)

    p_rows = []
    for distance in range(len(p_list)):
        bound = (71 - 5 * distance) // 12
        if bound < 0:
            continue
        row = centered_slice(p_list, distance, bound, center_a, center_b,
                             theta, 15, field)
        terms = nonzero_terms(row)
        p_rows.append({
            "distance": distance, "bound": bound,
            "nonzero_count": len(terms),
            "first_nonzero": [
                {"u_degree": degree, "value": str(field.to_sympy(value))}
                for degree, value in terms[:8]
            ],
        })

    q_rows = []
    for distance in range(len(q_list)):
        bound = (119 - 5 * distance) // 12
        if bound < 0:
            continue
        row = centered_slice(q_list, distance, bound, center_a, center_b,
                             theta, 25, field, q_support=True)
        terms = nonzero_terms(row)
        q_rows.append({
            "distance": distance, "bound": bound,
            "nonzero_count": len(terms),
            "first_nonzero": [
                {"u_degree": degree, "value": str(field.to_sympy(value))}
                for degree, value in terms[:8]
            ],
        })

    receipt = {
        "source_state": STATE.name,
        "source_state_sha256": state_hash,
        "rootwise_local_receipt": ROOTWISE.name,
        "root_parameter": "theta=CRootOf(x^10-3*x^5+1,1), theta^5=alpha=(3+sqrt(5))/2",
        "first_shear_cutoff": "15/4",
        "root_multiplicity": 6,
        "successor_corner": "(7/5,2)",
        "specialization": rootwise["exact_specialization"],
        "S_alpha": str(field.to_sympy(S_alpha)),
        "alpha_chart_initialization": {
            "leading_u6_coefficient": str(field.to_sympy(c0)),
            "center_a": str(field.to_sympy(center_a)),
            "center_b": str(field.to_sympy(center_b)),
            "S_alpha_zero": not bool(S_alpha),
        },
        "terminal_P_rows": p_rows,
        "terminal_Q_rows": q_rows,
        "conclusion": (
            "This is an exact pointwise replay of one depth-six local witness. "
            "It does not prove that the alpha-terminal equations eliminate the full branch."
        ),
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
