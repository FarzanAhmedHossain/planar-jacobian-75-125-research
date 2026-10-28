"""Exact alpha-local q=4,s=2 ideal from the depth-6 2+2 checkpoint.

The script specializes U1=h^2*S to U1=0 by p1=p2=0, keeps the t=1
q>=5/zero-U2 chart separate, imposes the alpha-local order equations,
the K12/K14 selector equations, and the division-free K15 residue. It then
optionally tests saturation by the two nonzero leading coefficients d*a.
No primary checkpoint is modified; no full-branch exclusion is asserted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from math import comb
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.rings import ring


ROOT = Path(__file__).resolve().parent
STATE = ROOT / "f2_75_125_sparse_ring_lift_22_terminal_roots6_q_replayed_state.json"
OUT = ROOT / "f2_75_125_depth6_alpha_q4s2_aux_ideal_agent_20261018.json"
P_RECEIPT = ROOT / "f2_75_125_sparse_terminal_q_22_terminal_roots6_verified.json"
Q_RECEIPT = ROOT / "f2_75_125_sparse_ring_lift_22_terminal_roots6_q_replayed.json"
SQRT5 = sp.sqrt(5)
ALPHA = (3 + SQRT5) / 2


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--groebner", action="store_true")
    parser.add_argument("--order", choices=("grevlex", "grlex", "lex"), default="grevlex")
    args = parser.parse_args()

    data = json.loads(STATE.read_text(encoding="utf-8"))
    state_sha256 = hashlib.sha256(STATE.read_bytes()).hexdigest()
    p_receipt = json.loads(P_RECEIPT.read_text(encoding="utf-8"))
    q_receipt = json.loads(Q_RECEIPT.read_text(encoding="utf-8"))
    if (p_receipt["branch"] != "2+2" or p_receipt["depth"] != 6
            or p_receipt["P_terminal_positions_checked_through5"] != 32
            or p_receipt["nonzero_positions"] != 0):
        raise AssertionError("the exact 32-row terminal-P receipt did not verify")
    terminal_p_slice_counts = [6, 6, 6, 5, 5, 4]
    if sum(terminal_p_slice_counts) != p_receipt["P_terminal_positions_checked_through5"]:
        raise AssertionError("terminal-P slice count does not match the verified total")
    if (q_receipt["output_sha256"] != state_sha256
            or q_receipt["terminal_Q_positions_checked"] != 64
            or q_receipt["terminal_Q_equations_added"] != 0):
        raise AssertionError("the input checkpoint does not match the exact Q replay receipt")
    field = QQ.algebraic_field(SQRT5)
    parameter_names = [f"p{i}" for i in range(56)]
    names = [*parameter_names, "sat_qs"]
    P, *gens = ring(names, field)
    old_names = data["generators"]
    new_index = {name: i for i, name in enumerate(names)}
    tau = field.from_sympy(ALPHA)
    tau_sym = ALPHA
    tau_poly = P.ground_new(tau)
    sqrt_field = field.from_sympy(SQRT5)

    def coeff(values):
        expr = sum((sp.Rational(v) * SQRT5 ** (len(values) - 1 - i)
                    for i, v in enumerate(values)), sp.S.Zero)
        return field.from_sympy(expr)

    def jet(terms, root, order=0, u1_zero=True):
        root_elem = field.from_sympy(sp.sympify(root))
        result = {}
        for sparse, values in terms:
            degree = 0
            mon = [0] * len(gens)
            discard = False
            for old_i, exponent in sparse:
                name = old_names[old_i]
                if name == "t":
                    degree = exponent
                else:
                    i = new_index[name]
                    if u1_zero and name in ("p1", "p2") and exponent:
                        discard = True
                        break
                    mon[i] = exponent
            if discard or degree < order:
                continue
            value = coeff(values) * field.convert(comb(degree, order)) * root_elem ** (degree - order)
            key = tuple(mon)
            result[key] = result.get(key, field.zero) + value
        return P.from_dict({m: c for m, c in result.items() if c})

    def constant(value):
        return P.ground_new(field.convert(value))

    def terms(poly):
        return len(poly)

    def degree(poly):
        return max((sum(m) for m in poly), default=-1)

    def variables(poly):
        return sorted({names[i] for m in poly for i, e in enumerate(m) if e})

    def signature(poly):
        digest = hashlib.sha256()
        for mon, value in sorted(poly.items()):
            digest.update(repr((mon, str(field.to_sympy(value)))).encode())
        return digest.hexdigest()

    def info(poly):
        return {"terms": terms(poly), "degree": degree(poly),
                "variables": variables(poly), "zero": not bool(poly),
                "sha256": signature(poly)}

    # Taylor coefficients at alpha; all equations are polynomial in p over K.
    hj = [jet(data["h"], tau_sym, j) for j in range(4)]
    c = next(iter(hj[2].values()))
    c1 = next(iter(hj[3].values())) if hj[3] else field.zero
    r0prime = tau - field.convert(1)
    U = {k: [jet(data["uk"][str(k)], tau_sym, j)
             for j in range(0, 6)] for k in range(2, 7)}

    # Exact root-order requirements for the selected alpha stratum:
    # ord(U2)=4, ord(U3)>=4, ord(U4)=2, ord(U5/t)>=2.
    order_requirements = []
    for j in range(4):
        order_requirements.append((f"U2_alpha_u{j}", U[2][j]))
    for j in range(4):
        order_requirements.append((f"U3_alpha_u{j}", U[3][j]))
    for j in range(2):
        order_requirements.append((f"U4_alpha_u{j}", U[4][j]))
    for j in range(2):
        order_requirements.append((f"U5_alpha_u{j}", U[5][j]))
    order_equations = [(tag, p) for tag, p in order_requirements if p]

    # Nonzero leading coefficients d (q=4) and a (s=2) are not divided by.
    d, a = U[2][4], U[4][2]
    f = U[6][0] / constant(tau)
    selector12 = d**2 - constant(3*c**3)*a
    selector14 = d**3 - constant(27*c**6)*f

    # Global quotient values at alpha, computed from Taylor jets. Since the
    # q=4/s=2 order equations hold and the t=1 orders are already deeper,
    # these are the values of B3=U3/h^2, L=(3C4-B2^2)/r0, T=U5/r0^2,
    # N=(27U6-t*B2^3)/r0. All divisions are by nonzero K scalars only.
    b2 = d / constant(c**2)
    b2_1 = (U[2][5] - constant(2*c*c1)*b2) / constant(c**2)
    b3 = U[3][4] / constant(c**2)
    c4 = a / constant(c)
    c4_1 = (U[4][3] - constant(c1)*c4) / constant(c)
    L = (3*c4_1 - 2*b2*b2_1) / constant(r0prime)
    T = U[5][2] / constant(c)
    N = (27*U[6][1] - b2**3 - constant(3*tau)*b2**2*b2_1) / constant(r0prime)
    k15 = tau_poly*b3*N**2 + constant(18*tau)*L*T*N + constant(243)*T**3

    equations = [p for _, p in order_equations]
    equations += [selector12, selector14, k15]
    tags = [tag for tag, _ in order_equations] + [
        "K12_selector_d2_minus_3c3a",
        "K14_selector_d3_minus_27c6f",
        "K15_division_free_residual",
    ]
    saturation = gens[-1]*d*a - 1

    # The only degree-one local order rows form a homogeneous linear system
    # over K=Q(sqrt(5)). Row-reduce it exactly before any Groebner attempt.
    # This eliminates only variables using nonzero field constants, never a
    # parameter or a leading coefficient. The original order rows remain in
    # `equations`; `reduced_exprs` is an equivalent presentation of their
    # ideal after solving these linear rows.
    linear_rows = [p for p in equations if degree(p) == 1]
    linear_vars = sorted(
        {sp.Symbol(name) for p in linear_rows for name in variables(p)},
        key=lambda s: names.index(str(s)),
    )
    linear_reduction = {"linear_row_count": len(linear_rows),
                        "linear_variable_count": len(linear_vars)}
    substitutions = {}
    if linear_rows:
        matrix, _ = sp.linear_eq_to_matrix([p.as_expr() for p in linear_rows], linear_vars)
        reduced_matrix, pivots = matrix.rref()
        for row_index, pivot_index in enumerate(pivots):
            substitutions[linear_vars[pivot_index]] = sp.expand(-sum(
                reduced_matrix[row_index, j] * linear_vars[j]
                for j in range(len(linear_vars)) if j != pivot_index
            ))
        linear_reduction.update({
            "rank": len(pivots),
            "pivot_variables": [str(linear_vars[i]) for i in pivots],
            "substitutions": {str(k): str(v) for k, v in substitutions.items()},
        })
    replacement_polys = {}
    for variable, expr in substitutions.items():
        linear = sp.Poly(expr, *linear_vars, domain=field)
        replacement = P.zero
        for mon, coefficient in linear.terms():
            if sum(mon) != 1:
                raise AssertionError(("unexpected nonlinear linear-row solution", expr))
            variable_index = mon.index(1)
            replacement += P.ground_new(field.from_sympy(coefficient)) * gens[
                new_index[str(linear_vars[variable_index])]
            ]
        replacement_polys[gens[new_index[str(variable)]]] = replacement
    reduced_equations = [p.compose(replacement_polys) for p in [*equations, saturation]]
    reduced_tags = [*tags, "saturation_sat_qs_d_a_minus_1"]
    reduced_equation_records = [
        {"tag": tag, "polynomial": poly}
        for tag, poly in zip(reduced_tags, reduced_equations) if poly
    ]
    reduced_generator_set_sha256 = hashlib.sha256("\n".join(
        f"{record['tag']}:{signature(record['polynomial'])}"
        for record in reduced_equation_records
    ).encode()).hexdigest()
    reduced_equations = [record["polynomial"] for record in reduced_equation_records]
    # Introduce exact quotient variables so K15 remains compact. The four
    # definitions are polynomial identities; no quotient variable is assumed
    # nonzero, and no division by B3, L, T, or N occurs.
    aux_names = ["aux_B3", "aux_L", "aux_T", "aux_N"]
    ext_names = [*aux_names, *names]
    R, *ext_gens = ring(ext_names, field)
    aux_b3, aux_l, aux_t, aux_n = ext_gens[:4]

    def lift(poly):
        return R.from_dict({(0, 0, 0, 0, *mon): coefficient
                            for mon, coefficient in poly.items()})

    def variables_ext(poly):
        return sorted({ext_names[i] for mon in poly
                       for i, exponent in enumerate(mon) if exponent})

    base_records = [record for record in reduced_equation_records
                    if record["tag"] != "K15_division_free_residual"]
    quotient_definitions = []
    for aux_name, aux_gen, quotient in zip(aux_names, ext_gens[:4], [b3, L, T, N]):
        quotient_lifted = lift(quotient.compose(replacement_polys))
        quotient_definitions.append({
            "tag": f"define_{aux_name}",
            "polynomial": aux_gen - quotient_lifted,
        })
    compact_k15 = (
        R.ground_new(tau) * aux_b3 * aux_n**2
        + R.ground_new(field.convert(18) * tau) * aux_l * aux_t * aux_n
        + R.ground_new(243) * aux_t**3
    )
    aux_equation_records = [
        {"tag": record["tag"], "polynomial": lift(record["polynomial"])}
        for record in base_records
    ] + quotient_definitions + [{
        "tag": "K15_compact_division_free",
        "polynomial": compact_k15,
    }]
    aux_equations = [record["polynomial"] for record in aux_equation_records]
    reduced_active_symbols = sorted(
        {sp.Symbol(name) for poly in aux_equations for name in variables_ext(poly)},
        key=lambda symbol: ext_names.index(str(symbol)),
    )
    reduced_generator_set_sha256 = hashlib.sha256("\n".join(
        f"{record['tag']}:{signature(record['polynomial'])}"
        for record in aux_equation_records
    ).encode()).hexdigest()
    linear_reduction.update({
        "post_substitution_variable_count": len(reduced_active_symbols),
        "post_substitution_variables": [str(symbol) for symbol in reduced_active_symbols],
        "post_substitution_equation_count": len(aux_equations),
        "post_substitution_generator_set_sha256": reduced_generator_set_sha256,
        "post_substitution_term_count": sum(len(poly) for poly in aux_equations),
        "post_substitution_max_polynomial_terms": max((len(poly) for poly in aux_equations), default=0),
        "post_substitution_generators": [
            {"tag": record["tag"], "terms": len(record["polynomial"]),
             "degree": degree(record["polynomial"]),
             "sha256": signature(record["polynomial"]),
             "polynomial": str(record["polynomial"].as_expr())}
            for record in aux_equation_records
        ],
        "auxiliary_quotients": {
            "definitions": [record["tag"] for record in quotient_definitions],
            "K15": "alpha*aux_B3*aux_N^2 + 18*alpha*aux_L*aux_T*aux_N + 243*aux_T^3",
            "no_quotient_zero_case_saturated": True,
        },
    })

    # The t=1 branch stays separate: on p1=p2=0 its U2 coefficients through
    # u^4 vanish, so it is q>=5 or U2=0, not the exact-q=4 alpha branch.
    at_one_u2 = [jet(data["uk"]["2"], sp.Integer(1), j) for j in range(5)]
    at_one = {f"U2_t1_u{j}": info(p) for j, p in enumerate(at_one_u2)}
    at_one.update({
        "retained_separate_from_alpha_q4_s2_ideal": True,
        "condition": "U2=0 or ord_1(U2)>=5",
        "U3_t1_u0_to_u3_zero": all(not jet(data["uk"]["3"], sp.Integer(1), j) for j in range(4)),
        "U4_t1_u0_to_u1_zero": all(not jet(data["uk"]["4"], sp.Integer(1), j) for j in range(2)),
        "U5_t1_u0_to_u1_zero": all(not jet(data["uk"]["5"], sp.Integer(1), j) for j in range(2)),
    })

    receipt = {
        "input_state": str(STATE.relative_to(ROOT.parent)),
        "input_state_sha256": state_sha256,
        "input_q_replay_receipt": str(Q_RECEIPT.relative_to(ROOT.parent)),
        "input_terminal_p_receipt": str(P_RECEIPT.relative_to(ROOT.parent)),
        "terminal_p_rows_through_distance_5": p_receipt["P_terminal_positions_checked_through5"],
        "terminal_p_slice_counts_through_distance_5": terminal_p_slice_counts,
        "terminal_p_nonzero_positions": p_receipt["nonzero_positions"],
        "terminal_q_positions_checked": q_receipt["terminal_Q_positions_checked"],
        "terminal_q_added_equations": q_receipt["terminal_Q_equations_added"],
        "parameters": 56,
        "saved_constraint_count": len(data.get("constraints", [])),
        "branch_specialization": "p1=p2=0, equivalently U1=h^2*S=0",
        "terminal_P_rows": "depth-6 checkpoint's resume guard covers all 32 forbidden positions through distance 5",
        "t1_separate_branch": at_one,
        "alpha_order_equations": [{"tag": tag, **info(poly)} for tag, poly in order_requirements],
        "alpha_K12_K14_premise": (
            "Use the global K12/K14 first-integral identities at alpha only on the selected "
            "rootwise profile ord(h)=2, ord(U2)=4, ord(U4)=2, ord(U3)>=4, "
            "ord(U5/t)>=2, with c,d,a nonzero. The q=4,s=2 audit checks the lower-pole "
            "residues and that correction modes are higher order; no second terminal-chart selector is added."
        ),
        "alpha_added_order_equation_count": len(order_equations),
        "alpha_leading_coefficients": {"d_q4": info(d), "a_s2": info(a)},
        "alpha_selectors": {"K12": info(selector12), "K14": info(selector14)},
        "alpha_quotients": {"B3": info(b3), "L": info(L), "T": info(T), "N": info(N)},
        "alpha_K15": info(k15),
        "saturation_equation": "sat_qs*d*a - 1 = 0",
        "zero_case_coverage": {
            "nonzero_factors_only": ["d=U2_alpha_u4", "a=U4_alpha_u2"],
            "saturation": "d*a != 0 on the selected q=4,s=2 stratum",
            "U3_B3_zero_allowed": True,
            "U5_T_zero_allowed": True,
            "N_zero_allowed": True,
            "K15_kept_division_free": "alpha*B3*N^2 + 18*alpha*L*T*N + 243*T^3 = 0",
            "if_N_zero": "the same equation gives T=0; no division by N is used",
            "if_T_zero": "the same equation gives B3*N^2=0; both resulting zero branches remain",
        },
        "saturated_generator_count_including_zero_orders": len(aux_equations),
        "constant_linear_reduction": linear_reduction,
        "groebner": {"attempted": bool(args.groebner), "order": args.order if args.groebner else None},
        "scope": "alpha-local U1=0, q=4,s=2,r>=4,v>=2 branch only; no full 2+2 exclusion",
    }

    if args.groebner:
        # Exact Groebner computation on the equivalent ideal after the
        # constant-coefficient row reduction above.
        reduced_exprs = [poly.as_expr() for poly in aux_equations]
        receipt["groebner"].update({"active_variable_count": len(reduced_active_symbols),
                                    "active_variables": [str(x) for x in reduced_active_symbols],
                                    "input_polynomial_count": len(reduced_exprs),
                                    "input_term_count": sum(len(p) for p in aux_equations),
                                    "linear_rows_eliminated_by_field_constants": len(pivots) if linear_rows else 0})
        G = sp.groebner(reduced_exprs, *reduced_active_symbols, order=args.order, domain=field, method="f5b")
        is_unit = len(G.polys) == 1 and G.polys[0].as_expr() == 1
        receipt["groebner"].update({"completed": True, "is_unit": is_unit,
                                    "basis_size": len(G.polys),
                                    "basis_terms": sum(len(poly.terms()) for poly in G.polys),
                                    "basis_sha256": hashlib.sha256(
                                        "\n".join(str(poly.as_expr()) for poly in G.polys).encode()
                                    ).hexdigest()})

    OUT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == "__main__":
    main()
