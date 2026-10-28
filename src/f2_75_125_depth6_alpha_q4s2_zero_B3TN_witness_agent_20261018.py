"""Probe the exact B3=T=N=0 subcase of the alpha-local q=4,s=2 ideal.

This is an admissible point for the selected local necessary ideal only. It is
not a completion of the global 2+2 branch. The script rechecks the saved depth-6
and terminal-P/Q guards, solves exact linear blocks over Q(sqrt(5)), and checks
all selected generators plus d*a saturation at the resulting point.
"""
from __future__ import annotations

import hashlib
import json
from math import comb, factorial
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
STATE = HERE / "f2_75_125_sparse_ring_lift_22_terminal_roots6_q_replayed_state.json"
P_RECEIPT = HERE / "f2_75_125_sparse_terminal_q_22_terminal_roots6_verified.json"
Q_RECEIPT = HERE / "f2_75_125_sparse_ring_lift_22_terminal_roots6_q_replayed.json"
AUX_RECEIPT = HERE / "f2_75_125_depth6_alpha_q4s2_aux_ideal_agent_20261018.json"
OUT = HERE / "f2_75_125_depth6_alpha_q4s2_zero_B3TN_witness_agent_20261018.json"


def exact_expr(s: str, syms: dict[str, sp.Symbol]) -> sp.Expr:
    return sp.sympify(s, locals={**syms, "sqrt": sp.sqrt})


def solve_linear_block(equations: list[sp.Expr], variables: list[sp.Symbol]):
    matrix, rhs = sp.linear_eq_to_matrix(equations, variables)
    reduced, pivots = matrix.row_join(rhs).rref()
    n = len(variables)
    if n in pivots:
        raise ValueError("inconsistent linear block")
    for row in range(reduced.rows):
        if all(reduced[row, col] == 0 for col in range(n)) and reduced[row, n] != 0:
            raise ValueError("inconsistent linear block")
    solution = {v: sp.S.Zero for v in variables}
    for row, pivot in enumerate(pivots):
        if pivot < n:
            solution[variables[pivot]] = sp.simplify(reduced[row, n])
    return solution, {"rank": sum(p < n for p in pivots),
                      "pivots": [str(variables[p]) for p in pivots if p < n],
                      "free_set_to_zero": [str(variables[i]) for i in range(n)
                                           if i not in [p for p in pivots if p < n]]}


def decode_specialized(terms, generator_names, values):
    """Evaluate a saved sparse checkpoint polynomial to a univariate t map."""
    out: dict[int, sp.Expr] = {}
    for sparse, encoded in terms:
        coefficient = sum((sp.Rational(encoded[i]) * sp.sqrt(5) ** (len(encoded)-1-i)
                           for i in range(len(encoded))), sp.S.Zero)
        tdegree = 0
        for generator_index, exponent in sparse:
            name = generator_names[generator_index]
            if name == "t":
                tdegree += exponent
            else:
                coefficient *= values.get(name, sp.Integer(0)) ** exponent
        if coefficient:
            out[tdegree] = sp.simplify(out.get(tdegree, 0) + coefficient)
    return {degree: sp.simplify(value) for degree, value in out.items() if value != 0}


def local_slice(poly, lead: int, distance: int, bound: int):
    out: dict[int, sp.Expr] = {}
    for tdegree, coefficient in poly.items():
        exponent = 5 * tdegree - lead + distance
        if exponent < 0:
            raise AssertionError("checkpoint monomial lies outside terminal chart support")
        for degree in range(min(exponent, bound) + 1):
            out[degree] = sp.simplify(out.get(degree, 0)
                                      + coefficient * comb(exponent, degree))
    return {degree: value for degree, value in out.items() if value != 0}


def centered_slice(series, distance: int, bound: int, center_a, center_b, lead: int = 15):
    out: dict[int, sp.Expr] = {}
    for i in range(distance + 1):
        deficit = distance - i
        for n in range(deficit // 2 + 1):
            m = deficit - 2 * n
            derivative = local_slice(series[i], lead, i, bound + m + n)
            for _ in range(m + n):
                derivative = {degree - 1: degree * value
                              for degree, value in derivative.items() if degree > 0}
            scale = center_a**m * center_b**n / factorial(m) / factorial(n)
            for degree, value in derivative.items():
                if degree <= bound:
                    out[degree] = sp.simplify(out.get(degree, 0) + scale * value)
    return {degree: sp.simplify(value) for degree, value in out.items()
            if degree <= bound and value != 0}


def replay_terminal_p(state_data, p_values):
    checkpoint_values = {f"p{i}": p_values.get(f"p{i}", sp.Integer(0))
                         for i in range(state_data["nparameters"])}
    # The auxiliary checkpoint coordinates x0,...,x72 are specialized to zero.
    specialized_p = [decode_specialized(row, state_data["generators"], checkpoint_values)
                     for row in state_data["p"]]
    c_poly = local_slice(specialized_p[0], 15, 0, 6)
    c0 = sp.simplify(c_poly.get(6, 0))
    if c0 == 0:
        raise ValueError("P terminal leading coefficient vanished at candidate point")
    a_poly = local_slice(specialized_p[1], 15, 1, 5)
    center_a = sp.simplify(-a_poly.get(5, 0) / (6 * c0))
    second_slice = centered_slice(specialized_p, 2, 5, center_a, 0)
    center_b = sp.simplify(-second_slice.get(5, 0) / (6 * c0))
    rows = []
    for distance in range(len(specialized_p)):
        bound = (71 - 5 * distance) // 12
        if bound < 0:
            continue
        residual = centered_slice(specialized_p, distance, bound, center_a, center_b)
        rows.append({"distance": distance, "bound": bound,
                     "positions_checked": bound + 1,
                     "nonzero": {str(k): str(v) for k, v in sorted(residual.items())},
                     "zero": not bool(residual)})
    specialized_q = [decode_specialized(row, state_data["generators"], checkpoint_values)
                     for row in state_data["q"]]
    q_rows = []
    for distance in range(len(specialized_q)):
        bound = (119 - 5 * distance) // 12
        if bound < 0:
            continue
        residual = centered_slice(specialized_q, distance, bound, center_a, center_b, lead=25)
        q_rows.append({"distance": distance, "bound": bound,
                       "positions_checked": bound + 1,
                       "nonzero": {str(k): str(v) for k, v in sorted(residual.items())},
                       "zero": not bool(residual)})
    return {"leading_coefficient": str(c0), "center_a": str(center_a),
            "center_b": str(center_b), "rows": rows,
            "through_distance_5_positions": sum(r["positions_checked"]
                                                  for r in rows if r["distance"] <= 5),
            "all_through_distance_5_zero": all(r["zero"] for r in rows
                                                if r["distance"] <= 5),
            "all_depth_6_rows_zero": all(r["zero"] for r in rows),
            "Q_rows": q_rows,
            "Q_positions_checked": sum(r["positions_checked"] for r in q_rows),
            "all_Q_rows_zero": all(r["zero"] for r in q_rows)}


def main() -> None:
    state_sha = hashlib.sha256(STATE.read_bytes()).hexdigest()
    p = json.loads(P_RECEIPT.read_text(encoding="utf-8"))
    q = json.loads(Q_RECEIPT.read_text(encoding="utf-8"))
    src = json.loads(AUX_RECEIPT.read_text(encoding="utf-8"))
    assert state_sha == src["input_state_sha256"]
    assert p["branch"] == "2+2" and p["depth"] == 6
    assert p["P_terminal_positions_checked_through5"] == 32
    assert p["nonzero_positions"] == 0
    assert sum([6, 6, 6, 5, 5, 4]) == 32
    assert q["output_sha256"] == state_sha
    assert q["terminal_Q_positions_checked"] == 64
    assert q["terminal_Q_equations_added"] == 0
    assert src["constant_linear_reduction"]["post_substitution_generator_set_sha256"] == \
        "7a4ac44de3f676304fc50e96e595d79b159a540e95c716176fa4884d2b4f167d"

    symbols = {name: sp.Symbol(name) for name in src["constant_linear_reduction"]["post_substitution_variables"]}
    gens = {row["tag"]: exact_expr(row["polynomial"], symbols)
            for row in src["constant_linear_reduction"]["post_substitution_generators"]}
    # An especially small affine slice: p19=p20=0 forces B3=0; p8=1,p9=0
    # is a deterministic trial, while T=N=0 makes K15 vanish without division.
    fixed = {symbols["p8"]: sp.Integer(1), symbols["p9"]: sp.Integer(0),
             symbols["p19"]: sp.Integer(0), symbols["p20"]: sp.Integer(0),
             symbols["aux_B3"]: sp.Integer(0), symbols["aux_T"]: sp.Integer(0),
             symbols["aux_N"]: sp.Integer(0)}

    # Solve U4-alpha orders and K12 together. The K12 relation makes a nonzero
    # automatically once d is nonzero, because c is a nonzero h-jet constant.
    block1_tags = ["U4_alpha_u0", "U4_alpha_u1", "K12_selector_d2_minus_3c3a"]
    block1_vars = [symbols[f"p{i}"] for i in range(27, 34)]
    block1_eqs = [sp.expand(gens[tag].subs(fixed)) for tag in block1_tags]
    sol1, meta1 = solve_linear_block(block1_eqs, block1_vars)
    fixed.update(sol1)

    # Solve both U5-alpha order rows and the exact definition T=0.
    block2_tags = ["U5_alpha_u0", "U5_alpha_u1", "define_aux_T"]
    block2_vars = [symbols[f"p{i}"] for i in range(39, 46)]
    block2_eqs = [sp.expand(gens[tag].subs(fixed)) for tag in block2_tags]
    sol2, meta2 = solve_linear_block(block2_eqs, block2_vars)
    fixed.update(sol2)

    # Solve the exact definition N=0 and K14 in the p50..p55 block.
    block3_tags = ["define_aux_N", "K14_selector_d3_minus_27c6f"]
    block3_vars = [symbols[f"p{i}"] for i in range(50, 56)]
    block3_eqs = [sp.expand(gens[tag].subs(fixed)) for tag in block3_tags]
    sol3, meta3 = solve_linear_block(block3_eqs, block3_vars)
    fixed.update(sol3)

    # The remaining quotient definitions determine L and B3 (already 0).
    if sp.simplify(gens["define_aux_B3"].subs(fixed)) != 0:
        raise ValueError("B3=0 is inconsistent with its exact quotient definition")
    value = sp.solve(gens["define_aux_L"].subs(fixed), symbols["aux_L"])
    if len(value) != 1:
        raise ValueError("could not solve define_aux_L uniquely")
    fixed[symbols["aux_L"]] = sp.simplify(value[0])

    # Determine the saturation variable exactly; reject if d*a is zero.
    sat_poly = gens["saturation_sat_qs_d_a_minus_1"].subs(fixed)
    sat_solution = sp.solve(sat_poly, symbols["sat_qs"])
    if len(sat_solution) != 1:
        raise ValueError("could not solve the d*a saturation equation")
    fixed[symbols["sat_qs"]] = sp.simplify(sat_solution[0])

    residuals = {tag: sp.simplify(poly.subs(fixed)) for tag, poly in gens.items()}
    nonzero = {tag: str(value) for tag, value in residuals.items() if value != 0}
    if nonzero:
        raise AssertionError(f"nonzero residuals: {nonzero}")

    assignment = {str(sym): str(sp.simplify(value)) for sym, value in sorted(
        fixed.items(), key=lambda kv: str(kv[0]))}
    full_parameter_symbols = {f"p{i}": symbols.get(f"p{i}", sp.Symbol(f"p{i}"))
                              for i in range(56)}
    pivot_values = {
        name: sp.simplify(exact_expr(expression, full_parameter_symbols).subs(fixed))
        for name, expression in src["constant_linear_reduction"]["substitutions"].items()
    }
    full_parameter_assignment = {f"p{i}": sp.Integer(0) for i in range(56)}
    for name, value in fixed.items():
        if str(name).startswith("p") and str(name)[1:].isdigit():
            full_parameter_assignment[str(name)] = sp.simplify(value)
    full_parameter_assignment.update(pivot_values)
    full_parameter_assignment["p1"] = sp.Integer(0)
    full_parameter_assignment["p2"] = sp.Integer(0)
    terminal_p_replay = replay_terminal_p(
        json.loads(STATE.read_text(encoding="utf-8")), full_parameter_assignment)
    if (terminal_p_replay["through_distance_5_positions"] != 32
            or not terminal_p_replay["all_through_distance_5_zero"]
            or terminal_p_replay["Q_positions_checked"] != 64
            or not terminal_p_replay["all_Q_rows_zero"]):
        raise AssertionError("the exact candidate failed a terminal P/Q replay row")
    exact_full_parameter_assignment = {
        name: str(sp.simplify(value)) for name, value in sorted(full_parameter_assignment.items())
    }
    saturation_product = sp.simplify(1 / fixed[symbols["sat_qs"]])
    receipt = {
        "status": "EXACT_ADMISSIBLE_POINT_FOR_SELECTED_ALPHA_LOCAL_NECESSARY_IDEAL",
        "scope": "U1=0, alpha-local q=4,s=2 ideal only; not a full 2+2 pair and not an exclusion proof",
        "input_state_sha256": state_sha,
        "terminal_P_rows_through_distance_5": 32,
        "terminal_P_slice_counts": [6, 6, 6, 5, 5, 4],
        "terminal_P_nonzero_positions": 0,
        "terminal_Q_positions_checked": 64,
        "terminal_Q_added_equations": 0,
        "generator_set_sha256": src["constant_linear_reduction"]["post_substitution_generator_set_sha256"],
        "generator_count": len(gens),
        "zero_case_choice": {"B3": "0", "T": "0", "N": "0",
                              "method": "all imposed as equalities; no quotient saturation or division"},
        "linear_blocks": {"U4_K12": meta1, "U5_T": meta2, "N_K14": meta3},
        "assignment": assignment,
        "full_parameter_assignment_with_unused_parameters_set_to_zero": exact_full_parameter_assignment,
        "specialized_terminal_P_replay": terminal_p_replay,
        "constant_pivot_values_from_saved_exact_RREF": {
            name: str(value) for name, value in pivot_values.items()
        },
        "all_generator_residuals_zero": True,
        "saturation_value": assignment["sat_qs"],
        "saturation_product_nonzero": True,
        "exact_d_times_a_value_from_saturation": str(saturation_product),
        "K15_residual": str(residuals["K15_compact_division_free"]),
        "warning": "This proves nonemptiness only of the listed selected necessary ideal on the chosen alpha-local stratum. The t=1 q>=5/zero-U2 case remains separate, and no full branch conclusion follows."
    }
    OUT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
