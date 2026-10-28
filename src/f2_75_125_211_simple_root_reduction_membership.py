"""Check new simple-root equations against exact d19 support Groebner bases."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ


SOURCE_HASH = "0f293c437b1cc18f379e3a8d0bb2e98647616e062fa321e91cfdc45ebe67926d"
DIVISIBILITY = Path(
    "research/f2_75_125_211_simple_root_divisibility_equations_20261025.json"
)
REPORTS = {
    "p3_mod_T": Path(
        "research/f2_75_125_sparse_ring_lift_211_roots19_u1zero_support3_gb.json"
    ),
    "p5_mod_S": Path(
        "research/f2_75_125_sparse_ring_lift_211_roots19_u1zero_support10_gb.json"
    ),
}
OUTPUT = Path(
    "research/f2_75_125_211_simple_root_reduction_membership_20261025.json"
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_equation(encoded, names, symbols, support_names):
    expression = sp.Integer(0)
    equation_support = set()
    for sparse, coefficient in encoded["terms"]:
        term = sp.Rational(coefficient)
        for index, exponent in sparse:
            name = names[int(index)]
            equation_support.add(name)
            term *= symbols[name] ** int(exponent)
        expression += term
    if not equation_support.issubset(support_names):
        raise ValueError(f"equation has unexpected support {equation_support}")
    return sp.Poly(expression, *(symbols[name] for name in support_names), domain=QQ)


def linear_rank_certificate(equations, variable_names, symbol_map, substitutions=None):
    substitutions = substitutions or {}
    variables = [symbol_map[name] for name in variable_names]
    rows = []
    for equation in equations:
        expression = (sp.Integer(0) if equation["normal_form"] is None else
                      sp.sympify(equation["normal_form"], locals=symbol_map))
        expression = sp.expand(expression.subs(substitutions))
        polynomial = sp.Poly(expression, *variables, domain=QQ)
        if polynomial.total_degree() > 1 or polynomial.coeff_monomial(1):
            raise ValueError("expected homogeneous linear normal forms")
        rows.append([
            polynomial.coeff_monomial(variable) for variable in variables
        ])

    selected_indices = []
    selected_rows = []
    current_rank = 0
    for index, row in enumerate(rows):
        candidate = sp.Matrix([*selected_rows, row])
        rank = candidate.rank()
        if rank > current_rank:
            selected_indices.append(index)
            selected_rows.append(row)
            current_rank = rank
    certificate = {
        "variables": variable_names,
        "rank": current_rank,
        "equation_count": len(rows),
        "row_indices_one_based": [index + 1 for index in selected_indices],
        "independent_rows": [[str(value) for value in row]
                             for row in selected_rows],
        "forces_all_variables_zero": current_rank == len(variable_names),
    }
    if current_rank == len(variable_names) and selected_rows:
        certificate["determinant"] = str(sp.Matrix(selected_rows).det())
    return certificate


def main():
    divisibility = json.loads(DIVISIBILITY.read_text(encoding="utf-8"))
    if divisibility["source_sha256"] != SOURCE_HASH:
        raise ValueError("divisibility equations are not bound to d19 source")
    names = divisibility["generator_names"]
    symbol_map = {name: sp.Symbol(name) for name in names}
    results = {}
    support_bases = {}
    support_variables = {}

    for band_key, report_path in REPORTS.items():
        report = json.loads(report_path.read_text(encoding="utf-8"))
        if report.get("source_sha256") != SOURCE_HASH:
            raise ValueError(f"{report_path} is bound to another source")
        support_names = report["selected_variables"]
        variables = [symbol_map[name] for name in support_names]
        basis = [
            sp.Poly(sp.sympify(text, locals=symbol_map), *variables, domain=QQ)
            for text in report["groebner_basis"]
        ]
        support_bases[band_key] = basis
        support_variables[band_key] = variables
        input_equations = divisibility["equations"][band_key][
            "coefficient_equations"
        ]
        equation_results = []
        for encoded in input_equations:
            polynomial = read_equation(
                encoded, names, symbol_map, support_names
            )
            _, remainder = sp.reduced(
                polynomial, basis, *variables, order="grevlex"
            )
            equation_results.append({
                "t_degree": encoded["t_degree"],
                "term_count": encoded["term_count"],
                "normal_form_is_zero": remainder.is_zero,
                "normal_form_term_count": len(remainder.terms()),
                "normal_form": None if remainder.is_zero else str(remainder.as_expr()),
            })
        results[band_key] = {
            "support_variables": support_names,
            "support_ideal_unit": report["unit_ideal"],
            "support_groebner_basis_length": len(basis),
            "equations": equation_results,
            "all_members_of_support_ideal": all(
                item["normal_form_is_zero"] for item in equation_results
            ),
        }

    p3_rank = linear_rank_certificate(
        results["p3_mod_T"]["equations"],
        results["p3_mod_T"]["support_variables"],
        symbol_map,
    )
    p5_rank_after_p3_zero = linear_rank_certificate(
        results["p5_mod_S"]["equations"],
        [name for name in results["p5_mod_S"]["support_variables"]
         if name not in {"p13", "p14", "p15"}],
        symbol_map,
        substitutions={symbol_map[name]: 0 for name in ("p13", "p14", "p15")},
    )

    p3_variables = support_variables["p3_mod_T"]
    p3_forms = [
        sp.Poly(
            sp.Integer(0) if equation["normal_form"] is None else
            sp.sympify(equation["normal_form"], locals=symbol_map),
            *p3_variables,
            domain=QQ,
        )
        for equation in results["p3_mod_T"]["equations"]
    ]
    augmented_support3 = sp.groebner(
        [poly.as_expr() for poly in support_bases["p3_mod_T"] + p3_forms],
        *p3_variables,
        order="grevlex",
        domain=QQ,
    )
    pure_powers = {}
    for variable in p3_variables:
        for polynomial in augmented_support3.polys:
            expression = sp.expand(polynomial.as_expr())
            if expression != 0 and expression.free_symbols.issubset({variable}):
                univariate = sp.Poly(expression, variable, domain=QQ)
                if len(univariate.terms()) == 1 and univariate.TC():
                    degree = univariate.degree()
                    if degree > 0:
                        pure_powers[str(variable)] = str(expression)
                        break

    line_parameter = sp.Symbol("r")
    line_substitutions = {
        symbol_map["p13"]: sp.Rational(452, 5) * line_parameter,
        symbol_map["p14"]: -14 * line_parameter,
        symbol_map["p15"]: line_parameter,
    }
    p3_on_support3_line = [
        sp.expand((sp.Integer(0) if equation["normal_form"] is None else
                   sp.sympify(equation["normal_form"], locals=symbol_map))
                  .subs(line_substitutions))
        for equation in results["p3_mod_T"]["equations"]
    ]
    p5_on_support3_line = [
        sp.expand((sp.Integer(0) if equation["normal_form"] is None else
                   sp.sympify(equation["normal_form"], locals=symbol_map))
                  .subs(line_substitutions))
        for equation in results["p5_mod_S"]["equations"]
    ]
    p34, p35 = symbol_map["p34"], symbol_map["p35"]
    p5_line_solution = sp.solve(p5_on_support3_line, (p34, p35), dict=True)
    if len(p5_line_solution) != 1:
        raise ValueError("expected a unique p34,p35 solution on the support3 line")

    OUTPUT.write_text(json.dumps({
        "scope": "exact ideal-membership check for added simple-root equations only",
        "d19_source_sha256": SOURCE_HASH,
        "divisibility_equation_sha256": digest(DIVISIBILITY),
        "support3_report_sha256": digest(REPORTS["p3_mod_T"]),
        "support10_report_sha256": digest(REPORTS["p5_mod_S"]),
        "results": results,
        "linear_consequences": {
            "after_reducing_p3_remainders_mod_support3_ideal": p3_rank,
            "after_also_setting_p13_p14_p15_zero_in_p5_remainders":
                p5_rank_after_p3_zero,
        },
        "augmented_support3_ideal": {
            "constructed_from": "exact support3 checkpoint Groebner basis plus the seven p3 mod T normal forms",
            "unit_ideal": any(poly.as_expr() == 1
                               for poly in augmented_support3.polys),
            "zero_dimensional": augmented_support3.is_zero_dimensional,
            "pure_power_in_each_variable": pure_powers,
            "forces_p13_p14_p15_zero_set_theoretically": (
                set(pure_powers) == {str(variable) for variable in p3_variables}
            ),
            "groebner_basis": [str(poly.as_expr()) for poly in augmented_support3.polys],
        },
        "restriction_to_known_support3_line": {
            "line": "(p13,p14,p15)=(452*r/5,-14*r,r)",
            "all_p3_remainders_vanish": all(value == 0
                                             for value in p3_on_support3_line),
            "p5_remainders_after_line": [str(value)
                                         for value in p5_on_support3_line],
            "solved_for": ["p34", "p35"],
            "solution": {str(key): str(value)
                         for key, value in p5_line_solution[0].items()},
        },
        "whole_2+1+1_branch_excluded_percent": 0.0,
    }, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value["all_members_of_support_ideal"]
                      for key, value in results.items()}, indent=2))
    for key, value in results.items():
        print(key, json.dumps(value["equations"], separators=(",", ":")))
    print("linear_consequences", json.dumps({
        "p3": p3_rank,
        "p5_after_p3_zero": p5_rank_after_p3_zero,
    }, separators=(",", ":")))
    print("augmented_support3", json.dumps({
        "unit_ideal": any(poly.as_expr() == 1
                           for poly in augmented_support3.polys),
        "zero_dimensional": augmented_support3.is_zero_dimensional,
        "pure_powers": pure_powers,
    }, separators=(",", ":")))
    print("restriction_to_support3_line", json.dumps({
        "all_p3_remainders_vanish": all(value == 0
                                         for value in p3_on_support3_line),
        "p5_remainders": [str(value) for value in p5_on_support3_line],
        "p34_p35_solution": {str(key): str(value)
                             for key, value in p5_line_solution[0].items()},
    }, separators=(",", ":")))
    print(f"report={OUTPUT}")


if __name__ == "__main__":
    main()
