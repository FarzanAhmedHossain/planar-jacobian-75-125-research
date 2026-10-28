"""Triangular exact graph reduction for the s=9, rho2 != 0 local rows.

The input rows are already exact rational polynomials from the support-tight
Magnus replay. This only parameterizes those necessary rows; it is not a branch
exclusion.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "magnus_75_125_s9_edge_support_case_split_20261021.json"
OUTPUT = ROOT / "magnus_75_125_s9_nonzero_rho2_graph_20261021.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_expr(source: str) -> sp.Expr:
    names = set(re.findall(r"\b[A-Za-z_]\w*\b", source))
    symbols = {name: sp.Symbol(name) for name in names}
    return sp.sympify(source, locals=symbols)


def main() -> None:
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    rows = {
        (item["band"], item["x_degree"]): parse_expr(item["factorized_equation"])
        for item in data["exact_branch"]["remaining_nonzero_rows_after_graph_substitution"]
    }
    expected_keys = {
        ("depth7", 5), ("depth7", 6), ("depth7", 7),
        ("depth8", 3), ("depth8", 4), ("depth8", 5), ("depth8", 6), ("depth8", 7),
    }
    if set(rows) != expected_keys:
        raise AssertionError("unexpected residual-row set")

    symbols = {name: sp.Symbol(name) for name in (
        "mu rho2 a4 u3 c4 d3 d4 d5 c5 c6 c7 c8"
    ).split()}
    replacements = {
        symbols["c4"]: symbols["a4"] * symbols["u3"],
    }
    rows = {key: sp.cancel(value.subs(replacements)) for key, value in rows.items()}

    # These pivots are safe on the already selected mu != 0, rho2 != 0 chart.
    pivot_order = [
        (("depth7", 5), "d3"),
        (("depth7", 6), "d4"),
        (("depth7", 7), "d5"),
        (("depth8", 4), "c5"),
        (("depth8", 5), "c6"),
        (("depth8", 6), "c7"),
        (("depth8", 7), "c8"),
    ]
    graph: dict[sp.Symbol, sp.Expr] = dict(replacements)
    pivots = []
    for key, name in pivot_order:
        variable = symbols[name]
        equation = sp.cancel(rows[key].subs(graph))
        poly = sp.Poly(equation, variable)
        if poly.degree() != 1:
            raise AssertionError(f"{key} is not linear in {name}: degree={poly.degree()}")
        coefficient = sp.factor(poly.coeff_monomial(variable))
        allowed = {symbols["mu"], symbols["rho2"]}
        if not coefficient.free_symbols.issubset(allowed) or coefficient == 0:
            raise AssertionError(f"unsafe pivot coefficient for {name}: {coefficient}")
        # A nonzero rational times mu*rho2^k is invertible on the guarded chart.
        coeff_poly = sp.Poly(coefficient, symbols["mu"], symbols["rho2"])
        if len(coeff_poly.terms()) != 1:
            raise AssertionError(f"pivot is not a single monomial: {coefficient}")
        solution = sp.factor(-poly.coeff_monomial(1) / coefficient)
        if sp.cancel(equation.subs(variable, solution)) != 0:
            raise AssertionError(f"pivot substitution did not kill {key}")
        graph[variable] = solution
        pivots.append({
            "row": f"{key[0]}_x{key[1]}",
            "pivot": name,
            "coefficient": str(coefficient),
            "solution": str(solution),
        })

    def apply_graph(expression: sp.Expr) -> sp.Expr:
        current = expression
        for _ in range(len(graph) + 2):
            updated = sp.expand(current.subs(graph, simultaneous=False))
            if updated == current:
                return updated
            current = updated
        raise AssertionError("pivot graph substitution did not stabilize")

    residual_rows = []
    for key in sorted(rows):
        reduced = sp.factor(apply_graph(rows[key]))
        if reduced != 0:
            residual_rows.append({"row": f"{key[0]}_x{key[1]}", "equation": str(reduced)})
    if residual_rows:
        raise AssertionError(f"unexpected residual equations after graph: {residual_rows}")

    for item in pivots:
        item["solution"] = str(sp.factor(apply_graph(symbols[item["pivot"]])))

    result = {
        "input_sha256": digest(INPUT),
        "guards": ["mu != 0", "rho2 != 0"],
        "inherited_relations": [
            "sigma0=0", "sigma1=2*a4*rho2/5", "c2=0", "u2=0", "c3=0",
            "d2=107*rho2*u3/103", "c4=a4*u3",
        ],
        "exact_linear_pivots": pivots,
        "residual_rows_after_graph": residual_rows,
        "scope": "the eight saved necessary rows on one local chart become a graph; no contradiction or full branch exclusion",
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
