"""Extract and verify global pivot charts for the two exact S(alpha) components.

Input is the archived Singular factorization log whose ideal-equality checks
were independently replayed.  The script verifies that the cubic component is
globally linear in p16 and the quartic component globally linear in p27 over
Q(sqrt(5)); these are charts for subsequent branch-specific terminal rows,
not exclusions of either component.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import sympy as sp


LOG = Path("research/f2_75_125_sparse_ring_lift_22_Salpha_slice_factorcheck_verified.log")
OUT = Path("research/f2_75_125_Salpha_factor_charts_20261018.json")


def main() -> None:
    raw = LOG.read_bytes()
    text = raw.decode("utf-8")
    expressions: dict[int, str] = {}
    for degree in (3, 4):
        match = re.search(rf"^\s*_\[{degree-1}\]=(.*)$", text, re.MULTILINE)
        if not match:
            raise ValueError(f"factor of weighted degree {degree} not found")
        expressions[degree] = match.group(1).strip()

    sqrt5 = sp.sqrt(5)
    a = sp.Symbol("a")
    variables = {
        name: sp.Symbol(name)
        for name in ("p1", "p7", "p8", "p9", "p16", "p17", "p18",
                     "p19", "p20", "p27", "p28", "p29", "p30", "p31",
                     "p32", "p33")
    }
    namespace = {"a": a, **variables}
    pivots = {3: "p16", 4: "p27"}
    records = {}
    factor_polynomials = {}
    solution_polynomials = {}
    for degree, source in expressions.items():
        expression = sp.sympify(source, locals=namespace).subs(a, sqrt5)
        expression = sp.expand(expression)
        factor_polynomials[degree] = expression
        pivot_name = pivots[degree]
        pivot = variables[pivot_name]
        derivative = sp.simplify(sp.diff(expression, pivot))
        if not derivative.is_number or derivative == 0:
            raise AssertionError(f"degree-{degree} pivot coefficient is not a nonzero constant")
        if derivative.has(*variables.values()):
            raise AssertionError(f"degree-{degree} pivot is not globally linear")
        solution = sp.cancel(-(expression - derivative * pivot) / derivative)
        solution_polynomials[degree] = solution
        if pivot in solution.free_symbols:
            raise AssertionError(f"degree-{degree} solution still contains its pivot")
        remainder = sp.simplify(sp.expand(expression.subs(pivot, solution)))
        if remainder != 0:
            raise AssertionError(
                f"degree-{degree} substitution does not annihilate its factor; "
                f"pivot coefficient={derivative}; remainder={str(remainder)[:800]}"
            )
        support = sorted(str(symbol) for symbol in expression.free_symbols)
        records[str(degree)] = {
            "weighted_degree": degree,
            "pivot": pivot_name,
            "pivot_coefficient_in_Qsqrt5": str(derivative),
            "pivot_coefficient_nonzero": True,
            "solution": str(solution),
            "factor_support": support,
            "factor_terms": len(sp.Poly(expression, *variables.values(), extension=sqrt5).terms()),
            "exact_substitution_remainder": "0",
            "component_dimension": 15,
        }

    p16 = variables["p16"]
    p27 = variables["p27"]
    intersection_factor = sp.expand(
        factor_polynomials[4].subs(p16, solution_polynomials[3])
    )
    intersection_pivot = sp.simplify(sp.diff(intersection_factor, p27))
    if intersection_pivot != sp.simplify(sp.diff(factor_polynomials[4], p27)):
        raise AssertionError("the F4 pivot changed on the F3 chart")
    intersection_solution = sp.cancel(
        -(intersection_factor - intersection_pivot * p27) / intersection_pivot
    )
    intersection_remainder = sp.simplify(
        sp.expand(intersection_factor.subs(p27, intersection_solution))
    )
    if intersection_remainder != 0:
        raise AssertionError("the intersection chart does not annihilate F4")

    receipt = {
        "input_log": str(LOG),
        "input_log_sha256": hashlib.sha256(raw).hexdigest(),
        "coefficient_field": "Q(sqrt(5))",
        "components": records,
        "intersection_chart": {
            "pivots_in_order": ["p16", "p27"],
            "second_pivot_coefficient_in_Qsqrt5": str(intersection_pivot),
            "dimension": 14,
            "exact_substitution_remainder": "0",
        },
        "scope": "Global pivot charts for the two components of the exact depth-seven S(alpha)-slice ideal; no component exclusion is claimed.",
    }
    OUT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
