"""Profile the exact residual ideal retained by a sparse F2 band lift."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from f2_75_125_sparse_terminal_q_postcheck import load_sparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    state = load_sparse(args.checkpoint)
    R, gens, constraints = state["R"], state["gens"], state["constraints"]

    variables = set()
    categories = Counter()
    degrees = Counter()
    constants = []
    linear_equations = []
    distinct = set()
    term_counts = []
    for tag, polynomial in constraints:
        categories[tag.split(":", 1)[0]] += 1
        distinct.add(polynomial)
        term_counts.append(len(polynomial))
        degree = max((sum(powers[1:]) for powers in polynomial), default=0)
        degrees[degree] += 1
        for powers in polynomial:
            variables.update(i for i, exponent in enumerate(powers[1:], start=1)
                             if exponent)
        if len(polynomial) == 1:
            (powers, coefficient), = polynomial.items()
            if not any(powers):
                constants.append((tag, str(coefficient)))
        if degree <= 1:
            linear_equations.append(tag)

    report = {
        "branch": state["data"]["branch"],
        "constraints": len(constraints),
        "distinct_polynomials": len(distinct),
        "variables_used": [str(gens[i]) for i in sorted(variables)],
        "variable_count": len(variables),
        "degree_distribution": dict(sorted(degrees.items())),
        "term_counts": {"min": min(term_counts, default=0),
                        "max": max(term_counts, default=0),
                        "total": sum(term_counts)},
        "nonzero_constant_constraints": constants,
        "linear_equations": linear_equations,
        "categories": dict(categories),
        "scope": "profile only; no ideal-membership or unit-ideal claim",
    }
    output = json.dumps(report, indent=2)
    if args.report:
        args.report.write_text(output + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
