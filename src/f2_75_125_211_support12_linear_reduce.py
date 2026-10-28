#!/usr/bin/env python3
"""Exact sparse Q-row reduction of support-12 polynomial generators.

This only compresses generators by constant rational linear combinations; it
does not compute a Groebner basis or exclude the 2+1+1 branch.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from flint import fmpq


INPUT = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_input.json")
EXPECTED_INPUT_SHA256 = "A9DF8405D43EAFFB98C724165C0C7A647C0952D7D8D1A8630CD26700898C6F56"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest().upper()


def add_scaled(target: dict, source: dict, scale: fmpq) -> None:
    if not scale:
        return
    for key, value in source.items():
        new = target.get(key, fmpq(0)) + scale * value
        if new:
            target[key] = new
        else:
            target.pop(key, None)


def clean(values: dict[int, fmpq]) -> dict[int, fmpq]:
    return {key: value for key, value in values.items() if value}


def main() -> None:
    input_hash = sha256(INPUT)
    if input_hash != EXPECTED_INPUT_SHA256:
        raise ValueError(f"support-12 input hash mismatch: {input_hash}")
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    names = tuple(data["selected_variables"])
    if len(names) != 12 or data["equation_count"] != len(data["equations"]):
        raise ValueError("support-12 extraction metadata mismatch")
    index = {name: i for i, name in enumerate(names)}

    pivots: dict[tuple[int, ...], tuple[dict, dict[int, fmpq]]] = {}
    independent: list[int] = []
    independent_pivots: list[list[int]] = []
    dependencies: dict[int, dict[int, fmpq]] = {}
    input_terms = 0

    for i, equation in enumerate(data["equations"]):
        row: dict[tuple[int, ...], fmpq] = {}
        for monomial, coefficient in equation["terms"]:
            exp = [0] * len(names)
            for name, power in monomial:
                if name not in index:
                    raise ValueError(f"unexpected variable {name} in {equation['tag']}")
                exp[index[name]] += int(power)
            key = tuple(exp)
            row[key] = row.get(key, fmpq(0)) + fmpq(str(coefficient))
            input_terms += 1
        row = {key: value for key, value in row.items() if value}
        combination: dict[int, fmpq] = {}
        while row:
            pivot = min(row)
            if pivot not in pivots:
                break
            coefficient = row[pivot]
            basis_row, basis_combination = pivots[pivot]
            add_scaled(row, basis_row, -coefficient)
            add_scaled(combination, basis_combination, coefficient)

        if not row:
            dependencies[i] = clean(combination)
            continue

        pivot = min(row)
        pivot_coefficient = row[pivot]
        combination[i] = fmpq(1)
        # Original row = residual + accumulated combination of old basis rows.
        # Therefore residual/pivot_coefficient is a combination of this input
        # row and the previously selected independent input rows.
        for key in list(combination):
            if key != i:
                combination[key] = -combination[key]
                if not combination[key]:
                    combination.pop(key)
        combination[i] = fmpq(1)
        basis_row = {key: value / pivot_coefficient for key, value in row.items()}
        basis_combination = {key: value / pivot_coefficient
                             for key, value in combination.items()}
        pivots[pivot] = (basis_row, basis_combination)
        independent.append(i)
        independent_pivots.append(list(pivot))

    tags = [equation["tag"] for equation in data["equations"]]
    result = {
        "input": str(INPUT),
        "input_sha256": input_hash,
        "branch": data["branch"],
        "selected_variables": list(names),
        "input_equations": len(data["equations"]),
        "input_terms": input_terms,
        "exact_Q_linear_rank": len(independent),
        "independent_equation_indices": independent,
        "independent_equation_tags": [tags[i] for i in independent],
        "independent_pivot_monomials": independent_pivots,
        "dependent_equation_relations": [
            {
                "equation_index": i,
                "equation_tag": tags[i],
                "equals": [[j, str(value)] for j, value in sorted(relation.items())],
            }
            for i, relation in sorted(dependencies.items())
        ],
        "basis_pivot_count": len(pivots),
        "scope": "Exact Q-linear generator compression of the support-12 subideal only; not a Groebner basis, characteristic-zero exclusion, or branch exclusion.",
    }
    output = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_linear_reduce.json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "input_sha256": input_hash,
        "input_equations": len(data["equations"]),
        "exact_Q_linear_rank": len(independent),
        "dependent_equations": len(dependencies),
        "output": str(output),
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
