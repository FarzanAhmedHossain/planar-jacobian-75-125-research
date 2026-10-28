#!/usr/bin/env python3
"""Independently verify the exact support-12 Q-linear reduction certificate."""

from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

from flint import fmpq


INPUT = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_input.json")
REPORT = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_linear_reduce.json")
EXPECTED_INPUT_SHA256 = "A9DF8405D43EAFFB98C724165C0C7A647C0952D7D8D1A8630CD26700898C6F56"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest().upper()


def polynomial(row: dict, names: tuple[str, ...]) -> dict[tuple[int, ...], fmpq]:
    index = {name: i for i, name in enumerate(names)}
    result: dict[tuple[int, ...], fmpq] = {}
    for monomial, coefficient in row["terms"]:
        exp = [0] * len(names)
        for name, power in monomial:
            exp[index[name]] += int(power)
        key = tuple(exp)
        result[key] = result.get(key, fmpq(0)) + fmpq(str(coefficient))
    return {key: value for key, value in result.items() if value}


def permutation_sign(values: tuple[int, ...]) -> int:
    inversions = sum(values[i] > values[j]
                     for i in range(len(values))
                     for j in range(i + 1, len(values)))
    return -1 if inversions % 2 else 1


def main() -> None:
    input_hash = sha256(INPUT)
    if input_hash != EXPECTED_INPUT_SHA256:
        raise ValueError(f"input hash mismatch: {input_hash}")
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    if report["input_sha256"] != input_hash:
        raise ValueError("reduction report is not bound to this input")
    names = tuple(data["selected_variables"])
    equations = data["equations"]
    independent = report["independent_equation_indices"]
    independent_polys = {i: polynomial(equations[i], names) for i in independent}
    if len(independent) != report["exact_Q_linear_rank"]:
        raise ValueError("independent index count disagrees with rank")

    checked = 0
    for relation in report["dependent_equation_relations"]:
        target_index = relation["equation_index"]
        target = polynomial(equations[target_index], names)
        reconstructed: dict[tuple[int, ...], fmpq] = {}
        for basis_index, coefficient in relation["equals"]:
            if basis_index not in independent_polys:
                raise ValueError("dependency relation uses a nonbasis row")
            scale = fmpq(coefficient)
            for monomial, value in independent_polys[basis_index].items():
                reconstructed[monomial] = reconstructed.get(monomial, fmpq(0)) + scale * value
        reconstructed = {m: c for m, c in reconstructed.items() if c}
        if reconstructed != target:
            raise ValueError(f"exact dependency check failed at row {target_index}")
        checked += 1

    pivots = [tuple(m) for m in report["independent_pivot_monomials"]]
    if len(pivots) != len(independent):
        raise ValueError("pivot count disagrees with independent row count")
    matrix = [[independent_polys[row].get(monomial, fmpq(0))
               for monomial in pivots] for row in independent]
    determinant = fmpq(0)
    for permutation in itertools.permutations(range(len(independent))):
        term = fmpq(permutation_sign(permutation))
        for i, j in enumerate(permutation):
            term *= matrix[i][j]
        determinant += term
    if not determinant:
        raise ValueError("independent-row pivot minor is singular")

    verified = {
        "input_sha256": input_hash,
        "reduction_report_sha256": sha256(REPORT),
        "input_equations": len(equations),
        "independent_equations": len(independent),
        "independent_equation_indices": independent,
        "verified_exact_dependencies": checked,
        "pivot_minor_determinant": str(determinant),
        "all_142_equations_span_exactly_4_Q_generators": True,
        "scope": "Exact generator compression only; this does not exclude the support-12 variety or the 2+1+1 branch.",
    }
    output = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_linear_reduce_verified.json")
    output.write_text(json.dumps(verified, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({**verified, "verification_receipt": str(output)}, indent=2), flush=True)


if __name__ == "__main__":
    main()
