#!/usr/bin/env python3
"""Export the verified support-12 QQ subsystem as a Singular input."""

from __future__ import annotations

import hashlib
import json
import argparse
from fractions import Fraction
from pathlib import Path


INPUT = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_input.json")
OUTPUT = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12.sing")
EXPECTED_INPUT_SHA256 = "A9DF8405D43EAFFB98C724165C0C7A647C0952D7D8D1A8630CD26700898C6F56"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=OUTPUT)
    ap.add_argument("--use-linear-reduction", action="store_true")
    ap.add_argument("--weighted-order", action="store_true")
    args = ap.parse_args()
    raw = INPUT.read_bytes()
    input_hash = hashlib.sha256(raw).hexdigest().upper()
    if input_hash != EXPECTED_INPUT_SHA256:
        raise ValueError(f"support-12 input hash mismatch: {input_hash}")
    data = json.loads(raw)
    names = tuple(data["selected_variables"])
    selected_indices = list(range(len(data["equations"])))
    verified = None
    if args.use_linear_reduction:
        reduction_path = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_linear_reduce.json")
        verified_path = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_linear_reduce_verified.json")
        reduction = json.loads(reduction_path.read_text(encoding="utf-8"))
        verified = json.loads(verified_path.read_text(encoding="utf-8"))
        if reduction.get("input_sha256") != input_hash or verified.get("input_sha256") != input_hash:
            raise ValueError("linear reduction is not bound to the support-12 input")
        if verified.get("reduction_report_sha256") != hashlib.sha256(reduction_path.read_bytes()).hexdigest().upper():
            raise ValueError("linear reduction verification receipt is stale")
        if not verified.get("all_142_equations_span_exactly_4_Q_generators"):
            raise ValueError("exact Q-linear reduction was not verified")
        selected_indices = list(verified["independent_equation_indices"])
    if args.weighted_order:
        weighted_path = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_weighted_homogeneity.json")
        weighted = json.loads(weighted_path.read_text(encoding="utf-8"))
        if weighted.get("input_sha256") != input_hash or not weighted.get(
                "all_142_input_equations_weight_homogeneous"):
            raise ValueError("weighted-homogeneity receipt is missing or stale")
        weights = tuple(weighted["primitive_positive_integer_weights"])
        if weights != (4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 2, 3):
            raise ValueError(f"unexpected weight vector: {weights}")
        order_text = f"wp({','.join(map(str, weights))})"
    else:
        order_text = "dp"
    normalized: dict[tuple, tuple[str, dict[tuple, Fraction]]] = {}
    input_term_count = 0
    for row_index in selected_indices:
        row = data["equations"][row_index]
        poly: dict[tuple, Fraction] = {}
        for monomial, coefficient in row["terms"]:
            powers = dict(monomial)
            exp = tuple(int(powers.get(name, 0)) for name in names)
            poly[exp] = poly.get(exp, Fraction(0)) + Fraction(coefficient)
            input_term_count += 1
        poly = {exp: c for exp, c in poly.items() if c}
        if not poly:
            continue
        scale = poly[sorted(poly)[0]]
        canonical = {exp: coeff / scale for exp, coeff in poly.items()}
        key = tuple((exp, canonical[exp]) for exp in sorted(canonical))
        normalized.setdefault(key, (row["tag"], canonical))

    def coefficient_text(value: Fraction) -> str:
        if value.denominator == 1:
            return str(value.numerator)
        return f"({value.numerator}/{value.denominator})"

    lines = [
        "// Exact support-12 subsystem; generated from hash-bound rational JSON.",
        f"// input_sha256={input_hash}",
        f"// source_checkpoint_sha256={data['source_sha256']}",
        f"// equations={data['equation_count']}; unique_up_to_QQ_scale={len(normalized)}",
        f"// selected_equation_indices={selected_indices}",
        f"ring S=0,({','.join(names)}),{order_text};",
    ]
    generator_names = []
    for i, (tag, poly) in enumerate(normalized.values(), 1):
        generator = f"f{i}"
        generator_names.append(generator)
        terms = []
        for exp in sorted(poly, reverse=True):
            monomial_parts = []
            for name, power in zip(names, exp):
                if power == 1:
                    monomial_parts.append(name)
                elif power > 1:
                    monomial_parts.append(f"{name}^{power}")
            monomial = "*".join(monomial_parts)
            coeff = coefficient_text(poly[exp])
            if monomial:
                terms.append(f"{coeff}*{monomial}")
            else:
                terms.append(coeff)
        lines.append(f"poly {generator}=" + "+\n".join(terms) + ";")
        lines.append(f"// tag {generator}: {tag}")
    lines.extend([
        "ideal I=" + ",".join(generator_names) + ";",
        "option(redSB);",
        "ideal G=slimgb(I);",
        "print(\"basis_size\");",
        "print(size(G));",
        "print(\"unit_remainder\");",
        "print(reduce(1,G));",
        f"write(\"{args.output.with_suffix('.basis.txt').as_posix()}\",G);",
        "quit;",
    ])
    args.output.write_bytes(("\n".join(lines) + "\n").encode("utf-8"))
    print(json.dumps({
        "input": str(INPUT),
        "input_sha256": input_hash,
        "branch": data["branch"],
        "input_equations": len(data["equations"]),
        "selected_equation_indices": selected_indices,
        "exact_linear_reduction_verified": verified is not None,
        "selected_equations": len(selected_indices),
        "unique_generators_up_to_QQ_scaling": len(normalized),
        "input_terms": input_term_count,
        "variables": names,
        "monomial_order": order_text,
        "weighted_order": order_text,
        "method": "Singular slimgb over QQ",
        "output": str(args.output),
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
