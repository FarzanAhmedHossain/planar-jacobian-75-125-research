#!/usr/bin/env python3
"""Check the d23 support-12 equations against the verified d18 input."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path


EXPECTED_D23_CHECKPOINT_SHA256 = "270f0eb737463c75d3bdfbcdfb3a303a2e7fa0af2da13f63a9e2d01ec5e04b59"
EXPECTED_D23_EQUATION_SIDECAR_SHA256 = "919624CFC7580899BF98BB28DCF70441127A269AA85DCEB9CCAE0E8670EB5588"
EXPECTED_D18_SUPPORT12_SHA256 = "A9DF8405D43EAFFB98C724165C0C7A647C0952D7D8D1A8630CD26700898C6F56"
EXPECTED_D18_REDUCTION_SHA256 = "01B1891AF54B9257F404840F113317FD11CDCCDA084ECC90AF7497892BC1AC73"
EXPECTED_D18_VERIFIED_SHA256 = "BE81B347A7A174809C6E962E921AE87B507AA25C363C34808AED12C0E9225125"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_equation(row: dict, rename: dict[str, str]) -> tuple:
    coefficients: dict[tuple, Fraction] = {}
    for monomial, raw_coefficient in row["terms"]:
        key = tuple(sorted((rename.get(name, name), int(power))
                           for name, power in monomial))
        coefficients[key] = coefficients.get(key, Fraction(0)) + Fraction(raw_coefficient)
    coefficients = {monomial: coefficient for monomial, coefficient in coefficients.items()
                    if coefficient}
    if not coefficients:
        return ()
    scale = coefficients[min(coefficients)]
    return tuple((monomial, str(coefficient / scale))
                 for monomial, coefficient in sorted(coefficients.items()))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("d23_equations", type=Path)
    parser.add_argument("d18_input", type=Path)
    parser.add_argument("d18_reduction", type=Path)
    parser.add_argument("d18_verified_reduction", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    d23_hash = sha256(args.d23_equations)
    d18_hash = sha256(args.d18_input).upper()
    if d23_hash.upper() != EXPECTED_D23_EQUATION_SIDECAR_SHA256:
        raise ValueError("unexpected d23 equation sidecar")
    if d18_hash != EXPECTED_D18_SUPPORT12_SHA256:
        raise ValueError("unexpected d18 support-12 input")
    d23 = json.loads(args.d23_equations.read_text(encoding="utf-8"))
    d18 = json.loads(args.d18_input.read_text(encoding="utf-8"))
    reduction = json.loads(args.d18_reduction.read_text(encoding="utf-8"))
    verified = json.loads(args.d18_verified_reduction.read_text(encoding="utf-8"))
    reduction_hash = sha256(args.d18_reduction).upper()
    verified_hash = sha256(args.d18_verified_reduction).upper()
    if reduction_hash != EXPECTED_D18_REDUCTION_SHA256:
        raise ValueError("unexpected d18 reduction report")
    if verified_hash != EXPECTED_D18_VERIFIED_SHA256:
        raise ValueError("unexpected d18 reduction receipt")

    if d23.get("source_sha256") != EXPECTED_D23_CHECKPOINT_SHA256:
        raise ValueError("d23 equation sidecar is not bound to the expected checkpoint")
    if d18.get("equation_count") != 142 or len(d18.get("equations", ())) != 142:
        raise ValueError("unexpected d18 equation count")
    if d23.get("equation_count") != 142 or len(d23.get("equations", ())) != 142:
        raise ValueError("unexpected d23 equation count")
    if (reduction.get("input_sha256") != d18_hash
            or verified.get("input_sha256") != d18_hash
            or verified.get("reduction_report_sha256") != reduction_hash
            or not verified.get(
                "all_142_equations_span_exactly_4_Q_generators")):
        raise ValueError("d18 Q-linear reduction receipt does not cover these rows")

    d18_variables = tuple(d18["selected_variables"])
    d23_variables = tuple(d23["variables"])
    if len(d18_variables) != 12 or len(d23_variables) != 12:
        raise ValueError("unexpected support dimension")
    if tuple(sorted(d18_variables)) != tuple(sorted(d23_variables)):
        # The known coordinate rename is d23 p15 -> d18 r.
        if set(d18_variables) != (set(d23_variables) - {"p15"}) | {"r"}:
            raise ValueError("support variables are not related by p15 -> r")

    d18_rows = Counter(canonical_equation(row, {}) for row in d18["equations"])
    d23_rows = Counter(canonical_equation(row, {"p15": "r"})
                        for row in d23["equations"])
    common = sum((d18_rows & d23_rows).values())
    result = {
        "d23_equation_sidecar": str(args.d23_equations),
        "d23_equation_sidecar_sha256": d23_hash,
        "d23_checkpoint_sha256": d23["source_sha256"],
        "d18_support12_input": str(args.d18_input),
        "d18_support12_input_sha256": d18_hash,
        "d18_reduction_report": str(args.d18_reduction),
        "d18_reduction_report_sha256": reduction_hash,
        "d18_verified_reduction": str(args.d18_verified_reduction),
        "d18_verified_reduction_sha256": verified_hash,
        "d23_equations": sum(d23_rows.values()),
        "d18_equations": sum(d18_rows.values()),
        "common_equations_after_scale_normalization": common,
        "same_equation_multiset_after_p15_to_r": d23_rows == d18_rows,
        "transferred_exact_Q_linear_rank": verified["independent_equations"],
        "transferred_exact_dependencies": verified["verified_exact_dependencies"],
        "scope": (
            "Exact transfer of the support-12 equation multiset and its verified Q-linear "
            "compression under p15 -> r. This is generator compression only; it does not "
            "exclude the support-12 variety, the U1=0 chart, or the 2+1+1 branch."
        ),
    }
    if not result["same_equation_multiset_after_p15_to_r"] or common != 142:
        raise ValueError("the d18 and d23 equations do not match exactly")
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
