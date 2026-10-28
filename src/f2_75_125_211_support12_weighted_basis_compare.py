#!/usr/bin/env python3
"""Compare the exact weighted QQ Groebner bases from Singular and FLINT."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from fractions import Fraction
from pathlib import Path


VARIABLES = (
    "p24", "p25", "p26", "p27", "p28", "p44", "p45", "p46",
    "p47", "p48", "p5", "r",
)
VARIABLE_RE = re.compile(r"([A-Za-z][A-Za-z0-9]*)(?:(?:\*\*(\d+))|(?:\^(\d+)))?")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_polynomial(source: str) -> dict[tuple[int, ...], Fraction]:
    """Parse the integer/rational coefficient syntax emitted by the two CASes."""
    result: dict[tuple[int, ...], Fraction] = {}
    for term in re.findall(r"[+-]?[^+-]+", source.replace(" ", "")):
        sign = -1 if term.startswith("-") else 1
        body = term.lstrip("+-")
        matches = list(VARIABLE_RE.finditer(body))
        exponents = [0] * len(VARIABLES)
        for match in matches:
            name = match.group(1)
            if name not in VARIABLES:
                raise ValueError(f"unknown variable {name!r} in term {term!r}")
            exponent = int(match.group(2) or match.group(3) or 1)
            exponents[VARIABLES.index(name)] += exponent
        coefficient_text = VARIABLE_RE.sub("", body).replace("*", "")
        coefficient = Fraction(coefficient_text) if coefficient_text else Fraction(1)
        monomial = tuple(exponents)
        result[monomial] = result.get(monomial, Fraction(0)) + sign * coefficient
    return {monomial: coefficient for monomial, coefficient in result.items() if coefficient}


def normalized(polynomial: dict[tuple[int, ...], Fraction]):
    if not polynomial:
        raise ValueError("zero polynomial in purported Groebner basis")
    leading_monomial = max(polynomial)
    leading_coefficient = polynomial[leading_monomial]
    return (
        {monomial: coefficient / leading_coefficient
         for monomial, coefficient in polynomial.items()},
        leading_monomial,
        leading_coefficient,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--research-dir", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    root = args.research_dir
    stem = "f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12"
    singular_path = root / f"{stem}_4gen_weighted.basis.txt"
    flint_path = root / f"{stem}_gb_4gen_weighted_flint.json"
    input_path = root / f"{stem}_input.json"
    reduction_path = root / f"{stem}_linear_reduce.json"
    receipt_path = root / f"{stem}_weighted_basis_compare.json"

    flint_data = json.loads(flint_path.read_text(encoding="utf-8"))
    singular_sources = singular_path.read_text(encoding="utf-8").split(",")
    singular_basis = [parse_polynomial(source) for source in singular_sources]
    flint_basis = [parse_polynomial(item["polynomial"]) for item in flint_data["basis"]]
    if len(singular_basis) != len(flint_basis):
        raise ValueError("CAS basis lengths differ")

    matches = []
    used_flint = set()
    for singular_index, singular_poly in enumerate(singular_basis):
        singular_norm, singular_lm, singular_lc = normalized(singular_poly)
        candidates = []
        for flint_index, flint_poly in enumerate(flint_basis):
            flint_norm, flint_lm, flint_lc = normalized(flint_poly)
            if singular_norm == flint_norm:
                candidates.append((flint_index, flint_lm, singular_lc / flint_lc))
        if len(candidates) != 1:
            raise ValueError(
                f"Singular basis element {singular_index} has {len(candidates)} exact matches"
            )
        flint_index, flint_lm, scalar = candidates[0]
        if flint_index in used_flint:
            raise ValueError(f"duplicate match to FLINT basis element {flint_index}")
        used_flint.add(flint_index)
        matches.append({
            "singular_index": singular_index,
            "flint_index": flint_index,
            "terms": len(singular_poly),
            "leading_monomial": list(singular_lm),
            "singular_equals_scalar_times_flint": str(scalar),
        })
    if used_flint != set(range(len(flint_basis))):
        raise ValueError("the matches are not a bijection")

    input_data = json.loads(input_path.read_text(encoding="utf-8"))
    reduction_data = json.loads(reduction_path.read_text(encoding="utf-8"))
    receipt = {
        "scope": "Exact polynomial-by-polynomial comparison up to nonzero QQ scalars; this verifies agreement of the two CAS outputs, not branch exclusion.",
        "field": "QQ",
        "monomial_order": flint_data["monomial_order"],
        "weights": flint_data["monomial_weights"],
        "basis_length": len(singular_basis),
        "singular_term_counts": [len(poly) for poly in singular_basis],
        "flint_term_counts": [len(poly) for poly in flint_basis],
        "matches": matches,
        "singular_basis_sha256": sha256(singular_path),
        "flint_report_sha256": sha256(flint_path),
        "support12_input_sha256": sha256(input_path),
        "linear_reduction_sha256": sha256(reduction_path),
        "flint_unit_ideal": flint_data["unit_ideal"],
        "origin_witness": flint_data["origin_witness"],
    }
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
