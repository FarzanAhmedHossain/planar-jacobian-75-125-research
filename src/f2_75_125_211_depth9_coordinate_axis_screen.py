#!/usr/bin/env python3
"""Exact common-zero screen on the 18 coordinate axes of the 2+1+1 ideal.

Restrict all 224 rational equations to each affine coordinate axis, compute
their exact gcd in Q[z], and remove the origin factor. This is a finite
coordinate-stratum check only; it is not a full branch exclusion.
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
from collections import defaultdict, Counter
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_SOURCE = HERE / "f2_75_125_211_depth9_nonzero_locus_sat_20261021.sing"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def parse_fraction(text: str) -> Fraction:
    text = text.strip().strip("()")
    return Fraction(text) if text else Fraction(1)


def parse_terms(expression: str):
    expression = expression.strip().rstrip(",").strip()
    pieces = []
    start = 0
    depth = 0
    for index, character in enumerate(expression):
        if character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
        elif character in "+-" and depth == 0 and index > start:
            pieces.append(expression[start:index].strip())
            start = index
    pieces.append(expression[start:].strip())

    for term in pieces:
        sign = 1
        if term.startswith("-"):
            sign = -1
            term = term[1:].strip()
        elif term.startswith("+"):
            term = term[1:].strip()
        factors = re.findall(r"p\d+(?:\^\d+)?", term)
        if not factors:
            coefficient = parse_fraction(term) * sign
            exponents = {}
        else:
            coefficient_text = term[: term.find(factors[0])].rstrip("*")
            coefficient = parse_fraction(coefficient_text) * sign
            exponents = {}
            for factor in factors:
                match = re.fullmatch(r"(p\d+)(?:\^(\d+))?", factor)
                assert match
                name, power = match.groups()
                exponents[name] = exponents.get(name, 0) + int(power or 1)
        yield coefficient, exponents


def trim(poly: list[Fraction]) -> list[Fraction]:
    while poly and poly[-1] == 0:
        poly.pop()
    return poly


def monic(poly: list[Fraction]) -> list[Fraction]:
    poly = trim(poly[:])
    if not poly:
        return []
    lead = poly[-1]
    return [coefficient / lead for coefficient in poly]


def divrem(a: list[Fraction], b: list[Fraction]):
    remainder = trim(a[:])
    b = trim(b[:])
    if not b:
        raise ZeroDivisionError
    quotient = [Fraction(0)] * max(0, len(remainder) - len(b) + 1)
    while len(remainder) >= len(b) and remainder:
        shift = len(remainder) - len(b)
        factor = remainder[-1] / b[-1]
        quotient[shift] += factor
        for i, coefficient in enumerate(b):
            remainder[shift + i] -= factor * coefficient
        trim(remainder)
    return trim(quotient), remainder


def gcd_poly(a: list[Fraction], b: list[Fraction]) -> list[Fraction]:
    a, b = trim(a[:]), trim(b[:])
    while b:
        _, remainder = divrem(a, b)
        a, b = b, remainder
    return monic(a)


def encode_poly(poly: list[Fraction]) -> list[str]:
    return [str(coefficient) for coefficient in poly]


def main() -> None:
    source = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else DEFAULT_SOURCE
    output = source.with_name(source.stem + "_coordinate_axis_screen.json")
    source_hash = sha256_file(source)
    variables: list[str] = []
    axis_rows = defaultdict(lambda: defaultdict(lambda: defaultdict(Fraction)))
    max_degrees = Counter()
    row_count = 0
    term_count = 0
    in_ideal = False

    with source.open(encoding="utf-8") as stream:
        for raw_line in stream:
            line = raw_line.strip()
            if not line:
                continue
            if not variables:
                match = re.search(r"ring r=0,\(([^)]*)\),dp;", line)
                if match:
                    variables = match.group(1).split(",")
            if line == "ideal I =":
                in_ideal = True
                continue
            if in_ideal and (
                line == ";" or line.startswith("print(") or line.startswith("ideal G=")
            ):
                break
            if not in_ideal:
                continue

            expression = line.rstrip(",;").strip()
            if not expression:
                continue
            row = row_count
            row_count += 1
            row_degrees = set()
            has_constant = False
            for coefficient, powers in parse_terms(expression):
                term_count += 1
                row_degrees.add(sum(powers.values()))
                if not powers:
                    has_constant = coefficient != 0
                    continue
                if len(powers) == 1:
                    name, exponent = next(iter(powers.items()))
                    axis_rows[name][row][exponent] += coefficient
            if has_constant:
                raise AssertionError(f"origin is not a zero at equation {row}")
            if not row_degrees:
                raise AssertionError(f"empty equation at row {row}")
            max_degrees[max(row_degrees)] += 1

    if not variables or not row_count or not term_count:
        raise AssertionError((len(variables), row_count, term_count))

    axis_results = {}
    for name in variables:
        common_gcd = None
        nonzero_restrictions = 0
        for row in range(row_count):
            by_degree = axis_rows[name].get(row, {})
            if not by_degree:
                continue
            polynomial = [by_degree.get(power, Fraction(0))
                          for power in range(max(by_degree) + 1)]
            trim(polynomial)
            if not polynomial:
                continue
            nonzero_restrictions += 1
            common_gcd = polynomial if common_gcd is None else gcd_poly(common_gcd, polynomial)
            if len(common_gcd) == 1:
                # The gcd is a unit, but continue to count all nonzero restrictions.
                continue

        if common_gcd is None:
            axis_results[name] = {
                "status": "entire_axis_survives_truncation",
                "nonzero_restricted_equations": 0,
            }
            continue

        gcd_value = monic(common_gcd)
        nonzero_root_factor = gcd_value[:]
        while len(nonzero_root_factor) > 1 and nonzero_root_factor[0] == 0:
            nonzero_root_factor = nonzero_root_factor[1:]
        survives = len(nonzero_root_factor) > 1
        axis_results[name] = {
            "status": "possible_nonzero_axis_points" if survives else "axis_excluded",
            "nonzero_restricted_equations": nonzero_restrictions,
            "common_gcd_ascending_coefficients_Q": encode_poly(gcd_value),
            "nonzero_root_factor_ascending_coefficients_Q": encode_poly(nonzero_root_factor),
        }

    result = {
        "source": source.name,
        "source_sha256": source_hash,
        "branch": "2+1+1",
        "ideal_equations": row_count,
        "monomials": term_count,
        "variables": variables,
        "ordinary_max_degree_profile": dict(sorted((str(k), v) for k, v in max_degrees.items())),
        "scope": "exact common zeros on the 18 affine coordinate axes only; not a full branch exclusion",
        "axis_results": axis_results,
        "axis_status_counts": {
            status: sum(value["status"] == status for value in axis_results.values())
            for status in sorted({value["status"] for value in axis_results.values()})
        },
    }
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "source_sha256": source_hash,
        "equations": row_count,
        "monomials": term_count,
        "ordinary_max_degree_profile": result["ordinary_max_degree_profile"],
        "axis_status_counts": result["axis_status_counts"],
        "axes": axis_results,
        "receipt": output.name,
    }, indent=2))


if __name__ == "__main__":
    main()
