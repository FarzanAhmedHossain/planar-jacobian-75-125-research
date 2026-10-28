#!/usr/bin/env python3
"""Exact gcd screen for one-variable equations in a saved 2+1+1 lift.

If the retained necessary equations include two univariate polynomials in a
single carried parameter with gcd 1 over Q, their common zero set is empty.
The output records an exact Bezout combination for replay.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from f2_75_125_quadratic_affine_profile_stream import (
    iter_json_array,
    top_metadata_and_constraints,
)


def trim(a: list[Fraction]) -> list[Fraction]:
    while len(a) > 1 and not a[-1]:
        a.pop()
    return a


def add(a: list[Fraction], b: list[Fraction]) -> list[Fraction]:
    out = [Fraction(0)] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += x
    return trim(out)


def neg(a: list[Fraction]) -> list[Fraction]:
    return [-x for x in a]


def sub(a: list[Fraction], b: list[Fraction]) -> list[Fraction]:
    return add(a, neg(b))


def mul(a: list[Fraction], b: list[Fraction]) -> list[Fraction]:
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def scale(a: list[Fraction], c: Fraction) -> list[Fraction]:
    return trim([c * x for x in a])


def divmod_poly(a: list[Fraction], b: list[Fraction]):
    if len(b) == 1 and not b[0]:
        raise ZeroDivisionError
    r = trim(a.copy())
    q = [Fraction(0)] * max(1, len(a) - len(b) + 1)
    while any(r) and len(r) >= len(b):
        degree = len(r) - len(b)
        coefficient = r[-1] / b[-1]
        q[degree] += coefficient
        for i, value in enumerate(b):
            r[degree + i] -= coefficient * value
        trim(r)
    return trim(q), trim(r)


def xgcd(a: list[Fraction], b: list[Fraction]):
    old_r, r = trim(a.copy()), trim(b.copy())
    old_s, s = [Fraction(1)], [Fraction(0)]
    old_t, t = [Fraction(0)], [Fraction(1)]
    while any(r):
        q, rem = divmod_poly(old_r, r)
        old_r, r = r, rem
        old_s, s = s, sub(old_s, mul(q, s))
        old_t, t = t, sub(old_t, mul(q, t))
    leading = old_r[-1]
    return (scale(old_r, 1 / leading), scale(old_s, 1 / leading),
            scale(old_t, 1 / leading))


def encode_poly(poly: list[Fraction]) -> list[str]:
    return [str(value) for value in trim(poly.copy())]


def source_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def analyze(source: Path, variable_name: str) -> dict:
    text = source.read_text(encoding="utf-8")
    metadata, constraint_pos = top_metadata_and_constraints(text)
    if metadata.get("format") != "f2-sparse-ring-v1":
        raise ValueError("unexpected checkpoint format")
    if metadata.get("branch") != "2+1+1" or metadata.get("quadratic"):
        raise ValueError("expected a rational 2+1+1 checkpoint")
    generators = metadata["generators"]
    target_index = generators.index(variable_name)
    equations: list[tuple[str, list[Fraction]]] = []
    for tag, terms in iter_json_array(text, constraint_pos):
        support = set()
        coefficients: dict[int, Fraction] = {}
        for sparse, values in terms:
            if len(values) != 1:
                raise ValueError("expected rational coefficients")
            coefficient = Fraction(values[0])
            power = 0
            for generator, exponent in sparse:
                support.add(generator)
                if generator == target_index:
                    power += exponent
            if support - {target_index}:
                break
            coefficients[power] = coefficients.get(power, Fraction(0)) + coefficient
        if support == {target_index}:
            equation = [Fraction(0)] * (max(coefficients, default=0) + 1)
            for exponent, coefficient in coefficients.items():
                equation[exponent] += coefficient
            equations.append((tag, trim(equation)))

    if not equations:
        gcd_poly = [Fraction(0)]
        bezout_coefficients: list[list[Fraction]] = []
    else:
        gcd_poly = equations[0][1]
        bezout_coefficients = [[Fraction(1)]]
        for _tag, equation in equations[1:]:
            gcd_poly, left, right = xgcd(gcd_poly, equation)
            bezout_coefficients = [mul(left, c) for c in bezout_coefficients]
            bezout_coefficients.append(right)
        if len(gcd_poly) == 1:
            normalization = 1 / gcd_poly[0]
            gcd_poly = [Fraction(1)]
            bezout_coefficients = [scale(c, normalization)
                                   for c in bezout_coefficients]

    identity = [Fraction(0)]
    for coefficient, (_tag, equation) in zip(bezout_coefficients, equations):
        identity = add(identity, mul(coefficient, equation))
    if equations and identity != gcd_poly:
        raise AssertionError("Bezout identity failed exact replay")

    return {
        "source": str(source),
        "source_sha256": source_sha256(source),
        "branch": metadata["branch"],
        "checkpoint_format": metadata["format"],
        "nparameters": int(metadata["nparameters"]),
        "variable": variable_name,
        "variable_generator_index": target_index,
        "univariate_constraint_count": len(equations),
        "univariate_constraints": [
            {"tag": tag, "coefficients_low_to_high": encode_poly(poly)}
            for tag, poly in equations
        ],
        "gcd_coefficients_low_to_high": encode_poly(gcd_poly),
        "bezout_coefficients_low_to_high": [
            encode_poly(poly) for poly in bezout_coefficients
        ],
        "bezout_identity_replayed": identity == gcd_poly,
        "status": (
            "univariate_inconsistent" if equations and len(gcd_poly) == 1
            else "nonconstant_common_factor" if equations
            else "no_univariate_constraints"
        ),
        "scope": (
            "An inconsistent gcd proves the retained necessary system empty on this "
            "checkpoint chart; branch completeness still depends on the pinned "
            "root-profile reduction and F2 normal-form bridge."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("variable", help="for example p79")
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    result = analyze(args.source, args.variable)
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        key: result[key] for key in (
            "source_sha256", "branch", "variable", "univariate_constraint_count",
            "gcd_coefficients_low_to_high", "bezout_identity_replayed", "status", "scope"
        )
    }, indent=2))
    print(f"report={args.report}")


if __name__ == "__main__":
    main()
