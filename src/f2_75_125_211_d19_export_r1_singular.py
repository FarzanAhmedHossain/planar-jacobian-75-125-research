#!/usr/bin/env python3
"""Export the exact d19 support ideal to Singular on the chart r=1."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp
from sympy.polys.rings import ring as polynomial_ring

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import f2_75_125_211_d19_support12_21_32_38_combined_groebner as builder

PRIME = 251
OUT_QQ = HERE / "f2_75_125_211_d19_r1_normalized_QQ_20261028.sing"
OUT_MOD = HERE / "f2_75_125_211_d19_r1_normalized_mod251_20261028.sing"
OUT_MANIFEST = HERE / "f2_75_125_211_d19_r1_normalized_singular_manifest_20261028.json"


class _CapturedSystem(Exception):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def capture_exact_system():
    captured = {}
    original = builder.ring_groebner

    def stop_after_construction(inputs, ring, **kwargs):
        captured["inputs"] = list(inputs)
        captured["ring"] = ring
        raise _CapturedSystem

    builder.QQ = sp.QQ
    builder.ring_groebner = stop_after_construction
    try:
        builder.main()
    except _CapturedSystem:
        pass
    finally:
        builder.ring_groebner = original
    if len(captured.get("inputs", ())) != 18:
        raise RuntimeError("failed to capture the checked 18-generator QQ system")
    return captured["inputs"]


def qstr(value) -> str:
    rational = sp.Rational(value)
    numerator, denominator = int(rational.p), int(rational.q)
    if denominator == 1:
        return str(numerator)
    return f"({numerator}/{denominator})"


def modq(value) -> int:
    rational = sp.Rational(value)
    denominator = int(rational.q) % PRIME
    if denominator == 0:
        raise ValueError("GF(251) divides an input coefficient denominator")
    return (int(rational.p) * pow(denominator, -1, PRIME)) % PRIME


def coefficient_text(value, modular: bool) -> str:
    if modular:
        return str(modq(value))
    return qstr(value)


def monomial_text(names, exponents) -> str:
    factors = []
    for name, exponent in zip(names, exponents):
        if exponent == 1:
            factors.append(name)
        elif exponent > 1:
            factors.append(f"{name}^{exponent}")
    return "*".join(factors) if factors else "1"


def polynomial_text(terms, names, modular: bool) -> str:
    parts = []
    for exponents, coefficient in terms:
        if modular:
            residue = modq(coefficient)
            if residue == 0:
                continue
            sign = "+" if residue <= PRIME // 2 else "-"
            magnitude = residue if residue <= PRIME // 2 else PRIME - residue
            coefficient_string = str(magnitude)
        else:
            rational = sp.Rational(coefficient)
            if rational == 0:
                continue
            sign = "+" if rational > 0 else "-"
            magnitude = abs(rational)
            coefficient_string = qstr(magnitude)

        monomial = monomial_text(names, exponents)
        if coefficient_string == "1" and monomial != "1":
            term = monomial
        else:
            term = coefficient_string if monomial == "1" else f"{coefficient_string}*{monomial}"
        parts.append((sign, term))

    if not parts:
        return "0"
    sign, term = parts[0]
    out = ("" if sign == "+" else "-") + term
    for sign, term in parts[1:]:
        out += sign + term
    return f"({out})"


def build_singular_source(inputs, all_names, modular: bool) -> tuple[str, list[int], list[str]]:
    r_index = all_names.index("r")
    names = tuple(name for i, name in enumerate(all_names) if i != r_index)
    source_polys = []
    term_counts = []
    constant_coefficients = []
    for poly in inputs:
        combined = {}
        for exponents, coefficient in poly.terms():
            reduced_exponents = tuple(power for i, power in enumerate(exponents)
                                      if i != r_index)
            if modular:
                combined[reduced_exponents] = (
                    combined.get(reduced_exponents, 0) + modq(coefficient)
                ) % PRIME
            else:
                combined[reduced_exponents] = (
                    combined.get(reduced_exponents, sp.Rational(0)) + sp.Rational(coefficient)
                )
        if modular:
            combined = {exponents: coefficient for exponents, coefficient in combined.items()
                        if coefficient % PRIME}
        else:
            combined = {exponents: coefficient for exponents, coefficient in combined.items()
                        if coefficient}
        term_counts.append(len(combined))
        zero_exp = (0,) * len(names)
        constant_coefficients.append(str(combined.get(zero_exp, 0)))
        source_polys.append(polynomial_text(list(combined.items()), names, modular))

    characteristic = str(PRIME) if modular else "0"
    ring_order = "dp"
    header = f"ring R={characteristic},({','.join(names)}),{ring_order};\n"
    ideal = "ideal I=" + ",\n".join(source_polys) + ";\n"
    analysis = (
        "option(redSB);\n"
        "ideal G=std(I);\n"
        "print(\"BASIS_SIZE=\"+string(size(G)));\n"
        "print(\"REDUCE_ONE=\"+string(reduce(1,G)));\n"
        "quit;\n"
    )
    return header + ideal + analysis, term_counts, constant_coefficients


def main() -> None:
    inputs = capture_exact_system()
    analysis_path = HERE / "f2_75_125_211_d19_support38_analysis_20261025.json"
    analysis = json.loads(analysis_path.read_text(encoding="utf-8"))
    names = tuple(analysis["selected_variables"])
    weights = tuple(int(value) for value in analysis["primitive_positive_integer_weights"])
    if len(names) != 28 or len(weights) != 28 or "r" not in names:
        raise AssertionError("source variable/weight metadata mismatch")

    qq_source, qq_counts, qq_constants = build_singular_source(inputs, names, False)
    mod_source, mod_counts, mod_constants = build_singular_source(inputs, names, True)
    OUT_QQ.write_text(qq_source, encoding="utf-8")
    OUT_MOD.write_text(mod_source, encoding="utf-8")

    source_paths = [builder.FAMILY21, builder.ANALYSIS21, builder.FAMILY32,
                    builder.ANALYSIS32, builder.FAMILY38, builder.ANALYSIS38,
                    builder.SUPPORT32_REDUNDANCY, builder.SUPPORT38_LOW_REDUCTION,
                    builder.SUPPORT12_GB, builder.SUPPORT12_INPUT, builder.SUPPORT12_REUSE]
    manifest = {
        "scope": "Exact substitution r=1 in the checked 18-generator d19 support-12/21/32/38 necessary ideal, exported to Singular over QQ and GF(251)",
        "branch": "2+1+1, U1=0; chart r!=0 normalized by weighted scaling",
        "variables_before_normalization": list(names),
        "weights_before_normalization": list(weights),
        "normalized_variables": [name for name in names if name != "r"],
        "generator_count": len(inputs),
        "qq_term_counts": qq_counts,
        "mod251_term_counts": mod_counts,
        "qq_constant_coefficients": qq_constants,
        "mod251_constant_coefficients": mod_constants,
        "prime": PRIME,
        "sources_sha256": {path.name: sha256_file(path) for path in source_paths},
        "qq_script": OUT_QQ.name,
        "qq_script_sha256": sha256_file(OUT_QQ),
        "mod251_script": OUT_MOD.name,
        "mod251_script_sha256": sha256_file(OUT_MOD),
        "warning": "A Singular unit-ideal result over GF(251) is only reconnaissance. Exact QQ verification and scope-matched handling of r=0 are required for branch exclusion.",
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "manifest": str(OUT_MANIFEST),
        "qq_script": str(OUT_QQ),
        "qq_sha256": manifest["qq_script_sha256"],
        "mod251_script": str(OUT_MOD),
        "mod251_sha256": manifest["mod251_script_sha256"],
        "qq_constants": qq_constants,
        "mod251_constants": mod_constants,
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
