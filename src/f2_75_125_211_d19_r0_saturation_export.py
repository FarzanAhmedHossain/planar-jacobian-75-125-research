#!/usr/bin/env python3
"""Export the exact d19 necessary ideal on the weighted-projective r=0 slice.

The resulting Singular inputs test whether the r=0 slice has any nonzero
solution by saturating against the irrelevant ideal.  A QQ unit ideal is an
exclusion only for this necessary subsystem and this slice; it is not, by
itself, an exclusion of the full 2+1+1 branch.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import f2_75_125_211_d19_export_r1_singular as exporter
import f2_75_125_211_d19_support12_21_32_38_combined_groebner as builder

PRIME = exporter.PRIME
OUT_QQ = HERE / "f2_75_125_211_d19_r0_saturated_QQ_20261028.sing"
OUT_MOD = HERE / "f2_75_125_211_d19_r0_saturated_mod251_20261028.sing"
OUT_MANIFEST = HERE / "f2_75_125_211_d19_r0_saturated_manifest_20261028.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def specialize_r0(inputs, names: tuple[str, ...]):
    """Substitute r=0 exactly, dropping terms with positive r-exponent."""
    r_index = names.index("r")
    variables = tuple(name for i, name in enumerate(names) if i != r_index)
    weights_data = json.loads(
        (HERE / "f2_75_125_211_d19_support38_analysis_20261025.json")
        .read_text(encoding="utf-8")
    )
    full_weights = tuple(
        int(value) for value in weights_data["primitive_positive_integer_weights"]
    )
    if len(full_weights) != len(names):
        raise AssertionError("weight vector does not align with the source variables")
    weights = tuple(value for i, value in enumerate(full_weights) if i != r_index)

    polynomials = []
    term_counts = []
    generator_degrees = []
    for source in inputs:
        reduced = []
        weighted_degrees = set()
        for exponents, coefficient in source.terms():
            if exponents[r_index] != 0:
                continue
            reduced_exponents = tuple(
                value for i, value in enumerate(exponents) if i != r_index
            )
            weighted_degrees.add(
                sum(weight * exponent for weight, exponent in zip(weights, reduced_exponents))
            )
            reduced.append((reduced_exponents, coefficient))
        if len(weighted_degrees) > 1:
            raise AssertionError("r=0 specialization broke weighted homogeneity")
        generator_degrees.append(next(iter(weighted_degrees), None))
        term_counts.append(len(reduced))
        polynomials.append(reduced)
    return variables, weights, polynomials, term_counts, generator_degrees


def build_singular(polynomials, names, weights, modular: bool) -> str:
    prime = PRIME if modular else 0
    ring_order = f"wp({','.join(map(str, weights))})"
    header = f"ring R={prime},({','.join(names)}),{ring_order};\n"
    source_polynomials = [
        exporter.polynomial_text(terms, names, modular) for terms in polynomials
    ]
    ideal_text = "ideal I=" + ",\n".join(source_polynomials) + ";\n"
    irrelevant_text = "ideal irrelevant=" + ",".join(names) + ";\n"
    body = (
        "option(redSB);\n"
        "ideal G=std(I);\n"
        "list satResult=sat(G,irrelevant);\n"
        "ideal nonzeroLocus=std(satResult[1]);\n"
        'print("SAT_BASIS_SIZE="+string(size(nonzeroLocus)));\n'
        'print("SAT_REDUCE_ONE="+string(reduce(1,nonzeroLocus)));\n'
        "quit;\n"
    )
    return header + ideal_text + irrelevant_text + body


def main() -> None:
    inputs = exporter.capture_exact_system()
    analysis = json.loads(
        (HERE / "f2_75_125_211_d19_support38_analysis_20261025.json")
        .read_text(encoding="utf-8")
    )
    names = tuple(analysis["selected_variables"])
    if len(inputs) != 18 or len(names) != 28 or "r" not in names:
        raise AssertionError("source d19 system metadata mismatch")

    variables, weights, polynomials, term_counts, generator_degrees = specialize_r0(
        inputs, names
    )
    # Keep separate exact-QQ and modular rings.  Reusing the same QQ source
    # terms is safe because conversion to GF(251) is checked denominator by
    # denominator by exporter.modq.
    qq_source = build_singular(polynomials, variables, weights, modular=False)
    mod_source = build_singular(polynomials, variables, weights, modular=True)
    OUT_QQ.write_text(qq_source, encoding="utf-8")
    OUT_MOD.write_text(mod_source, encoding="utf-8")

    source_paths = [
        builder.FAMILY21, builder.ANALYSIS21, builder.FAMILY32,
        builder.ANALYSIS32, builder.FAMILY38, builder.ANALYSIS38,
        builder.SUPPORT32_REDUNDANCY, builder.SUPPORT38_LOW_REDUCTION,
        builder.SUPPORT12_GB, builder.SUPPORT12_INPUT, builder.SUPPORT12_REUSE,
    ]
    manifest = {
        "scope": (
            "Exact r=0 weighted-projective slice of the 18-generator d19 "
            "support-12/21/32/38 necessary ideal; saturation by all remaining "
            "variables tests for nonzero solutions."
        ),
        "branch": "2+1+1, U1=0; r=0 slice",
        "variables_before_specialization": list(names),
        "weights_before_specialization": list(analysis["primitive_positive_integer_weights"]),
        "variables_after_specialization": list(variables),
        "weights_after_specialization": list(weights),
        "generator_count": len(inputs),
        "generator_term_counts_after_r0": term_counts,
        "generator_weight_degrees_after_r0": generator_degrees,
        "prime": PRIME,
        "sources_sha256": {path.name: sha256_file(path) for path in source_paths},
        "qq_script": OUT_QQ.name,
        "qq_script_sha256": sha256_file(OUT_QQ),
        "mod251_script": OUT_MOD.name,
        "mod251_script_sha256": sha256_file(OUT_MOD),
        "warning": (
            "A QQ unit result excludes nonzero r=0 solutions only for this "
            "necessary subsystem; the full support/source-attachment bridge "
            "and other branch coverage remain required."
        ),
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "manifest": str(OUT_MANIFEST),
        "qq_script": str(OUT_QQ),
        "mod251_script": str(OUT_MOD),
        "variables_after_specialization": len(variables),
        "generators": len(inputs),
        "term_counts": term_counts,
        "generator_weight_degrees": generator_degrees,
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
