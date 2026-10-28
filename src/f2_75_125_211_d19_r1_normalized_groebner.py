#!/usr/bin/env python3
"""Finite-field Groebner exploration of the d19 subsystem on r=1.

The substitution r=1 is a weighted-projective chart for the r != 0 locus of
the weighted-homogeneous necessary system.  A modular result is reconnaissance
only; any unit-ideal claim must be replayed over QQ.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import time

import sympy as sp
from sympy.polys.domains import GF
from sympy.polys.rings import ring as polynomial_ring
from sympy.polys.groebnertools import groebner as ring_groebner

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import f2_75_125_211_d19_support12_21_32_38_combined_groebner as builder

PRIME = 251
OUTPUT = HERE / "f2_75_125_211_d19_r1_normalized_mod251_20261028.json"


class _CapturedSystem(Exception):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def capture_system():
    captured = {}
    original = builder.ring_groebner

    def stop_after_construction(inputs, ring, **kwargs):
        captured["inputs"] = list(inputs)
        captured["ring"] = ring
        raise _CapturedSystem

    builder.QQ = GF(PRIME)
    builder.ring_groebner = stop_after_construction
    try:
        builder.main()
    except _CapturedSystem:
        pass
    finally:
        builder.ring_groebner = original
    if len(captured.get("inputs", ())) != 18:
        raise RuntimeError("failed to capture the checked 18-generator system")
    return captured["inputs"]


def main() -> None:
    start = time.monotonic()
    inputs = capture_system()
    analysis_path = HERE / "f2_75_125_211_d19_support38_analysis_20261025.json"
    analysis = json.loads(analysis_path.read_text(encoding="utf-8"))
    names = tuple(analysis["selected_variables"])
    weights = tuple(int(value) for value in analysis["primitive_positive_integer_weights"])
    r_index = names.index("r")
    reduced_names = names[:r_index] + names[r_index + 1:]
    reduced_weights = weights[:r_index] + weights[r_index + 1:]
    symbols = sp.symbols(" ".join(reduced_names))

    def weighted_grevlex(monomial):
        return (sum(w * e for w, e in zip(reduced_weights, monomial)),
                tuple(reversed([-e for e in monomial])))

    target_ring, *_ = polynomial_ring(symbols, GF(PRIME), order=weighted_grevlex)
    normalized_inputs = []
    input_term_counts = []
    for poly in inputs:
        combined = {}
        for exponents, coefficient in poly.terms():
            reduced_exponents = tuple(value for i, value in enumerate(exponents)
                                      if i != r_index)
            combined[reduced_exponents] = (
                combined.get(reduced_exponents, 0) + int(coefficient)
            ) % PRIME
        normalized = target_ring.from_dict({exponents: coefficient
                                             for exponents, coefficient in combined.items()
                                             if coefficient % PRIME})
        normalized_inputs.append(normalized)
        input_term_counts.append(len(normalized.terms()))

    print(json.dumps({"stage": "normalized_groebner_start",
                      "prime": PRIME, "chart": "r=1",
                      "variables": len(reduced_names),
                      "generators": len(normalized_inputs),
                      "generator_terms": input_term_counts}, indent=2), flush=True)
    basis = ring_groebner(normalized_inputs, target_ring, method="f5b")
    unit_ideal = any(poly == target_ring.one for poly in basis)
    basis_records = []
    for poly in basis:
        terms = poly.terms()
        basis_records.append({
            "polynomial": str(poly),
            "term_count": len(terms),
            "leading_monomial": list(terms[0][0]) if terms else None,
        })
    result = {
        "scope": "GF(251) Groebner exploration of the support-12/21/32/38 d19 necessary ideal on weighted chart r=1; not a QQ certificate or full branch exclusion",
        "branch": "2+1+1, U1=0",
        "normalization": "r=1, representing the r!=0 weighted-projective locus over an algebraic closure",
        "prime": PRIME,
        "original_variables": list(names),
        "original_weights": list(weights),
        "normalized_variables": list(reduced_names),
        "input_generators": len(normalized_inputs),
        "input_term_counts": input_term_counts,
        "groebner_method": "F5B",
        "basis_length": len(basis_records),
        "unit_ideal": unit_ideal,
        "elapsed_seconds": time.monotonic() - start,
        "source_sha256": {
            path.name: sha256_file(path) for path in (
                builder.FAMILY21, builder.ANALYSIS21, builder.FAMILY32,
                builder.ANALYSIS32, builder.FAMILY38, builder.ANALYSIS38,
                builder.SUPPORT32_REDUNDANCY, builder.SUPPORT38_LOW_REDUCTION,
                builder.SUPPORT12_GB, builder.SUPPORT12_INPUT, builder.SUPPORT12_REUSE,
            )
        },
        "basis": basis_records,
        "warning": "A modular unit ideal only identifies a lead. Reconstruct and verify an exact rational identity on r=1 before concluding that all necessary-system solutions have r=0.",
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "complete", "basis_length": len(basis_records),
                      "unit_ideal": unit_ideal,
                      "elapsed_seconds": result["elapsed_seconds"],
                      "report": str(OUTPUT)}, indent=2), flush=True)


if __name__ == "__main__":
    main()
