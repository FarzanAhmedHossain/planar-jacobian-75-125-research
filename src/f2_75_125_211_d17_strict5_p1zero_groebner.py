#!/usr/bin/env python3
"""Exact p1=0 slice of the strict five-parameter d17 2+1+1 subsystem.

This is a necessary-subsystem diagnostic. It streams the immutable d17
checkpoint, applies the already-certified p79=0 specialization, retains only
complete equations supported on {p1,p7,p15,p16,p17}, then sets p1=0 and
computes an exact Groebner basis over QQ in the remaining four variables.
"""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

from f2_75_125_211_streaming_checkpoint_scan import JsonStream
from sympy.polys.domains import QQ
from sympy.polys.groebnertools import groebner as ring_groebner
from sympy.polys.rings import ring as polynomial_ring


EXPECTED_SOURCE_SHA256 = (
    "c205df23dc70fb51df2142c482c3145bc7fc25c6be95f28131818c9c06ca3812"
)
TARGET = ("p1", "p7", "p15", "p16", "p17")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    source = Path("research/f2_75_125_sparse_ring_lift_211_roots17_dd_20261024_d17_p79zero_snapshot.json")
    report_path = Path("research/f2_75_125_sparse_ring_lift_211_roots17_dd_20261024_d17_strict5_p1zero_gb.json")
    source_hash = sha256(source)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise ValueError(f"d17 checkpoint SHA-256 mismatch: {source_hash}")

    reader = JsonStream(source)
    try:
        metadata, _ = reader.top_level_constraints()
        if metadata.get("format") != "f2-sparse-ring-v1" or metadata.get("branch") != "2+1+1":
            raise ValueError("expected the rational d17 2+1+1 checkpoint")
        generators = metadata["generators"]
        target_ids = {generators.index(name): i for i, name in enumerate(TARGET)}
        p79_id = generators.index("p79")
        equations: dict[tuple, tuple[str, object]] = {}
        scanned_rows = 0
        selected_rows = 0
        selected_terms = 0
        p1zero_terms = 0

        for tag, sparse_terms in reader.array_values():
            scanned_rows += 1
            terms: dict[tuple[int, ...], Fraction] = defaultdict(Fraction)
            outside = False
            for sparse, values in sparse_terms:
                if len(values) != 1:
                    raise ValueError("expected rational checkpoint coefficients")
                monomial = [0] * len(TARGET)
                for generator_id, exponent in sparse:
                    exponent = int(exponent)
                    if not exponent:
                        continue
                    if generator_id == p79_id:
                        break  # valid specialization p79=0 kills this monomial
                    target_index = target_ids.get(generator_id)
                    if target_index is None:
                        outside = True
                        break
                    monomial[target_index] = exponent
                else:
                    terms[tuple(monomial)] += Fraction(values[0])
                    continue
                if outside:
                    break
            if outside:
                continue
            terms = {m: c for m, c in terms.items() if c}
            if not terms:
                continue
            support = {i for i in range(len(TARGET)) if any(m[i] for m in terms)}
            if len(support) != len(TARGET):
                continue
            selected_rows += 1
            selected_terms += len(terms)

            # Apply p1=0 by retaining exactly the terms with p1 exponent zero.
            sliced: dict[tuple[int, ...], Fraction] = defaultdict(Fraction)
            for monomial, coefficient in terms.items():
                if monomial[0] == 0:
                    sliced[monomial[1:]] += coefficient
            sliced = {m: c for m, c in sliced.items() if c}
            if not sliced:
                continue
            p1zero_terms += len(sliced)
            # Nonzero scalar duplicates are removed after constructing QQ-ring elements.
            key = tuple(sorted(sliced.items()))
            equations.setdefault(key, (tag, sliced))
    finally:
        reader.close()

    names = TARGET[1:]
    R, *gens = polynomial_ring(names, QQ, order="grevlex")
    unique = {}
    for tag, terms in equations.values():
        encoded = {
            monomial: QQ.convert(coefficient.numerator) / QQ.convert(coefficient.denominator)
            for monomial, coefficient in terms.items()
        }
        polynomial = R.from_dict(encoded)
        if polynomial:
            monic = polynomial.monic()
            unique.setdefault(tuple(sorted(monic.terms())), (tag, monic))
    inputs = [item[1] for item in unique.values()]
    basis = ring_groebner(inputs, R, method="f5b") if inputs else []
    unit = any(poly == R.one for poly in basis)
    result = {
        "source": str(source),
        "source_sha256": source_hash,
        "source_branch": metadata["branch"],
        "forced_specialization": {"p79": "0", "p1": "0"},
        "retention_semantics": (
            "After p79=0, retain only whole equations supported in "
            "QQ[p1,p7,p15,p16,p17], then quotient by p1=0."
        ),
        "scanned_checkpoint_rows": scanned_rows,
        "selected_complete_five_variable_rows": selected_rows,
        "selected_sparse_terms": selected_terms,
        "nonzero_terms_after_p1_zero": p1zero_terms,
        "unique_equations_up_to_rational_scale": len(inputs),
        "variables": list(names),
        "groebner_method": "F5B",
        "monomial_order": "grevlex",
        "groebner_basis_length": len(basis),
        "groebner_basis": [str(poly) for poly in basis],
        "unit_ideal": unit,
        "scope": (
            "An exact unit ideal proves this p1=0 slice of the necessary "
            "five-variable subsystem is empty. A nonunit basis does not "
            "exclude the slice or the full 2+1+1 branch."
        ),
    }
    report_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "source_sha256", "scanned_checkpoint_rows",
        "selected_complete_five_variable_rows", "selected_sparse_terms",
        "nonzero_terms_after_p1_zero", "unique_equations_up_to_rational_scale",
        "variables", "groebner_basis_length", "unit_ideal", "scope"
    )}, indent=2), flush=True)
    print(f"report={report_path}", flush=True)


if __name__ == "__main__":
    main()
