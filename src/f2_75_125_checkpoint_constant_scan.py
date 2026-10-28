#!/usr/bin/env python3
"""Stream an F2 lift checkpoint and detect nonzero scalar constraints.

A nonzero constant among the retained necessary equations is an immediate
certificate that the represented branch is empty.  The converse is not true:
absence of a scalar equation does not prove the constraint ideal is nonempty.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re

import ijson


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--expected-sha256")
    args = parser.parse_args()

    source = args.checkpoint
    source_hash = sha256(source)
    if args.expected_sha256 and source_hash.lower() != args.expected_sha256.lower():
        raise SystemExit(f"checkpoint SHA mismatch: {source_hash}")

    with source.open("rb") as stream:
        branch = next(ijson.items(stream, "branch"))
    with source.open("rb") as stream:
        quadratic = next(ijson.items(stream, "quadratic"))
    with source.open("rb") as stream:
        generators = list(ijson.items(stream, "generators.item"))

    constants = []
    nonzero_constant_terms = []
    variables_in_constraints = set()
    constraint_count = 0
    total_terms = 0
    with source.open("rb") as stream:
        for tag, terms in ijson.items(stream, "constraints.item"):
            constraint_count += 1
            total_terms += len(terms)
            for sparse, values in terms:
                if sparse:
                    for index, exponent in sparse:
                        index = int(index)
                        if int(exponent):
                            if index >= len(generators):
                                raise SystemExit(f"constraint variable index {index} is out of range")
                            variables_in_constraints.add(generators[index])
                    continue
                coefficients = [Fraction(str(value)) for value in values]
                if any(coefficients):
                    nonzero_constant_terms.append({
                        "tag": tag,
                        "coefficient_basis": [str(value) for value in coefficients],
                    })
            if len(terms) != 1:
                continue
            sparse, values = terms[0]
            if sparse:
                continue
            coefficients = [Fraction(str(value)) for value in values]
            if any(coefficients):
                constants.append({
                    "tag": tag,
                    "coefficient_basis": [str(value) for value in coefficients],
                })

    result = {
        "status": "NONZERO_CONSTANT_CONSTRAINT_FOUND" if constants else "NO_NONZERO_CONSTANT_CONSTRAINT",
        "branch": branch,
        "quadratic_coefficient_field": bool(quadratic),
        "source": source.name,
        "source_sha256": source_hash,
        "constraint_count": constraint_count,
        "total_constraint_terms": total_terms,
        "nonzero_constant_constraints": constants,
        "nonzero_constant_term_count": len(nonzero_constant_terms),
        "nonzero_constant_term_constraints": nonzero_constant_terms,
        "variables_in_constraints": sorted(variables_in_constraints),
        "nonparameter_variables_in_constraints": sorted(
            name for name in variables_in_constraints
            if not re.fullmatch(r"p\d+", name)
        ),
        "all_ring_variables_zero_satisfies_retained_equations": not nonzero_constant_terms,
        "all_free_parameter_variables_zero_satisfies_retained_equations": (
            not nonzero_constant_terms and all(
                re.fullmatch(r"p\d+", name) for name in variables_in_constraints
            )
        ),
        "interpretation": (
            "A listed nonzero scalar equation proves the represented necessary system empty."
            if constants else
            "No single-equation scalar contradiction was found. The all-zero specialization result concerns only retained finite equations, not the existence of a Keller pair."
        ),
    }
    receipt = args.receipt or source.with_name(source.stem + "_constant_scan.json")
    receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"receipt={receipt}")


if __name__ == "__main__":
    main()
