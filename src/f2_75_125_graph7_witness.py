"""Exact point witness for the 330-equation F3 depth-nine constraint system."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import f2_75_125_quadratic_affine_profile as core
from f2_75_125_quadratic_affine_profile_stream import (
    iter_json_array,
    top_metadata_and_constraints,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def monomial_value(monomial, assignment):
    value = core.ONE
    for index, exponent in monomial:
        for _ in range(exponent):
            value = core.mul(value, assignment[index])
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("linear_receipt", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")

    source_hash = sha256(args.checkpoint)
    linear = json.loads(args.linear_receipt.read_text(encoding="utf-8"))
    if linear.get("source_sha256") != source_hash:
        raise ValueError("linear receipt does not match source checkpoint")
    text = args.checkpoint.read_text(encoding="utf-8")
    metadata, pos = top_metadata_and_constraints(text)
    names = metadata["generators"]
    monomials = [
        tuple((int(index), int(exponent)) for index, exponent in monomial)
        for monomial in linear["monomials"]
    ]
    rref = [
        [core.from_json(value) for value in row]
        for row in linear["rref_field_elements"]
    ]
    pivots = linear["pivot_columns"]
    if pivots != [0, 1, 2, 3] or len(rref) != 4:
        raise ValueError("unexpected RREF shape")
    assignment = [(index + 2, 0) for index in range(len(names))]
    pivot_indices = [monomials[column][0][0] for column in pivots]
    for row_index, pivot in enumerate(pivots):
        value = core.ZERO
        for column, monomial in enumerate(monomials):
            if column in pivots:
                continue
            value = core.add(
                value,
                core.mul(rref[row_index][column], monomial_value(monomial, assignment)),
            )
        assignment[pivot_indices[row_index]] = core.sub(core.ZERO, value)

    nonzero_rows = []
    for tag, saved in iter_json_array(text, pos):
        value = core.ZERO
        for sparse, coefficient_values in saved:
            monomial = tuple((int(index), int(exponent))
                             for index, exponent in sparse if exponent)
            value = core.add(
                value,
                core.mul(core.from_json(coefficient_values),
                         monomial_value(monomial, assignment)),
            )
        if not core.iszero(value):
            nonzero_rows.append({"tag": tag, "value": core.show(value)})
    if nonzero_rows:
        raise AssertionError(f"point failed {len(nonzero_rows)} constraints")

    result = {
        "source_checkpoint": str(args.checkpoint),
        "source_checkpoint_sha256": source_hash,
        "linear_receipt": str(args.linear_receipt),
        "component": "F3=0, 2+2 necessary constraint system",
        "field": "Q(sqrt(5))",
        "free_value_rule": "every generator is assigned its position index plus 2 in Q before solving the four pivots",
        "pivot_generators": [names[index] for index in pivot_indices],
        "constraint_count": 330,
        "nonzero_constraint_count": 0,
        "all_constraints_exactly_satisfied": True,
        "generator_values": {
            name: core.show(value) for name, value in zip(names, assignment)
        },
        "zero_generator_values": [
            name for name, value in zip(names, assignment) if core.iszero(value)
        ],
        "scope": (
            "This is an exact common zero of the saved finite necessary-equation "
            "system through distance nine. It is not a Jacobian pair and does not "
            "verify any inequalities or source-to-carrier conditions omitted from "
            "that checkpoint."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "component", "constraint_count", "nonzero_constraint_count",
        "all_constraints_exactly_satisfied", "pivot_generators",
        "zero_generator_values", "scope",
    )}, indent=2))


if __name__ == "__main__":
    main()
