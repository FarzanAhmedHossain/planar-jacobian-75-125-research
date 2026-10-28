"""Exact linear elimination of the six smallest post-graph F3 equations.

After the p1,p2 graph substitution, six support-six equations are cubic but
linear in (p5,p6,p7,p8,p0^3,p0*p3).  Row-reducing them over Q(sqrt(5)) can
expose the next exact parameter relations without a large Groebner basis.
"""

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


TAGS = {
    "terminal-alpha-P-distance3:u^2",
    "terminal-alpha-P-distance3:u^3",
    "terminal-alpha-P-distance3:u^4",
    "terminal-alpha-Q-distance3:u^6",
    "terminal-alpha-Q-distance3:u^7",
    "terminal-alpha-Q-distance3:u^8",
}
VARIABLES = ("p5", "p6", "p7", "p8", "U=p0^3", "V=p0*p3")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def show_row(row):
    return [core.show(value) for value in row]


def rref(rows):
    matrix = [list(row) for row in rows]
    if not matrix:
        return [], []
    ncols = len(matrix[0])
    pivot_cols = []
    pivot_row = 0
    for column in range(ncols):
        selected = next((i for i in range(pivot_row, len(matrix))
                         if not core.iszero(matrix[i][column])), None)
        if selected is None:
            continue
        matrix[pivot_row], matrix[selected] = matrix[selected], matrix[pivot_row]
        inverse = core.inv(matrix[pivot_row][column])
        matrix[pivot_row] = [core.mul(value, inverse)
                             for value in matrix[pivot_row]]
        for i in range(len(matrix)):
            if i == pivot_row:
                continue
            factor = matrix[i][column]
            if core.iszero(factor):
                continue
            matrix[i] = [core.sub(value, core.mul(factor, pivot_value))
                         for value, pivot_value in zip(matrix[i], matrix[pivot_row])]
        pivot_cols.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    nonzero_rows = [row for row in matrix if any(not core.iszero(x) for x in row)]
    return nonzero_rows, pivot_cols


def profile(checkpoint: Path, stream_profile: Path) -> dict:
    checkpoint_hash = sha256(checkpoint)
    profile_data = json.loads(stream_profile.read_text(encoding="utf-8"))
    if profile_data.get("source_sha256") != checkpoint_hash:
        raise ValueError("support profile does not match the sparse checkpoint")
    text = checkpoint.read_text(encoding="utf-8")
    metadata, pos = top_metadata_and_constraints(text)
    names = metadata["generators"]
    expected_monomials = [
        tuple(sorted(((names.index("p5"), 1),))),
        tuple(sorted(((names.index("p6"), 1),))),
        tuple(sorted(((names.index("p7"), 1),))),
        tuple(sorted(((names.index("p8"), 1),))),
        tuple(sorted(((names.index("p0"), 3),))),
        tuple(sorted(((names.index("p0"), 1), (names.index("p3"), 1)))),
    ]
    target = {}
    for tag, saved in iter_json_array(text, pos):
        if tag not in TAGS:
            continue
        term_map = {}
        for sparse, values in saved:
            monomial = tuple(sorted((int(index), int(power))
                                    for index, power in sparse if power))
            term_map[monomial] = core.from_json(values)
        if set(term_map) != set(expected_monomials):
            raise ValueError(f"unexpected monomial support in {tag}")
        lead = term_map[expected_monomials[4]]
        if core.iszero(lead):
            raise ValueError(f"missing p0^3 coefficient in {tag}")
        inverse = core.inv(lead)
        row = tuple(core.mul(term_map[m], inverse) for m in expected_monomials)
        target[tag] = row
    if set(target) != TAGS:
        raise ValueError("not all six audited cubic equations were found")

    unique_rows = []
    for tag in sorted(target):
        row = target[tag]
        if row not in unique_rows:
            unique_rows.append(row)
    reduced, pivots = rref(unique_rows)
    if any(core.iszero(row[-2]) and core.iszero(row[-1])
           and all(core.iszero(x) for x in row[:4]) for row in reduced):
        raise AssertionError("unexpected zero row in reduced system")

    coefficient_digest = hashlib.sha256()
    for tag in sorted(target):
        coefficient_digest.update(tag.encode("utf-8"))
        coefficient_digest.update(repr(target[tag]).encode("utf-8"))
    return {
        "checkpoint": str(checkpoint),
        "checkpoint_sha256": checkpoint_hash,
        "stream_profile": str(stream_profile),
        "component": "F3=0",
        "field": "Q(sqrt(5))",
        "source_rows": sorted(TAGS),
        "source_row_count": len(target),
        "unique_normalized_equations": len(unique_rows),
        "variables": list(VARIABLES),
        "equation_normalization": "each row divided by its nonzero p0^3 coefficient",
        "normalized_rows": {
            tag: show_row(target[tag]) for tag in sorted(target)
        },
        "normalized_row_digest": coefficient_digest.hexdigest(),
        "rref": [show_row(row) for row in reduced],
        "pivot_columns": pivots,
        "rank": len(pivots),
        "scope": (
            "Exact linear elimination of the selected cubic necessary equations. "
            "U and V stand for p0^3 and p0*p3; translating their relations back "
            "to p0,p3 and using any nonzero chart conditions is a separate step."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("stream_profile", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    result = profile(args.checkpoint, args.stream_profile)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "component", "source_row_count", "unique_normalized_equations",
        "variables", "normalized_row_digest", "rref", "pivot_columns",
        "rank", "scope",
    )}, indent=2))


if __name__ == "__main__":
    main()
