"""Exact row reduction of the eight F3 distance-five equations."""

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


TARGET_TAGS = {
    *(f"terminal-alpha-P-distance5:u^{i}" for i in range(4)),
    *(f"terminal-alpha-Q-distance5:u^{i}" for i in range(4, 8)),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rref(rows):
    matrix = [list(row) for row in rows]
    if not matrix:
        return [], []
    pivots = []
    pivot_row = 0
    for column in range(len(matrix[0])):
        selected = next((i for i in range(pivot_row, len(matrix))
                         if not core.iszero(matrix[i][column])), None)
        if selected is None:
            continue
        matrix[pivot_row], matrix[selected] = matrix[selected], matrix[pivot_row]
        inverse = core.inv(matrix[pivot_row][column])
        matrix[pivot_row] = [core.mul(x, inverse) for x in matrix[pivot_row]]
        for i in range(len(matrix)):
            if i == pivot_row or core.iszero(matrix[i][column]):
                continue
            factor = matrix[i][column]
            matrix[i] = [core.sub(x, core.mul(factor, y))
                         for x, y in zip(matrix[i], matrix[pivot_row])]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return [row for row in matrix if any(not core.iszero(x) for x in row)], pivots


def encode_coefficient(value):
    rational, radical = value
    if not radical:
        return [str(rational)]
    return [str(radical), str(rational)]


def profile(checkpoint: Path) -> dict:
    source_hash = sha256(checkpoint)
    text = checkpoint.read_text(encoding="utf-8")
    metadata, pos = top_metadata_and_constraints(text)
    names = metadata["generators"]
    monomials = [
        ((names.index(f"p{i}"), 1),) for i in range(23, 17, -1)
    ] + [
        tuple(sorted(((names.index("p3"), 1), (names.index("p8"), 1)))),
        tuple(sorted(((names.index("p0"), 1), (names.index("p16"), 1)))),
        tuple(sorted(((names.index("p0"), 1), (names.index("p15"), 1)))),
        tuple(sorted(((names.index("p0"), 1), (names.index("p3"), 2)))),
        tuple(sorted(((names.index("p0"), 2), (names.index("p8"), 1)))),
        tuple(sorted(((names.index("p0"), 3), (names.index("p3"), 1)))),
        ((names.index("p0"), 5),),
    ]
    if len(monomials) != 13 or len(set(monomials)) != 13:
        raise AssertionError("distance-five support basis is malformed")
    rows = {}
    p05_column = monomials.index(((names.index("p0"), 5),))
    for tag, saved in iter_json_array(text, pos):
        if tag not in TARGET_TAGS:
            continue
        term_map = {
            tuple(sorted((int(index), int(power)) for index, power in sparse if power)):
            core.from_json(values)
            for sparse, values in saved
        }
        if set(term_map) != set(monomials):
            raise ValueError(f"unexpected monomial support in {tag}")
        row = [term_map[monomial] for monomial in monomials]
        if core.iszero(row[p05_column]):
            raise ValueError(f"p0^5 coefficient vanishes in {tag}")
        inverse = core.inv(row[p05_column])
        rows[tag] = [core.mul(value, inverse) for value in row]
    if set(rows) != TARGET_TAGS:
        raise ValueError("not all eight distance-five rows were found")

    reduced, pivots = rref([rows[tag] for tag in sorted(rows)])
    digest = hashlib.sha256()
    for tag in sorted(rows):
        digest.update(tag.encode("utf-8"))
        digest.update(repr(rows[tag]).encode("utf-8"))
    return {
        "source_checkpoint": str(checkpoint),
        "source_sha256": source_hash,
        "component": "F3=0",
        "field": "Q(sqrt(5))",
        "source_rows": sorted(TARGET_TAGS),
        "source_row_count": len(rows),
        "normalized_row_digest": digest.hexdigest(),
        "variables": [
            "p23", "p22", "p21", "p20", "p19", "p18",
            "p3*p8", "p0*p16", "p0*p15", "p0*p3^2",
            "p0^2*p8", "p0^3*p3", "p0^5",
        ],
        "monomials": [
            [[index, exponent] for index, exponent in monomial]
            for monomial in monomials
        ],
        "normalization": "each row divided by its nonzero p0^5 coefficient",
        "normalized_rows": {
            tag: [core.show(value) for value in rows[tag]]
            for tag in sorted(rows)
        },
        "rref": [[core.show(value) for value in row] for row in reduced],
        "rref_field_elements": [
            [encode_coefficient(value) for value in row] for row in reduced
        ],
        "pivot_columns": pivots,
        "rank": len(pivots),
        "scope": (
            "Exact linear elimination in six parameters and seven formal base "
            "monomials. Translating the resulting graph back into p0,p3,p8,p15,p16 "
            "and applying it to all constraints is still required."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    result = profile(args.checkpoint)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "component", "source_row_count", "normalized_row_digest", "variables",
        "rref", "pivot_columns", "rank",
    )}, indent=2))


if __name__ == "__main__":
    main()
