"""Exact linear reduction of the ten smallest F3 distance-four equations."""

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
    *(f"terminal-alpha-P-distance4:u^{i}" for i in range(5)),
    *(f"terminal-alpha-Q-distance4:u^{i}" for i in range(4, 9)),
}
VARIABLES = (
    "p10", "p11", "p12", "p13", "p14", "p15", "p16",
    "U=p0^4", "V=p0^2*p3", "W=p0*p8", "Z=p3^2",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rref(rows):
    matrix = [list(row) for row in rows]
    ncols = len(matrix[0]) if matrix else 0
    pivots = []
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
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    nonzero_rows = [row for row in matrix if any(not core.iszero(x) for x in row)]
    return nonzero_rows, pivots


def profile(checkpoint: Path, stream_profile: Path) -> dict:
    checkpoint_hash = sha256(checkpoint)
    audited = json.loads(stream_profile.read_text(encoding="utf-8"))
    if audited.get("source_sha256") != checkpoint_hash:
        raise ValueError("profile and checkpoint hashes disagree")
    text = checkpoint.read_text(encoding="utf-8")
    metadata, pos = top_metadata_and_constraints(text)
    names = metadata["generators"]
    expected_monomials = [
        *((names.index(f"p{i}"), 1) for i in range(10, 17)),
        ((names.index("p0"), 4),),
        tuple(sorted(((names.index("p0"), 2), (names.index("p3"), 1)))),
        tuple(sorted(((names.index("p0"), 1), (names.index("p8"), 1)))),
        ((names.index("p3"), 2),),
    ]
    expected_monomials = [tuple(item) if isinstance(item, tuple) else item
                          for item in expected_monomials]
    # The first seven entries above are single (generator, exponent) pairs.
    expected_monomials = [
        ((item[0], item[1]),) if len(item) == 2 and isinstance(item[0], int)
        else tuple(item)
        for item in expected_monomials
    ]

    raw_rows = {}
    for tag, saved in iter_json_array(text, pos):
        if tag not in TARGET_TAGS:
            continue
        term_map = {}
        for sparse, values in saved:
            monomial = tuple(sorted((int(index), int(power))
                                    for index, power in sparse if power))
            term_map[monomial] = core.from_json(values)
        if set(term_map) != set(expected_monomials) or len(saved) != 11:
            raise ValueError(f"unexpected support in {tag}")
        lead = term_map[expected_monomials[7]]
        if core.iszero(lead):
            raise ValueError(f"p0^4 coefficient vanishes in {tag}")
        inverse = core.inv(lead)
        raw_rows[tag] = tuple(core.mul(term_map[m], inverse)
                              for m in expected_monomials)
    if set(raw_rows) != TARGET_TAGS:
        raise ValueError("not all ten distance-four equations were found")

    unique_rows = []
    for tag in sorted(raw_rows):
        if raw_rows[tag] not in unique_rows:
            unique_rows.append(raw_rows[tag])
    reduced, pivots = rref(unique_rows)
    digest = hashlib.sha256()
    for tag in sorted(raw_rows):
        digest.update(tag.encode("utf-8"))
        digest.update(repr(raw_rows[tag]).encode("utf-8"))
    return {
        "checkpoint": str(checkpoint),
        "checkpoint_sha256": checkpoint_hash,
        "stream_profile": str(stream_profile),
        "component": "F3=0",
        "field": "Q(sqrt(5))",
        "source_rows": sorted(raw_rows),
        "source_row_count": len(raw_rows),
        "unique_normalized_equations": len(unique_rows),
        "variables": list(VARIABLES),
        "normalization": "each row divided by its nonzero p0^4 coefficient",
        "normalized_rows": {
            tag: [core.show(value) for value in raw_rows[tag]]
            for tag in sorted(raw_rows)
        },
        "normalized_row_digest": digest.hexdigest(),
        "rref": [[core.show(value) for value in row] for row in reduced],
        "pivot_columns": pivots,
        "rank": len(pivots),
        "scope": (
            "Exact linear elimination with formal U,V,W,Z; substituting their "
            "definitions and any nonzero chart conditions remains necessary."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("profile", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    result = profile(args.checkpoint, args.profile)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "component", "source_row_count", "unique_normalized_equations",
        "variables", "normalized_row_digest", "rref", "pivot_columns", "rank",
    )}, indent=2))


if __name__ == "__main__":
    main()
