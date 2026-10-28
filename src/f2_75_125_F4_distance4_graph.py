"""Exact F4 distance-four row graph and substitution through all constraints."""

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
    *(f"terminal-alpha-P-distance4:u^{i}" for i in range(1, 5)),
    *(f"terminal-alpha-Q-distance4:u^{i}" for i in range(5, 9)),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rref(rows):
    matrix = [list(row) for row in rows]
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


def multiply_poly(left, right):
    result = {}
    for monomial_a, coefficient_a in left.items():
        for monomial_b, coefficient_b in right.items():
            monomial = tuple(a + b for a, b in zip(monomial_a, monomial_b))
            value = core.mul(coefficient_a, coefficient_b)
            updated = core.add(result.get(monomial, core.ZERO), value)
            if core.iszero(updated):
                result.pop(monomial, None)
            else:
                result[monomial] = updated
    return result


def apply_one(saved, size, eliminated, powers, cache):
    lookup = {index: i for i, index in enumerate(eliminated)}
    result = {}
    for sparse, coefficient_values in saved:
        exponents = [0] * len(eliminated)
        original = [0] * size
        for index, exponent in sparse:
            index = int(index)
            exponent = int(exponent)
            if index in lookup:
                exponents[lookup[index]] = exponent
            else:
                original[index] = exponent
        key = tuple(exponents)
        if key not in cache:
            expanded = {(0,) * size: core.ONE}
            for i, exponent in enumerate(key):
                if exponent:
                    expanded = multiply_poly(expanded, powers[i][exponent])
            cache[key] = expanded
        coefficient = core.from_json(coefficient_values)
        for shift, shift_coefficient in cache[key].items():
            monomial = tuple(a + b for a, b in zip(original, shift))
            value = core.mul(coefficient, shift_coefficient)
            updated = core.add(result.get(monomial, core.ZERO), value)
            if core.iszero(updated):
                result.pop(monomial, None)
            else:
                result[monomial] = updated
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("source_receipt", type=Path)
    parser.add_argument("stream_profile", type=Path)
    parser.add_argument("linear_receipt", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("graph_receipt", type=Path)
    args = parser.parse_args()
    for path in (args.linear_receipt, args.output, args.graph_receipt):
        if path.exists():
            raise FileExistsError(f"refusing to overwrite {path}")

    source_hash = sha256(args.checkpoint)
    provenance = json.loads(args.source_receipt.read_text(encoding="utf-8"))
    stream_profile = json.loads(args.stream_profile.read_text(encoding="utf-8"))
    if provenance.get("component") != "F4=0":
        raise ValueError("source receipt is not the F4 component")
    if provenance.get("output_checkpoint_sha256") != source_hash:
        raise ValueError("source receipt does not certify this checkpoint")
    if stream_profile.get("source_sha256") != source_hash:
        raise ValueError("stream profile does not match checkpoint")

    text = args.checkpoint.read_text(encoding="utf-8")
    metadata, pos = top_metadata_and_constraints(text)
    names = metadata["generators"]
    monomials = [
        ((names.index(f"p{i}"), 1),) for i in range(16, 10, -1)
    ] + [
        ((names.index("p3"), 2),),
        tuple(sorted(((names.index("p0"), 1), (names.index("p4"), 1)))),
        tuple(sorted(((names.index("p0"), 2), (names.index("p3"), 1)))),
        ((names.index("p0"), 4),),
    ]
    if len(monomials) != 10 or len(set(monomials)) != 10:
        raise AssertionError("F4 distance-four support basis is malformed")
    normalizer = monomials.index(((names.index("p0"), 4),))
    rows = {}
    for tag, saved in iter_json_array(text, pos):
        if tag not in TARGET_TAGS:
            continue
        term_map = {
            tuple(sorted((int(i), int(e)) for i, e in sparse if e)):
            core.from_json(values)
            for sparse, values in saved
        }
        if set(term_map) != set(monomials):
            raise ValueError(f"unexpected F4 distance-four support in {tag}")
        row = [term_map[m] for m in monomials]
        if core.iszero(row[normalizer]):
            raise ValueError(f"p0^4 coefficient vanishes in {tag}")
        inverse = core.inv(row[normalizer])
        rows[tag] = [core.mul(value, inverse) for value in row]
    if set(rows) != TARGET_TAGS:
        raise ValueError("not all eight F4 distance-four rows were found")

    reduced, pivots = rref([rows[tag] for tag in sorted(rows)])
    if pivots != [0, 1, 2, 3] or len(reduced) != 4:
        raise ValueError(f"unexpected F4 distance-four pivot structure {pivots}")
    digest = hashlib.sha256()
    for tag in sorted(rows):
        digest.update(tag.encode("utf-8"))
        digest.update(repr(rows[tag]).encode("utf-8"))
    linear = {
        "source_checkpoint": str(args.checkpoint),
        "source_checkpoint_sha256": source_hash,
        "source_receipt": str(args.source_receipt),
        "stream_profile": str(args.stream_profile),
        "component": "F4=0",
        "branch": "2+2",
        "field": "Q(sqrt(5))",
        "source_rows": sorted(TARGET_TAGS),
        "source_row_count": len(rows),
        "normalized_row_digest": digest.hexdigest(),
        "variables": ["p16", "p15", "p14", "p13", "p12", "p11",
                      "p3^2", "p0*p4", "p0^2*p3", "p0^4"],
        "monomials": [
            [[index, exponent] for index, exponent in monomial]
            for monomial in monomials
        ],
        "normalization": "each row divided by the nonzero p0^4 coefficient",
        "normalized_rows": {tag: [core.show(x) for x in rows[tag]]
                            for tag in sorted(rows)},
        "rref": [[core.show(x) for x in row] for row in reduced],
        "rref_field_elements": [[encode_coefficient(x) for x in row]
                                 for row in reduced],
        "pivot_columns": pivots,
        "rank": len(pivots),
        "scope": "Exact linear graph from the F4 terminal-alpha distance-four rows.",
    }

    eliminated = tuple(names.index(f"p{i}") for i in range(16, 12, -1))
    zero = (0,) * len(names)
    replacements = []
    for row in reduced:
        polynomial = {}
        for column, monomial in enumerate(monomials):
            if column in pivots or core.iszero(row[column]):
                continue
            powers_vector = [0] * len(names)
            for index, exponent in monomial:
                powers_vector[index] = exponent
            polynomial[tuple(powers_vector)] = core.sub(core.ZERO, row[column])
        if zero in polynomial or not polynomial:
            raise AssertionError("unexpected constant or zero F4 graph row")
        replacements.append(polynomial)

    maxima = [0] * len(eliminated)
    for _tag, saved in iter_json_array(text, pos):
        for sparse, _coefficient in saved:
            for index, exponent in sparse:
                index = int(index)
                if index in eliminated:
                    maxima[eliminated.index(index)] = max(
                        maxima[eliminated.index(index)], int(exponent))
    powers = []
    for replacement, maximum in zip(replacements, maxima):
        values = [{zero: core.ONE}]
        for _ in range(maximum):
            values.append(multiply_poly(values[-1], replacement))
        powers.append(values)

    output_data = {
        "format": "f2-sparse-ring-v1",
        "branch": metadata["branch"],
        "quadratic": True,
        "generators": names,
        "nparameters": int(metadata["nparameters"]),
        "compaction": {
            "source_checkpoint_sha256": source_hash,
            "component": "F4=0",
            "linear_graph_digest": linear["normalized_row_digest"],
            "method": "exact substitution of the rank-four distance-four graph",
            "eliminated_parameters": ["p16", "p15", "p14", "p13"],
            "scope": "constraint-only necessary system restricted to the verified exact graph",
        },
        "constraints": [],
    }
    cache = {}
    constraint_count = source_terms = output_terms = zero_count = 0
    zero_tags = set()
    with args.output.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps({k: v for k, v in output_data.items()
                                 if k != "constraints"}, separators=(",", ":"))[:-1])
        stream.write(',"constraints":[')
        first = True
        for tag, saved in iter_json_array(text, pos):
            constraint_count += 1
            source_terms += len(saved)
            transformed = apply_one(saved, len(names), eliminated, powers, cache)
            output_terms += len(transformed)
            if not transformed:
                zero_count += 1
                zero_tags.add(tag)
            if tag in TARGET_TAGS and transformed:
                raise AssertionError(f"F4 graph failed to kill {tag}")
            if any(any(monomial[index] for index in eliminated)
                   for monomial in transformed):
                raise AssertionError("an eliminated F4 parameter reappeared")
            encoded = [
                [[[i, power] for i, power in enumerate(monomial) if power],
                 encode_coefficient(coefficient)]
                for monomial, coefficient in sorted(transformed.items())
            ]
            if not first:
                stream.write(",")
            json.dump([tag, encoded], stream, separators=(",", ":"))
            first = False
        stream.write("]}\n")
    if not TARGET_TAGS.issubset(zero_tags):
        raise AssertionError("not all F4 distance-four rows vanished")

    graph_receipt = {
        "source_checkpoint": str(args.checkpoint),
        "source_checkpoint_sha256": source_hash,
        "source_receipt": str(args.source_receipt),
        "source_stream_profile": str(args.stream_profile),
        "linear_receipt": str(args.linear_receipt),
        "component": "F4=0",
        "branch": "2+2",
        "field": "Q(sqrt(5))",
        "source_constraint_count": constraint_count,
        "source_constraint_terms": source_terms,
        "transformed_constraint_terms": output_terms,
        "zero_constraints_after_graph": zero_count,
        "zero_constraint_tags_after_graph": sorted(zero_tags),
        "maximum_exponents_eliminated": {
            name: exponent for name, exponent in zip(("p16", "p15", "p14", "p13"), maxima)
        },
        "polynomial_power_cache_entries": len(cache),
        "output_checkpoint": str(args.output),
        "output_checkpoint_sha256": sha256(args.output),
        "scope": "Exact necessary-system graph restriction, not an exclusion until the residual ideal is certified unit.",
    }
    args.linear_receipt.write_text(json.dumps(linear, indent=2) + "\n", encoding="utf-8")
    args.graph_receipt.write_text(json.dumps(graph_receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: graph_receipt[key] for key in (
        "component", "source_constraint_count", "source_constraint_terms",
        "transformed_constraint_terms", "zero_constraints_after_graph",
        "maximum_exponents_eliminated", "polynomial_power_cache_entries",
        "output_checkpoint_sha256", "scope",
    )}, indent=2))


if __name__ == "__main__":
    main()
