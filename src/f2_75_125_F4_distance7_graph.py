"""Exact F4 distance-seven graph and substitution through all constraints."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import f2_75_125_quadratic_affine_profile as core
import f2_75_125_F4_distance4_graph as exact
from f2_75_125_quadratic_affine_profile_stream import (
    iter_json_array,
    top_metadata_and_constraints,
)


TARGET_TAGS = {
    *(f"terminal-alpha-P-distance7:u^{i}" for i in range(4)),
    *(f"terminal-alpha-Q-distance7:u^{i}" for i in range(4, 8)),
}


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

    source_hash = exact.sha256(args.checkpoint)
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
    new_monomials = [((names.index(f"p{i}"), 1),) for i in range(35, 29, -1)]
    row_maps = {}
    union = set(new_monomials)
    for tag, saved in iter_json_array(text, pos):
        if tag not in TARGET_TAGS:
            continue
        term_map = {
            tuple(sorted((int(i), int(e)) for i, e in sparse if e)):
            core.from_json(values)
            for sparse, values in saved
        }
        if not term_map:
            raise ValueError(f"unexpected zero F4 target row {tag}")
        row_maps[tag] = term_map
        union.update(term_map)
    if set(row_maps) != TARGET_TAGS:
        raise ValueError("not all eight F4 distance-seven rows were found")
    monomials = new_monomials + sorted(union - set(new_monomials))
    normalizer = ((names.index("p0"), 7),)
    if normalizer not in monomials:
        raise ValueError("p0^7 normalization monomial is absent")
    normalizer_column = monomials.index(normalizer)
    rows = {}
    for tag, term_map in row_maps.items():
        row = [term_map.get(m, core.ZERO) for m in monomials]
        if core.iszero(row[normalizer_column]):
            raise ValueError(f"p0^7 coefficient vanishes in {tag}")
        inverse = core.inv(row[normalizer_column])
        rows[tag] = [core.mul(value, inverse) for value in row]

    reduced, pivots = exact.rref([rows[tag] for tag in sorted(rows)])
    if pivots != [0, 1, 2, 3] or len(reduced) != 4:
        raise ValueError(f"unexpected F4 distance-seven pivot structure {pivots}")
    digest = hashlib.sha256()
    for tag in sorted(rows):
        digest.update(tag.encode("utf-8"))
        digest.update(repr(rows[tag]).encode("utf-8"))
    labels = [
        names[monomial[0][0]] if len(monomial) == 1 else
        "*".join(f"{names[i]}^{exponent}" for i, exponent in monomial)
        for monomial in monomials
    ]
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
        "variables": labels,
        "monomials": [[[i, e] for i, e in m] for m in monomials],
        "normalization": "each row divided by its nonzero p0^7 coefficient",
        "normalized_rows": {tag: [core.show(x) for x in rows[tag]]
                            for tag in sorted(rows)},
        "rref": [[core.show(x) for x in row] for row in reduced],
        "rref_field_elements": [[exact.encode_coefficient(x) for x in row]
                                 for row in reduced],
        "pivot_columns": pivots,
        "rank": len(pivots),
        "scope": "Exact linear graph from the F4 terminal-alpha distance-seven rows.",
    }

    zero = (0,) * len(names)
    eliminated = tuple(names.index(f"p{i}") for i in range(35, 31, -1))
    replacements = []
    for row in reduced:
        polynomial = {}
        for column, monomial in enumerate(monomials):
            if column in pivots or core.iszero(row[column]):
                continue
            powers = [0] * len(names)
            for index, exponent in monomial:
                powers[index] = exponent
            polynomial[tuple(powers)] = core.sub(core.ZERO, row[column])
        if zero in polynomial or not polynomial:
            raise AssertionError("unexpected constant or zero F4 distance-seven graph row")
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
            values.append(exact.multiply_poly(values[-1], replacement))
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
            "method": "exact substitution of the rank-four distance-seven graph",
            "eliminated_parameters": ["p35", "p34", "p33", "p32"],
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
            transformed = exact.apply_one(saved, len(names), eliminated, powers, cache)
            output_terms += len(transformed)
            if not transformed:
                zero_count += 1
                zero_tags.add(tag)
            if tag in TARGET_TAGS and transformed:
                raise AssertionError(f"F4 distance-seven graph failed to kill {tag}")
            if any(any(monomial[index] for index in eliminated)
                   for monomial in transformed):
                raise AssertionError("an eliminated F4 parameter reappeared")
            encoded = [
                [[[i, power] for i, power in enumerate(monomial) if power],
                 exact.encode_coefficient(coefficient)]
                for monomial, coefficient in sorted(transformed.items())
            ]
            if not first:
                stream.write(",")
            json.dump([tag, encoded], stream, separators=(",", ":"))
            first = False
        stream.write("]}\n")
    if not TARGET_TAGS.issubset(zero_tags):
        raise AssertionError("not all F4 distance-seven rows vanished")

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
            name: exponent for name, exponent in zip(("p35", "p34", "p33", "p32"), maxima)
        },
        "polynomial_power_cache_entries": len(cache),
        "output_checkpoint": str(args.output),
        "output_checkpoint_sha256": exact.sha256(args.output),
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
