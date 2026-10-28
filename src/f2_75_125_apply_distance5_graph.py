"""Apply the exact distance-five F3 linear graph to all reduced constraints."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import f2_75_125_quadratic_affine_profile as core
from f2_75_125_graph_distance5_linear_eliminate import (
    TARGET_TAGS,
    profile as audit_graph,
)
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


def encode_coefficient(value):
    rational, radical = value
    if not radical:
        return [str(rational)]
    return [str(radical), str(rational)]


def graph_from_receipt(names, receipt):
    monomials = [
        tuple(sorted((int(index), int(power)) for index, power in monomial))
        for monomial in receipt["monomials"]
    ]
    rref = [
        [core.from_json(value) for value in row]
        for row in receipt["rref_field_elements"]
    ]
    pivots = receipt["pivot_columns"]
    if pivots != [0, 1, 2, 3] or len(rref) != 4:
        raise ValueError("unexpected distance-five pivot structure")
    eliminated = tuple(names.index(f"p{i}") for i in range(23, 19, -1))
    zero = (0,) * len(names)
    replacements = {}
    for row_index, pivot in enumerate(pivots):
        pivot_monomial = monomials[pivot]
        if len(pivot_monomial) != 1 or pivot_monomial[0][1] != 1:
            raise ValueError("pivot is not a parameter variable")
        pivot_index = pivot_monomial[0][0]
        if pivot_index != eliminated[row_index]:
            raise ValueError("unexpected eliminated parameter order")
        poly = {}
        for column, monomial in enumerate(monomials):
            if column in pivots:
                continue
            coefficient = core.sub(core.ZERO, rref[row_index][column])
            if core.iszero(coefficient):
                continue
            powers = [0] * len(names)
            for index, exponent in monomial:
                powers[index] = exponent
            key = tuple(powers)
            updated = core.add(poly.get(key, core.ZERO), coefficient)
            if core.iszero(updated):
                poly.pop(key, None)
            else:
                poly[key] = updated
        if zero in poly or not poly:
            raise ValueError("unexpected constant or zero graph row")
        replacements[names[pivot_index]] = poly
    return eliminated, replacements, zero


def apply_one(saved, size, eliminated, powers, cache):
    result = {}
    eliminated_lookup = {index: i for i, index in enumerate(eliminated)}
    for sparse, coefficient_values in saved:
        exponents = [0] * len(eliminated)
        original = [0] * size
        for index, exponent in sparse:
            if index in eliminated_lookup:
                exponents[eliminated_lookup[index]] = exponent
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
    parser.add_argument("stream_profile", type=Path)
    parser.add_argument("linear_elimination", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("receipt", type=Path)
    args = parser.parse_args()
    for path in (args.output, args.receipt):
        if path.exists():
            raise FileExistsError(f"refusing to overwrite {path}")

    source_hash = sha256(args.checkpoint)
    source_profile = json.loads(args.stream_profile.read_text(encoding="utf-8"))
    linear_receipt = json.loads(args.linear_elimination.read_text(encoding="utf-8"))
    if source_profile.get("source_sha256") != source_hash:
        raise ValueError("stream profile does not match source constraints")
    if linear_receipt.get("source_sha256") != source_hash:
        raise ValueError("linear-elimination receipt does not match source constraints")
    replay = audit_graph(args.checkpoint)
    for key in ("normalized_row_digest", "rref", "rref_field_elements",
                "pivot_columns", "rank"):
        if replay.get(key) != linear_receipt.get(key):
            raise ValueError(f"linear-elimination receipt failed replay at {key}")
    if replay["rank"] != 4:
        raise ValueError("expected the audited rank-four distance-five graph")

    text = args.checkpoint.read_text(encoding="utf-8")
    metadata, constraint_pos = top_metadata_and_constraints(text)
    names = metadata["generators"]
    eliminated, replacement, zero = graph_from_receipt(names, linear_receipt)
    maxima = [0] * len(eliminated)
    for _tag, saved in iter_json_array(text, constraint_pos):
        for sparse, _coefficient in saved:
            for index, exponent in sparse:
                if index in eliminated:
                    maxima[eliminated.index(index)] = max(
                        maxima[eliminated.index(index)], exponent)
    powers = []
    for name, maximum in zip(replacement, maxima):
        values = [{zero: core.ONE}]
        for _ in range(maximum):
            values.append(multiply_poly(values[-1], replacement[name]))
        powers.append(values)

    output_data = {
        "format": "f2-sparse-ring-v1",
        "branch": metadata["branch"],
        "quadratic": True,
        "generators": names,
        "nparameters": int(metadata["nparameters"]),
        "compaction": {
            "source_checkpoint_sha256": source_hash,
            "component": "F3=0",
            "linear_graph_digest": linear_receipt["normalized_row_digest"],
            "method": "exact substitution of the rank-four distance-five linear graph",
            "eliminated_parameters": [f"p{i}" for i in range(23, 19, -1)],
            "scope": "constraint-only necessary system restricted to the verified exact graph",
        },
        "constraints": [],
    }
    cache = {}
    constraint_count = source_term_count = output_term_count = zero_count = 0
    zero_tags = set()
    with args.output.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps({k: v for k, v in output_data.items()
                                 if k != "constraints"}, separators=(",", ":"))[:-1])
        stream.write(',"constraints":[')
        first = True
        for tag, saved in iter_json_array(text, constraint_pos):
            constraint_count += 1
            source_term_count += len(saved)
            transformed = apply_one(saved, len(names), eliminated, powers, cache)
            output_term_count += len(transformed)
            if not transformed:
                zero_count += 1
                zero_tags.add(tag)
            if tag in TARGET_TAGS and transformed:
                raise AssertionError(f"distance-five graph failed to kill {tag}")
            if any(any(monomial[index] for index in eliminated)
                   for monomial in transformed):
                raise AssertionError(f"an eliminated generator reappeared in {tag}")
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
        raise AssertionError("not all distance-five rows vanished")
    output_hash = sha256(args.output)
    receipt = {
        "source_checkpoint": str(args.checkpoint),
        "source_checkpoint_sha256": source_hash,
        "source_stream_profile": str(args.stream_profile),
        "source_linear_elimination": str(args.linear_elimination),
        "component": "F3=0",
        "branch": metadata["branch"],
        "quadratic_field": "Q(sqrt(5))",
        "source_constraint_count": constraint_count,
        "source_constraint_terms": source_term_count,
        "transformed_constraint_terms": output_term_count,
        "zero_constraints_after_graph": zero_count,
        "zero_constraint_tags_after_graph": sorted(zero_tags),
        "distance_five_rows_vanishing": sorted(TARGET_TAGS),
        "maximum_exponents_eliminated": {
            f"p{i}": exponent for i, exponent in zip(range(23, 19, -1), maxima)
        },
        "polynomial_power_cache_entries": len(cache),
        "output_checkpoint": str(args.output),
        "output_checkpoint_sha256": output_hash,
        "scope": (
            "Exact necessary-system graph restriction. This is not an exclusion "
            "unless the transformed ideal is proved unit and the reduction is "
            "attached to the actual hypothetical pair."
        ),
    }
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: receipt[key] for key in (
        "component", "source_constraint_count", "source_constraint_terms",
        "transformed_constraint_terms", "zero_constraints_after_graph",
        "distance_five_rows_vanishing", "maximum_exponents_eliminated",
        "polynomial_power_cache_entries", "output_checkpoint_sha256", "scope",
    )}, indent=2))


if __name__ == "__main__":
    main()
