#!/usr/bin/env python3
"""Stream a sparse FLINT checkpoint through exact zero specializations.

The selected parameter generators are set to zero, their monomials are
discarded, zero constraint rows are removed, and the remaining parameter
names are compacted.  No polynomial expansion or coefficient arithmetic is
needed because every specialization is zero.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
from io import StringIO
from pathlib import Path

import ijson

ENCODER = json.JSONEncoder(separators=(",", ":"))


def consume_value(events, first):
    prefix, event, value = first
    if event == "start_array":
        result = []
        while True:
            current = next(events)
            if current[0] == prefix and current[1] == "end_array":
                return result
            result.append(consume_value(events, current))
    if event == "start_map":
        result = {}
        key = None
        while True:
            current = next(events)
            if current[0] == prefix and current[1] == "end_map":
                return result
            if current[0] == prefix and current[1] == "map_key":
                key = current[2]
            else:
                if key is None:
                    raise ValueError("malformed JSON map event stream")
                result[key] = consume_value(events, current)
                key = None
    if event in ("end_array", "end_map", "map_key"):
        raise ValueError(f"unexpected JSON event {event!r} at {prefix!r}")
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def compact_generators(names: list[str], zero_names: set[str]) -> list[str]:
    result = []
    for name in names:
        if name in zero_names:
            continue
        match = re.fullmatch(r"p(\d+)", name)
        if match:
            index = int(match.group(1))
            shift = sum(int(value[1:]) < index for value in zero_names)
            name = f"p{index-shift}"
        result.append(name)
    return result


def transform_polynomial(events, first, stream, zero_indices: set[int]) -> bool:
    prefix, event, _ = first
    if event == "null":
        stream.write("null")
        return False
    if event != "start_array":
        raise ValueError(f"expected a polynomial at {prefix!r}")
    stream.write("[")
    first_term = True
    while True:
        current = next(events)
        if current[0] == prefix and current[1] == "end_array":
            break
        sparse, values = consume_value(events, current)
        if any(index in zero_indices and exponent for index, exponent in sparse):
            continue
        shifted = []
        for index, exponent in sparse:
            if index in zero_indices:
                continue
            shift = sum(removed < index for removed in zero_indices)
            shifted.append([index - shift, exponent])
        if not first_term:
            stream.write(",")
        first_term = False
        stream.write(ENCODER.encode([shifted, values]))
    stream.write("]")
    return not first_term


def transform_polynomial_array(events, first, stream, zero_indices):
    prefix, event, _ = first
    if event != "start_array":
        raise ValueError(f"expected polynomial array at {prefix!r}")
    stream.write("[")
    first_item = True
    while True:
        current = next(events)
        if current[0] == prefix and current[1] == "end_array":
            break
        if not first_item:
            stream.write(",")
        first_item = False
        transform_polynomial(events, current, stream, zero_indices)
    stream.write("]")


def transform_constraints(events, first, stream, zero_indices):
    prefix, event, _ = first
    if event != "start_array":
        raise ValueError("malformed constraints array")
    stream.write("[")
    first_item = True
    count_in = count_out = 0
    while True:
        current = next(events)
        if current[0] == prefix and current[1] == "end_array":
            break
        if current[1] != "start_array":
            raise ValueError("malformed constraint row")
        count_in += 1
        row_prefix = current[0]
        tag = consume_value(events, next(events))
        polynomial_event = next(events)
        buffer = StringIO()
        nonzero = transform_polynomial(events, polynomial_event, buffer, zero_indices)
        end_row = next(events)
        if end_row[0] != row_prefix or end_row[1] != "end_array":
            raise ValueError("malformed constraint row end")
        if not nonzero:
            continue
        if not first_item:
            stream.write(",")
        first_item = False
        count_out += 1
        stream.write("[")
        stream.write(ENCODER.encode(tag))
        stream.write(",")
        stream.write(buffer.getvalue())
        stream.write("]")
    stream.write("]")
    return count_in, count_out


def transform_uk(events, first, stream, zero_indices):
    prefix, event, _ = first
    if event != "start_map":
        raise ValueError("malformed uk map")
    stream.write("{")
    first_item = True
    while True:
        current = next(events)
        if current[0] == prefix and current[1] == "end_map":
            break
        if current[1] != "map_key":
            raise ValueError("malformed uk map key")
        key = current[2]
        polynomial_event = next(events)
        if not first_item:
            stream.write(",")
        first_item = False
        stream.write(ENCODER.encode(key))
        stream.write(":")
        transform_polynomial(events, polynomial_event, stream, zero_indices)
    stream.write("}")


def specialize(source: Path, destination: Path, zero_names: set[str]) -> dict:
    source_hash = sha256(source)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    events = iter(ijson.parse(source.open("rb")))
    first = next(events)
    if first[1] != "start_map":
        raise ValueError("checkpoint root is not an object")

    old_names = None
    zero_indices = None
    old_nparameters = None
    new_nparameters = None
    old_compaction = None
    constraint_counts = None
    first_field = True
    saw_compaction = False

    def field(stream, name):
        nonlocal first_field
        if not first_field:
            stream.write(",")
        first_field = False
        stream.write(ENCODER.encode(name))
        stream.write(":")

    with temporary.open("w", encoding="utf-8", newline="") as stream:
        stream.write("{")
        while True:
            current = next(events)
            if current[0] == "" and current[1] == "end_map":
                break
            if current[0] != "" or current[1] != "map_key":
                raise ValueError(f"unexpected top-level JSON event {current!r}")
            key = current[2]
            value_event = next(events)
            field(stream, key)

            if key == "generators":
                old_names = consume_value(events, value_event)
                if any(name not in old_names for name in zero_names):
                    raise ValueError("a requested zero parameter is absent")
                zero_indices = {old_names.index(name) for name in zero_names}
                old_p_indices = {int(name[1:]) for name in zero_names}
                new_names = compact_generators(old_names, zero_names)
                stream.write(ENCODER.encode(new_names))
            elif key == "nparameters":
                old_nparameters = int(consume_value(events, value_event))
                if old_names is None or zero_indices is None:
                    raise ValueError("generator metadata must precede nparameters")
                allocated_zeros = sum(
                    int(name[1:]) < old_nparameters
                    for name in zero_names
                )
                if allocated_zeros != len(zero_names):
                    raise ValueError("all specialized parameters must be allocated")
                new_nparameters = old_nparameters - allocated_zeros
                stream.write(str(new_nparameters))
            elif key in {"h", "v", "S", "B", "kappa", "kappa2", "kappa3"}:
                if zero_indices is None:
                    raise ValueError("generator metadata must precede polynomials")
                transform_polynomial(events, value_event, stream, zero_indices)
            elif key == "constraints":
                if zero_indices is None:
                    raise ValueError("generator metadata must precede constraints")
                constraint_counts = transform_constraints(
                    events, value_event, stream, zero_indices
                )
            elif key == "uk":
                transform_uk(events, value_event, stream, zero_indices)
            elif key in {"p", "q", "w"}:
                transform_polynomial_array(events, value_event, stream, zero_indices)
            elif key == "compaction":
                old_compaction = consume_value(events, value_event)
                saw_compaction = True
                compacted = {
                    "source": str(source),
                    "source_sha256": source_hash,
                    "zero_assignments": {name: "0" for name in sorted(zero_names)},
                    "removed_parameters": sorted(zero_names),
                    "previous": old_compaction,
                }
                stream.write(ENCODER.encode(compacted))
            else:
                value = consume_value(events, value_event)
                stream.write(ENCODER.encode(value))

        if not saw_compaction:
            field(stream, "compaction")
            stream.write(ENCODER.encode({
                "source": str(source),
                "source_sha256": source_hash,
                "zero_assignments": {name: "0" for name in sorted(zero_names)},
                "removed_parameters": sorted(zero_names),
                "previous": old_compaction,
            }))
        stream.write("}")
        stream.flush()
        os.fsync(stream.fileno())

    if old_names is None or zero_indices is None or new_nparameters is None:
        raise ValueError("checkpoint metadata is incomplete")
    os.replace(temporary, destination)
    return {
        "source": str(source),
        "source_sha256": source_hash,
        "destination": str(destination),
        "destination_sha256": sha256(destination),
        "zero_assignments": {name: "0" for name in sorted(zero_names)},
        "source_parameters": old_nparameters,
        "destination_parameters": new_nparameters,
        "constraints_before_after": constraint_counts,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--zero", action="append", required=True,
                        help="force a parameter to zero, e.g. --zero p79")
    args = parser.parse_args()
    result = specialize(args.source, args.destination, set(args.zero))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
