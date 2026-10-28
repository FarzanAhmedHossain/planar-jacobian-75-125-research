#!/usr/bin/env python3
"""Stream selected terminal-Q constraints from an F2 sparse-ring checkpoint.

This deliberately stops after the checkpoint's ``constraints`` array and
materializes only matching polynomials.  It is intended for inspecting the
active parameter support of a few retained residuals without loading the
checkpoint (which can exceed 1 GB) into memory.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "f2_75_125_flint_deps"))
import ijson  # noqa: E402


TARGET_TAGS = {
    "terminal-Q-distance17:0:t^0",
    "terminal-Q-distance17:0:t^1",
    "terminal-Q-distance17:0:t^2",
    "terminal-Q-distance17:postsolve:t^0",
    "terminal-Q-distance17:postsolve:t^1",
    "terminal-Q-distance17:postsolve:t^2",
}


def skip_value(events, first):
    """Consume one JSON value from an ijson event stream with O(1) memory."""
    _, event, _ = first
    if event not in ("start_array", "start_map"):
        return
    depth = 1
    while depth:
        _, event, _ = next(events)
        if event in ("start_array", "start_map"):
            depth += 1
        elif event in ("end_array", "end_map"):
            depth -= 1


def read_value(events, first):
    """Build one selected encoded polynomial value."""
    prefix, event, value = first
    if event == "start_array":
        result = []
        while True:
            current = next(events)
            if current[0] == prefix and current[1] == "end_array":
                return result
            result.append(read_value(events, current))
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
                    raise ValueError("malformed map in checkpoint")
                result[key] = read_value(events, current)
                key = None
    if event in ("end_array", "end_map", "map_key"):
        raise ValueError(f"unexpected JSON event {event!r} at {prefix!r}")
    return value


def summarize(tag, polynomial, generators):
    active = set()
    max_total_degree = 0
    coefficient_digests = hashlib.sha256()
    for sparse, values in polynomial:
        degree = 0
        for index, power in sparse:
            active.add(index)
            degree += power
        max_total_degree = max(max_total_degree, degree)
        coefficient_digests.update(",".join(values).encode("ascii"))
        coefficient_digests.update(b"\n")
    return {
        "tag": tag,
        "terms": len(polynomial),
        "active_generator_indices": sorted(active),
        "active_generators": [generators[i] for i in sorted(active)],
        "max_total_degree": max_total_degree,
        "coefficient_sequence_sha256": coefficient_digests.hexdigest(),
    }


def extract(source: Path, destination: Path) -> None:
    generators = None
    found = {}
    with source.open("rb") as source_stream, destination.open(
        "w", encoding="utf-8", newline="\n"
    ) as output:
        output.write("{\"source_checkpoint\":")
        json.dump(source.name, output)
        output.write(",\"source_checkpoint_sha256\":"
                     "\"5b705296d05d6405645217f3d2413ebc963bc7529195bc4d5795967b92ded995\",")
        output.write("\"constraints\":[")
        first_output = True

        events = iter(ijson.parse(source_stream))
        first = next(events)
        if first[1] != "start_map":
            raise ValueError("checkpoint root is not an object")

        while True:
            current = next(events)
            if current[0] == "" and current[1] == "end_map":
                break
            if current[0] != "" or current[1] != "map_key":
                raise ValueError(f"unexpected top-level event {current!r}")
            key = current[2]
            first_value = next(events)

            if key == "generators":
                generators = read_value(events, first_value)
            elif key == "constraints":
                if first_value[1] != "start_array":
                    raise ValueError("constraints is not an array")
                while True:
                    item = next(events)
                    if item[0] == "constraints" and item[1] == "end_array":
                        break
                    if item[1] != "start_array":
                        raise ValueError("malformed constraint entry")
                    tag_event = next(events)
                    tag = tag_event[2]
                    polynomial_event = next(events)
                    if tag in TARGET_TAGS:
                        polynomial = read_value(events, polynomial_event)
                        if generators is None:
                            raise ValueError("generator list was not read")
                        summary = summarize(tag, polynomial, generators)
                        found[tag] = summary
                        if not first_output:
                            output.write(",")
                        first_output = False
                        output.write("{\"summary\":")
                        json.dump(summary, output, separators=(",", ":"))
                        output.write(",\"encoded_polynomial\":")
                        json.dump(polynomial, output, separators=(",", ":"))
                        output.write("}")
                    else:
                        skip_value(events, polynomial_event)
                    end_item = next(events)
                    if end_item[1] != "end_array":
                        raise ValueError("malformed constraint tuple")
                break
            else:
                skip_value(events, first_value)

        output.write("],\"found_tags\":")
        json.dump(sorted(found), output, separators=(",", ":"))
        output.write("}\n")

    missing = sorted(TARGET_TAGS - set(found))
    if missing:
        raise ValueError(f"missing target constraints: {missing}")
    print(json.dumps({"source": str(source), "destination": str(destination),
                      "found": [found[tag] for tag in sorted(found)]},
                     indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: extract_terminalq_residuals.py INPUT.json OUTPUT.json")
    extract(Path(sys.argv[1]), Path(sys.argv[2]))
