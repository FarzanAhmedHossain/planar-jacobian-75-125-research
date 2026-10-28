#!/usr/bin/env python3
"""Streaming checkpoint loader for exact F2 sparse-lift continuations.

The original checkpoint loader materializes a large JSON state before
decoding it.  This runner streams one encoded polynomial at a time into the
selected sparse-ring backend, then delegates the recurrence and receipt to
the source lift module unchanged.  It is intended for memory-constrained
continuations; the checkpoint format and mathematics are unchanged.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import ijson
import sympy as sp
from sympy.polys.domains import QQ

import f2_75_125_sparse_ring_lift as lift

JSON_ENCODER = json.JSONEncoder(separators=(",", ":"))


def consume_value(events, first):
    """Build one JSON value from an ijson event stream."""
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


def decode_polynomial(events, first, old_names, field, ring):
    """Decode one sparse polynomial without retaining its JSON term list."""
    prefix, event, _ = first
    if event == "null":
        return None
    if event != "start_array":
        raise ValueError(f"expected encoded polynomial at {prefix!r}")
    accumulated = ring.zero
    terms = {}
    chunk_limit = 32768
    build_chunk = getattr(ring, "from_prevalidated_dict", ring.from_dict)
    while True:
        current = next(events)
        if current[0] == prefix and current[1] == "end_array":
            break
        sparse, values = consume_value(events, current)
        powers = [0] * len(old_names)
        for old_index, exponent in sparse:
            powers[int(old_index)] = int(exponent)
        if field == QQ:
            coefficient = field.dtype(values[0])
        else:
            ground_coefficients = [field.dom.dtype(value) for value in values]
            coefficient = field.dtype(ground_coefficients, field.mod, field.dom)
        terms[tuple(powers)] = coefficient
        if len(terms) >= chunk_limit:
            accumulated += build_chunk(terms)
            terms = {}
    if terms:
        accumulated += build_chunk(terms)
    return accumulated


def streaming_from_checkpoint(cls, path, *, verify_saved_terminal_rows=True):
    """Resume a saved lift while keeping only one JSON polynomial in memory."""
    path = Path(path)
    events = iter(ijson.parse(path.open("rb")))
    first = next(events)
    if first[1] != "start_map":
        raise ValueError("checkpoint root is not an object")

    meta = {}
    obj = cls.__new__(cls)
    obj.constraints = []
    obj.p, obj.q, obj.w = [], [], []
    obj.uk = {}
    ring = None
    field = None
    old_names = None

    polynomial_keys = {"h", "v", "S", "B", "kappa", "kappa2", "kappa3"}
    while True:
        current = next(events)
        prefix, event, value = current
        if prefix == "" and event == "end_map":
            break
        if prefix != "" or event != "map_key":
            raise ValueError(f"unexpected top-level JSON event {current!r}")
        key = value
        first_value = next(events)

        if key in ("branch", "quadratic", "generators", "nparameters",
                   "receipt", "compaction"):
            meta[key] = consume_value(events, first_value)
            if key == "generators":
                old_names = meta[key]
            if key == "quadratic":
                field = QQ.algebraic_field(sp.sqrt(5)) if meta[key] else QQ
            if old_names is not None and field is not None and ring is None:
                obj.tag = meta.get("branch")
                obj.field = field
                obj.ring, *generators = lift.ring(old_names, field)
                indices = {name: i for i, name in enumerate(old_names)}
                if "t" not in indices:
                    raise ValueError("checkpoint has no t generator")
                obj.t = generators[indices["t"]]
                parameter_names = sorted(
                    (name for name in old_names
                     if name.startswith("p") and name[1:].isdigit()),
                    key=lambda name: int(name[1:]),
                )
                if parameter_names != [f"p{i}" for i in range(len(parameter_names))]:
                    raise ValueError("checkpoint parameters are not contiguous from p0")
                unknown_names = sorted(
                    (name for name in old_names
                     if name.startswith("x") and name[1:].isdigit()),
                    key=lambda name: int(name[1:]),
                )
                if unknown_names != [f"x{i}" for i in range(73)]:
                    raise ValueError("checkpoint must retain x0 through x72")
                obj.parameters = [generators[indices[name]] for name in parameter_names]
                obj.unknowns = [generators[indices[name]] for name in unknown_names]
                obj.new_start = 1 + len(obj.parameters)
                obj.affine_ring = (
                    obj.ring.prefix(obj.new_start)
                    if getattr(obj.ring, "backend_name", None) == "flint-mpoly"
                    else None
                )
                obj.compaction = meta.get("compaction")
                ring = obj.ring

        elif key in polynomial_keys:
            if ring is None:
                raise ValueError("polynomial data appeared before ring metadata")
            setattr(obj, {"S": "S"}.get(key, key),
                    decode_polynomial(events, first_value, old_names, field, ring))
        elif key in ("constraints", "p", "q", "w"):
            if ring is None:
                raise ValueError("array polynomial data appeared before ring metadata")
            array_prefix = key
            if first_value[1] != "start_array":
                raise ValueError(f"expected array for {key}")
            destination = obj.constraints if key == "constraints" else getattr(obj, key)
            while True:
                item = next(events)
                if item[0] == array_prefix and item[1] == "end_array":
                    break
                if key == "constraints":
                    if item[1] != "start_array":
                        raise ValueError("malformed constraint entry")
                    tag_event = next(events)
                    tag = tag_event[2]
                    polynomial_event = next(events)
                    polynomial = decode_polynomial(
                        events, polynomial_event, old_names, field, ring
                    )
                    end_item = next(events)
                    if end_item[1] != "end_array":
                        raise ValueError("malformed constraint tuple")
                    destination.append((tag, polynomial))
                else:
                    destination.append(
                        decode_polynomial(events, item, old_names, field, ring)
                    )

        elif key == "uk":
            if ring is None or first_value[1] != "start_map":
                raise ValueError("malformed uk map or missing ring")
            while True:
                item = next(events)
                if item[0] == "uk" and item[1] == "end_map":
                    break
                if item[0] != "uk" or item[1] != "map_key":
                    raise ValueError("malformed uk map entry")
                degree = int(item[2])
                polynomial_event = next(events)
                obj.uk[degree] = decode_polynomial(
                    events, polynomial_event, old_names, field, ring
                )

        else:
            meta[key] = consume_value(events, first_value)

    if ring is None or "branch" not in meta:
        raise ValueError("checkpoint is missing ring metadata")
    # Older checkpoints may omit correction scalars that have not yet been
    # allocated at the saved depth.  Match SparseLift.from_checkpoint: those
    # absent coefficients represent zero until their resonance distance.
    for name in ("kappa", "kappa2", "kappa3"):
        if not hasattr(obj, name):
            setattr(obj, name, obj.ring.zero)
    obj.nparameters = meta["nparameters"]
    obj.receipt = meta["receipt"]
    obj.compaction = meta.get("compaction")
    obj.u = obj.t * obj.h
    obj.started = time.monotonic()
    obj.local_C = obj.local_coefficient(obj.local_slice(0, 6), 6)
    zero_powers = (0,) * len(obj.ring.gens)
    c0 = obj.local_C[zero_powers]
    obj.center_a = -obj.local_coefficient(obj.local_slice(1, 5), 5) / (6 * c0)
    obj.center_b = obj.ring.zero
    second = obj.centered_slice(2, 5)
    obj.center_b = -obj.local_coefficient(second, 5) / (6 * c0)
    if verify_saved_terminal_rows:
        for distance in range(len(obj.p)):
            bound = (71 - 5 * distance) // 12
            if bound >= 0 and obj.centered_slice(distance, bound):
                raise AssertionError(f"saved P terminal rows fail at distance {distance}")
    print(f"stream_resumed={path} branch={obj.tag} depth={len(obj.p)-1} "
          f"parameters={obj.nparameters} constraints={len(obj.constraints)}", flush=True)
    return obj


def write_json_value(stream, value):
    # Encode each metadata value or polynomial term into one buffered write.
    # The previous iterencode loop issued several Python-level writes per
    # sparse term, which dominates checkpoint time for billion-byte states.
    stream.write(JSON_ENCODER.encode(value))


def write_polynomial(stream, polynomial, field):
    """Stream an FLINT polynomial by native terms, avoiding a full to_dict."""
    if polynomial is None:
        stream.write("null")
        return
    stream.write("[")
    if hasattr(polynomial, "data") and hasattr(polynomial.data, "terms"):
        # The FLINT adapter's items() walks indexed term storage lazily.
        # fmpq_mpoly.terms() materializes both monomial and coefficient lists;
        # on the depth-16 checkpoint this drove peak private memory above 9 GB.
        terms = (polynomial.items() if hasattr(polynomial, "items")
                 else polynomial.data.terms())
    else:
        terms = iter(sorted(polynomial.items()))
    first = True
    for powers, coefficient in terms:
        if not first:
            stream.write(",")
        first = False
        sparse = [[i, int(power)] for i, power in enumerate(powers) if power]
        values = [coefficient] if field == QQ else coefficient.to_list()
        write_json_value(stream, [sparse, [str(value) for value in values]])
    stream.write("]")


def streaming_save_checkpoint(self, path):
    """Write the ordinary checkpoint schema without building a JSON mirror."""
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        stream.write("{")
        scalar_fields = (
            ("format", "f2-sparse-ring-v1"),
            ("branch", self.tag),
            ("quadratic", self.field != QQ),
            ("polynomial_backend", getattr(self.ring, "backend_name", "sympy")),
            ("generators", [str(symbol) for symbol in self.ring.symbols]),
            ("nparameters", self.nparameters),
            ("receipt", self.receipt),
        )
        first = True
        for key, value in scalar_fields:
            if not first:
                stream.write(",")
            first = False
            write_json_value(stream, key)
            stream.write(":")
            write_json_value(stream, value)

        def field_name(key):
            nonlocal first
            if not first:
                stream.write(",")
            first = False
            write_json_value(stream, key)
            stream.write(":")

        field_name("constraints")
        stream.write("[")
        for index, (tag, polynomial) in enumerate(self.constraints):
            if index:
                stream.write(",")
            stream.write("[")
            write_json_value(stream, tag)
            stream.write(",")
            write_polynomial(stream, polynomial, self.field)
            stream.write("]")
        stream.write("]")

        for key in ("h", "v", "S", "B", "kappa", "kappa2", "kappa3"):
            field_name(key)
            write_polynomial(stream, getattr(self, key), self.field)

        field_name("uk")
        stream.write("{")
        for index, degree in enumerate(sorted(self.uk)):
            if index:
                stream.write(",")
            write_json_value(stream, str(degree))
            stream.write(":")
            write_polynomial(stream, self.uk[degree], self.field)
        stream.write("}")

        for key in ("p", "q", "w"):
            field_name(key)
            stream.write("[")
            for index, polynomial in enumerate(getattr(self, key)):
                if index:
                    stream.write(",")
                write_polynomial(stream, polynomial, self.field)
            stream.write("]")

        if getattr(self, "compaction", None) is not None:
            field_name("compaction")
            write_json_value(stream, self.compaction)
        stream.write("}\n")
    temporary.replace(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    self.last_checkpoint_sha256 = digest.hexdigest()
    print(f"stream_checkpoint={path} bytes={path.stat().st_size} "
          f"sha256={self.last_checkpoint_sha256}", flush=True)


def main():
    lift.SparseLift.from_checkpoint = classmethod(streaming_from_checkpoint)
    lift.SparseLift.save_checkpoint = streaming_save_checkpoint
    lift.main()


if __name__ == "__main__":
    main()
