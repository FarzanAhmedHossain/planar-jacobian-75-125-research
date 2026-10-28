#!/usr/bin/env python3
"""Stream an exact parameter relabel over a sparse F2 lift checkpoint.

This avoids materializing a multi-gigabyte JSON document.  Every sparse
monomial is remapped by generator name and checked through an inverse map;
coefficient encodings are copied as parsed JSON values without arithmetic.
The resulting checkpoint preserves the input's relaxed/retained equation
system and is suitable for the ordinary streaming loader.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

from f2_75_125_211_streaming_checkpoint_scan import JsonStream


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        while chunk := source.read(8 << 20):
            digest.update(chunk)
    return digest.hexdigest()


class StreamWriter:
    def __init__(self, path: Path):
        self.path = path
        self.handle = path.open("wb")
        self.digest = hashlib.sha256()
        self.parts: list[str] = []
        self.characters = 0
        self.byte_count = 0

    def emit(self, value: str) -> None:
        self.parts.append(value)
        self.characters += len(value)
        if self.characters >= (1 << 20):
            self.flush()

    def flush(self) -> None:
        if not self.parts:
            return
        data = "".join(self.parts).encode("ascii")
        self.handle.write(data)
        self.digest.update(data)
        self.byte_count += len(data)
        self.parts.clear()
        self.characters = 0

    def close(self) -> None:
        self.flush()
        self.handle.flush()
        os.fsync(self.handle.fileno())
        self.handle.close()


def compact(source: Path, support_path: Path, output: Path, receipt_path: Path,
            reserve_new_parameters: int) -> dict:
    source = source.resolve()
    support_path = support_path.resolve()
    output = output.resolve()
    receipt_path = receipt_path.resolve()
    temporary = output.with_suffix(output.suffix + ".partial")
    receipt_temporary = receipt_path.with_suffix(receipt_path.suffix + ".partial")
    if source == output or source == temporary:
        raise ValueError("output must not replace the source checkpoint")
    if reserve_new_parameters < 0:
        raise ValueError("reserve must be nonnegative")
    for path in (output, receipt_path, temporary, receipt_temporary):
        if path.exists():
            raise FileExistsError(f"refusing to overwrite {path}")

    support = json.loads(support_path.read_text(encoding="utf-8"))
    input_sha = sha256_file(source)
    if support.get("source_sha256", "").lower() != input_sha:
        raise ValueError("support scan does not match the current source SHA-256")
    if support.get("constraint_rows_scanned") != 0:
        raise ValueError("this compact path expects the verified empty-constraint relaxed checkpoint")

    reader = JsonStream(source)
    writer = StreamWriter(temporary)
    try:
        metadata, _ = reader.top_level_constraints()
        if metadata.get("format") != "f2-sparse-ring-v1":
            raise ValueError("unexpected sparse checkpoint format")
        if metadata.get("branch") != support.get("branch"):
            raise ValueError("support branch does not match checkpoint")

        old_names = metadata.get("generators")
        if not isinstance(old_names, list) or not old_names:
            raise ValueError("checkpoint generator list is missing")
        old_index = {name: index for index, name in enumerate(old_names)}
        if len(old_index) != len(old_names):
            raise ValueError("checkpoint has duplicate generator names")
        old_parameters = sorted(
            (name for name in old_names
             if name.startswith("p") and name[1:].isdigit()),
            key=lambda name: int(name[1:]),
        )
        if old_parameters != [f"p{i}" for i in range(len(old_parameters))]:
            raise ValueError("source parameter slots are not contiguous")
        if old_names != ["t", *old_parameters, *(f"x{i}" for i in range(73))]:
            raise ValueError("source generators are not ordered as t, p-slots, x0..x72")
        if len(old_parameters) != support.get("parameter_slot_count"):
            raise ValueError("support scan parameter-slot count does not match source")

        active_names = sorted(
            support.get("active_parameters", []),
            key=lambda name: int(name[1:]) if name.startswith("p") else -1,
        )
        if (not active_names or len(set(active_names)) != len(active_names)
                or any(name not in old_index for name in active_names)
                or any(not name.startswith("p") or not name[1:].isdigit()
                       for name in active_names)):
            raise ValueError("support scan has an invalid active-parameter list")
        if len(active_names) != support.get("active_parameter_count"):
            raise ValueError("support active-parameter count does not match its list")
        old_nparameters = int(metadata["nparameters"])
        if any(int(name[1:]) >= old_nparameters for name in active_names):
            raise ValueError("an active parameter lies beyond the allocated prefix")

        active_map = {name: f"p{i}" for i, name in enumerate(active_names)}
        capacity = len(active_names) + reserve_new_parameters
        new_names = ["t", *(f"p{i}" for i in range(capacity)),
                     *(f"x{i}" for i in range(73))]
        old_to_new_index: dict[int, int] = {old_index["t"]: 0}
        new_to_old_index: dict[int, int] = {0: old_index["t"]}
        for old_name, new_name in active_map.items():
            old_i, new_i = old_index[old_name], 1 + int(new_name[1:])
            old_to_new_index[old_i] = new_i
            new_to_old_index[new_i] = old_i
        for i in range(73):
            old_i = old_index[f"x{i}"]
            new_i = 1 + capacity + i
            old_to_new_index[old_i] = new_i
            new_to_old_index[new_i] = old_i

        metadata["generators"] = new_names
        metadata["nparameters"] = len(active_names)

        writer.emit("{")
        first_field = True

        def emit_field(key: str, value_writer) -> None:
            nonlocal first_field
            if not first_field:
                writer.emit(",")
            first_field = False
            writer.emit(json.dumps(key, separators=(",", ":")))
            writer.emit(":")
            value_writer()

        def emit_json(value) -> None:
            writer.emit(json.dumps(value, separators=(",", ":"),
                                   ensure_ascii=True, allow_nan=False))

        for key, value in metadata.items():
            emit_field(key, lambda value=value: emit_json(value))

        # The input is a prior-constraint-relaxed necessary-system superset.
        # Certify the constraints are truly empty before preserving that scope.
        if reader.peek() != "]":
            raise ValueError("source constraints are not empty")
        reader.take("]")
        emit_field("constraints", lambda: writer.emit("[]"))
        reader.take(",")

        counts = {"polynomials": 0, "terms": 0, "inverse_checks": 0}
        next_progress = 1_000_000

        def write_polynomial(field_name: str) -> None:
            nonlocal next_progress
            reader.take("[")
            writer.emit("[")
            first_term = True
            while reader.peek() != "]":
                if not first_term:
                    reader.take(",")
                    writer.emit(",")
                first_term = False
                reader.take("[")
                sparse = reader.value()
                reader.take(",")
                coefficients = reader.value()
                reader.take("]")

                mapped = []
                for old_generator_index, exponent in sparse:
                    if exponent == 0:
                        raise ValueError("sparse term unexpectedly stores a zero exponent")
                    try:
                        mapped.append([old_to_new_index[old_generator_index], exponent])
                    except KeyError as exc:
                        old_name = old_names[old_generator_index]
                        raise ValueError(
                            f"support scan omitted active generator {old_name} in {field_name}"
                        ) from exc
                mapped.sort(key=lambda pair: pair[0])
                back = sorted(
                    ([new_to_old_index[new_index], exponent]
                     for new_index, exponent in mapped),
                    key=lambda pair: pair[0],
                )
                if back != sorted(sparse, key=lambda pair: pair[0]):
                    raise AssertionError("monomial inverse relabel check failed")

                writer.emit("[")
                emit_json(mapped)
                writer.emit(",")
                emit_json(coefficients)
                writer.emit("]")
                counts["terms"] += 1
                counts["inverse_checks"] += 1
                if counts["terms"] >= next_progress:
                    print(
                        f"stream_compaction_terms={counts['terms']} field={field_name} "
                        f"input_offset={reader.absolute_position()}",
                        flush=True,
                    )
                    next_progress += 1_000_000
            reader.take("]")
            writer.emit("]")
            counts["polynomials"] += 1

        def write_poly_array(field_name: str) -> None:
            reader.take("[")
            writer.emit("[")
            first = True
            index = 0
            while reader.peek() != "]":
                if not first:
                    reader.take(",")
                    writer.emit(",")
                first = False
                write_polynomial(f"{field_name}[{index}]")
                index += 1
            reader.take("]")
            writer.emit("]")

        def write_poly_map(field_name: str) -> None:
            reader.take("{")
            writer.emit("{")
            first = True
            while reader.peek() != "}":
                if not first:
                    reader.take(",")
                    writer.emit(",")
                first = False
                key = reader.value()
                reader.take(":")
                emit_json(key)
                writer.emit(":")
                write_polynomial(f"{field_name}[{key}]")
            reader.take("}")
            writer.emit("}")

        single_polynomials = {"h", "v", "S", "B", "kappa", "kappa2", "kappa3"}
        sequence_polynomials = {"p", "q", "w"}
        seen = set()
        compaction_manifest = {
            "format": "exact-parameter-relabel-v1",
            "source_checkpoint": str(source),
            "source_checkpoint_sha256": input_sha,
            "old_nparameters": old_nparameters,
            "new_nparameters": len(active_names),
            "active_old_parameters": active_names,
            "old_to_new_parameter_names": active_map,
            "reserved_new_parameters": reserve_new_parameters,
            "generator_count_before": len(old_names),
            "generator_count_after": len(new_names),
            "exact_relabeling_verified": True,
        }
        previous_compaction = None
        while reader.peek() != "}":
            if seen:
                reader.take(",")
            key = reader.value()
            reader.take(":")
            if key == "compaction":
                previous_compaction = reader.value()
                continue
            if key in single_polynomials:
                if key == "B" and reader.peek() == "n":
                    emit_field(key, lambda: emit_json(reader.value()))
                else:
                    emit_field(key, lambda key=key: write_polynomial(key))
            elif key in sequence_polynomials:
                emit_field(key, lambda key=key: write_poly_array(key))
            elif key == "uk":
                emit_field(key, lambda: write_poly_map("uk"))
            else:
                raise ValueError(f"unexpected checkpoint field {key!r}")
            seen.add(key)
        reader.take("}")
        if reader.peek():
            raise ValueError("unexpected trailing data after source checkpoint")
        required = {"h", "v", "S", "kappa", "uk", "p", "q", "w"}
        if not required.issubset(seen):
            raise ValueError(f"checkpoint misses required fields: {sorted(required-seen)}")
        if counts["terms"] != support.get("serialized_terms_scanned"):
            raise ValueError("streamed term count differs from active-support audit")
        if counts["polynomials"] != support.get("polynomial_arrays_scanned"):
            raise ValueError("streamed polynomial count differs from active-support audit")
        if counts["inverse_checks"] != counts["terms"]:
            raise AssertionError("not every serialized monomial passed its inverse check")

        compaction_manifest["term_count_before_and_after"] = counts["terms"]
        compaction_manifest["polynomial_count"] = counts["polynomials"]
        compaction_manifest["monomial_inverse_checks"] = counts["inverse_checks"]
        compaction_manifest["previous_compaction"] = previous_compaction
        compaction_manifest["exact_relabeling_verified"] = (
            counts["inverse_checks"] == counts["terms"]
        )
        writer.emit(",")
        writer.emit(json.dumps("compaction"))
        writer.emit(":")
        emit_json(compaction_manifest)
        writer.emit("}\n")
        writer.close()
    except BaseException:
        reader.close()
        if not writer.handle.closed:
            writer.handle.close()
        raise
    else:
        reader.close()

    output_sha = writer.digest.hexdigest()
    output_bytes = writer.byte_count
    os.replace(temporary, output)
    receipt = {
        **compaction_manifest,
        "output_checkpoint": str(output),
        "output_checkpoint_sha256": output_sha,
        "output_bytes": output_bytes,
        "active_old_generator_indices": sorted(
            old_index[name] for name in active_names
        ),
        "input_constraints": 0,
        "scope": (
            "Exact coordinate relabel of the source's relaxed necessary-system "
            "superset; no equation was added or removed. Compatibility does not "
            "exclude any branch."
        ),
    }
    with receipt_temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(receipt, handle, indent=2, ensure_ascii=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(receipt_temporary, receipt_path)
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("support", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--reserve", type=int, default=21,
                        help="fresh parameter slots following the active relabelled parameters")
    args = parser.parse_args()
    print(json.dumps(compact(args.source, args.support, args.output,
                             args.receipt, args.reserve), indent=2), flush=True)


if __name__ == "__main__":
    main()
