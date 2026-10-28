"""Temporarily detach and exactly reattach a sparse checkpoint's constraints.

This supports memory-bounded continuation of the coefficient recurrence.  It
is safe only when the detached constraint polynomials depend on already-active
parameters and the continuation introduces fresh parameters; the streaming
profile must certify that precondition.  The merge copies the original JSON
constraint array byte-for-byte before adding the new constraints.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import mmap
import shutil
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def skip_ws(text: str, pos: int) -> int:
    while pos < len(text) and text[pos].isspace():
        pos += 1
    return pos


def skip_json_value(text: str, pos: int, decoder: json.JSONDecoder) -> int:
    """Return the end of one JSON value without decoding nested arrays."""
    if text[pos] not in "[{":
        _value, end = decoder.raw_decode(text, pos)
        return end
    closing = {"[": "]", "{": "}"}
    stack = [closing[text[pos]]]
    in_string = False
    escaped = False
    i = pos + 1
    while stack:
        if i >= len(text):
            raise ValueError("unterminated JSON composite")
        char = text[i]
        if in_string:
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == '"':
                in_string = False
        elif char == '"':
            in_string = True
        elif char in "[{":
            stack.append(closing[char])
        elif char in "]}":
            if not stack or char != stack[-1]:
                raise ValueError("unbalanced JSON composite")
            stack.pop()
        i += 1
    return i


def strip_constraints(source: Path, stripped: Path, detached: Path, receipt: Path) -> None:
    for path in (stripped, detached, receipt):
        if path.exists():
            raise FileExistsError(path)
    text = source.read_text(encoding="utf-8")
    decoder = json.JSONDecoder()
    pos = skip_ws(text, 0)
    if text[pos] != "{":
        raise ValueError("checkpoint root is not an object")
    pos += 1
    state = {}
    raw_constraints = None
    while True:
        pos = skip_ws(text, pos)
        if text[pos] == "}":
            break
        key, pos = decoder.raw_decode(text, pos)
        pos = skip_ws(text, pos)
        if text[pos] != ":":
            raise ValueError("malformed top-level JSON object")
        value_pos = skip_ws(text, pos + 1)
        end = skip_json_value(text, value_pos, decoder)
        if key == "constraints":
            if text[value_pos] != "[":
                raise ValueError("constraints field is not an array")
            raw_constraints = text[value_pos:end].encode("utf-8")
            state[key] = []
        else:
            state[key] = json.loads(text[value_pos:end])
        pos = skip_ws(text, end)
        if text[pos] == ",":
            pos += 1
        elif text[pos] == "}":
            break
        else:
            raise ValueError("malformed top-level separator")

    if raw_constraints is None:
        raise ValueError("checkpoint has no constraints field")
    if state.get("format") != "f2-sparse-ring-v1":
        raise ValueError("unexpected sparse checkpoint format")
    if state.get("quadratic"):
        raise ValueError("this helper is restricted to rational checkpoints")
    if not raw_constraints.startswith(b"[") or not raw_constraints.endswith(b"]"):
        raise ValueError("detached constraints are not a JSON array")

    stripped.write_text(
        json.dumps(state, separators=(",", ":")) + "\n", encoding="utf-8"
    )
    detached.write_bytes(raw_constraints)
    receipt.write_text(json.dumps({
        "source_checkpoint": str(source),
        "source_sha256": sha256_file(source),
        "branch": state["branch"],
        "active_parameter_count": state["nparameters"],
        "detached_constraint_array_bytes": len(raw_constraints),
        "detached_constraint_array_sha256": hashlib.sha256(raw_constraints).hexdigest(),
        "stripped_checkpoint": str(stripped),
        "stripped_checkpoint_sha256": sha256_file(stripped),
        "detached_array": str(detached),
        "scope": "Recurrence-only checkpoint; original equations are detached and must be reattached before interpreting the result.",
    }, indent=2) + "\n", encoding="utf-8")
    print(f"stripped={stripped} bytes={stripped.stat().st_size}")
    print(f"detached={detached} bytes={detached.stat().st_size}")
    print(f"receipt={receipt}")


def array_end(mapped: mmap.mmap, start: int) -> int:
    if mapped[start:start + 1] != b"[":
        raise ValueError("expected a JSON array")
    depth = 0
    in_string = False
    escaped = False
    for i in range(start, len(mapped)):
        char = mapped[i]
        if in_string:
            if escaped:
                escaped = False
            elif char == 0x5C:
                escaped = True
            elif char == 0x22:
                in_string = False
        elif char == 0x22:
            in_string = True
        elif char in (0x5B, 0x7B):
            depth += 1
        elif char in (0x5D, 0x7D):
            depth -= 1
            if depth == 0:
                return i + 1
            if depth < 0:
                break
    raise ValueError("unterminated JSON array")


def copy_n(source, destination, count: int) -> None:
    remaining = count
    while remaining:
        block = source.read(min(1024 * 1024, remaining))
        if not block:
            raise EOFError("unexpected end of source while copying JSON fragment")
        destination.write(block)
        remaining -= len(block)


def merge_constraints(new_only: Path, detached: Path, merged: Path, receipt: Path) -> None:
    for path in (merged, receipt):
        if path.exists():
            raise FileExistsError(path)
    old_size = detached.stat().st_size
    if old_size < 2:
        raise ValueError("detached constraint array is empty or malformed")
    with detached.open("rb") as old:
        if old.read(1) != b"[":
            raise ValueError("detached array has no opening bracket")
        old.seek(old_size - 1)
        if old.read(1) != b"]":
            raise ValueError("detached array has no closing bracket")

    with new_only.open("rb") as new, mmap.mmap(new.fileno(), 0, access=mmap.ACCESS_READ) as mapped:
        if b'"format":"f2-sparse-ring-v1"' not in mapped[:256]:
            raise ValueError("new checkpoint has an unexpected format")
        if b'"branch":"2+1+1"' not in mapped[:512]:
            raise ValueError("new checkpoint has an unexpected branch")
        key_pos = mapped.find(b'"constraints":[')
        if key_pos < 0:
            raise ValueError("new checkpoint has no compact constraints field")
        start = key_pos + len(b'"constraints":')
        end = array_end(mapped, start)
        new_has_items = end - start > 2
        old_has_items = old_size > 2
        with merged.open("wb") as out:
            new.seek(0)
            copy_n(new, out, start)
            out.write(b"[")
            if old_has_items:
                old.seek(1)
                copy_n(old, out, old_size - 2)
            if old_has_items and new_has_items:
                out.write(b",")
            if new_has_items:
                new.seek(start + 1)
                copy_n(new, out, end - start - 2)
            out.write(b"]")
            new.seek(end)
            shutil.copyfileobj(new, out, length=1024 * 1024)
        old_body_bytes = old_size - 2
        prefix_digest = hashlib.sha256()
        with merged.open("rb") as check:
            check.seek(start + 1)
            remaining = old_body_bytes
            while remaining:
                block = check.read(min(1024 * 1024, remaining))
                if not block:
                    raise EOFError("merged file ended inside copied constraints")
                prefix_digest.update(block)
                remaining -= len(block)
        with detached.open("rb") as old:
            old.seek(1)
            old_hash = hashlib.sha256()
            remaining = old_body_bytes
            while remaining:
                block = old.read(min(1024 * 1024, remaining))
                if not block:
                    raise EOFError("detached constraint body ended early")
                old_hash.update(block)
                remaining -= len(block)
            old_digest = old_hash.hexdigest()
        if prefix_digest.hexdigest() != old_digest:
            raise AssertionError("original constraints were not copied byte-for-byte")

    receipt.write_text(json.dumps({
        "new_only_checkpoint": str(new_only),
        "new_only_sha256": sha256_file(new_only),
        "detached_constraint_array": str(detached),
        "detached_constraint_array_sha256": sha256_file(detached),
        "merged_checkpoint": str(merged),
        "merged_checkpoint_sha256": sha256_file(merged),
        "merged_bytes": merged.stat().st_size,
        "original_constraint_body_prefix_verified": True,
        "scope": "Original constraints are preserved verbatim and concatenated with recurrence continuation constraints; this is not an exclusion certificate.",
    }, indent=2) + "\n", encoding="utf-8")
    print(f"merged={merged} bytes={merged.stat().st_size}")
    print("original_constraint_body_prefix_verified=True")
    print(f"receipt={receipt}")


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    strip = sub.add_parser("strip")
    strip.add_argument("source", type=Path)
    strip.add_argument("stripped", type=Path)
    strip.add_argument("detached", type=Path)
    strip.add_argument("receipt", type=Path)
    merge = sub.add_parser("merge")
    merge.add_argument("new_only", type=Path)
    merge.add_argument("detached", type=Path)
    merge.add_argument("merged", type=Path)
    merge.add_argument("receipt", type=Path)
    args = parser.parse_args()
    if args.command == "strip":
        strip_constraints(args.source, args.stripped, args.detached, args.receipt)
    else:
        merge_constraints(args.new_only, args.detached, args.merged, args.receipt)


if __name__ == "__main__":
    main()
