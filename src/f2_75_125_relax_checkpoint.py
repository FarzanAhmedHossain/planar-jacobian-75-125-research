"""Make a sound relaxed copy of a sparse F2 checkpoint.

The exact lift recurrence can continue without replaying its old residual
equations: retaining the checkpoint's already-solved state polynomials while
dropping those equations enlarges the parameter space.  Thus a contradiction
in the continued relaxed system is still a valid exclusion, while a surviving
system proves nothing about the original branch.

This utility preserves every byte of the state after the constraints array,
replaces only that array by [], and writes an audit sidecar.  It is specific to
the checkpoint writer's top-level order, where constraints is followed by h.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import mmap
import os
from pathlib import Path


CONSTRAINTS_KEY = b',"constraints":['
CONSTRAINTS_END = b'],"h":'
CHUNK = 64 << 20


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--source-sha256", required=True)
    args = parser.parse_args()

    source = args.source.resolve()
    destination = args.destination.resolve()
    sidecar = destination.with_suffix(destination.suffix + ".relaxation.json")
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    expected = args.source_sha256.lower()
    if len(expected) != 64 or any(c not in "0123456789abcdef" for c in expected):
        raise SystemExit("--source-sha256 must be a 64-digit hexadecimal SHA-256")
    if not source.is_file():
        raise SystemExit(f"source does not exist: {source}")
    if destination.exists() or sidecar.exists() or temporary.exists():
        raise SystemExit("refusing to overwrite an existing destination or sidecar")

    source_bytes = source.stat().st_size
    if source_bytes == 0:
        raise SystemExit("source checkpoint is empty")

    output_hash = hashlib.sha256()
    with source.open("rb") as src, mmap.mmap(src.fileno(), 0, access=mmap.ACCESS_READ) as data:
        key_pos = data.find(CONSTRAINTS_KEY)
        if key_pos < 0:
            raise SystemExit("top-level constraints array marker was not found")
        metadata_prefix = data[:key_pos]
        try:
            metadata = json.loads(metadata_prefix + b"}")
        except json.JSONDecodeError as exc:
            raise SystemExit(f"checkpoint metadata prefix is invalid JSON: {exc}")
        if metadata.get("format") != "f2-sparse-ring-v1":
            raise SystemExit("unexpected checkpoint format")
        if metadata.get("branch") != "2+1+1":
            raise SystemExit("this relaxation helper is restricted to the 2+1+1 branch")
        if not isinstance(metadata.get("receipt"), list) or not metadata["receipt"]:
            raise SystemExit("expected a nonempty checkpoint receipt list")

        array_start = key_pos + len(CONSTRAINTS_KEY) - 1
        if data[array_start:array_start + 1] != b"[":
            raise SystemExit("constraints marker did not end at an array opening")
        end_pos = data.find(CONSTRAINTS_END, array_start + 1)
        if end_pos < 0 or data[end_pos + len(CONSTRAINTS_END):end_pos + len(CONSTRAINTS_END) + 1] != b"[":
            raise SystemExit("could not validate the constraints-to-h state boundary")

        omitted_bytes = end_pos + 1 - array_start
        metadata["receipt"].append({
            "step": "resume-relaxation",
            "source_sha256": expected,
            "prior_constraints": "omitted",
            "omitted_constraint_array_bytes": omitted_bytes,
            "soundness": (
                "state polynomials retained; dropping saved equations enlarges "
                "the solution set, so only a contradiction in the relaxed "
                "continuation can exclude the original branch"
            ),
        })
        header = (json.dumps(metadata, ensure_ascii=True, separators=(",", ":"))[:-1]
                  + ',"constraints":[]').encode("ascii")

        destination.parent.mkdir(parents=True, exist_ok=True)
        try:
            with temporary.open("xb") as out:
                out.write(header)
                output_hash.update(header)
                # Drop the original closing bracket, retaining the comma and
                # the untouched h/v/S/... state following the constraints.
                position = end_pos + 1
                while position < source_bytes:
                    stop = min(position + CHUNK, source_bytes)
                    block = data[position:stop]
                    out.write(block)
                    output_hash.update(block)
                    position = stop
                out.flush()
                os.fsync(out.fileno())
            os.replace(temporary, destination)
        except BaseException:
            if temporary.exists():
                temporary.unlink()
            raise

    output_bytes = destination.stat().st_size
    audit = {
        "operation": "relax_saved_constraints",
        "source": str(source),
        "source_sha256": expected,
        "source_bytes": source_bytes,
        "destination": str(destination),
        "destination_sha256": output_hash.hexdigest(),
        "destination_bytes": output_bytes,
        "constraints_array_bytes_omitted": omitted_bytes,
        "state_after_constraints_preserved_byte_for_byte": True,
        "interpretation": (
            "The destination encodes a necessary-system superset. A proven "
            "contradiction there excludes the source branch; compatibility "
            "does not establish a candidate or an exclusion."
        ),
    }
    sidecar.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2), flush=True)


if __name__ == "__main__":
    main()
