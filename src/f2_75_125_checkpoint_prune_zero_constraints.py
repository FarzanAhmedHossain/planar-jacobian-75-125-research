"""Remove identically zero equations from an exact sparse-lift checkpoint."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def prune(source: Path, output: Path, receipt_path: Path) -> dict:
    raw = source.read_bytes()
    state = json.loads(raw)
    if state.get("format") != "f2-sparse-ring-v1":
        raise ValueError("unexpected checkpoint format")
    if output.exists() or receipt_path.exists():
        raise FileExistsError("refusing to overwrite an existing artifact")

    constraints = state.get("constraints", [])
    old_count = len(constraints)
    removed = sum(not polynomial for _tag, polynomial in constraints)
    state["constraints"] = [(tag, polynomial) for tag, polynomial in constraints
                            if polynomial]
    if len(state["constraints"]) != old_count - removed:
        raise AssertionError("zero-row pruning count mismatch")

    state["compaction"] = dict(state.get("compaction") or {})
    state["compaction"]["additional_zero_constraint_pruning"] = {
        "source_checkpoint": str(source),
        "source_checkpoint_sha256": hashlib.sha256(raw).hexdigest(),
        "zero_constraints_removed": removed,
        "constraints_before": old_count,
        "constraints_after": len(state["constraints"]),
        "scope": "Removed only empty polynomial equations; zero locus unchanged.",
    }
    encoded = (json.dumps(state, separators=(",", ":")) + "\n").encode("utf-8")
    output.write_bytes(encoded)
    receipt = {
        **state["compaction"]["additional_zero_constraint_pruning"],
        "output_checkpoint": str(output),
        "output_checkpoint_sha256": hashlib.sha256(encoded).hexdigest(),
        "output_bytes": len(encoded),
    }
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("receipt", type=Path)
    args = parser.parse_args()
    print(json.dumps(prune(args.source, args.output, args.receipt), indent=2))


if __name__ == "__main__":
    main()
