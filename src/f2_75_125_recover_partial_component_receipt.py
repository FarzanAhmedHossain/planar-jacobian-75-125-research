"""Reconstruct a provenance sidecar for an exact partial F3/F4 checkpoint.

The component lift writes its sparse checkpoint after each distance but only
writes the JSON receipt on normal completion.  If an interrupted run left a
valid checkpoint without that final receipt, this program loads the sparse
state, rechecks its component provenance against the last completed receipt,
recomputes its exact digest, and writes a clearly marked continuation sidecar.
It is not an exclusion certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from f2_75_125_Salpha_component_lift_20261018 import FACTOR_LOG
from f2_75_125_sparse_ring_lift import SparseLift


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def recover(state_path: Path, template_path: Path, out_path: Path,
            component: int) -> dict:
    if out_path.exists():
        raise FileExistsError(f"refusing to overwrite {out_path}")
    template = json.loads(template_path.read_text(encoding="utf-8"))
    expected_component = f"F{component}=0"
    expected_tag = f"Salpha-component-F{component}"
    factor_hash = sha256(FACTOR_LOG)
    if (template.get("component") != expected_component
            or template.get("component_tag") != expected_tag
            or template.get("factor_log_sha256") != factor_hash):
        raise ValueError("template is not the verified receipt for this component")

    lift = SparseLift.from_checkpoint(state_path, verify_saved_terminal_rows=False)
    if lift.tag != "2+2" or not lift.compaction:
        raise ValueError("expected a compacted quadratic-field F2 2+2 state")
    original_hash = lift.compaction.get("source_checkpoint_sha256")
    if original_hash != template.get("original_source_checkpoint_sha256"):
        raise ValueError("partial state and template have different source provenance")
    tagged = [poly for tag, poly in lift.constraints if tag == expected_tag]
    if len(tagged) != 1:
        raise ValueError("partial state must contain exactly one component tag")
    if (template.get("component_parameter_elimination") is not None
            and tagged[0] != lift.ring.zero):
        raise ValueError("component parameter was reportedly eliminated but its row is not zero")

    numeric_steps = [entry["step"] for entry in lift.receipt
                     if isinstance(entry.get("step"), int)]
    if not numeric_steps:
        raise ValueError("partial state receipt has no numeric distance")
    through_distance = max(numeric_steps)
    receipt = dict(template)
    receipt.update({
        "source_checkpoint": str(state_path),
        "source_checkpoint_sha256": sha256(state_path),
        "factor_log_sha256": factor_hash,
        "original_source_checkpoint_sha256": original_hash,
        "through_distance": through_distance,
        "parameter_count": lift.nparameters,
        "constraint_count": len(lift.constraints),
        "constraint_term_count": sum(len(poly) for _, poly in lift.constraints),
        "blocks": lift.receipt,
        "digest": lift.digest(),
        "recovered_from_partial_checkpoint": True,
        "scope": ("recovered exact partial checkpoint sidecar for continuation; "
                  "not an exclusion certificate"),
    })
    out_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return {
        "state": str(state_path),
        "sidecar": str(out_path),
        "component": expected_component,
        "through_distance": through_distance,
        "parameters": lift.nparameters,
        "constraints": len(lift.constraints),
        "constraint_terms": receipt["constraint_term_count"],
        "state_sha256": receipt["source_checkpoint_sha256"],
        "digest": receipt["digest"],
        "component_row_is_zero": not tagged[0],
        "scope": receipt["scope"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--component", type=int, choices=(3, 4), required=True)
    args = parser.parse_args()
    print(json.dumps(recover(args.state, args.template, args.out, args.component),
                     indent=2), flush=True)


if __name__ == "__main__":
    main()
