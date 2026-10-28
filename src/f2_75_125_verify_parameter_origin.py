"""Verify whether a sparse F2 checkpoint vanishes at the zero-tail point.

This checks every serialized monomial exactly at the structural level: each
nonzero term must contain at least one carried parameter. It is not a test of
the full Keller equation; a zero-tail solution only demonstrates that the
stored necessary subsystem omits the inhomogeneous constant-Jacobian row.
"""

from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path

from f2_75_125_quadratic_affine_profile_stream import (
    iter_json_array,
    top_metadata_and_constraints,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("receipt", type=Path)
    args = parser.parse_args()

    raw = args.checkpoint.read_bytes()
    source_hash = hashlib.sha256(raw).hexdigest()
    script_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    text = raw.decode("utf-8")
    metadata, position = top_metadata_and_constraints(text)
    names = metadata["generators"]
    parameter_names = [f"p{i}" for i in range(int(metadata["nparameters"]))]
    parameter_indices = {names.index(name) for name in parameter_names}

    constraint_count = 0
    term_count = 0
    parameter_degrees: dict[int, int] = {}
    terms_surviving_at_zero: list[dict[str, object]] = []
    for tag, saved in iter_json_array(text, position):
        constraint_count += 1
        for sparse, coefficient in saved:
            term_count += 1
            parameter_degree = sum(
                power for index, power in sparse if index in parameter_indices
            )
            parameter_degrees[parameter_degree] = parameter_degrees.get(parameter_degree, 0) + 1
            if parameter_degree == 0:
                terms_surviving_at_zero.append({
                    "tag": tag,
                    "sparse_monomial": sparse,
                    "coefficient": coefficient,
                })

    receipt = {
        "checkpoint": str(args.checkpoint),
        "checkpoint_sha256": source_hash,
        "verifier_sha256": script_hash,
        "recorded_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "branch": metadata.get("branch"),
        "nparameters": int(metadata["nparameters"]),
        "constraint_count": constraint_count,
        "serialized_term_count": term_count,
        "parameter_degree_histogram": {
            str(key): parameter_degrees[key] for key in sorted(parameter_degrees)
        },
        "terms_surviving_at_all_parameters_zero": len(terms_surviving_at_zero),
        "surviving_term_samples": terms_surviving_at_zero[:10],
        "zero_parameter_point_solves_entire_saved_subsystem": not terms_surviving_at_zero,
        "scope": (
            "Exact check of the saved necessary-equation checkpoint only. If true, "
            "it does not give a Keller pair: the full inhomogeneous constant-Jacobian "
            "equation and branch nonvanishing data are outside this checkpoint."
        ),
    }
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
