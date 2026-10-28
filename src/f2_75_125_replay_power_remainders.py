"""Replay exact corrected-power remainder equations from a sparse lift checkpoint.

Older checkpoints predate serialization of the retained-constraint list.  This
replay reconstructs Q-distance coefficients from the saved P/Q slices and
restores any exact divisibility remainders before a deeper continuation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from f2_75_125_sparse_ring_lift import SparseLift


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()

    lift = SparseLift.from_checkpoint(args.checkpoint)
    original_p = list(lift.p)
    original_q = list(lift.q)
    depth = len(original_p) - 1
    if len(original_q) != len(original_p):
        raise ValueError("saved P and Q slices have different depths")

    recovered = []
    # The chart defines p_0,p_1,q_0,q_1 as base data; new_q is used only
    # for the solved slices starting at distance 2.
    for distance in range(2, depth + 1):
        lift.p = original_p[:distance]
        lift.q = original_q[:distance]
        before = len(lift.constraints)
        quotient = lift.new_q(distance, original_p[distance])
        quotient_matches = quotient == original_q[distance]
        added_tags = [tag for tag, _ in lift.constraints[before:]]
        recovered.append({
            "distance": distance,
            "saved_Q_slice_recovered": quotient_matches,
            "new_remainder_constraints": len(added_tags),
            "tags": added_tags,
        })
        if not quotient_matches:
            raise ArithmeticError(
                f"replayed Q slice at distance {distance} does not match checkpoint"
            )

    lift.p = original_p
    lift.q = original_q

    source_before = len(lift.constraints)
    for distance in range(depth + 1):
        lift.impose_q_source_support(distance)
    source_added = len(lift.constraints) - source_before

    q_positions = 0
    q_terminal_added = 0
    q_by_distance = []
    for distance in range(depth + 1):
        bound = (119 - 5 * distance) // 12
        if bound < 0:
            continue
        chart = lift.centered_q_slice(distance, bound)
        nonzero = 0
        for degree in range(bound + 1):
            q_positions += 1
            value = lift.local_coefficient(chart, degree)
            if value:
                nonzero += 1
                q_terminal_added += 1
                lift.constraints.append(
                    (f"terminal-Q-distance{distance}:u^{degree}", value)
                )
        q_by_distance.append({
            "distance": distance,
            "positions": bound + 1,
            "nonzero_equations": nonzero,
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    lift.save_checkpoint(args.output)
    state_hash = hashlib.sha256(args.output.read_bytes()).hexdigest()
    report = {
        "branch": lift.tag,
        "depth": depth,
        "input_state": str(args.checkpoint),
        "output_state": str(args.output),
        "output_sha256": state_hash,
        "saved_Q_slices_recovered": all(
            row["saved_Q_slice_recovered"] for row in recovered
        ),
        "base_Q_slices_preserved_without_replay": [0, 1],
        "corrected_power_replay": recovered,
        "source_Q_equations_added": source_added,
        "terminal_Q_positions_checked": q_positions,
        "terminal_Q_equations_added": q_terminal_added,
        "terminal_Q_by_distance": q_by_distance,
        "retained_constraint_count": len(lift.constraints),
        "retained_constraint_tags": [tag for tag, _ in lift.constraints],
        "scope": (
            "exact replay of saved corrected-power relations, Q source support, "
            "and Q terminal halfspace conditions; necessary constraints only"
        ),
    }
    args.receipt.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
