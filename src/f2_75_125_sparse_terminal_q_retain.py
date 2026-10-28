"""Append all available Q-side terminal halfspace equations to a lift state."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from f2_75_125_sparse_ring_lift import SparseLift


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()

    lift = SparseLift.from_checkpoint(args.checkpoint)
    depth = len(lift.q) - 1
    for distance in range(depth + 1):
        lift.impose_q_source_support(distance)
    total_positions = 0
    nonzero_positions = 0
    by_distance = []
    for distance in range(depth + 1):
        bound = (119 - 5 * distance) // 12
        if bound < 0:
            continue
        chart = lift.centered_q_slice(distance, bound)
        count = 0
        for degree in range(bound + 1):
            total_positions += 1
            value = lift.local_coefficient(chart, degree)
            if value:
                count += 1
                nonzero_positions += 1
                lift.constraints.append((f"terminal-Q-distance{distance}:u^{degree}", value))
        by_distance.append({"distance": distance, "bound": bound,
                            "positions": bound + 1, "nonzero": count})

    lift.save_checkpoint(args.output)
    report = {"branch": lift.tag, "depth": depth,
              "candidate_positions": total_positions,
              "nonzero_equations_retained": nonzero_positions,
              "by_distance": by_distance,
              "output_state": str(args.output),
              "scope": "exact necessary Q terminal equations appended to the parameter ideal"}
    args.receipt.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
