#!/usr/bin/env python3
"""Measure selected parameter powers in the exact d18 checkpoint."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

from f2_75_125_211_streaming_checkpoint_scan import JsonStream


SOURCE = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_state.json")
EXPECTED_SHA256 = "e2e35061e883437806d6486d8f61fdc4cc6b963c49655366a6795db210d874e7"
TARGETS = ("p13", "p14", "p15")
SELECTED = ("p5", "p34", "p35", "p36", "p37", "p38", "p39")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    source_hash = sha256(SOURCE)
    if source_hash != EXPECTED_SHA256:
        raise ValueError(f"unexpected checkpoint hash: {source_hash}")
    reader = JsonStream(SOURCE)
    try:
        metadata, _ = reader.top_level_constraints()
        generators = metadata["generators"]
        selected_ids = {generators.index(name): name for name in SELECTED}
        maxima = Counter()
        histograms = {name: Counter() for name in SELECTED}
        joint_histogram = Counter()
        joint_maximum = 0
        row_count = term_count = 0
        for _, terms in reader.array_values():
            row_count += 1
            for sparse, _ in terms:
                term_count += 1
                selected_powers = {"p34": 0, "p35": 0}
                for generator, exponent in sparse:
                    name = selected_ids.get(generator)
                    if name is None:
                        continue
                    exponent = int(exponent)
                    maxima[name] = max(maxima[name], exponent)
                    histograms[name][exponent] += 1
                    if name in selected_powers:
                        selected_powers[name] = exponent
                joint = selected_powers["p34"] + selected_powers["p35"]
                joint_histogram[joint] += 1
                joint_maximum = max(joint_maximum, joint)
    finally:
        reader.close()
    report = {
        "source": str(SOURCE),
        "source_sha256": source_hash,
        "branch": metadata["branch"],
        "checkpoint_rows": row_count,
        "sparse_terms": term_count,
        "selected_parameters": list(SELECTED),
        "maximum_exponent_by_parameter": dict(maxima),
        "maximum_joint_exponent_p34_p35": joint_maximum,
        "joint_exponent_occurrence_counts_p34_p35": {
            str(power): count for power, count in sorted(joint_histogram.items())
        },
        "exponent_occurrence_counts": {
            name: {str(power): count for power, count in sorted(hist.items())}
            for name, hist in histograms.items()
        },
        "scope": "Power-size diagnostic only; no equations or terms were projected or changed.",
    }
    output = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support8_parameter_exponents.json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "source_sha256": source_hash,
        "checkpoint_rows": row_count,
        "sparse_terms": term_count,
        "maximum_exponent_by_parameter": report["maximum_exponent_by_parameter"],
        "maximum_joint_exponent_p34_p35": report["maximum_joint_exponent_p34_p35"],
        "report": str(output),
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
