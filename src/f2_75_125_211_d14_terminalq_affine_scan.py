#!/usr/bin/env python3
"""Stream terminal-Q constraints and check their degree in p44,...,p50."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

import ijson


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "f2_75_125_sparse_ring_lift_211_roots14_p43zero_d14_terminal_state.json"
OUTPUT = ROOT / "f2_75_125_211_d14_terminalq_affine_scan.json"
EXPECTED_SHA256 = "0412d7637384badd19ff22cc7ce4813636cdabd51897ec13508e0668b7e44272"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    source_hash = sha256(SOURCE)
    if source_hash != EXPECTED_SHA256:
        raise SystemExit(f"checkpoint SHA mismatch: {source_hash}")

    with SOURCE.open("rb") as stream:
        generators = list(ijson.items(stream, "generators.item"))
    new_indices = {
        index for index, name in enumerate(generators)
        if name in {f"p{i}" for i in range(44, 51)}
    }
    if len(new_indices) != 7:
        raise SystemExit(f"expected p44,...,p50, found indices {new_indices}")

    rows = []
    with SOURCE.open("rb") as stream:
        for tag, terms in ijson.items(stream, "constraints.item"):
            if "terminal-Q-distance14" not in tag and "Q-source-support-distance14" not in tag:
                continue
            histogram: Counter[int] = Counter()
            support: set[str] = set()
            for sparse, _values in terms:
                degree = 0
                for index, exponent in sparse:
                    name = generators[int(index)]
                    if int(index) in new_indices:
                        degree += int(exponent)
                    if name.startswith("p"):
                        support.add(name)
                histogram[degree] += 1
            rows.append({
                "tag": tag,
                "term_count": sum(histogram.values()),
                "support_parameters": sorted(support, key=lambda s: int(s[1:])),
                "degree_in_p44_p50": {str(k): v for k, v in sorted(histogram.items())},
                "max_degree_in_p44_p50": max(histogram, default=0),
            })

    result = {
        "source": SOURCE.name,
        "source_sha256": source_hash,
        "new_parameters": [f"p{i}" for i in range(44, 51)],
        "constraint_count": len(rows),
        "constraints": rows,
        "interpretation": "Exact term scan of the saved checkpoint. Affine dependence is only a necessary algebraic reduction; this audit does not prove the full branch empty.",
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "source_sha256": source_hash,
        "constraint_count": len(rows),
        "constraints": rows,
        "output": OUTPUT.name,
    }, indent=2))


if __name__ == "__main__":
    main()
