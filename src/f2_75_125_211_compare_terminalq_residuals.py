#!/usr/bin/env python3
"""Compare extracted sparse terminal-Q polynomials by exact term data."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "f2_75_125_flint_deps"))
import ijson  # noqa: E402


def main(path: Path) -> None:
    digests = {}
    with path.open("rb") as stream:
        for entry in ijson.items(stream, "constraints.item"):
            summary = entry["summary"]
            tag = summary["tag"]
            digest = hashlib.sha256()
            polynomial = entry["encoded_polynomial"]
            for term in polynomial:
                digest.update(json.dumps(term, separators=(",", ":")).encode("ascii"))
                digest.update(b"\n")
            digests[tag] = (len(polynomial), digest.hexdigest())
            print(f"{tag}\t{len(polynomial)}\t{digest.hexdigest()}")

    pairs = (
        ("terminal-Q-distance17:0:t^0", "terminal-Q-distance17:postsolve:t^0"),
        ("terminal-Q-distance17:0:t^1", "terminal-Q-distance17:postsolve:t^1"),
        ("terminal-Q-distance17:0:t^2", "terminal-Q-distance17:postsolve:t^2"),
    )
    for left, right in pairs:
        if left not in digests or right not in digests:
            raise ValueError(f"missing comparison pair {left!r}, {right!r}")
        print(f"EXACT_DUPLICATE\t{left}\t{right}\t{digests[left] == digests[right]}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: compare_terminalq_residuals.py EXTRACTED.json")
    main(Path(sys.argv[1]))
