#!/usr/bin/env python3
"""Export the four checked support-12 equations on the weighted r=1 chart."""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
from pathlib import Path
import sys

import sympy as sp
from sympy.polys.rings import ring

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import f2_75_125_211_d19_export_r1_singular as full_export

PRIME = 251
STEM = "f2_75_125_211_d19_r1_support12"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    sink = io.StringIO()
    with contextlib.redirect_stdout(sink):
        polys = full_export.capture_exact_system()
    source_names = tuple(json.loads((HERE / "f2_75_125_211_d19_support38_analysis_20261025.json").read_text())[
        "selected_variables"
    ])
    # The last four generators are the checked support-12 subsystem.
    source_polys = polys[-4:]
    active = tuple(i for i in range(len(source_names))
                   if any(any(m[i] for m, c in p.terms() if c) for p in source_polys))
    names = tuple(source_names[i] for i in active if source_names[i] != "r")
    active = tuple(i for i in active if source_names[i] != "r")
    ring_names = (*names, "r")
    R, *xs = ring(ring_names, sp.QQ)
    reduced = []
    for p in source_polys:
        # Set r=1 by dropping its exponent while keeping every term.
        q = R.from_dict({tuple(powers[i] for i in active) + (0,): sp.Rational(coeff)
                         for powers, coeff in p.terms()})
        reduced.append(q)

    outputs = {}
    for modular in (False, True):
        source, counts, constants = full_export.build_singular_source(reduced, ring_names, modular)
        suffix = "mod251" if modular else "QQ"
        path = HERE / f"{STEM}_{suffix}_20261028.sing"
        path.write_text(source, encoding="utf-8")
        outputs[suffix] = {"path": path.name, "sha256": sha256(path),
                           "term_counts": counts, "constants": constants}

    manifest_path = HERE / f"{STEM}_manifest_20261028.json"
    manifest = {
        "scope": "four checked support-12 d19 necessary equations, weighted chart r=1",
        "branch": "2+1+1, U1=0; r=p15 nonzero chart only",
        "variables": names,
        "weights": [int(json.loads((HERE / "f2_75_125_211_d19_support38_analysis_20261025.json").read_text())[
            "primitive_positive_integer_weights"][source_names.index(n)]) for n in names],
        "generator_count": len(reduced),
        "outputs": outputs,
        "full_manifest_sha256": sha256(HERE / "f2_75_125_211_d19_r1_normalized_singular_manifest_20261028.json"),
        "warning": "Subset result only; exact QQ unit would exclude only this normalized necessary subsystem's chart if all scope links are verified.",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"manifest": manifest_path.name, "variables": names,
                      "outputs": outputs}, indent=2), flush=True)


if __name__ == "__main__":
    main()
