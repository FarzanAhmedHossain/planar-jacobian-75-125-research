#!/usr/bin/env python3
"""Export a six-generator exact subset on the weighted r=1 chart."""

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

INDICES = (0, 1, 14, 15, 16, 17)
STEM = "f2_75_125_211_d19_r1_support12_support21"
LABELS = (
    "support21:corrected-power-distance13:t^10",
    "support21:corrected-power-distance13:t^11",
    "support21:19:1:t^38",
    "support21:19:1:t^39",
    "support32:corrected-power-distance13:t^10",
    "support32:corrected-power-distance13:t^11",
    "support38:corrected-power-distance14:t^10",
    "support38:corrected-power-distance14:t^11",
    "support38:corrected-power-distance14:t^12",
    "support38:corrected-power-distance14:t^13",
    "support38:corrected-power-distance14:t^14",
    "support38:14:0:t^5",
    "support38:corrected-power-distance18:t^9",
    "support38:18:0:t^4",
    "support12-G0", "support12-G1", "support12-G2", "support12-G3",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    with contextlib.redirect_stdout(io.StringIO()):
        full = full_export.capture_exact_system()
    analysis = json.loads((HERE / "f2_75_125_211_d19_support38_analysis_20261025.json").read_text())
    source_names = tuple(analysis["selected_variables"])
    source_weights = tuple(map(int, analysis["primitive_positive_integer_weights"]))
    selected = tuple(full[i] for i in INDICES)
    active = tuple(i for i, n in enumerate(source_names)
                   if n != "r" and any(any(m[i] for m, c in p.terms() if c) for p in selected))
    names = tuple(source_names[i] for i in active)
    ring_names = (*names, "r")
    R, *_ = ring(ring_names, sp.QQ)
    r_index = source_names.index("r")
    reduced = []
    for p in selected:
        q = R.from_dict({tuple(m[i] for i in active) + (0,): sp.Rational(c)
                         for m, c in p.terms()})
        reduced.append(q)

    files = {}
    for modular in (False, True):
        source, counts, constants = full_export.build_singular_source(reduced, ring_names, modular)
        # Use the natural positive weighted monomial order on the normalized variables.
        weights = [source_weights[source_names.index(n)] for n in names]
        first, rest = source.split("\n", 1)
        assert first.endswith(",dp;")
        first = first[:-3] + f"wp({','.join(map(str, weights))});"
        source = first + "\n" + rest
        suffix = "mod251" if modular else "QQ"
        path = HERE / f"{STEM}_{suffix}_wp_20261028.sing"
        path.write_text(source, encoding="utf-8")
        files[suffix] = {"file": path.name, "sha256": sha256(path),
                         "term_counts": counts, "constants": constants}

    manifest = {
        "scope": "Exact six-generator d19 necessary subset on weighted chart r=1",
        "branch": "2+1+1, U1=0, r=p15 nonzero chart only",
        "source_generator_indices_zero_based": list(INDICES),
        "source_labels": [LABELS[i] for i in INDICES],
        "normalized_variables": names,
        "normalized_weights": [source_weights[source_names.index(n)] for n in names],
        "files": files,
        "base_manifest_sha256": sha256(HERE / "f2_75_125_211_d19_r1_normalized_singular_manifest_20261028.json"),
        "warning": "These are necessary equations only; even a unit result excludes only this chart of this subsystem until every source and branch implication is verified.",
    }
    out = HERE / f"{STEM}_manifest_20261028.json"
    out.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"manifest": out.name, "variables": names, "files": files}, indent=2), flush=True)


if __name__ == "__main__":
    main()
