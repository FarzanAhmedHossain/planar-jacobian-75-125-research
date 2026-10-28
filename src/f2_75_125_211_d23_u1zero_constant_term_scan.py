#!/usr/bin/env python3
"""Check exact constant coefficients in the complete d23 checkpoint."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from f2_75_125_211_streaming_checkpoint_scan import JsonStream


EXPECTED_SOURCE_SHA256 = "270f0eb737463c75d3bdfbcdfb3a303a2e7fa0af2da13f63a9e2d01ec5e04b59"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    source_hash = sha256(args.source)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise ValueError(f"unexpected d23 checkpoint: {source_hash}")
    reader = JsonStream(args.source)
    rows = sparse_terms = 0
    nonzero_constant_rows = []
    try:
        metadata, _ = reader.top_level_constraints()
        if (metadata.get("format") != "f2-sparse-ring-v1"
                or metadata.get("branch") != "2+1+1"
                or metadata.get("quadratic")
                or int(metadata.get("nparameters", -1)) != 211):
            raise ValueError("unexpected checkpoint scope")
        for tag, terms in reader.array_values():
            rows += 1
            constant = Fraction(0)
            for sparse, values in terms:
                sparse_terms += 1
                if len(values) != 1:
                    raise ValueError(f"non-rational coefficient in {tag}")
                if all(int(exponent) == 0 for _generator, exponent in sparse):
                    constant += Fraction(values[0])
            if constant:
                nonzero_constant_rows.append({"tag": tag, "constant": str(constant)})
            if rows % 200 == 0:
                print(f"rows={rows} sparse_terms={sparse_terms} "
                      f"nonzero_constant_rows={len(nonzero_constant_rows)}", flush=True)
    finally:
        reader.close()

    result = {
        "source": str(args.source),
        "source_sha256": source_hash,
        "branch": metadata["branch"],
        "parameters": metadata["nparameters"],
        "rows_scanned": rows,
        "sparse_terms_scanned": sparse_terms,
        "nonzero_constant_rows": nonzero_constant_rows,
        "all_checkpoint_equations_vanish_at_parameter_origin": not nonzero_constant_rows,
        "substitution_p34_p35_has_no_constant_terms": True,
        "scope": (
            "The exact L=M=0 formulas for p34 and p35 are polynomials with zero "
            "constant term. Therefore this scan decides whether every saved necessary "
            "equation still vanishes at the origin after that substitution. An origin "
            "witness proves the equality ideal nonunit, but does not address open-chart "
            "nonvanishing conditions or exclude the branch."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    print(f"report={args.output}", flush=True)


if __name__ == "__main__":
    main()
