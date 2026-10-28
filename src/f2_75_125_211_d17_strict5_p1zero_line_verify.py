#!/usr/bin/env python3
"""Replay a two-parameter solution family in the strict d17 five-var slice.

For p1=0, the exact F5B basis gives p16=-14*p17 and
p15=(452/5)*p17, while p7 is free. This verifier checks the resulting
parameterization directly in every original complete five-variable equation
selected from the immutable d17 checkpoint, without relying on parsing the
saved Groebner basis.
"""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

from f2_75_125_211_streaming_checkpoint_scan import JsonStream


EXPECTED_SOURCE_SHA256 = (
    "c205df23dc70fb51df2142c482c3145bc7fc25c6be95f28131818c9c06ca3812"
)
TARGET = {"p1", "p7", "p15", "p16", "p17"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    source = Path("research/f2_75_125_sparse_ring_lift_211_roots17_dd_20261024_d17_p79zero_snapshot.json")
    report_path = Path("research/f2_75_125_sparse_ring_lift_211_roots17_dd_20261024_d17_strict5_p1zero_line_verify.json")
    source_hash = sha256(source)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise ValueError(f"d17 checkpoint SHA-256 mismatch: {source_hash}")

    reader = JsonStream(source)
    try:
        metadata, _ = reader.top_level_constraints()
        if metadata.get("format") != "f2-sparse-ring-v1" or metadata.get("branch") != "2+1+1":
            raise ValueError("expected the rational d17 2+1+1 checkpoint")
        generators = metadata["generators"]
        p79_id = generators.index("p79")
        name_by_id = {i: name for i, name in enumerate(generators)}
        selected_rows = selected_terms = 0
        residual_rows = []
        for tag, sparse_terms in reader.array_values():
            polynomial = []
            outside = False
            for sparse, values in sparse_terms:
                if len(values) != 1:
                    raise ValueError("expected rational checkpoint coefficients")
                if any(index == p79_id and exponent for index, exponent in sparse):
                    continue
                term = {name_by_id[index]: int(exponent)
                        for index, exponent in sparse if exponent}
                if any(name not in TARGET for name in term):
                    outside = True
                    break
                polynomial.append((term, Fraction(values[0])))
            if outside or not polynomial:
                continue
            support = {name for term, _coefficient in polynomial for name in term}
            if support != TARGET:
                continue
            selected_rows += 1
            selected_terms += len(polynomial)

            # p1=0, p7=s, p17=r, p16=-14r, p15=(452/5)r.
            residual: dict[tuple[int, int], Fraction] = defaultdict(Fraction)
            for term, coefficient in polynomial:
                if term.get("p1", 0):
                    continue
                e7 = term.get("p7", 0)
                e15 = term.get("p15", 0)
                e16 = term.get("p16", 0)
                e17 = term.get("p17", 0)
                r_degree = e15 + e16 + e17
                coefficient *= Fraction(452, 5) ** e15
                coefficient *= Fraction(-14, 1) ** e16
                residual[(e7, r_degree)] += coefficient
            residual = {m: c for m, c in residual.items() if c}
            if residual:
                residual_rows.append({
                    "tag": tag,
                    "nonzero_monomials": len(residual),
                    "first_terms": [[list(m), str(c)] for m, c in list(residual.items())[:4]],
                })

        if selected_rows != 112:
            raise AssertionError(f"expected 112 complete five-variable equations, found {selected_rows}")
        if residual_rows:
            raise AssertionError(f"candidate line fails {len(residual_rows)} selected equations")
        result = {
            "source": str(source),
            "source_sha256": source_hash,
            "branch": metadata["branch"],
            "p79_specialization": "0",
            "selected_complete_equations": selected_rows,
            "source_terms_in_selected_equations": selected_terms,
            "checked_parameterization": {
                "p1": "0",
                "p7": "s (free)",
                "p15": "(452/5)*r",
                "p16": "-14*r",
                "p17": "r (free)",
            },
            "all_equations_vanish_identically_in_QQ[s,r]": True,
            "nonzero_residual_rows": residual_rows,
            "scope": (
                "This proves the strict five-variable necessary subsystem admits a "
                "two-parameter family on p1=0. It does not solve the remaining "
                "checkpoint equations or construct a Keller pair."
            ),
        }
        report_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(result, indent=2), flush=True)
        print(f"report={report_path}", flush=True)
    finally:
        reader.close()


if __name__ == "__main__":
    main()
