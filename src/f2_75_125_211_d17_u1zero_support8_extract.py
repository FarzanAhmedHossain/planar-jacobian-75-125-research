#!/usr/bin/env python3
"""Extract the smallest-support exact subsystem after the U1=0 reductions."""

from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

from f2_75_125_211_streaming_checkpoint_scan import JsonStream
from f2_75_125_211_d17_u1zero_strict5_reduction_scan import check_basis_claim


EXPECTED_SOURCE_SHA256 = (
    "c205df23dc70fb51df2142c482c3145bc7fc25c6be95f28131818c9c06ca3812"
)
KILLED_NAMES = {"p0", "p1", "p79"}
REPLACED_NAMES = {"p15", "p16", "p17"}
EXPECTED_MIN_SUPPORT_ROWS = 143


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    source = Path("research/f2_75_125_sparse_ring_lift_211_roots17_dd_20261024_d17_p79zero_snapshot.json")
    basis_path = Path("research/f2_75_125_sparse_ring_lift_211_roots17_dd_20261024_d17_strict5_p1zero_gb.json")
    output = Path("research/f2_75_125_sparse_ring_lift_211_roots17_dd_20261024_d17_u1zero_support8_inputs.json")
    source_hash = sha256(source)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise ValueError(f"d17 checkpoint SHA-256 mismatch: {source_hash}")
    basis_hash, basis_claim = check_basis_claim(basis_path)

    reader = JsonStream(source)
    try:
        metadata, _ = reader.top_level_constraints()
        if metadata.get("format") != "f2-sparse-ring-v1" or metadata.get("branch") != "2+1+1":
            raise ValueError("expected the rational d17 2+1+1 checkpoint")
        generators = metadata["generators"]
        nparameters = int(metadata["nparameters"])
        parameter_ids = {generators.index(f"p{i}"): f"p{i}" for i in range(nparameters)}
        killed_ids = {generators.index(name) for name in KILLED_NAMES}
        replaced_ids = {generators.index(name): name for name in REPLACED_NAMES}

        groups: dict[tuple[str, ...], list] = defaultdict(list)
        support_counts: Counter[tuple[str, ...]] = Counter()
        row_count = input_terms = auxiliary_terms = 0
        for tag, sparse_terms in reader.array_values():
            row_count += 1
            polynomial: dict[tuple[tuple[str, int], ...], Fraction] = defaultdict(Fraction)
            has_auxiliary = False
            for sparse, values in sparse_terms:
                input_terms += 1
                if len(values) != 1:
                    raise ValueError("expected rational checkpoint coefficients")
                if any(index in killed_ids and exponent for index, exponent in sparse):
                    continue
                powers: dict[str, int] = defaultdict(int)
                coefficient = Fraction(values[0])
                for generator_id, exponent in sparse:
                    exponent = int(exponent)
                    if not exponent or generator_id in killed_ids:
                        continue
                    if generator_id in replaced_ids:
                        name = replaced_ids[generator_id]
                        powers["r"] += exponent
                        if name == "p15":
                            coefficient *= Fraction(452, 5) ** exponent
                        elif name == "p16":
                            coefficient *= Fraction(-14, 1) ** exponent
                    elif generator_id in parameter_ids:
                        powers[parameter_ids[generator_id]] += exponent
                    else:
                        has_auxiliary = True
                        auxiliary_terms += 1
                monomial = tuple(sorted((name, exponent) for name, exponent in powers.items() if exponent))
                polynomial[monomial] += coefficient
            polynomial = {m: c for m, c in polynomial.items() if c}
            if not polynomial or has_auxiliary:
                continue
            support = tuple(sorted({name for monomial in polynomial for name, _ in monomial}))
            support_counts[support] += 1
            if len(support) != 8:
                continue
            groups[support].append({
                "tag": tag,
                "terms": [[[[name, exponent] for name, exponent in monomial], str(coefficient)]
                          for monomial, coefficient in sorted(polynomial.items())],
            })

        if auxiliary_terms:
            raise ValueError(f"unexpected auxiliary variable occurrences: {auxiliary_terms}")
        if sum(len(rows) for rows in groups.values()) != EXPECTED_MIN_SUPPORT_ROWS:
            raise AssertionError(
                f"expected {EXPECTED_MIN_SUPPORT_ROWS} size-8 equations, got "
                f"{sum(len(rows) for rows in groups.values())}"
            )
        payload = {
            "source": str(source),
            "source_sha256": source_hash,
            "strict5_basis_report": str(basis_path),
            "strict5_basis_report_sha256": basis_hash,
            "branch": metadata["branch"],
            "specialization": {"p79": "0", "p0": "0", "p1": "0"},
            "additional_set_theoretic_reduction": {
                "source": basis_claim,
                "p15": "(452/5)*r",
                "p16": "-14*r",
                "p17": "r",
            },
            "checkpoint_equations_scanned": row_count,
            "checkpoint_terms_scanned": input_terms,
            "smallest_nonzero_support_size": 8,
            "smallest_support_equation_count": EXPECTED_MIN_SUPPORT_ROWS,
            "smallest_support_groups": {
                ",".join(support): {"variables": list(support), "equations": rows}
                for support, rows in sorted(groups.items())
            },
            "all_support_counts_at_size_8": {
                ",".join(support): count
                for support, count in sorted(support_counts.items()) if len(support) == 8
            },
            "scope": (
                "These are complete equations supported on the smallest remaining "
                "parameter sets after the exact U1=0 restrictions. A unit ideal in "
                "one group excludes that U1=0 slice; a nonunit basis does not."
            ),
        }
        output.write_text(json.dumps(payload, separators=(",", ":")) + "\n", encoding="utf-8")
        payload["output_sha256"] = sha256(output)
        payload["output_bytes"] = output.stat().st_size
        print(json.dumps({key: payload[key] for key in (
            "source_sha256", "strict5_basis_report_sha256",
            "checkpoint_equations_scanned", "checkpoint_terms_scanned",
            "smallest_nonzero_support_size", "smallest_support_equation_count",
            "smallest_support_groups", "all_support_counts_at_size_8",
            "output_sha256", "output_bytes", "scope"
        )}, indent=2), flush=True)
        print(f"inputs={output}", flush=True)
    finally:
        reader.close()


if __name__ == "__main__":
    main()
