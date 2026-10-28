"""Extract the two original support-21 equations reducing to support 17."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


SOURCE = Path(
    "research/f2_75_125_sparse_ring_lift_211_roots19_u1zero_streamed_state.json"
)
SOURCE_HASH = "0f293c437b1cc18f379e3a8d0bb2e98647616e062fa321e91cfdc45ebe67926d"
REDUCTION = Path(
    "research/f2_75_125_211_simple_root_reduction_membership_20261025.json"
)
OUTPUT = Path(
    "research/f2_75_125_211_reduced_support17_equations_20261025.json"
)
REDUCED_SUPPORT = {
    "p5", "r", "p24", "p25", "p26", "p27", "p28", "p29",
    "p36", "p37", "p38", "p39", "p44", "p45", "p46", "p47", "p48",
}
CANDIDATE_SUPPORT = REDUCED_SUPPORT | {"p13", "p14", "p15", "p34", "p35"}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    reduction = json.loads(REDUCTION.read_text(encoding="utf-8"))
    if reduction["d19_source_sha256"] != SOURCE_HASH:
        raise ValueError("divisibility reduction receipt belongs to another source")
    solution = reduction["restriction_to_known_support3_line"]["solution"]
    if set(solution) != {"p34", "p35"}:
        raise ValueError("unexpected p34,p35 solution variables")

    research = Path(__file__).resolve().parent
    sys.path.insert(0, str(research))
    from f2_75_125_flint_mpoly_backend import FlintRing
    from f2_75_125_211_streaming_checkpoint_scan import JsonStream
    from flint import fmpq

    reader = JsonStream(SOURCE)
    try:
        metadata, _ = reader.top_level_constraints()
        if metadata.get("branch") != "2+1+1" or metadata.get("quadratic"):
            raise ValueError("expected rational F2 2+1+1 checkpoint")
        names = metadata["generators"]
        if int(metadata["nparameters"]) != 174:
            raise ValueError("unexpected active parameter count")
        ids = {name: index for index, name in enumerate(names)}
        ring = FlintRing(names)
        g = {name: ring.gens[index] for index, name in enumerate(names)}
        r = g["p174"]
        substitutions = {
            "p13": ring.ground_new(fmpq(452, 5))*r,
            "p14": -14*r,
            "p15": r,
            "p34": (
                45*g["p36"] + fmpq(1654, 5)*g["p37"]
                + fmpq(7491, 5)*g["p38"] - fmpq(1449178, 25)*g["p39"]
                + fmpq(612164, 125)*g["p5"]*r
            ),
            "p35": (
                -10*g["p36"] - 55*g["p37"] - fmpq(1102, 5)*g["p38"]
                + fmpq(40222, 5)*g["p39"] - fmpq(251692, 375)*g["p5"]*r
            ),
        }

        rows_scanned = candidate_rows_composed = 0
        selected = []
        output_support_counts = {}
        for tag, sparse_terms in reader.array_values():
            rows_scanned += 1
            if rows_scanned % 200 == 0:
                print(f"rows={rows_scanned} candidates={candidate_rows_composed} "
                      f"selected={len(selected)}", flush=True)
            original_support = set()
            for sparse, _values in sparse_terms:
                original_support.update(
                    names[int(index)] for index, exponent in sparse if exponent
                )
            if not original_support.issubset(CANDIDATE_SUPPORT):
                continue
            candidate_rows_composed += 1
            terms = {}
            for sparse, values in sparse_terms:
                if len(values) != 1:
                    raise ValueError(f"non-rational coefficient in {tag}")
                powers = [0]*len(names)
                for index, exponent in sparse:
                    powers[int(index)] = int(exponent)
                key = tuple(powers)
                terms[key] = terms.get(key, fmpq(0)) + fmpq(values[0])
            reduced = ring.from_prevalidated_dict({
                powers: coefficient for powers, coefficient in terms.items()
                if coefficient
            }).compose(substitutions)
            encoded = []
            reduced_support = set()
            for powers, coefficient in reduced.items():
                monomial = []
                for index, exponent in enumerate(powers):
                    if exponent:
                        name = "r" if names[index] == "p174" else names[index]
                        reduced_support.add(name)
                        monomial.append([name, int(exponent)])
                encoded.append([monomial, str(coefficient)])
            support_key = tuple(sorted(reduced_support))
            output_support_counts[support_key] = output_support_counts.get(support_key, 0)+1
            if reduced_support == REDUCED_SUPPORT:
                selected.append({
                    "tag": tag,
                    "original_support": sorted(original_support),
                    "reduced_support": sorted(reduced_support),
                    "term_count": len(encoded),
                    "terms": encoded,
                })
        if len(selected) != 2:
            raise ValueError(f"expected 2 source rows, found {len(selected)}")
    finally:
        reader.close()

    source_hash_now = sha256(SOURCE)
    if source_hash_now != SOURCE_HASH:
        raise ValueError(f"source checkpoint hash mismatch: {source_hash_now}")
    result = {
        "scope": "exact two-row reduced support-17 subsystem; no branch exclusion",
        "source_sha256": SOURCE_HASH,
        "source_sha256_rechecked": source_hash_now,
        "reduction_receipt_sha256": sha256(REDUCTION),
        "rows_scanned": rows_scanned,
        "candidate_rows_composed": candidate_rows_composed,
        "candidate_original_support_domain": sorted(CANDIDATE_SUPPORT),
        "output_support_histogram_among_candidates": {
            ",".join(key): value
            for key, value in sorted(output_support_counts.items())
        },
        "reduced_support": sorted(REDUCED_SUPPORT),
        "equation_count": len(selected),
        "equations": selected,
        "substitution": {
            "p13": "452*r/5", "p14": "-14*r", "p15": "r",
            "p34": solution["p34"], "p35": solution["p35"],
        },
        "whole_2+1+1_branch_excluded_percent": 0.0,
    }
    OUTPUT.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({
        "rows_scanned": rows_scanned,
        "equations": [{"tag": row["tag"], "term_count": row["term_count"]}
                      for row in selected],
        "report": str(OUTPUT),
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
