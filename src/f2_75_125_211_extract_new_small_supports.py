"""Save the genuinely new d19 equations with reduced support <=17."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path
import sys
import time


SOURCE = Path(
    "research/f2_75_125_sparse_ring_lift_211_roots19_u1zero_streamed_state.json"
)
SOURCE_HASH = "0f293c437b1cc18f379e3a8d0bb2e98647616e062fa321e91cfdc45ebe67926d"
REDUCTION = Path(
    "research/f2_75_125_211_simple_root_reduction_membership_20261025.json"
)
EXPECTED_HISTOGRAM = {12: 142, 13: 1, 17: 2}
OUTPUT = Path(
    "research/f2_75_125_211_new_small_support_equations_20261025.json"
)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    if sha256(SOURCE) != SOURCE_HASH:
        raise ValueError("d19 checkpoint hash mismatch")
    reduction = json.loads(REDUCTION.read_text(encoding="utf-8"))
    if reduction["d19_source_sha256"] != SOURCE_HASH:
        raise ValueError("divisibility receipt is not bound to d19 checkpoint")
    expected_solution = reduction["restriction_to_known_support3_line"]["solution"]

    research = Path(__file__).resolve().parent
    sys.path.insert(0, str(research))
    from f2_75_125_flint_mpoly_backend import FlintRing
    from f2_75_125_211_streaming_checkpoint_scan import JsonStream
    from flint import fmpq

    reader = JsonStream(SOURCE)
    started = time.monotonic()
    try:
        metadata, _ = reader.top_level_constraints()
        if metadata.get("branch") != "2+1+1" or metadata.get("quadratic"):
            raise ValueError("expected rational F2 2+1+1 checkpoint")
        names = metadata["generators"]
        nparameters = int(metadata["nparameters"])
        if nparameters != 174:
            raise ValueError(f"unexpected parameter count {nparameters}")
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
        support_histogram = Counter()
        saved_equations = []
        rows_scanned = input_terms = zero_rows = 0

        for tag, sparse_terms in reader.array_values():
            rows_scanned += 1
            if rows_scanned % 200 == 0:
                print(f"rows={rows_scanned} small_supports="
                      f"{dict(sorted(support_histogram.items()))}", flush=True)
            terms = {}
            for sparse, values in sparse_terms:
                input_terms += 1
                if len(values) != 1:
                    raise ValueError(f"non-rational coefficient in {tag}")
                powers = [0]*len(names)
                for index, exponent in sparse:
                    powers[int(index)] = int(exponent)
                key = tuple(powers)
                terms[key] = terms.get(key, fmpq(0))+fmpq(values[0])
            polynomial = ring.from_prevalidated_dict({
                powers: coefficient for powers, coefficient in terms.items()
                if coefficient
            }).compose(substitutions)
            if not len(polynomial):
                zero_rows += 1
                continue

            support = set()
            encoded = []
            too_large = False
            for powers, coefficient in polynomial.items():
                monomial = []
                for index, exponent in enumerate(powers):
                    if not exponent:
                        continue
                    if 1 <= index <= nparameters:
                        name = names[index]
                    elif names[index] == "p174":
                        name = "r"
                    else:
                        raise ValueError(f"unexpected variable {names[index]} in {tag}")
                    support.add(name)
                    monomial.append([name, int(exponent)])
                if len(support) > 17:
                    too_large = True
                    break
                encoded.append([monomial, str(coefficient)])
            if too_large:
                continue
            size = len(support)
            support_histogram[size] += 1
            if size in (13, 17):
                encoded.sort(key=lambda term: tuple((v, e) for v, e in term[0]))
                saved_equations.append({
                    "tag": tag,
                    "support_size": size,
                    "support": sorted(support),
                    "term_count": len(encoded),
                    "terms": encoded,
                })
        reader.close()

    finally:
        try:
            reader.close()
        except Exception:
            pass

    actual_histogram = {int(k): int(v) for k, v in support_histogram.items()}
    if actual_histogram != EXPECTED_HISTOGRAM:
        raise ValueError(f"small-support histogram mismatch: {actual_histogram}")
    if len(saved_equations) != 3:
        raise ValueError(f"expected 3 new equations, found {len(saved_equations)}")

    source_hash_now = sha256(SOURCE)
    if source_hash_now != SOURCE_HASH:
        raise ValueError(f"source checkpoint changed: {source_hash_now}")
    result = {
        "scope": "exact newly reduced support-13/support-17 equations; support-12 reused from verified d18 subsystem",
        "source_sha256": SOURCE_HASH,
        "d19_divisibility_receipt_sha256": sha256(REDUCTION),
        "source_sha256_rechecked": source_hash_now,
        "substitution": {
            "p13": "452*r/5", "p14": "-14*r", "p15": "r",
            "p34": expected_solution["p34"],
            "p35": expected_solution["p35"],
        },
        "rows_scanned": rows_scanned,
        "input_sparse_terms": input_terms,
        "zero_rows_after_substitution": zero_rows,
        "reduced_support_histogram_through_size17": {
            str(k): v for k, v in sorted(actual_histogram.items())
        },
        "equations": saved_equations,
        "whole_2+1+1_branch_excluded_percent": 0.0,
        "elapsed_seconds": time.monotonic()-started,
    }
    OUTPUT.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({
        "rows_scanned": rows_scanned,
        "small_support_histogram": result["reduced_support_histogram_through_size17"],
        "new_equations": [{"tag": row["tag"], "support_size": row["support_size"],
                           "support": row["support"], "terms": row["term_count"]}
                          for row in saved_equations],
        "elapsed_seconds": result["elapsed_seconds"],
        "report": str(OUTPUT),
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
