"""Extract the complete original support-32 family after d19 substitutions."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import time

from flint import fmpq


ROOT = Path("research")
SOURCE = ROOT / "f2_75_125_sparse_ring_lift_211_roots19_u1zero_streamed_state.json"
SCAN = ROOT / "f2_75_125_sparse_ring_lift_211_u1zero_streamed_scan_all_supports.json"
REDUCTION = ROOT / "f2_75_125_211_simple_root_reduction_membership_20261025.json"
OUTPUT = ROOT / "f2_75_125_211_d19_support32_family_input_20261025.json"
SOURCE_SHA256 = "0f293c437b1cc18f379e3a8d0bb2e98647616e062fa321e91cfdc45ebe67926d"
TARGET_SIZE = 32
EXPECTED_ROWS = 117


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main(target_size: int = TARGET_SIZE,
         output: Path = OUTPUT) -> None:
    if sha256(SOURCE) != SOURCE_SHA256:
        raise ValueError("d19 checkpoint hash mismatch")
    scan = json.loads(SCAN.read_text(encoding="utf-8"))
    reduction = json.loads(REDUCTION.read_text(encoding="utf-8"))
    if scan.get("source_sha256") != SOURCE_SHA256:
        raise ValueError("support scan is not bound to d19")
    if reduction.get("d19_source_sha256") != SOURCE_SHA256:
        raise ValueError("substitution receipt is not bound to d19")

    research = Path(__file__).resolve().parent
    sys.path.insert(0, str(research))
    from f2_75_125_flint_mpoly_backend import FlintRing
    from f2_75_125_211_streaming_checkpoint_scan import JsonStream

    reader = JsonStream(SOURCE)
    started = time.monotonic()
    saved = []
    rows_scanned = input_terms = candidate_rows = zero_rows = 0
    reduced_support_histogram: dict[int, int] = {}
    try:
        metadata, _ = reader.top_level_constraints()
        if metadata.get("branch") != "2+1+1" or metadata.get("quadratic"):
            raise ValueError("expected rational F2 2+1+1 checkpoint")
        names = metadata["generators"]
        if int(metadata["nparameters"]) != 174:
            raise ValueError("unexpected parameter count")

        group = scan["most_repeated_supports"][str(target_size)]
        if len(group) != 1:
            raise ValueError(f"support-{target_size} scan metadata mismatch")
        expected_rows = int(group[0]["constraints"])
        # The support report stores zero-based p-indices (p_i), while the
        # checkpoint generator list has a leading auxiliary `t` at position 0.
        original_support = {names[int(index)+1] for index in group[0]["variables"]}
        if len(original_support) != target_size:
            raise ValueError(f"support-{target_size} variable list has wrong size")

        ring = FlintRing(names)
        g = {name: ring.gens[index] for index, name in enumerate(names)}
        r = g["p174"]
        substitutions = {
            "p13": ring.ground_new(fmpq(452, 5))*r,
            "p14": -14*r,
            "p15": r,
            "p34": (
                45*g["p36"] + fmpq(1654, 5)*g["p37"]
                + fmpq(7491, 5)*g["p38"]
                - fmpq(1449178, 25)*g["p39"]
                + fmpq(612164, 125)*g["p5"]*r
            ),
            "p35": (
                -10*g["p36"] - 55*g["p37"] - fmpq(1102, 5)*g["p38"]
                + fmpq(40222, 5)*g["p39"]
                - fmpq(251692, 375)*g["p5"]*r
            ),
        }

        for tag, sparse_terms in reader.array_values():
            rows_scanned += 1
            if rows_scanned % 200 == 0:
                print(json.dumps({
                    "rows_scanned": rows_scanned,
                    "candidate_rows": candidate_rows,
                    "saved_rows": len(saved),
                    "input_terms": input_terms,
                }), flush=True)
            row_support: set[str] = set()
            terms = {}
            for sparse, values in sparse_terms:
                input_terms += 1
                if len(values) != 1:
                    raise ValueError(f"non-rational coefficient in {tag}")
                exponents = [0] * len(names)
                for index, exponent in sparse:
                    index, exponent = int(index), int(exponent)
                    exponents[index] = exponent
                    if exponent:
                        row_support.add(names[index])
                key = tuple(exponents)
                terms[key] = terms.get(key, fmpq(0)) + fmpq(values[0])
            if row_support != original_support:
                continue
            candidate_rows += 1
            polynomial = ring.from_prevalidated_dict({
                exponents: coefficient
                for exponents, coefficient in terms.items()
                if coefficient
            }).compose(substitutions)
            if not len(polynomial):
                zero_rows += 1
                continue
            reduced_support = set()
            for exponents, _coefficient in polynomial.items():
                reduced_support.update(
                    "r" if names[index] == "p174" else names[index]
                    for index, exponent in enumerate(exponents) if exponent
                )
            reduced_support_histogram[len(reduced_support)] = (
                reduced_support_histogram.get(len(reduced_support), 0) + 1
            )
            encoded = []
            for exponents, coefficient in polynomial.items():
                monomial = []
                for index, exponent in enumerate(exponents):
                    if exponent:
                        name = "r" if names[index] == "p174" else names[index]
                        monomial.append([name, int(exponent)])
                encoded.append([monomial, str(coefficient)])
            encoded.sort(key=lambda term: tuple((name, exponent)
                                                 for name, exponent in term[0]))
            saved.append({
                "tag": tag,
                "original_support": sorted(original_support),
                "reduced_support": sorted(reduced_support),
                "term_count": len(encoded),
                "terms": encoded,
            })
        reader.close()
    finally:
        try:
            reader.close()
        except Exception:
            pass

    if rows_scanned != 1421 or candidate_rows != expected_rows:
        raise ValueError(
            f"unexpected checkpoint/support count: rows={rows_scanned}, "
            f"support{target_size}={candidate_rows}"
        )
    if sha256(SOURCE) != SOURCE_SHA256:
        raise ValueError("d19 checkpoint changed during extraction")
    result = {
        "scope": (
            f"Complete exact original support-{target_size} family after the verified d19 "
            "substitutions. Source support is tested before substitution; "
            "reduced supports may be smaller."
        ),
        "source_sha256": SOURCE_SHA256,
        "support_scan_sha256": sha256(SCAN),
        "substitution_receipt_sha256": sha256(REDUCTION),
        "branch": metadata["branch"],
        "original_support": sorted(original_support),
        "checkpoint_rows_scanned": rows_scanned,
        "input_sparse_terms_scanned": input_terms,
        "candidate_rows": candidate_rows,
        "zero_rows_after_substitution": zero_rows,
        "reduced_support_size_histogram": {
            str(k): v for k, v in sorted(reduced_support_histogram.items())
        },
        "equation_count": len(saved),
        "equations": saved,
        "elapsed_seconds": time.monotonic()-started,
        "scope_limit": (
            "This is a necessary equation family only. A nonunit ideal does "
            "not exclude the full F2 branch."
        ),
    }
    output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({
        "rows_scanned": rows_scanned,
        "candidate_rows": candidate_rows,
        "saved_rows": len(saved),
        "zero_rows": zero_rows,
        "reduced_support_size_histogram": result["reduced_support_size_histogram"],
        "elapsed_seconds": result["elapsed_seconds"],
        "output": str(output),
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
