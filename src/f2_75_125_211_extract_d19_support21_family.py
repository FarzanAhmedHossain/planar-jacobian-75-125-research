"""Extract one exact support-21 family after the certified d19 substitutions."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path
import sys
import time

from flint import fmpq


ROOT = Path("research")
SOURCE = ROOT / "f2_75_125_sparse_ring_lift_211_roots19_u1zero_streamed_state.json"
REDUCTION = ROOT / "f2_75_125_211_simple_root_reduction_membership_20261025.json"
PROBE = ROOT / "f2_75_125_211_large_support_divisibility_probe_20261025.json"
OUTPUT = ROOT / "f2_75_125_211_d19_support21_family_input_20261025.json"
SOURCE_SHA256 = "0f293c437b1cc18f379e3a8d0bb2e98647616e062fa321e91cfdc45ebe67926d"
TARGET_SUPPORT = {
    "p24", "p25", "p26", "p27", "p28", "p36", "p37", "p38", "p39",
    "p44", "p45", "p46", "p47", "p48", "p5", "p53", "p54", "p55",
    "p56", "p57", "r",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if sha256(SOURCE) != SOURCE_SHA256:
        raise ValueError("d19 checkpoint hash mismatch")
    reduction = json.loads(REDUCTION.read_text(encoding="utf-8"))
    probe = json.loads(PROBE.read_text(encoding="utf-8"))
    if reduction.get("d19_source_sha256") != SOURCE_SHA256:
        raise ValueError("divisor substitution receipt is not bound to d19")
    if set(probe.get("reduced_support", [])) != TARGET_SUPPORT:
        raise ValueError("probe support differs from the hard-coded family")

    research = Path(__file__).resolve().parent
    sys.path.insert(0, str(research))
    from f2_75_125_flint_mpoly_backend import FlintRing
    from f2_75_125_211_streaming_checkpoint_scan import JsonStream

    reader = JsonStream(SOURCE)
    started = time.monotonic()
    saved = []
    rows_scanned = input_terms = zero_rows = 0
    support_histogram: Counter[int] = Counter()
    try:
        metadata, _ = reader.top_level_constraints()
        if metadata.get("branch") != "2+1+1" or metadata.get("quadratic"):
            raise ValueError("expected rational F2 2+1+1 checkpoint")
        names = metadata["generators"]
        if int(metadata["nparameters"]) != 174:
            raise ValueError("unexpected parameter count")
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
                print(f"rows={rows_scanned}; family_rows={len(saved)}; "
                      f"input_terms={input_terms}", flush=True)
            terms = {}
            for sparse, values in sparse_terms:
                input_terms += 1
                if len(values) != 1:
                    raise ValueError(f"non-rational coefficient in {tag}")
                exponents = [0] * len(names)
                for index, exponent in sparse:
                    exponents[int(index)] = int(exponent)
                key = tuple(exponents)
                terms[key] = terms.get(key, fmpq(0)) + fmpq(values[0])
            polynomial = ring.from_prevalidated_dict({
                exponents: coefficient
                for exponents, coefficient in terms.items()
                if coefficient
            }).compose(substitutions)
            if not len(polynomial):
                zero_rows += 1
                continue
            support_native = set()
            for exponents, _coefficient in polynomial.items():
                support_native.update(
                    names[index] for index, exponent in enumerate(exponents)
                    if exponent
                )
            support = {"r" if name == "p174" else name for name in support_native}
            support_histogram[len(support)] += 1
            if support != TARGET_SUPPORT:
                continue

            encoded = []
            for exponents, coefficient in polynomial.items():
                monomial = []
                for index, exponent in enumerate(exponents):
                    if not exponent:
                        continue
                    name = "r" if names[index] == "p174" else names[index]
                    monomial.append([name, int(exponent)])
                encoded.append([monomial, str(coefficient)])
            encoded.sort(key=lambda term: tuple((name, exponent)
                                                 for name, exponent in term[0]))
            saved.append({
                "tag": tag,
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

    if rows_scanned != 1421:
        raise ValueError(f"unexpected checkpoint row count {rows_scanned}")
    if not saved:
        raise ValueError("the target support-21 family was not found")
    if sha256(SOURCE) != SOURCE_SHA256:
        raise ValueError("d19 checkpoint changed during extraction")

    result = {
        "scope": (
            "Complete exact family of transformed d19 equations whose full "
            "parameter support equals the selected support-21 set. Equations "
            "with any other support are not included."
        ),
        "source_sha256": SOURCE_SHA256,
        "divisibility_receipt_sha256": sha256(REDUCTION),
        "probe_receipt_sha256": sha256(PROBE),
        "branch": metadata["branch"],
        "target_support": sorted(TARGET_SUPPORT),
        "checkpoint_rows_scanned": rows_scanned,
        "input_sparse_terms_scanned": input_terms,
        "zero_rows_after_specialization": zero_rows,
        "support_size_histogram": {
            str(k): v for k, v in sorted(support_histogram.items())
        },
        "equation_count": len(saved),
        "equations": saved,
        "elapsed_seconds": time.monotonic()-started,
        "scope_limit": (
            "A unit ideal from these necessary equations excludes the chart; a "
            "nonunit ideal does not exclude the full 2+1+1 branch."
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "rows_scanned": rows_scanned,
        "family_equations": len(saved),
        "family_term_count_min": min(len(row["terms"]) for row in saved),
        "family_term_count_max": max(len(row["terms"]) for row in saved),
        "elapsed_seconds": result["elapsed_seconds"],
        "output": str(OUTPUT),
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
