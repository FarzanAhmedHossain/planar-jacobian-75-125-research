"""Verify whether the reduced d19 support-12 subsystem is the known d18 one."""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path
import sys
import time


SOURCE = Path(
    "research/f2_75_125_sparse_ring_lift_211_roots19_u1zero_streamed_state.json"
)
SOURCE_HASH = "0f293c437b1cc18f379e3a8d0bb2e98647616e062fa321e91cfdc45ebe67926d"
OLD_INPUT = Path(
    "research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_input.json"
)
REDUCTION = Path(
    "research/f2_75_125_211_simple_root_reduction_membership_20261025.json"
)
REDUCED_SCAN = Path(
    "research/f2_75_125_211_support3_line_divisor_full_scan_20261025.json"
)
OUTPUT = Path(
    "research/f2_75_125_211_verify_d19_support12_reuse_20261025.json"
)
TARGET_SUPPORT = {
    "p24", "p25", "p26", "p27", "p28",
    "p44", "p45", "p46", "p47", "p48", "p5", "r",
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    old = json.loads(OLD_INPUT.read_text(encoding="utf-8"))
    if old["selected_variables"] != sorted(TARGET_SUPPORT):
        raise ValueError("unexpected d18 support-12 variable list")
    expected = {}
    for row in old["equations"]:
        equation = {}
        for monomial, coefficient in row["terms"]:
            powers = tuple(sorted((name, int(exponent))
                                  for name, exponent in monomial if exponent))
            value = Fraction(coefficient)
            equation[powers] = equation.get(powers, Fraction(0)) + value
        expected[row["tag"]] = {
            powers: coefficient for powers, coefficient in equation.items()
            if coefficient
        }

    reduction = json.loads(REDUCTION.read_text(encoding="utf-8"))
    reduced_scan = json.loads(REDUCED_SCAN.read_text(encoding="utf-8"))
    if reduction["d19_source_sha256"] != SOURCE_HASH:
        raise ValueError("d19 reduction receipt does not match the source checkpoint")
    if reduced_scan["source_sha256"] != SOURCE_HASH:
        raise ValueError("d19 transformed scan does not match the source checkpoint")
    if reduced_scan["most_repeated_supports"]["12"][0]["constraints"] != len(expected):
        raise ValueError("reduced support-12 count differs from d18 subsystem size")

    research = Path(__file__).resolve().parent
    sys.path.insert(0, str(research))
    from f2_75_125_flint_mpoly_backend import FlintRing
    from f2_75_125_211_streaming_checkpoint_scan import JsonStream
    from flint import fmpq

    reader = JsonStream(SOURCE)
    start = time.monotonic()
    try:
        metadata, _ = reader.top_level_constraints()
        if metadata.get("branch") != "2+1+1" or metadata.get("quadratic"):
            raise ValueError("unexpected checkpoint branch or field")
        names = metadata["generators"]
        if int(metadata["nparameters"]) != 174:
            raise ValueError("unexpected active parameter count")
        ring = FlintRing(names)
        generators = {name: ring.gens[index]
                      for index, name in enumerate(names)}
        r = generators["p174"]
        substitutions = {
            "p13": ring.ground_new(fmpq(452, 5))*r,
            "p14": -14*r,
            "p15": r,
            "p34": (
                45*generators["p36"] + fmpq(1654, 5)*generators["p37"]
                + fmpq(7491, 5)*generators["p38"]
                - fmpq(1449178, 25)*generators["p39"]
                + fmpq(612164, 125)*generators["p5"]*r
            ),
            "p35": (
                -10*generators["p36"] - 55*generators["p37"]
                - fmpq(1102, 5)*generators["p38"]
                + fmpq(40222, 5)*generators["p39"]
                - fmpq(251692, 375)*generators["p5"]*r
            ),
        }
        found = {}
        rows_scanned = terms_scanned = 0
        for tag, sparse_terms in reader.array_values():
            rows_scanned += 1
            if rows_scanned % 200 == 0:
                print(f"rows={rows_scanned} matched={len(found)}", flush=True)
            if tag not in expected:
                continue
            terms = {}
            for sparse, values in sparse_terms:
                terms_scanned += 1
                if len(values) != 1:
                    raise ValueError(f"non-rational coefficient in {tag}")
                powers = [0]*len(names)
                for generator, exponent in sparse:
                    powers[int(generator)] = int(exponent)
                key = tuple(powers)
                terms[key] = terms.get(key, fmpq(0)) + fmpq(values[0])
            polynomial = ring.from_prevalidated_dict({
                powers: coefficient for powers, coefficient in terms.items()
                if coefficient
            }).compose(substitutions)
            actual = {}
            support = set()
            for powers, coefficient in polynomial.items():
                monomial = []
                for index, exponent in enumerate(powers):
                    if exponent:
                        name = "r" if names[index] == "p174" else names[index]
                        support.add(name)
                        monomial.append((name, int(exponent)))
                key = tuple(sorted(monomial))
                actual[key] = actual.get(key, Fraction(0)) + Fraction(str(coefficient))
            actual = {powers: value for powers, value in actual.items() if value}
            found[tag] = {
                "equal_to_d18_equation": actual == expected[tag],
                "support_matches": support == TARGET_SUPPORT,
                "expected_terms": len(expected[tag]),
                "actual_terms": len(actual),
            }
    finally:
        reader.close()

    missing = sorted(set(expected)-set(found))
    mismatches = {tag: value for tag, value in found.items()
                  if not value["equal_to_d18_equation"] or not value["support_matches"]}
    result = {
        "scope": "exact tagged-equation comparison under the d19 divisor substitution",
        "d19_source_sha256": SOURCE_HASH,
        "d18_support12_input_sha256": sha256(OLD_INPUT),
        "d18_support12_equation_count": len(expected),
        "d19_rows_scanned": rows_scanned,
        "d19_sparse_terms_decoded_for_matching_rows": terms_scanned,
        "matched_d18_tags": len(found),
        "missing_d18_tags": missing,
        "mismatch_count": len(mismatches),
        "mismatches": mismatches,
        "d19_reduced_support12_total_count": reduced_scan[
            "most_repeated_supports"]["12"][0]["constraints"],
        "same_142_equations_and_full_support12_group": (
            not missing and not mismatches and
            reduced_scan["most_repeated_supports"]["12"][0]["constraints"]
            == len(expected)
        ),
        "d18_support12_basis_unit": False,
        "interpretation": (
            "If exact equality is confirmed, the transformed d19 support-12 group "
            "is the previously computed nonunit d18 subsystem; it gives no new "
            "exclusion by itself."
        ),
        "elapsed_seconds": time.monotonic()-start,
    }
    OUTPUT.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "d19_rows_scanned", "matched_d18_tags", "missing_d18_tags",
        "mismatch_count", "d19_reduced_support12_total_count",
        "same_142_equations_and_full_support12_group", "elapsed_seconds",
    )}, indent=2), flush=True)
    print(f"report={OUTPUT}")


if __name__ == "__main__":
    main()
