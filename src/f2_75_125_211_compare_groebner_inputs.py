#!/usr/bin/env python3
"""Diagnostic comparison of five-variable projections after p79=0.

These projections discard monomials containing other parameters and therefore
are not necessary equations unless those parameters are separately set to zero.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

from f2_75_125_211_streaming_checkpoint_scan import JsonStream


TARGET = {"p1", "p7", "p15", "p16", "p17"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def row_digests(terms: list, generators: list[str], p79_index: int):
    combined: dict[tuple[tuple[str, int], ...], Fraction] = defaultdict(Fraction)
    for sparse, values in terms:
        if len(values) != 1:
            raise ValueError("expected rational checkpoint coefficients")
        if any(index == p79_index and exponent for index, exponent in sparse):
            continue  # exact specialization p79=0
        monomial: list[tuple[str, int]] = []
        outside = False
        for index, exponent in sparse:
            if not exponent:
                continue
            name = generators[index]
            if name in TARGET:
                monomial.append((name, int(exponent)))
            else:
                outside = True
                break
        if outside:
            # Match the Gröbner driver's projection onto the five target
            # variables: discard this monomial, not the whole equation.
            continue
        combined[tuple(sorted(monomial))] += Fraction(values[0])

    combined = {monomial: coefficient for monomial, coefficient in combined.items()
                if coefficient}
    support = {name for monomial in combined for name, _power in monomial}
    if support != TARGET or not combined:
        return None
    decoded = sorted(combined.items())

    def digest(rows) -> str:
        payload = [[list(monomial), str(coefficient)]
                   for monomial, coefficient in rows]
        encoded = json.dumps(payload, separators=(",", ":"), ensure_ascii=False)
        return hashlib.sha256(encoded.encode("utf-8")).hexdigest()

    exact = digest(decoded)
    scale = decoded[0][1]
    normalized = [(monomial, coefficient / scale)
                  for monomial, coefficient in decoded]
    return exact, digest(normalized)


def scan(path: Path, expected_sha256: str) -> dict:
    actual_sha256 = sha256_file(path)
    if actual_sha256 != expected_sha256.lower():
        raise ValueError(f"SHA-256 mismatch for {path}: {actual_sha256}")
    reader = JsonStream(path)
    try:
        metadata, _ = reader.top_level_constraints()
        if (metadata.get("format") != "f2-sparse-ring-v1"
                or metadata.get("branch") != "2+1+1"
                or metadata.get("quadratic")):
            raise ValueError("expected a rational 2+1+1 checkpoint")
        generators = metadata["generators"]
        p79_index = generators.index("p79")
        if not TARGET.issubset(generators):
            raise ValueError("checkpoint lacks a target parameter")

        exact_rows: Counter[str] = Counter()
        scalar_rows: Counter[str] = Counter()
        exact_tags: dict[str, list[str]] = defaultdict(list)
        scalar_tags: dict[str, list[str]] = defaultdict(list)
        selected = constraint_count = 0
        for row in reader.array_values():
            if not isinstance(row, list) or len(row) != 2:
                raise ValueError("malformed constraint row")
            tag, terms = row
            constraint_count += 1
            digests = row_digests(terms, generators, p79_index)
            if digests is None:
                continue
            exact_digest, scalar_digest = digests
            selected += 1
            exact_rows[exact_digest] += 1
            scalar_rows[scalar_digest] += 1
            if len(exact_tags[exact_digest]) < 4:
                exact_tags[exact_digest].append(tag)
            if len(scalar_tags[scalar_digest]) < 4:
                scalar_tags[scalar_digest].append(tag)
        return {
            "source": str(path),
            "source_sha256": actual_sha256,
            "nparameters": int(metadata["nparameters"]),
            "constraints": constraint_count,
            "selected_rows_after_p79_zero": selected,
            "exact_row_multiset": exact_rows,
            "scalar_normalized_row_multiset": scalar_rows,
            "exact_tags": exact_tags,
            "scalar_tags": scalar_tags,
        }
    finally:
        reader.close()


def unmatched_tags(left: Counter, left_tags: dict, right: Counter, limit: int = 8) -> list[str]:
    missing = left - right
    result = []
    for digest, count in missing.items():
        result.extend(left_tags[digest][:count])
        if len(result) >= limit:
            break
    return result[:limit]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("left", type=Path)
    parser.add_argument("left_sha256")
    parser.add_argument("right", type=Path)
    parser.add_argument("right_sha256")
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    left = scan(args.left, args.left_sha256)
    right = scan(args.right, args.right_sha256)

    def common(a: Counter, b: Counter) -> int:
        return sum((a & b).values())

    report = {
        "target_support": sorted(TARGET),
        "specialization": {"p79": "0"},
        "left": {key: left[key] for key in (
            "source", "source_sha256", "nparameters", "constraints",
            "selected_rows_after_p79_zero")},
        "right": {key: right[key] for key in (
            "source", "source_sha256", "nparameters", "constraints",
            "selected_rows_after_p79_zero")},
        "exact_polynomial_multisets_equal": (
            left["exact_row_multiset"] == right["exact_row_multiset"]
        ),
        "scalar_normalized_equation_multisets_equal": (
            left["scalar_normalized_row_multiset"]
            == right["scalar_normalized_row_multiset"]
        ),
        "exact_common_rows": common(left["exact_row_multiset"], right["exact_row_multiset"]),
        "scalar_normalized_common_rows": common(
            left["scalar_normalized_row_multiset"],
            right["scalar_normalized_row_multiset"]),
        "left_only_exact_tags": unmatched_tags(
            left["exact_row_multiset"], left["exact_tags"], right["exact_row_multiset"]),
        "right_only_exact_tags": unmatched_tags(
            right["exact_row_multiset"], right["exact_tags"], left["exact_row_multiset"]),
        "left_only_scalar_equation_tags": unmatched_tags(
            left["scalar_normalized_row_multiset"], left["scalar_tags"],
            right["scalar_normalized_row_multiset"]),
        "right_only_scalar_equation_tags": unmatched_tags(
            right["scalar_normalized_row_multiset"], right["scalar_tags"],
            left["scalar_normalized_row_multiset"]),
        "scope": (
            "Exact comparison of projected five-variable polynomials after p79=0. "
            "Discarding monomials in outside parameters does not preserve the zero "
            "locus, so this is diagnostic only and is not a valid necessary subsystem "
            "without further zero-specializations."
        ),
    }
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"report={args.report}")


if __name__ == "__main__":
    main()
