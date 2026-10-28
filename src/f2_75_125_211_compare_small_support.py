#!/usr/bin/env python3
"""Compare a selected sparse-constraint support across two checkpoints.

Rows are decoded one at a time. The report compares both literal equation
multisets and equation multisets after normalizing each polynomial by a
nonzero rational scalar, which preserves the zero locus of each equation.
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


def polynomial_digests(terms: list, generators: list[str]) -> tuple[str, str] | None:
    decoded = []
    for sparse, values in terms:
        if len(values) != 1:
            raise ValueError("expected rational checkpoint coefficients")
        monomial = tuple(sorted(
            (generators[index], int(exponent))
            for index, exponent in sparse if exponent
        ))
        coefficient = Fraction(values[0])
        if coefficient:
            decoded.append((monomial, coefficient))
    if not decoded:
        return None
    decoded.sort(key=lambda term: term[0])

    def digest_for(rows) -> str:
        payload = [[list(monomial), str(coefficient)]
                   for monomial, coefficient in rows]
        encoded = json.dumps(payload, separators=(",", ":"), ensure_ascii=False)
        return hashlib.sha256(encoded.encode("utf-8")).hexdigest()

    exact = digest_for(decoded)
    scale = decoded[0][1]
    normalized = [(monomial, coefficient / scale)
                  for monomial, coefficient in decoded]
    return exact, digest_for(normalized)


def scan(path: Path, expected_sha256: str) -> dict:
    actual_sha256 = sha256_file(path)
    if actual_sha256 != expected_sha256.lower():
        raise ValueError(f"SHA-256 mismatch for {path}: {actual_sha256}")

    reader = JsonStream(path)
    try:
        metadata, _ = reader.top_level_constraints()
        if metadata.get("format") != "f2-sparse-ring-v1":
            raise ValueError("unexpected checkpoint format")
        if metadata.get("branch") != "2+1+1" or metadata.get("quadratic"):
            raise ValueError("expected rational 2+1+1 checkpoint")
        generators = metadata.get("generators")
        if not isinstance(generators, list) or not TARGET.issubset(generators):
            raise ValueError("checkpoint lacks the requested parameter names")

        exact_rows: Counter[str] = Counter()
        scalar_rows: Counter[str] = Counter()
        exact_tags: dict[str, list[str]] = defaultdict(list)
        scalar_tags: dict[str, list[str]] = defaultdict(list)
        selected = 0
        constraint_count = 0

        for row in reader.array_values():
            if not isinstance(row, list) or len(row) != 2:
                raise ValueError("malformed constraint row")
            tag, terms = row
            constraint_count += 1
            support: set[str] = set()
            outside = False
            for sparse, _values in terms:
                for index, exponent in sparse:
                    if not exponent:
                        continue
                    name = generators[index]
                    if name == "t":
                        continue
                    if name in TARGET:
                        support.add(name)
                    else:
                        outside = True
            if support != TARGET or outside:
                continue

            result = polynomial_digests(terms, generators)
            if result is None:
                continue
            exact_digest, scalar_digest = result
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
            "target_support": sorted(TARGET),
            "selected_rows": selected,
            "exact_row_multiset": exact_rows,
            "scalar_normalized_row_multiset": scalar_rows,
            "exact_tags": exact_tags,
            "scalar_tags": scalar_tags,
        }
    finally:
        reader.close()


def common_count(left: Counter, right: Counter) -> int:
    return sum((left & right).values())


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
    exact_common = common_count(left["exact_row_multiset"], right["exact_row_multiset"])
    scalar_common = common_count(left["scalar_normalized_row_multiset"],
                                 right["scalar_normalized_row_multiset"])
    report = {
        "target_support": sorted(TARGET),
        "left": {key: left[key] for key in
                 ("source", "source_sha256", "nparameters", "constraints", "selected_rows")},
        "right": {key: right[key] for key in
                  ("source", "source_sha256", "nparameters", "constraints", "selected_rows")},
        "exact_polynomial_multisets_equal": (
            left["exact_row_multiset"] == right["exact_row_multiset"]
        ),
        "scalar_normalized_equation_multisets_equal": (
            left["scalar_normalized_row_multiset"]
            == right["scalar_normalized_row_multiset"]
        ),
        "exact_common_rows": exact_common,
        "scalar_normalized_common_rows": scalar_common,
        "left_only_exact_tags": unmatched_tags(
            left["exact_row_multiset"], left["exact_tags"],
            right["exact_row_multiset"]
        ),
        "right_only_exact_tags": unmatched_tags(
            right["exact_row_multiset"], right["exact_tags"],
            left["exact_row_multiset"]
        ),
        "left_only_scalar_equation_tags": unmatched_tags(
            left["scalar_normalized_row_multiset"], left["scalar_tags"],
            right["scalar_normalized_row_multiset"]
        ),
        "right_only_scalar_equation_tags": unmatched_tags(
            right["scalar_normalized_row_multiset"], right["scalar_tags"],
            left["scalar_normalized_row_multiset"]
        ),
        "scope": (
            "Exact comparison of saved equations with exactly this parameter support; "
            "it does not prove inconsistency or exclude the full branch."
        ),
    }
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"report={args.report}")


if __name__ == "__main__":
    main()
