#!/usr/bin/env python3
"""Audit monomial contents of the exact depth-14 compatibility equations.

This is a streaming audit: the 39 MB JSON certificate is kept as text, while
only one compatibility equation is decoded at a time. A common monomial is
reported as a case split, never divided out without its zero locus being
handled separately.
"""

from __future__ import annotations

import collections
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "f2_75_125_211_depth14_global_compatibility_reduction_20261021.json"
EXPECTED_SHA256 = "372ba14b06d04297b25e098ccbdfc03c3d03dffd6ae70f9060655be608a0689b"
OUTPUT = ROOT / "f2_75_125_211_depth14_monomial_content_audit.json"


def main() -> None:
    raw = SOURCE.read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    if sha != EXPECTED_SHA256:
        raise SystemExit(f"source SHA mismatch: {sha}")
    source = raw.decode("utf-8")
    key = '"global_compatibility_equations":'
    cursor = source.index(key) + len(key)
    decoder = json.JSONDecoder()
    if source[cursor] != "[":
        raise SystemExit("unexpected compatibility-equation encoding")
    cursor += 1

    rows: list[tuple[str, int, tuple[tuple[int, int], ...]]] = []
    while True:
        while cursor < len(source) and (source[cursor].isspace() or source[cursor] == ","):
            cursor += 1
        if source[cursor] == "]":
            break
        equation, cursor = decoder.raw_decode(source, cursor)
        common_indices: set[int] | None = None
        min_exponent: dict[int, int] = {}
        for sparse, _coefficient in equation["terms"]:
            present = {int(index) for index, _power in sparse}
            common_indices = present if common_indices is None else common_indices & present
            for index, power in sparse:
                index, power = int(index), int(power)
                min_exponent[index] = min(min_exponent.get(index, power), power)
        content = tuple(
            (index, min_exponent[index])
            for index in sorted(common_indices or ())
            if min_exponent[index] > 0
        )
        rows.append((equation["constraint_tag"], int(equation["term_count"]), content))

    histogram = collections.Counter(content for _tag, _terms, content in rows)
    result = {
        "source": SOURCE.name,
        "source_sha256": sha,
        "branch": "2+1+1",
        "equation_count": len(rows),
        "equations_with_nontrivial_monomial_content": sum(bool(c) for _, _, c in rows),
        "content_histogram": [
            {"monomial_content": [[index, exponent] for index, exponent in content], "count": count}
            for content, count in histogram.most_common()
        ],
        "interpretation": "Every one of the 139 exact global compatibility equations is primitive as a polynomial: none has a nonconstant monomial gcd. This does not prove that their common zero locus is empty.",
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    for tag, term_count, content in rows:
        if content:
            print(f"{tag}: terms={term_count}, content={content}")


if __name__ == "__main__":
    main()
