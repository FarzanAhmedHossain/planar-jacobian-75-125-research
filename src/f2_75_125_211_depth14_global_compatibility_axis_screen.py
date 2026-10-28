#!/usr/bin/env python3
"""Exact axis screen of the 139 global depth-14 compatibility equations.

These equations are weighted-homogeneous of degree 14 in 42 old parameters.
On a coordinate axis, each equation restricts either to zero or to a scalar
multiple of one variable's pure power. Thus one pass identifies exactly
which coordinate axes are excluded by this necessary ideal. This is not a
full-branch exclusion.
"""

from __future__ import annotations

import hashlib
import json
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_depth14_global_compatibility_reduction_20261021.json"
OUTPUT = HERE / "f2_75_125_211_depth14_global_compatibility_axis_screen.json"
CHUNK = 1 << 20
VARIABLES = [*(f"p{i}" for i in range(36)), *(f"p{i}" for i in range(37, 43))]
WEIGHTS_BY_INDEX = {
    **{0: 1, 1: 2, 2: 3, 3: 3, 4: 3},
    **{i: 4 for i in range(5, 13)},
    **{i: 5 for i in range(13, 20)},
    **{i: 6 for i in range(20, 25)},
    **{i: 7 for i in range(25, 30)},
    **{i: 8 for i in range(30, 36)},
    **{i: 9 for i in range(37, 43)},
}
WEIGHTS = [WEIGHTS_BY_INDEX[int(name[1:])] for name in VARIABLES]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def source_hash_and_array_offset(path: Path, key: str) -> tuple[str, int]:
    digest = hashlib.sha256()
    needle = json.dumps(key).encode("ascii")
    overlap = b""
    base = 0
    offset = None
    with path.open("rb") as stream:
        while chunk := stream.read(CHUNK):
            digest.update(chunk)
            data = overlap + chunk
            pos = 0
            while True:
                found = data.find(needle, pos)
                if found < 0:
                    break
                suffix = data[found + len(needle) :]
                j = 0
                while j < len(suffix) and suffix[j] in b" \t\r\n":
                    j += 1
                if j < len(suffix) and suffix[j] == ord(":"):
                    j += 1
                    while j < len(suffix) and suffix[j] in b" \t\r\n":
                        j += 1
                    if j < len(suffix) and suffix[j] == ord("["):
                        offset = base - len(overlap) + found + len(needle) + j + 1
                        break
                pos = found + 1
            if offset is not None:
                while chunk := stream.read(CHUNK):
                    digest.update(chunk)
                break
            overlap = data[-(len(needle) + 32) :]
            base += len(chunk)
    if offset is None:
        raise ValueError(f"could not locate array {key!r}")
    return digest.hexdigest(), offset


def array_items(path: Path, offset: int):
    decoder = json.JSONDecoder()
    with path.open("rb") as stream:
        stream.seek(offset)
        buffer = ""
        while True:
            buffer = buffer.lstrip()
            if buffer.startswith(","):
                buffer = buffer[1:]
                continue
            if buffer.startswith("]"):
                return
            if not buffer:
                raw = stream.read(CHUNK)
                if not raw:
                    raise ValueError("truncated compatibility-equation array")
                buffer = raw.decode("utf-8")
                continue
            try:
                item, end = decoder.raw_decode(buffer)
            except json.JSONDecodeError:
                raw = stream.read(CHUNK)
                if not raw:
                    raise ValueError("truncated compatibility-equation item")
                buffer += raw.decode("utf-8")
                continue
            yield item
            buffer = buffer[end:]


def main() -> None:
    source = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else SOURCE
    output = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else OUTPUT
    source_hash, array_offset = source_hash_and_array_offset(
        source, "global_compatibility_equations"
    )
    hit_rows: dict[str, list[dict[str, object]]] = {name: [] for name in VARIABLES}
    equations = term_count = 0
    support_histogram: Counter[int] = Counter()

    for row in array_items(source, array_offset):
            equations += 1
            if row["constraint_row_zero_based"] < 0:
                raise ValueError("invalid source row index")
            seen_powers: dict[str, Fraction] = {}
            for sparse, coefficient_text in row["terms"]:
                term_count += 1
                coefficient = Fraction(coefficient_text)
                weighted_degree = 0
                support = []
                exponent_by_variable = {}
                for index, power in sparse:
                    if index < 0 or index >= len(VARIABLES):
                        raise ValueError(f"unexpected variable index {index}")
                    weighted_degree += WEIGHTS[index] * power
                    if power:
                        support.append(index)
                        exponent_by_variable[index] = power
                if weighted_degree != 14:
                    raise AssertionError(f"row {row['constraint_tag']} has weighted degree {weighted_degree}")
                support_histogram[len(support)] += 1
                if len(support) == 1:
                    index = support[0]
                    name = VARIABLES[index]
                    power = exponent_by_variable[index]
                    seen_powers[name] = seen_powers.get(name, Fraction(0)) + coefficient

            for name, coefficient in seen_powers.items():
                if coefficient:
                    hit_rows[name].append({
                        "constraint_row_zero_based": row["constraint_row_zero_based"],
                        "constraint_tag": row["constraint_tag"],
                        "coefficient": str(coefficient),
                    })

            if equations % 25 == 0:
                print(f"axis_screen equations={equations} terms={term_count}", flush=True)

    if equations != 139:
        raise AssertionError(f"expected 139 global compatibility equations, found {equations}")

    axes = {}
    for name in VARIABLES:
        axes[name] = {
            "status": "axis_excluded" if hit_rows[name] else "entire_axis_survives_truncation",
            "pure_power_rows": hit_rows[name],
        }

    result = {
        "source": source.name,
        "source_sha256": source_hash,
        "branch": "2+1+1",
        "scope": "exact coordinate-axis restrictions of 139 weighted-degree-14 necessary compatibility equations; not full-branch exclusion",
        "equations": equations,
        "sparse_terms": term_count,
        "weighted_degree": 14,
        "variables": VARIABLES,
        "weights": WEIGHTS,
        "term_support_size_histogram": dict(sorted(support_histogram.items())),
        "axis_count": len(axes),
        "axes": axes,
        "axis_status_counts": {
            status: sum(entry["status"] == status for entry in axes.values())
            for status in sorted({entry["status"] for entry in axes.values()})
        },
    }
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    summary = {key: result[key] for key in (
        "source_sha256", "equations", "sparse_terms", "axis_count",
        "axis_status_counts", "term_support_size_histogram",
    )}
    summary["excluded_axes"] = [name for name, row in axes.items() if row["status"] == "axis_excluded"]
    summary["surviving_axes"] = [name for name, row in axes.items() if row["status"] != "axis_excluded"]
    summary["receipt"] = output.name
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
