#!/usr/bin/env python3
"""Export the exact independent 80-row basis of the 149 terminal-Q equations.

The source/span hashes are checked before extraction.  Since the rows are
homogeneous of the same weighted degree and the span computation is over Q,
the selected rows generate exactly the same Q-vector space as all 149
forcing polynomials.  This is an input-size reduction only, not an ideal
emptiness or branch-exclusion certificate.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_terminalq_global_compatibility_reduction_20261021.json"
SPAN = HERE / "f2_75_125_211_terminalq_global_compatibility_exact_span_20261021.json"
OUTPUT = HERE / "f2_75_125_211_terminalq_exact_independent_generators_20261021.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def iter_json_array_items(path: Path, key: str):
    """Stream one top-level JSON array without loading the document at once."""
    marker = f'"{key}"'
    decoder = json.JSONDecoder()
    chunk_size = 1 << 20

    with path.open("r", encoding="utf-8") as stream:
        overlap = ""
        buffer = ""
        while True:
            chunk = stream.read(chunk_size)
            if not chunk:
                raise ValueError(f"top-level key {key!r} was not found")
            candidate = overlap + chunk
            marker_at = candidate.find(marker)
            if marker_at >= 0:
                buffer = candidate[marker_at + len(marker):]
                break
            overlap = candidate[-(len(marker) + 8):]

        index = 0
        while True:
            while index < len(buffer) and buffer[index].isspace():
                index += 1
            if index < len(buffer):
                if buffer[index] != ":":
                    raise ValueError(f"malformed value after top-level key {key!r}")
                index += 1
                break
            chunk = stream.read(chunk_size)
            if not chunk:
                raise ValueError(f"missing array after top-level key {key!r}")
            buffer = chunk
            index = 0

        while True:
            while index < len(buffer) and buffer[index].isspace():
                index += 1
            if index == len(buffer):
                chunk = stream.read(chunk_size)
                if not chunk:
                    raise ValueError(f"unterminated array for top-level key {key!r}")
                buffer = chunk
                index = 0
                continue
            if buffer[index] != "[":
                raise ValueError(f"top-level value {key!r} is not an array")
            index += 1
            break

        while True:
            while index < len(buffer) and (buffer[index].isspace() or buffer[index] == ","):
                index += 1
            if index == len(buffer):
                chunk = stream.read(chunk_size)
                if not chunk:
                    raise ValueError(f"unterminated array for top-level key {key!r}")
                buffer = chunk
                index = 0
                continue
            if buffer[index] == "]":
                return
            try:
                item, end = decoder.raw_decode(buffer, index)
            except json.JSONDecodeError:
                chunk = stream.read(chunk_size)
                if not chunk:
                    raise
                buffer = buffer[index:] + chunk
                index = 0
                continue
            yield item
            buffer = buffer[end:]
            index = 0


def main() -> None:
    span = json.loads(SPAN.read_text(encoding="utf-8"))
    source_hash = sha256(SOURCE)
    if span["source_sha256"] != source_hash:
        raise SystemExit("span receipt does not match the source equations")
    if span["equation_count"] != 149 or span["exact_vector_space_rank"] != 80:
        raise SystemExit("unexpected equation count or exact span rank")

    selected = set(span["independent_equation_rows"])
    if len(selected) != 80:
        raise SystemExit("independent row index list is not a rank-80 set")

    rows = []
    for row in iter_json_array_items(SOURCE, "compatibility_equations"):
        if row["row"] in selected:
            rows.append(row)

    rows.sort(key=lambda row: row["row"])
    if [row["row"] for row in rows] != sorted(selected):
        raise SystemExit("could not recover every independent source row")

    result = {
        "branch": span["branch"],
        "scope": "80 independent source generators spanning exactly the 149 terminal-Q forcing polynomials over Q; specialized p43=0 terminal chart only",
        "coefficient_field": "Q",
        "source": SOURCE.name,
        "source_sha256": source_hash,
        "span_receipt": SPAN.name,
        "span_receipt_sha256": sha256(SPAN),
        "exact_vector_space_rank": 80,
        "input_equation_count": 149,
        "exported_generator_count": len(rows),
        "generators": rows,
        "interpretation": "Exact generator-span reduction only; this does not prove that the generated ideal is empty or exclude the full 2+1+1 branch.",
    }
    OUTPUT.write_text(json.dumps(result, separators=(",", ":")) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "branch", "source_sha256", "exact_vector_space_rank",
        "input_equation_count", "exported_generator_count", "interpretation",
    )}, indent=2), flush=True)
    print(f"output={OUTPUT.name} sha256={sha256(OUTPUT)} bytes={OUTPUT.stat().st_size}", flush=True)


if __name__ == "__main__":
    main()
