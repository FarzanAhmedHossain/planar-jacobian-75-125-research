"""Stream a sparse F2 checkpoint and find parameter slots used in its state.

This scan decodes sparse monomial exponent lists only. It does not construct
polynomials or alter the checkpoint. A parameter can be dropped from a later
ring compaction only if it is absent from both every retained equation and
every saved state polynomial; the resulting support report is not an
exclusion certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from f2_75_125_211_streaming_checkpoint_scan import JsonStream


SINGLE_POLYNOMIAL_FIELDS = {"h", "v", "S", "B", "kappa", "kappa2", "kappa3"}
POLYNOMIAL_ARRAY_FIELDS = {"uk", "p", "q", "w"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def scan_polynomial(reader: JsonStream, counts: Counter[int], label: str,
                    progress: dict[str, int]) -> int:
    """Scan one serialized polynomial, retaining only generator occurrence counts."""
    reader.take("[")
    first = True
    terms = 0
    while reader.peek() != "]":
        if not first:
            reader.take(",")
        first = False
        reader.take("[")
        sparse = reader.value()
        reader.take(",")
        reader.value()  # coefficient encoding; support scanning needs no arithmetic
        reader.take("]")
        terms += 1
        for generator, exponent in sparse:
            if exponent:
                counts[int(generator)] += 1
        if terms % 1_000_000 == 0:
            print(f"field={label} terms={terms}", flush=True)
    reader.take("]")
    progress["polynomial_count"] += 1
    progress["term_count"] += terms
    print(f"field={label} complete terms={terms} "
          f"polynomials={progress['polynomial_count']} "
          f"cumulative_terms={progress['term_count']}", flush=True)
    return terms


def scan_polynomial_array(reader: JsonStream, counts: Counter[int], label: str,
                          progress: dict[str, int]) -> int:
    reader.take("[")
    first = True
    terms = 0
    while reader.peek() != "]":
        if not first:
            reader.take(",")
        first = False
        terms += scan_polynomial(reader, counts, label, progress)
    reader.take("]")
    return terms


def scan_constraints(reader: JsonStream, counts: Counter[int],
                     progress: dict[str, int]) -> int:
    rows = 0
    first = True
    while reader.peek() != "]":
        if not first:
            reader.take(",")
        first = False
        reader.take("[")
        tag = reader.value()
        reader.take(",")
        scan_polynomial(reader, counts, f"constraint:{tag}", progress)
        reader.take("]")
        rows += 1
    reader.take("]")
    return rows


def scan_uk_object(reader: JsonStream, counts: Counter[int],
                   progress: dict[str, int]) -> int:
    reader.take("{")
    first = True
    terms = 0
    while reader.peek() != "}":
        if not first:
            reader.take(",")
        first = False
        key = reader.value()
        reader.take(":")
        terms += scan_polynomial(reader, counts, f"uk:{key}", progress)
    reader.take("}")
    return terms


def scan(source: Path, expected_sha256: str) -> dict:
    actual_sha256 = sha256(source)
    if actual_sha256.lower() != expected_sha256.lower():
        raise ValueError(f"checkpoint hash mismatch: {actual_sha256}")

    reader = JsonStream(source)
    counts: Counter[int] = Counter()
    progress = {"polynomial_count": 0, "term_count": 0}
    try:
        metadata, _ = reader.top_level_constraints()
        if metadata.get("format") != "f2-sparse-ring-v1":
            raise ValueError("unexpected checkpoint format")
        if metadata.get("branch") != "2+1+1" or metadata.get("quadratic"):
            raise ValueError("expected a rational 2+1+1 checkpoint")
        generators = metadata["generators"]
        p_indices = {index for index, name in enumerate(generators)
                     if name.startswith("p") and name[1:].isdigit()}
        nparameters = int(metadata["nparameters"])

        constraint_rows = scan_constraints(reader, counts, progress)
        if constraint_rows != int(metadata.get("retained_constraint_count", 0)):
            # Some checkpoint metadata omits the count; only enforce when present.
            if "retained_constraint_count" in metadata:
                raise AssertionError("constraint row count disagrees with metadata")

        reader.take(",")
        fields_seen = []
        while reader.peek() != "}":
            key = reader.value()
            reader.take(":")
            fields_seen.append(key)
            if key in SINGLE_POLYNOMIAL_FIELDS:
                scan_polynomial(reader, counts, key, progress)
            elif key == "uk":
                scan_uk_object(reader, counts, progress)
            elif key in {"p", "q", "w"}:
                scan_polynomial_array(reader, counts, key, progress)
            else:
                reader.value()
            if reader.peek() == ",":
                reader.take(",")
            elif reader.peek() != "}":
                raise ValueError(f"unexpected top-level separator after {key}")
        reader.take("}")
    finally:
        reader.close()

    active_parameters = sorted(index for index in p_indices if counts[index])
    inactive_parameters = sorted(p_indices - set(active_parameters))
    report = {
        "source": str(source.resolve()),
        "source_sha256": actual_sha256,
        "branch": metadata["branch"],
        "metadata_nparameters": nparameters,
        "parameter_slot_count": len(p_indices),
        "constraint_rows_scanned": constraint_rows,
        "polynomial_arrays_scanned": progress["polynomial_count"],
        "serialized_terms_scanned": progress["term_count"],
        "fields_scanned": fields_seen,
        "active_parameters": [generators[index] for index in active_parameters],
        "inactive_parameters": [generators[index] for index in inactive_parameters],
        "active_parameter_count": len(active_parameters),
        "inactive_parameter_count": len(inactive_parameters),
        "scope": (
            "Exact support-only scan of every equation and saved state polynomial. "
            "An inactive slot is absent from the entire checkpoint and may be "
            "removed by a semantics-preserving generator compaction; this report "
            "does not exclude any branch."
        ),
    }
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("expected_sha256")
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    if args.report.exists():
        raise FileExistsError(args.report)
    report = scan(args.checkpoint, args.expected_sha256)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in (
        "source_sha256", "metadata_nparameters", "parameter_slot_count",
        "constraint_rows_scanned", "polynomial_arrays_scanned",
        "serialized_terms_scanned", "active_parameter_count",
        "inactive_parameter_count", "inactive_parameters", "scope",
    )}, indent=2), flush=True)
    print(f"report={args.report}", flush=True)


if __name__ == "__main__":
    main()
