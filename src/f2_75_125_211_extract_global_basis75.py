#!/usr/bin/env python3
"""Extract the exact independent Q-linear basis of the 139 compatibility rows."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import ijson


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_depth14_global_compatibility_reduction_20261021.json"
RANK_RECEIPT = HERE / "f2_75_125_211_global_compatibility_exact_span.json"
OUT = HERE / "f2_75_125_211_global_compatibility_basis75.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    source_hash = sha256_file(SOURCE)
    rank_hash = sha256_file(RANK_RECEIPT)
    rank = json.loads(RANK_RECEIPT.read_text(encoding="utf-8"))
    if rank["source_sha256"] != source_hash or rank["coefficient_field"] != "Q":
        raise AssertionError("exact span receipt does not match the source over Q")
    selected = [int(index) for index in rank["independent_source_rows"]]
    if len(selected) != rank["exact_vector_space_rank"] or len(selected) != 75:
        raise AssertionError("unexpected exact independent-row count")
    selected_set = set(selected)

    with SOURCE.open("rb") as stream:
        old_parameters = list(ijson.items(stream, "old_parameters.item"))
    with SOURCE.open("rb") as stream:
        weights = list(ijson.items(stream, "old_parameter_weights.item"))
    rows = []
    with SOURCE.open("rb") as stream:
        for equation in ijson.items(stream, "global_compatibility_equations.item"):
            if int(equation["constraint_row_zero_based"]) in selected_set:
                rows.append(equation)
    rows.sort(key=lambda equation: int(equation["constraint_row_zero_based"]))
    if [int(row["constraint_row_zero_based"]) for row in rows] != sorted(selected):
        raise AssertionError("failed to extract every independent source row")

    output = {
        "branch": "2+1+1",
        "coefficient_field": "Q",
        "scope": (
            "The 75 rows selected by exact Q-linear row reduction of the 139 global "
            "weighted-degree-14 compatibility equations. Their Q-span equals the full "
            "139-row span according to the paired exact-span receipt. This is not a "
            "Groebner basis or a branch-exclusion certificate."
        ),
        "source": SOURCE.name,
        "source_sha256": source_hash,
        "exact_span_receipt": RANK_RECEIPT.name,
        "exact_span_receipt_sha256": rank_hash,
        "old_parameters": old_parameters,
        "old_parameter_weights": [int(value) for value in weights],
        "original_equation_count": int(rank["equation_count"]),
        "exact_vector_space_rank": int(rank["exact_vector_space_rank"]),
        "dependent_equation_count": len(rank["dependent_source_rows"]),
        "basis_source_rows": rows,
    }
    OUT.write_text(json.dumps(output, separators=(",", ":")) + "\n", encoding="utf-8")
    print(json.dumps({
        "source_sha256": source_hash,
        "exact_span_receipt_sha256": rank_hash,
        "original_equation_count": output["original_equation_count"],
        "exact_vector_space_rank": output["exact_vector_space_rank"],
        "dependent_equation_count": output["dependent_equation_count"],
        "basis_rows_extracted": len(rows),
        "output_bytes": OUT.stat().st_size,
        "output_sha256": sha256_file(OUT),
        "output": str(OUT),
    }, indent=2))


if __name__ == "__main__":
    main()
