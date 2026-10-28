"""Stream all d19 necessary equations supported on a chosen variable union.

The selected union is the union of the verified support-12/21/32/38 slice and
support-15 slice.  We keep only whole checkpoint equations whose complete
parameter support lies in that union, then add the exact support-3 radical
line.  No monomials are projected away.  A Singular unit ideal would exclude
this necessary-condition intersection only.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from f2_75_125_211_streaming_checkpoint_scan import JsonStream


HERE = Path(__file__).resolve().parent
DEFAULT_SOURCE = HERE / "f2_75_125_sparse_ring_lift_211_roots19_u1zero_streamed_state.json"
SUPPORT15 = HERE / "f2_75_125_sparse_ring_lift_211_roots19_u1zero_support15_gb.json"
SUPPORT38 = HERE / "f2_75_125_211_d19_support38_analysis_20261025.json"
SUPPORT3_LINE = HERE / "f2_75_125_211_support3_radical_line_20261025.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fmt_coeff(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def emit_poly(handle, index: int, terms, generator_names, selected_names) -> int:
    """Write a complete sparse polynomial as Singular syntax; return terms."""
    handle.write(f"poly f{index}=")
    first = True
    count = 0
    for sparse, values in terms:
        if len(values) != 1:
            raise ValueError("expected rational checkpoint coefficients")
        coefficient = Fraction(values[0])
        if not coefficient:
            continue
        factors = []
        for generator_index, exponent in sparse:
            exponent = int(exponent)
            if not exponent:
                continue
            name = generator_names[int(generator_index)]
            if name not in selected_names:
                raise AssertionError("selected row unexpectedly contains an outside variable")
            factors.append(name if exponent == 1 else f"{name}^{exponent}")
        negative = coefficient < 0
        absolute = abs(coefficient)
        if first:
            if negative:
                handle.write("-")
        else:
            handle.write("-" if negative else "+")
        if absolute != 1 or not factors:
            handle.write(fmt_coeff(absolute))
            if factors:
                handle.write("*")
        handle.write("*".join(factors) if factors else "")
        first = False
        count += 1
    if first:
        handle.write("0")
    handle.write(";\n")
    return count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument(
        "--output", type=Path,
        default=HERE / "f2_75_125_211_d19_union_line_slice_20261025.sing",
    )
    parser.add_argument(
        "--receipt", type=Path,
        default=HERE / "f2_75_125_211_d19_union_line_slice_20261025.json",
    )
    parser.add_argument(
        "--basis-output",
        default=str(HERE / "f2_75_125_211_d19_union_line_slice_20261025_basis.txt"),
    )
    args = parser.parse_args()

    support15 = json.loads(SUPPORT15.read_text(encoding="utf-8"))
    support38 = json.loads(SUPPORT38.read_text(encoding="utf-8"))
    line = json.loads(SUPPORT3_LINE.read_text(encoding="utf-8"))
    expected_hashes = {support15["source_sha256"], line["source_sha256"]}
    if len(expected_hashes) != 1:
        raise ValueError("support-15 and support-3 reports have different source hashes")
    expected_hash = next(iter(expected_hashes))
    if support15.get("branch") != "2+1+1" or support15.get("unit_ideal"):
        raise ValueError("expected the exact nonunit support-15 report")
    if (line.get("branch") != "2+1+1, U1=0, d19 necessary system"
            or line.get("geometric_zero_set") != "(p13,p14,p15) = (452*s/5,-14*s,s), s in algebraic closure"
            or not line.get("line_containment_verified")
            or not line.get("basis_restrictions_verified")):
        raise ValueError("support-3 radical line report does not verify the expected line")

    names = tuple(dict.fromkeys((
        *support38["selected_variables"],
        *support15["selected_variables"],
    )))
    if not {"p13", "p14", "p15", "p5"} <= set(names):
        raise ValueError("selected variable union omits support-3/15 variables")
    selected = set(names)

    reader = JsonStream(args.source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    included_tags = []
    included_rows = included_terms = skipped_rows = 0
    try:
        metadata, _ = reader.top_level_constraints()
        if (metadata.get("format") != "f2-sparse-ring-v1"
                or metadata.get("branch") != "2+1+1"
                or metadata.get("quadratic")):
            raise ValueError("expected a rational 2+1+1 sparse checkpoint")
        generator_names = tuple(metadata["generators"])
        selected_indices = {generator_names.index(name) for name in names}
        with args.output.open("w", encoding="ascii", newline="\n") as handle:
            handle.write("// Complete necessary equations supported on the selected union.\n")
            handle.write(f"// source_sha256_expected={expected_hash}\n")
            handle.write(f"// support15_report={SUPPORT15.name}\n")
            handle.write(f"// support38_analysis={SUPPORT38.name}\n")
            handle.write(f"// support3_line_report={SUPPORT3_LINE.name}\n")
            handle.write(f"ring R=0,({','.join(names)}),dp;\n")
            for tag, terms in reader.array_values():
                support = set()
                for sparse, _values in terms:
                    for generator_index, exponent in sparse:
                        if int(exponent):
                            support.add(int(generator_index))
                if not support <= selected_indices:
                    skipped_rows += 1
                    continue
                term_count = emit_poly(handle, included_rows, terms, generator_names, selected)
                included_tags.append(tag)
                included_terms += term_count
                included_rows += 1
                if included_rows % 25 == 0:
                    print(f"included_rows={included_rows} skipped_rows={skipped_rows}", flush=True)
            # These linear equations define the exact geometric zero set of
            # the support-3 basis, as verified in SUPPORT3_LINE.
            handle.write(f"poly f{included_rows}=p14+14*p15;\n")
            handle.write(f"poly f{included_rows+1}=5*p13-452*p15;\n")
            total = included_rows + 2
            handle.write("ideal I=" + ",".join(f"f{i}" for i in range(total)) + ";\n")
            handle.write("option(redSB);\n")
            handle.write("ideal G=slimgb(I);\n")
            handle.write('print("EXACT_D19_UNION_LINE_SLICE_RESULT");\n')
            handle.write('print("basis_size"); print(size(G));\n')
            handle.write('print("one_remainder"); print(reduce(1,G));\n')
            handle.write(f'write("{args.basis_output}",string(G));\n')
            handle.write("quit;\n")
    finally:
        reader.close()

    actual_hash = sha256(args.source)
    receipt = {
        "source": str(args.source),
        "source_sha256": actual_hash,
        "expected_source_sha256": expected_hash,
        "source_hash_matches_reports": actual_hash == expected_hash,
        "branch": metadata["branch"],
        "variables": list(names),
        "complete_constraints_selected": included_rows,
        "constraints_skipped_for_outside_variables": skipped_rows,
        "selected_sparse_terms": included_terms,
        "added_radical_line_equations": ["p14+14*p15", "5*p13-452*p15"],
        "input_generators": included_rows + 2,
        "input": str(args.output),
        "input_sha256": sha256(args.output),
        "scope": (
            "A Singular unit ideal proves no geometric zero of these complete necessary equations "
            "on the support-3 line; it does not by itself exclude the other support-3 components "
            "or the full 2+1+1 branch."
        ),
    }
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == "__main__":
    main()
