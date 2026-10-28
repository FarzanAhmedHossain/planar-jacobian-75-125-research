"""Export the d19 overlap ideal cut by the exact support-three radical line.

The support-12/21/32/38 builder validates and constructs its 18 exact
necessary generators.  We add a verified four-polynomial Groebner basis for
the support-15 slice and the two linear equations defining the exact support-3
zero line.  A unit ideal would exclude this intersection only; a nonunit result
does not exclude the whole 2+1+1 branch.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import f2_75_125_211_d19_support12_21_32_38_combined_groebner as builder


SOURCE = HERE / "f2_75_125_sparse_ring_lift_211_roots19_u1zero_streamed_state.json"
SUPPORT15 = HERE / "f2_75_125_sparse_ring_lift_211_roots19_u1zero_support15_gb.json"
SUPPORT3_LINE = HERE / "f2_75_125_211_support3_radical_line_20261025.json"


class ExportComplete(Exception):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default=str(HERE / "f2_75_125_211_d19_support12_21_32_38_plus_support15_line_20261025.sing"),
    )
    parser.add_argument(
        "--basis-output",
        default=str(HERE / "f2_75_125_211_d19_support12_21_32_38_plus_support15_line_20261025_basis.txt"),
    )
    args = parser.parse_args()

    support15 = json.loads(SUPPORT15.read_text(encoding="utf-8"))
    support3 = json.loads(SUPPORT3_LINE.read_text(encoding="utf-8"))
    expected_source_hash = support15["source_sha256"]
    if expected_source_hash != support3["source_sha256"]:
        raise ValueError("support-15 and support-3 reports do not share a source checkpoint")
    if support15.get("branch") != "2+1+1" or support15.get("unit_ideal"):
        raise ValueError("expected the nonunit exact support-15 report")
    if support15.get("groebner_basis_length") != 4 or len(support15.get("groebner_basis", [])) != 4:
        raise ValueError("support-15 Groebner basis report is incomplete")
    if not support3.get("line_containment_verified") or not support3.get("basis_restrictions_verified"):
        raise ValueError("support-3 radical line report has not passed its exact checks")

    # Checkpoint hash validation is deferred until the 1 GB checkpoint has
    # finished being read by the active lift process.  The reports themselves
    # carry the exact source hash, recorded in the emitted Singular header.
    names15 = tuple(support15["selected_variables"])
    for expression in support15["groebner_basis"]:
        used = set(re.findall(r"\bp\d+\b", expression))
        if not used <= set(names15):
            raise ValueError(f"support-15 basis uses undeclared variables: {sorted(used-set(names15))}")

    output = Path(args.output)
    basis_output = Path(args.basis_output)

    def export_exact_inputs(inputs, _ring, _method):
        analysis = json.loads(builder.ANALYSIS38.read_text(encoding="utf-8"))
        names28 = tuple(analysis["selected_variables"])
        names = tuple(dict.fromkeys((*names28, *names15)))
        if len(inputs) != 18 or len(names28) != 28:
            raise ValueError("verified support-12/21/32/38 system dimensions changed")
        # All generators are exact rational polynomials.  Singular's dp order
        # avoids asserting a weight compatibility between the separate slices.
        lines = [
            "// Exact 2+1+1 d19 necessary ideal plus support-3 radical line.",
            "// Unit ideal excludes this intersection of necessary-condition slices only.",
            f"// source_sha256={expected_source_hash}",
            f"// checkpoint_sha256_current={sha256(SOURCE)}",
            f"// support15_report_sha256={sha256(SUPPORT15)}",
            f"// support3_line_report_sha256={sha256(SUPPORT3_LINE)}",
            f"ring R=0,({','.join(names)}),dp;",
        ]
        all_expressions = [str(poly).replace("**", "^") for poly in inputs]
        all_expressions.extend(expr.replace("**", "^") for expr in support15["groebner_basis"])
        all_expressions.extend(("p14+14*p15", "5*p13-452*p15"))
        for index, expression in enumerate(all_expressions):
            lines.append(f"poly f{index}={expression};")
        lines.append("ideal I=" + ",".join(f"f{index}" for index in range(len(all_expressions))) + ";")
        lines.append("option(redSB);")
        lines.append("ideal G=slimgb(I);")
        lines.append('print("EXACT_SUPPORT_INTERSECTION_RESULT");')
        lines.append('print("basis_size"); print(size(G));')
        lines.append('print("one_remainder"); print(reduce(1,G));')
        lines.append(f'write("{basis_output.as_posix()}",string(G));')
        lines.append("quit;")
        output.parent.mkdir(parents=True, exist_ok=True)
        basis_output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("\n".join(lines) + "\n", encoding="ascii")
        print(f"singular_input={output}", flush=True)
        print(f"basis_output={basis_output}", flush=True)
        print(f"generators={len(all_expressions)} variables={len(names)}", flush=True)
        print(f"input_sha256={hashlib.sha256(output.read_bytes()).hexdigest()}", flush=True)
        raise ExportComplete

    builder.ring_groebner = export_exact_inputs
    try:
        builder.main()
    except ExportComplete:
        return


if __name__ == "__main__":
    main()
