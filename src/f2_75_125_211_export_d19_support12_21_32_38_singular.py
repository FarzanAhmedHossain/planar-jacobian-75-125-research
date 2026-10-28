"""Export the verified 2+1+1 d19 support ideal to Singular syntax.

The existing builder performs all source-hash, exact row-reduction, and
weighted-homogeneity checks before this exporter intercepts the Groebner call.
No algebra is changed: the exported input consists of those same 18 rational
generators in the 28-variable support ring.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import f2_75_125_211_d19_support12_21_32_38_combined_groebner as builder


class ExportComplete(Exception):
    pass


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default=str(HERE / "f2_75_125_211_d19_support12_21_32_38_20261025.sing"),
    )
    parser.add_argument(
        "--basis-output",
        default=str(HERE / "f2_75_125_211_d19_support12_21_32_38_20261025_basis.txt"),
    )
    args = parser.parse_args()
    output = Path(args.output)
    basis_output = Path(args.basis_output)
    output.parent.mkdir(parents=True, exist_ok=True)
    basis_output.parent.mkdir(parents=True, exist_ok=True)

    def export_exact_inputs(inputs, _ring, _method):
        analysis = json.loads(builder.ANALYSIS38.read_text(encoding="utf-8"))
        names = analysis["selected_variables"]
        weights = [int(value) for value in analysis["primitive_positive_integer_weights"]]
        if len(inputs) != 18 or len(names) != 28 or len(weights) != 28:
            raise ValueError("verified support system dimensions changed")

        sources = (
            builder.FAMILY21,
            builder.ANALYSIS21,
            builder.FAMILY32,
            builder.ANALYSIS32,
            builder.FAMILY38,
            builder.ANALYSIS38,
            builder.SUPPORT32_REDUNDANCY,
            builder.SUPPORT38_LOW_REDUCTION,
            builder.SUPPORT12_GB,
            builder.SUPPORT12_INPUT,
            builder.SUPPORT12_REUSE,
        )
        lines = [
            "// Exact 2+1+1 d19 support-12/21/32/38 necessary ideal.",
            "// Coefficients are in QQ; weights define a global weighted order.",
        ]
        for source in sources:
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            lines.append(f"// {source.name} sha256={digest}")
        lines.append(f"ring R=0,({','.join(names)}),wp({','.join(map(str, weights))});")
        for index, polynomial in enumerate(inputs):
            expression = str(polynomial).replace("**", "^")
            lines.append(f"poly f{index}={expression};")
        lines.append("ideal I=" + ",".join(f"f{index}" for index in range(len(inputs))) + ";")
        lines.append("option(redSB);")
        lines.append("ideal G=slimgb(I);")
        lines.append('print("EXACT_SUPPORT_IDEAL_RESULT");')
        lines.append('print("basis_size"); print(size(G));')
        lines.append('print("one_remainder"); print(reduce(1,G));')
        lines.append(f'write("{basis_output.as_posix()}",string(G));')
        lines.append("quit;")
        output.write_text("\n".join(lines) + "\n", encoding="ascii")
        print(f"singular_input={output}", flush=True)
        print(f"basis_output={basis_output}", flush=True)
        print(f"input_sha256={hashlib.sha256(output.read_bytes()).hexdigest()}", flush=True)
        raise ExportComplete

    builder.ring_groebner = export_exact_inputs
    try:
        builder.main()
    except ExportComplete:
        return


if __name__ == "__main__":
    main()
