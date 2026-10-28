"""Export a saved sparse-ring lift checkpoint as an exact Singular ideal.

The output ideal consists of every retained polynomial constraint in the
checkpoint.  For quadratic checkpoints the coefficient field is Q(sqrt(5));
otherwise it is Q.  This is a reproducibility helper, not by itself an
exclusion certificate: Singular must report the unit ideal, and the resulting
lift representation should be independently checked before claiming one.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def rational(text: str) -> str:
    """Normalize a Python Fraction string into a Singular rational literal."""
    if "/" in text:
        numerator, denominator = text.split("/", 1)
        return f"({numerator}/{denominator})"
    return f"({text})"


def coefficient(values: list[str], quadratic: bool) -> str:
    if not quadratic:
        if len(values) != 1:
            raise ValueError(f"expected one rational coefficient, got {values}")
        return rational(values[0])
    # AlgebraicField.to_list trims the zero sqrt(5) coefficient, so rational
    # elements in a quadratic checkpoint can be stored as a one-entry list.
    if len(values) == 1:
        return rational(values[0])
    if len(values) != 2:
        raise ValueError(f"expected Q(sqrt(5)) basis coefficients, got {values}")
    # The sparse checkpoint stores coefficients in descending powers of sqrt(5).
    return f"({rational(values[0])}*a+{rational(values[1])})"


def export(checkpoint: Path, output: Path,
           zero_variable: str | None = None,
           algorithm: str = "std") -> dict[str, object]:
    state = json.loads(checkpoint.read_text(encoding="utf-8"))
    if state.get("format") != "f2-sparse-ring-v1":
        raise ValueError("unexpected checkpoint format")
    constraints = state.get("constraints", [])
    generators = state["generators"]
    substituted_index: int | None = None
    if zero_variable is not None:
        if zero_variable not in generators:
            raise ValueError(f"unknown generator {zero_variable}")
        substituted_index = generators.index(zero_variable)
        # Require an explicit retained equation c*zero_variable=0 before
        # specializing.  This makes the simplification a certified
        # consequence of the checkpoint ideal, not a guessed parameter value.
        proves_zero = False
        for _tag, poly in constraints:
            if len(poly) == 1:
                sparse, values = poly[0]
                if sparse == [[substituted_index, 1]]:
                    if any(value != "0" for value in values):
                        proves_zero = True
                        break
        if not proves_zero:
            raise ValueError(
                f"no retained scalar-multiple equation for {zero_variable}=0")

        specialized = []
        for tag, poly in constraints:
            kept = [(sparse, values) for sparse, values in poly
                    if all(index != substituted_index for index, _power in sparse)]
            if kept:
                specialized.append((tag, kept))
        constraints = specialized

    active = sorted({index for _tag, poly in constraints
                     for sparse, _values in poly for index, _power in sparse})
    variables = [generators[index] for index in active]
    if not variables:
        variables = ["dummy"]

    def polynomial(terms: list[list[object]]) -> str:
        pieces: list[str] = []
        for sparse, values in terms:
            term = coefficient(values, bool(state["quadratic"]))
            for index, power in sparse:
                term += f"*{generators[index]}" + (f"^{power}" if power != 1 else "")
            pieces.append(term)
        return "+".join(pieces) if pieces else "0"

    lines = ["// Generated from an exact sparse-ring checkpoint."]
    if state["quadratic"]:
        lines.append(f"ring r=(0,a),({','.join(variables)}),dp;")
        lines.append("minpoly=a2-5;")
    else:
        lines.append(f"ring r=0,({','.join(variables)}),dp;")
    lines.append("ideal I =")
    if constraints:
        for i, (_tag, poly) in enumerate(constraints):
            suffix = "," if i + 1 < len(constraints) else ";"
            lines.append("  " + polynomial(poly) + suffix)
    else:
        lines.append("  0;")
    if algorithm not in {"std", "slimgb"}:
        raise ValueError("algorithm must be std or slimgb")
    lines.extend([
        'print("GENERATORS");',
        "print(size(I));",
        'print("ACTIVE_VARIABLES");',
        f'print("{len(active)}");',
        'print("COMPUTING_STANDARD_BASIS");',
        f"ideal G={algorithm}(I);",
        'print("BASIS_SIZE");',
        "print(size(G));",
        'if (size(G)==1 && G[1]==1) {',
        '  print("UNIT_IDEAL");',
        "  matrix L=lift(I,1);",
        '  print("LIFT_CERTIFICATE");',
        "  print(L);",
        "  poly C=0;",
        "  for (int j=1; j<=size(L); j++) { C=C+L[j,1]*I[j]; }",
        '  print("LIFTED_COMBINATION_CHECK");',
        "  print(C);",
        '  if (C==1) { print("CERTIFICATE_VERIFIED"); }',
        '  else { print("CERTIFICATE_CHECK_FAILED"); }',
        "} else {",
        '  print("NOT_UNIT_IDEAL");',
        "}",
        "",
    ])
    output.write_text("\n".join(lines), encoding="utf-8")
    return {
        "branch": state["branch"],
        "quadratic": bool(state["quadratic"]),
        "constraint_count": len(constraints),
        "term_count": sum(len(poly) for _tag, poly in constraints),
        "active_variable_count": len(active),
        "active_variables": variables,
        "certified_zero_substitution": zero_variable,
        "algorithm": algorithm,
        "singular_file": str(output),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--set-zero", metavar="GENERATOR",
                        help="specialize a generator already forced to zero")
    parser.add_argument("--algorithm", choices=("std", "slimgb"), default="std")
    args = parser.parse_args()
    print(json.dumps(export(args.checkpoint, args.output, args.set_zero,
                            args.algorithm), indent=2))


if __name__ == "__main__":
    main()
