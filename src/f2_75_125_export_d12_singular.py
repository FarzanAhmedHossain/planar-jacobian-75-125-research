"""Export and pin the D12 23-parameter subsystem for exact Singular replay."""

from __future__ import annotations

import hashlib
from pathlib import Path


INPUT = Path(__file__).with_name(
    "f2_75_125_d24_u1zero_d12_window23_verified_constraints.tsv"
)
OUTPUT = Path(__file__).with_name("f2_75_125_d24_u1zero_d12_window23.sing")
EXPECTED_INPUT_SHA256 = (
    "ade20ac189079fba316e98958191681e0c60ebaa2188a4e04230b65e32d248e5"
)
CHECKPOINT_SHA256 = (
    "2369f36ce13d0d4336053ba0a7a546aa4407d04d77ee544aa4e313dbc88f2c79"
)
GENERATOR_INDICES = (
    6, 16, 25, 26, 27, 28, 29,
    35, 36, 37, 38, 39, 40,
    45, 46, 47, 48, 49,
    54, 55, 56, 57, 58,
)
EXPECTED_ROWS = 428


def term(encoded: str) -> tuple[str, bool]:
    coefficient, monomial = encoded.split(":", 1)
    if monomial == "1":
        return f"({coefficient})", True
    factors = []
    for factor in monomial.split(","):
        index_text, exponent_text = factor.split("^", 1)
        index = int(index_text)
        if index not in GENERATOR_INDICES:
            raise ValueError(f"unexpected generator index {index}")
        factors.append(f"v{index}^{int(exponent_text)}")
    return f"({coefficient})*" + "*".join(factors), False


def main() -> None:
    digest = hashlib.sha256(INPUT.read_bytes()).hexdigest()
    if digest != EXPECTED_INPUT_SHA256:
        raise ValueError(f"input digest mismatch: {digest}")
    rows = []
    constant_equations = []
    for line_number, line in enumerate(INPUT.read_text(encoding="utf-8").splitlines(), 1):
        tag, body = line.split("\t", 1)
        summands = [term(item) for item in body.split(";")]
        if any(is_constant for _expression, is_constant in summands):
            constant_equations.append(tag)
        rows.append(" + ".join(expression for expression, _ in summands))
    if len(rows) != EXPECTED_ROWS:
        raise ValueError(f"expected {EXPECTED_ROWS} equations, found {len(rows)}")

    variables = ",".join(f"v{index}" for index in GENERATOR_INDICES)
    lines = [
        f"// D24 checkpoint SHA256 (bound by extraction receipt): {CHECKPOINT_SHA256}",
        f"// D12 TSV SHA256: {digest}",
        f"// Exact necessary subsystem: {len(rows)} equations in {len(GENERATOR_INDICES)} parameters.",
        "// Constant-term equation tags: " + (", ".join(constant_equations) or "none"),
        f"ring r=0,({variables}),dp;",
    ]
    lines.extend(f"poly f{i} = {polynomial};" for i, polynomial in enumerate(rows, 1))
    lines.append("ideal I = " + ",".join(f"f{i}" for i in range(1, len(rows) + 1)) + ";")
    lines.extend(
        [
            'print("INPUT_SIZE");',
            "print(size(I));",
            'print("SLIMGB_START");',
            "ideal G = slimgb(I);",
            'print("SLIMGB_COMPLETE");',
            'print("BASIS_SIZE");',
            "print(size(G));",
            'print("REDUCE_ONE");',
            "print(reduce(1,G));",
            'print("STD_START");',
            "ideal H = std(I);",
            'print("STD_COMPLETE");',
            'print("STD_BASIS_SIZE");',
            "print(size(H));",
            'print("STD_REDUCE_ONE");',
            "print(reduce(1,H));",
            "quit;",
        ]
    )
    OUTPUT.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\r\n")
    print(
        f"wrote={OUTPUT} equations={len(rows)} variables={len(GENERATOR_INDICES)} "
        f"constant_equations={len(constant_equations)} "
        f"sha256={hashlib.sha256(OUTPUT.read_bytes()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
