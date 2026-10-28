"""Export the pinned D24 18-parameter subsystem to Singular syntax."""

from __future__ import annotations

import hashlib
from pathlib import Path


INPUT = Path(__file__).with_name("f2_75_125_d24_u1zero_window18_constraints.tsv")
OUTPUT = Path(__file__).with_name("f2_75_125_d24_u1zero_window18.sing")
EXPECTED_INPUT_SHA256 = (
    "515046d7e36dd3c1d39972d7071b174a4b8356a6d8339a1e6126b9f36763be01"
)
CHECKPOINT_SHA256 = (
    "2369f36ce13d0d4336053ba0a7a546aa4407d04d77ee544aa4e313dbc88f2c79"
)
GENERATOR_INDICES = (6, 16, 25, 26, 27, 28, 29, 35, 36, 37, 38, 39, 40, 45, 46, 47, 48, 49)


def singular_term(encoded: str) -> str:
    coefficient, monomial = encoded.split(":", 1)
    if monomial == "1":
        return f"({coefficient})"
    factors = []
    for factor in monomial.split(","):
        index_text, exponent_text = factor.split("^", 1)
        index = int(index_text)
        if index not in GENERATOR_INDICES:
            raise ValueError(f"unexpected generator index {index}")
        factors.append(f"v{index}^{int(exponent_text)}")
    return f"({coefficient})*" + "*".join(factors)


def main() -> None:
    digest = hashlib.sha256(INPUT.read_bytes()).hexdigest()
    if digest != EXPECTED_INPUT_SHA256:
        raise ValueError(f"input digest mismatch: {digest}")
    rows = []
    for line in INPUT.read_text(encoding="utf-8").splitlines():
        _tag, body = line.split("\t", 1)
        rows.append(" + ".join(singular_term(term) for term in body.split(";")))
    if len(rows) != 286:
        raise ValueError(f"expected 286 equations, found {len(rows)}")

    variables = ",".join(f"v{index}" for index in GENERATOR_INDICES)
    lines = [
        f"// D24 checkpoint SHA256: {CHECKPOINT_SHA256}",
        f"// TSV SHA256: {digest}",
        f"// Exact necessary subsystem: {len(rows)} equations in 18 parameters.",
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
    OUTPUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote={OUTPUT} equations={len(rows)} sha256={hashlib.sha256(OUTPUT.read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
