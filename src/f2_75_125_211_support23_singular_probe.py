#!/usr/bin/env python3
"""Prepare a modular Singular standard-basis probe for the 35 support-23 rows.

The input equations are exact over Q. This script reduces their coefficients
modulo a prime and writes a reproducible Singular input using the original
weighted grading. The modular result is reconnaissance only.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_depth14_global_compatibility_reduction_20261021.json"
OUTPUT = HERE / "f2_75_125_211_support23_singular_probe_20261021.sing"
LOG = HERE / "f2_75_125_211_support23_singular_probe_20261021.log"
RECEIPT = HERE / "f2_75_125_211_support23_singular_probe_20261021.json"
PRIME = 32749
ROWS = set(range(86, 121))
VARIABLE_INDICES = [*range(5), *range(6, 13), *range(14, 25)]
WEIGHTS = [1, 2, 3, 3, 3, *([4] * 7), *([5] * 6), *([6] * 5)]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def mod_fraction(value: str) -> int:
    q = Fraction(value)
    denominator = q.denominator % PRIME
    if denominator == 0:
        raise ValueError(f"prime {PRIME} divides a coefficient denominator")
    return (q.numerator % PRIME) * pow(denominator, -1, PRIME) % PRIME


def main() -> None:
    data = json.loads(SOURCE.read_text(encoding="utf-8"))
    equations = [
        row for row in data["global_compatibility_equations"]
        if row["constraint_row_zero_based"] in ROWS
    ]
    expected_support = [f"p{i}" for i in VARIABLE_INDICES]
    if len(equations) != 35 or WEIGHTS != [1, 2, 3, 3, 3, *([4] * 7), *([5] * 6), *([6] * 5)]:
        raise AssertionError("unexpected 35-row input or weighted variable census")

    variable_names = [f"p{i}" for i in VARIABLE_INDICES]
    variable_position = {index: position for position, index in enumerate(VARIABLE_INDICES)}
    rendered = []
    total_terms = 0
    for row in equations:
        if row["support_variables"] != expected_support:
            raise AssertionError(f"support mismatch at {row['constraint_tag']}")
        terms = []
        for sparse, coefficient in row["terms"]:
            coeff = mod_fraction(coefficient)
            if not coeff:
                continue
            monomial_parts = []
            weighted_degree = 0
            for variable_index, exponent in sparse:
                if variable_index not in variable_position:
                    raise AssertionError(f"unexpected variable index {variable_index}")
                pos = variable_position[variable_index]
                weighted_degree += WEIGHTS[pos] * exponent
                name = variable_names[pos]
                monomial_parts.append(name if exponent == 1 else f"{name}^{exponent}")
            if weighted_degree != 14:
                raise AssertionError(f"wrong weighted degree at {row['constraint_tag']}")
            monomial = "*".join(monomial_parts) or "1"
            terms.append(f"{coeff}*{monomial}")
        if not terms:
            raise AssertionError(f"equation vanished modulo {PRIME}: {row['constraint_tag']}")
        rendered.append("(" + "+".join(terms) + ")")
        total_terms += len(terms)

    script_lines = [
        f"ring r={PRIME},({','.join(variable_names)}),wp({','.join(map(str, WEIGHTS))});",
        "option(redSB);",
        "option(prot);",
        f"print(\"START prime={PRIME} rows={len(rendered)} terms={total_terms}\");",
        "ideal I=",
        ",\n".join(rendered) + ";",
        "timer=1;",
        "ideal G=slimgb(I);",
        "print(\"STANDARD_BASIS_DONE\");",
        "print(\"BASIS_SIZE\"); print(size(G));",
        "print(\"DIMENSION\"); print(dim(G));",
        "if(dim(G)==0)",
        "{",
        "  print(\"NILPOTENCE_START\");",
        "  for(int i=1;i<=nvars(r);i++)",
        "  {",
        "    int found=0;",
        "    for(int e=1;e<=32;e++)",
        "    {",
        "      poly target=var(i)^e;",
        "      if(reduce(target,G)==0) { print(string(var(i))+\"^\"+string(e)); found=1; break; }",
        "    }",
        "    if(found==0) { print(string(var(i))+\"^not_found_through_32\"); }",
        "  }",
        "}",
        "quit;",
    ]
    OUTPUT.write_text("\n".join(script_lines) + "\n", encoding="utf-8")
    receipt = {
        "source": SOURCE.name,
        "source_sha256": sha256_file(SOURCE),
        "branch": "2+1+1",
        "scope": "35 support-23 exact compatibility equations reduced modulo a prime; diagnostic only",
        "prime": PRIME,
        "variables": variable_names,
        "weights": WEIGHTS,
        "selected_source_rows": sorted(ROWS),
        "equation_count": len(rendered),
        "term_count_after_modular_reduction": total_terms,
        "weighted_degree": 14,
        "ring_order": f"wp({','.join(map(str, WEIGHTS))})",
        "singular_input": OUTPUT.name,
        "singular_log": LOG.name,
        "interpretation": "A modular standard-basis result is reconnaissance, not a characteristic-zero or branch-exclusion proof.",
    }
    RECEIPT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
