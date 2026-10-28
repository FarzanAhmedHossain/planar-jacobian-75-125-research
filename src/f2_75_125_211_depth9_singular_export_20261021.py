#!/usr/bin/env python3
"""Export the exact 2+1+1 distance-8/9 residual ideal to Singular.

This is a necessary truncated system in 18 parameters.  The generated
Singular script computes its characteristic-zero standard basis and tests
whether saturation by the maximal ideal removes every nonzero zero.  Even if
only the origin remains, the full F2 branch still needs to be coupled to it.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_sparse_ring_lift_211_terminal_roots9_full_state.json"
OUTPUT = HERE / "f2_75_125_211_depth9_nonzero_locus_sat_20261021.sing"
VARIABLES = (
    "p1", "p7", "p15", "p16", "p17", "p24", "p25", "p26", "p27",
    "p28", "p29", "p30", "p36", "p37", "p38", "p39", "p40", "p41",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def format_fraction(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"({value.numerator}/{value.denominator})"


def format_equation(encoded_terms, generator_names, variable_positions):
    pieces: list[str] = []
    total_degrees: set[int] = set()
    for sparse, values in encoded_terms:
        if len(values) != 1:
            raise ValueError("expected rational coefficient in the 2+1+1 source")
        coefficient = Fraction(values[0])
        exponents: dict[str, int] = {}
        total_degree = 0
        for generator_index, exponent in sparse:
            name = generator_names[generator_index]
            if name not in variable_positions:
                raise ValueError(f"equation contains unexpected variable {name}")
            exponents[name] = exponent
            total_degree += exponent
        total_degrees.add(total_degree)
        monomial = "*".join(
            name if exponents[name] == 1 else f"{name}^{exponents[name]}"
            for name in VARIABLES if exponents.get(name, 0)
        )
        magnitude = abs(coefficient)
        coefficient_text = "" if magnitude == 1 and monomial else format_fraction(magnitude)
        term = coefficient_text + ("*" if coefficient_text and monomial else "") + monomial
        if not monomial:
            term = format_fraction(magnitude)
        if not pieces:
            pieces.append(("-" if coefficient < 0 else "") + term)
        else:
            pieces.append((" - " if coefficient < 0 else " + ") + term)
    if not pieces:
        raise ValueError("zero equation unexpectedly present")
    return "".join(pieces), total_degrees


def main() -> None:
    data = json.loads(SOURCE.read_text(encoding="utf-8"))
    if data["branch"] != "2+1+1" or data["quadratic"]:
        raise ValueError("unexpected source branch or coefficient field")
    generator_names = data["generators"]
    positions = {name: i for i, name in enumerate(VARIABLES)}
    expressions: list[str] = []
    max_degrees: Counter[int] = Counter()
    term_count = 0
    tags: list[str] = []
    origin_solution = True
    for tag, terms in data["constraints"]:
        expression, row_degrees = format_equation(terms, generator_names, positions)
        max_degrees[max(row_degrees)] += 1
        origin_solution = origin_solution and 0 not in row_degrees
        expressions.append(expression)
        tags.append(tag)
        term_count += len(terms)

    if (len(expressions) != 224 or term_count != 33353
            or max_degrees != Counter({8: 112, 9: 112}) or not origin_solution):
        raise AssertionError("the saved distance-8/9 system changed")

    ring_variables = ",".join(VARIABLES)
    ideal_terms = ",\n  ".join(expressions)
    irrelevant_generators = ",".join(VARIABLES)
    source_hash = sha256_file(SOURCE)
    script = f'''// Generated from exact rational source; branch 2+1+1, distances 8 and 9.
// Source SHA-256: {source_hash}
// This truncated necessary system is not a full-branch exclusion.
LIB "elim.lib";
ring r=0,({ring_variables}),dp;
option(redSB);
ideal I =
  {ideal_terms}
;
print("SOURCE_SHA256={source_hash}");
print("INPUT_GENERATORS={len(expressions)}");
print("INPUT_TERMS={term_count}");
timer=1;
ideal G=slimgb(I);
print("STDBASIS_DONE");
print("STDBASIS_SIZE="+string(size(G)));
ideal irrelevant={irrelevant_generators};
list satResult=sat(G,irrelevant);
ideal nonzeroLocus=std(satResult[1]);
print("NONZERO_LOCUS_SAT_SIZE="+string(size(nonzeroLocus)));
if (size(nonzeroLocus)==1 && nonzeroLocus[1]==1)
{{
  print("ONLY_ORIGIN");
}}
else
{{
  print("NONZERO_LOCUS_NOT_EMPTY_OR_NOT_PROVED");
  print(nonzeroLocus);
}}
quit;
'''
    OUTPUT.write_text(script, encoding="utf-8")
    print(json.dumps({
        "source": str(SOURCE),
        "source_sha256": source_hash,
        "output": str(OUTPUT),
        "output_sha256": sha256_file(OUTPUT),
        "branch": data["branch"],
        "variables": list(VARIABLES),
        "constraints": len(expressions),
        "terms": term_count,
        "max_ordinary_degrees": dict(sorted(max_degrees.items())),
        "origin_satisfies_truncation": origin_solution,
        "constraint_tags": tags,
        "scope": "exact truncated necessary system only; not an exclusion of the F2 branch",
    }, indent=2))


if __name__ == "__main__":
    main()
