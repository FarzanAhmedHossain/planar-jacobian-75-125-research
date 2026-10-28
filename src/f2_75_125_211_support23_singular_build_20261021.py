"""Build a provenance-checked modular Singular probe for support-23 equations.

The four selected equations are the exact-Q basis for the 35 support-23 rows.
The generated Singular computation is only modular reconnaissance; it does not
exclude the full branch.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

import ijson


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_depth14_global_compatibility_reduction_20261021.json"
CERTIFICATE = HERE / "f2_75_125_211_support23_exact_rank4_certificate.json"
SINGULAR = HERE / "f2_75_125_211_support23_mod32749_20261021.sing"
MANIFEST = HERE / "f2_75_125_211_support23_mod32749_20261021.json"
PRIME = 32749
ROWS = (86, 87, 88, 89)
VARIABLES = tuple([*range(5), *range(6, 13), *range(14, 25)])
NAMES = tuple(f"p{i}" for i in VARIABLES)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def reduce_mod(value: str) -> int:
    q = Fraction(value)
    denominator = q.denominator % PRIME
    if denominator == 0:
        raise ValueError("prime divides an input coefficient denominator")
    return (q.numerator % PRIME) * pow(denominator, -1, PRIME) % PRIME


def main() -> None:
    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    source_hash = sha256(SOURCE)
    certificate_hash = sha256(CERTIFICATE)
    if certificate["source_sha256"] != source_hash:
        raise AssertionError("rank-four certificate belongs to another source")
    if certificate["exact_rank"] != 4:
        raise AssertionError("expected the exact rank-four support basis")
    if tuple(row["row"] for row in certificate["selected_source_rows"]) != ROWS:
        raise AssertionError("selected rows differ from the exact rank-four basis")

    equations = {}
    with SOURCE.open("rb") as stream:
        for row in ijson.items(stream, "global_compatibility_equations.item"):
            row_index = row["constraint_row_zero_based"]
            if row_index not in ROWS:
                continue
            if row["support_variables"] != list(NAMES):
                raise AssertionError(f"unexpected support in source row {row_index}")
            terms = []
            for sparse, coefficient in row["terms"]:
                value = reduce_mod(coefficient)
                if not value:
                    continue
                monomial = []
                for variable, power in sparse:
                    if variable not in VARIABLES:
                        raise AssertionError(f"unexpected variable p{variable}")
                    monomial.append(f"p{variable}" if power == 1 else f"p{variable}^{power}")
                factors = "*".join(monomial) or "1"
                terms.append((value, factors))
            equations[row_index] = terms

    if tuple(sorted(equations)) != ROWS:
        raise AssertionError("did not extract all four basis rows")

    polynomials = []
    for row in ROWS:
        terms = []
        for value, monomial in equations[row]:
            terms.append(monomial if value == 1 else f"{value}*{monomial}")
        polynomials.append("(" + "+".join(terms) + ")")

    lines = [
        f"ring r={PRIME},({','.join(NAMES)}),dp;",
        "option(redSB);",
        "timer=1;",
        "ideal I=" + ",\n".join(polynomials) + ";",
        "ideal G=std(I);",
        'print("basis_size="+string(size(G)));',
        'print("quotient_dimension="+string(dim(G)));',
        'print("remainder_p0_25="+string(reduce(p0^25,G)));',
        'print("remainder_p14_5="+string(reduce(p14^5,G)));',
        'print("remainder_p20_4="+string(reduce(p20^4,G)));',
        'print("remainder_p21_4="+string(reduce(p21^4,G)));',
        'print("remainder_p22_4="+string(reduce(p22^4,G)));',
        'print("remainder_p23_4="+string(reduce(p23^4,G)));',
        'print("remainder_p24_4="+string(reduce(p24^4,G)));',
        "quit;",
        "",
    ]
    script = "\n".join(lines)
    SINGULAR.write_text(script, encoding="utf-8")
    manifest = {
        "scope": "modular Singular Groebner probe of four exact basis equations for the 35 support-23 compatibility rows",
        "branch": "2+1+1",
        "field": f"GF({PRIME})",
        "source": SOURCE.name,
        "source_sha256": source_hash,
        "rank4_certificate": CERTIFICATE.name,
        "rank4_certificate_sha256": certificate_hash,
        "selected_rows_zero_based": list(ROWS),
        "variables": list(NAMES),
        "nonzero_term_counts": {str(row): len(equations[row]) for row in ROWS},
        "singular_script": SINGULAR.name,
        "singular_script_sha256": hashlib.sha256(script.encode("utf-8")).hexdigest(),
        "warning": "Modular Gröbner output is reconnaissance only; exact Q replay and attachment to the branch are required for any exclusion.",
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
