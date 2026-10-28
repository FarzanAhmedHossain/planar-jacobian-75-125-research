"""Factor the ten smallest distance-four equations after both F3 graphs."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import sympy as sp

from f2_75_125_quadratic_affine_profile_stream import (
    iter_json_array,
    top_metadata_and_constraints,
)


TARGET_TAGS = {
    *(f"terminal-alpha-P-distance4:u^{i}" for i in range(5)),
    *(f"terminal-alpha-Q-distance4:u^{i}" for i in range(4, 9)),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def decode_coefficient(values):
    root = sp.sqrt(5)
    return sp.Add(*(sp.Rational(value) * root ** (len(values) - 1 - i)
                    for i, value in enumerate(values)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("profile", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--groebner", action="store_true")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    checkpoint_hash = sha256(args.checkpoint)
    profile = json.loads(args.profile.read_text(encoding="utf-8"))
    if profile.get("source_sha256") != checkpoint_hash:
        raise ValueError("profile and checkpoint hashes disagree")
    profile_tags = {row["tag"] for row in profile["constraint_profiles"]}
    if not TARGET_TAGS <= profile_tags:
        raise ValueError("not all expected distance-four rows are in the profile")

    text = args.checkpoint.read_text(encoding="utf-8")
    metadata, constraint_pos = top_metadata_and_constraints(text)
    names = metadata["generators"]
    parameter_indices = {
        names.index(f"p{i}"): i for i in range(int(metadata["nparameters"]))
    }
    extracted = {}
    for tag, saved in iter_json_array(text, constraint_pos):
        if tag not in TARGET_TAGS:
            continue
        support = sorted({parameter_indices[index]
                          for sparse, _values in saved
                          for index, power in sparse
                          if power and index in parameter_indices})
        if len(support) != 10 or len(saved) != 11:
            raise ValueError(f"unexpected support or term count in {tag}")
        symbols = tuple(sp.Symbol(f"p{i}") for i in support)
        local = {index: i for i, index in enumerate(support)}
        terms = {}
        for sparse, values in saved:
            monomial = [0] * len(support)
            for index, power in sparse:
                if index == 0 and power:
                    raise ValueError(f"unexpected t power in {tag}")
                if index in parameter_indices and power:
                    monomial[local[parameter_indices[index]]] = power
            terms[tuple(monomial)] = decode_coefficient(values)
        polynomial = sp.Poly.from_dict(
            terms, *symbols, domain=sp.QQ.algebraic_field(sp.sqrt(5))
        )
        content, factors = sp.factor_list(polynomial.as_expr(), extension=sp.sqrt(5))
        extracted[tag] = {
            "support": [f"p{i}" for i in support],
            "term_count": len(polynomial.terms()),
            "degree": polynomial.total_degree(),
            "content": str(content),
            "factors": [{"factor": str(factor), "multiplicity": multiplicity}
                        for factor, multiplicity in factors],
            "polynomial": polynomial,
            "symbols": symbols,
        }

    if set(extracted) != TARGET_TAGS:
        raise ValueError("the selected distance-four row set was incomplete")
    factor_multiplicity = Counter()
    for row in extracted.values():
        for factor in row["factors"]:
            expr = sp.sympify(factor["factor"], locals={
                f"p{i}": sp.Symbol(f"p{i}") for i in range(64)
            })
            normalized = sp.Poly(expr, extension=sp.sqrt(5)).monic().as_expr()
            factor_multiplicity[str(normalized)] += 1

    result = {
        "checkpoint": str(args.checkpoint),
        "checkpoint_sha256": checkpoint_hash,
        "profile": str(args.profile),
        "component": "F3=0",
        "selected_row_count": len(extracted),
        "rows": {tag: {key: value for key, value in row.items()
                       if key not in ("polynomial", "symbols")}
                 for tag, row in sorted(extracted.items())},
        "common_monic_factors_by_row_count": dict(factor_multiplicity),
        "scope": "factorization of ten selected necessary equations; not an exclusion",
    }
    if args.groebner:
        groups = {}
        for row in extracted.values():
            key = tuple(row["support"])
            groups.setdefault(key, []).append(row)
        result["groebner"] = []
        for support, rows in groups.items():
            symbols = rows[0]["symbols"]
            basis = sp.groebner([row["polynomial"].as_expr() for row in rows],
                                *symbols, extension=sp.sqrt(5))
            result["groebner"].append({
                "support": list(support),
                "row_count": len(rows),
                "basis": [str(poly) for poly in basis],
                "unit": basis.contains(sp.S.One),
            })
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "component", "selected_row_count", "common_monic_factors_by_row_count",
        "groebner", "scope",
    ) if key in result}, indent=2))


if __name__ == "__main__":
    main()
