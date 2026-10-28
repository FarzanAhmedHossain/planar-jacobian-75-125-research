"""Linearize tiny-support checkpoint equations by their monomial features.

For a selected small-support subsystem, the equations may be linear in a
short list of monomials even when they are nonlinear in the parameters.  An
exact RREF over Q(sqrt(5)) then exposes consequences without a multivariate
Groebner computation.  This does not itself prove that a surviving branch is
empty; it is an exact necessary-equation reduction.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import f2_75_125_quadratic_affine_profile as core
from f2_75_125_quadratic_affine_profile_stream import (
    iter_json_array,
    top_metadata_and_constraints,
)


def monomial_name(support: list[int], powers: tuple[int, ...]) -> str:
    factors = [f"p{var}^{power}" if power != 1 else f"p{var}"
               for var, power in zip(support, powers) if power]
    return "*".join(factors) if factors else "1"


def audit(path: Path, max_support: int) -> dict:
    text = path.read_text(encoding="utf-8")
    metadata, position = top_metadata_and_constraints(text)
    if metadata.get("format") != "f2-sparse-ring-v1":
        raise ValueError("unexpected checkpoint format")
    names = metadata["generators"]
    nparameters = int(metadata["nparameters"])
    p_indices = {names.index(f"p{i}"): i for i in range(nparameters)}

    selected = []
    for tag, polynomial in iter_json_array(text, position):
        support = sorted({p_indices[index]
                          for sparse, _ in polynomial
                          for index, power in sparse
                          if power and index in p_indices})
        if polynomial and 0 < len(support) <= max_support:
            selected.append((tag, polynomial, support))

    groups: dict[tuple[int, ...], list[tuple[str, list]]] = defaultdict(list)
    for tag, polynomial, support in selected:
        groups[tuple(support)].append((tag, polynomial))

    report_groups = []
    for support_tuple, equations in sorted(groups.items()):
        support = list(support_tuple)
        exponent_vectors: set[tuple[int, ...]] = set()
        decoded = []
        for tag, polynomial in equations:
            terms = defaultdict(lambda: core.ZERO)
            for sparse, coefficient_data in polynomial:
                powers = [0] * len(support)
                for index, power in sparse:
                    if not power:
                        continue
                    if index == 0:
                        raise ValueError(f"formal t-variable in {tag}")
                    if index not in p_indices:
                        raise ValueError(f"auxiliary generator in {tag}")
                    variable = p_indices[index]
                    powers[support.index(variable)] = power
                exponent = tuple(powers)
                terms[exponent] = core.add(
                    terms[exponent], core.from_json(coefficient_data)
                )
            terms = {power: coefficient for power, coefficient in terms.items()
                     if not core.iszero(coefficient)}
            exponent_vectors.update(terms)
            decoded.append((tag, terms))

        columns = sorted(exponent_vectors)
        rows = []
        for tag, terms in decoded:
            rows.append([core.ZERO] + [terms.get(power, core.ZERO)
                                       for power in columns])
        rank, pivots, inconsistency = core.rank_affine(rows)
        if inconsistency is not None:
            raise AssertionError("homogeneous monomial-feature equations became inconsistent")
        report_groups.append({
            "support": support,
            "equation_count": len(equations),
            "feature_count": len(columns),
            "feature_rank": rank,
            "features": [monomial_name(support, power) for power in columns],
            "rref": [
                {
                    "pivot_feature": monomial_name(
                        support, columns[row["pivot_column"]]
                    ),
                    "coefficients": row["row"][1:],
                }
                for row in pivots
            ],
            "equation_tags": [tag for tag, _ in equations],
        })

    return {
        "source_checkpoint": str(path),
        "source_sha256": __import__("hashlib").sha256(path.read_bytes()).hexdigest(),
        "branch": metadata["branch"],
        "parameter_count": nparameters,
        "max_support": max_support,
        "small_support_groups": report_groups,
        "scope": "Exact linear consequences in monomial features; no exclusion claim.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path, nargs="+")
    parser.add_argument("--max-support", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    results = [audit(path, args.max_support) for path in args.checkpoint]
    args.output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    for result in results:
        print(f"{result['source_checkpoint']}: branch={result['branch']} ")
        for group in result["small_support_groups"]:
            print(json.dumps({
                "support": group["support"],
                "equations": group["equation_count"],
                "features": group["features"],
                "rank": group["feature_rank"],
                "pivots": [row["pivot_feature"] for row in group["rref"]],
            }))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
