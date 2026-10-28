"""Compact exact support-incidence audit of sparse F2 checkpoints.

The audit reads sparse checkpoint constraints without converting coefficients
to symbolic expressions.  It reports connected components of the variable-
constraint incidence hypergraph and the smallest-support rows.  This is
diagnostic only: a disconnected incidence graph can expose independent ideals,
but connectivity or small support does not itself prove exclusion.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from f2_75_125_quadratic_affine_profile_stream import (
    iter_json_array,
    top_metadata_and_constraints,
)


class UnionFind:
    def __init__(self, size: int):
        self.parent = list(range(size))
        self.rank = [0] * size

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit(path: Path, max_reported_support: int) -> dict:
    text = path.read_text(encoding="utf-8")
    metadata, position = top_metadata_and_constraints(text)
    if metadata.get("format") != "f2-sparse-ring-v1":
        raise ValueError("unexpected checkpoint format")

    names = metadata["generators"]
    count = int(metadata["nparameters"])
    parameter_index = {names.index(f"p{i}"): i for i in range(count)}
    uf = UnionFind(count)
    active: set[int] = set()
    rows: list[dict] = []
    support_histogram: Counter[int] = Counter()
    zero_rows = 0
    constant_rows: list[str] = []

    for tag, polynomial in iter_json_array(text, position):
        if not polynomial:
            zero_rows += 1
            continue
        support: set[int] = set()
        for sparse, _coeff in polynomial:
            for generator, power in sparse:
                if power and generator in parameter_index:
                    support.add(parameter_index[generator])
        support_histogram[len(support)] += 1
        if not support:
            constant_rows.append(tag)
        if support:
            ordered = sorted(support)
            active.update(ordered)
            for variable in ordered[1:]:
                uf.union(ordered[0], variable)
        if len(support) <= max_reported_support:
            rows.append({"tag": tag, "support": sorted(support),
                         "term_count": len(polynomial)})

    component_vars: dict[int, list[int]] = defaultdict(list)
    for variable in sorted(active):
        component_vars[uf.find(variable)].append(variable)
    components = []
    for root, variables in sorted(component_vars.items(), key=lambda item: min(item[1])):
        vset = set(variables)
        component_rows = [row for row in rows if vset.intersection(row["support"])]
        components.append({
            "variables": variables,
            "variable_count": len(variables),
            "reported_small_rows": component_rows,
        })

    return {
        "source": str(path),
        "source_sha256": sha256_file(path),
        "branch": metadata["branch"],
        "quadratic_coefficients": bool(metadata.get("quadratic")),
        "parameter_count": count,
        "constraint_count": sum(support_histogram.values()) + zero_rows,
        "nonzero_constraint_count": sum(support_histogram.values()),
        "zero_constraint_count": zero_rows,
        "support_size_histogram": {str(k): v for k, v in sorted(support_histogram.items())},
        "constant_constraint_tags": constant_rows,
        "active_variable_count": len(active),
        "inactive_variables": sorted(set(range(count)) - active),
        "incidence_component_count": len(components),
        "components": components,
        "scope": "Exact support graph only; not an ideal or branch exclusion certificate.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path, nargs="+")
    parser.add_argument("--max-reported-support", type=int, default=8)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    results = [audit(path, args.max_reported_support) for path in args.checkpoint]
    args.output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    for result in results:
        print(json.dumps({key: result[key] for key in (
            "source", "branch", "parameter_count", "constraint_count",
            "zero_constraint_count", "support_size_histogram",
            "constant_constraint_tags", "active_variable_count",
            "inactive_variables", "incidence_component_count",
        )}))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
