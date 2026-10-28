"""Stream exact variable-support incidence from a detached sparse constraint array.

This reports combinatorial support structure only; it does not claim an
algebraic exclusion.  It is intentionally low-memory so it can run alongside
the exact branch continuations.
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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def incidence(source: Path, nparameters: int) -> dict:
    text = source.read_text(encoding="utf-8")
    if text.lstrip().startswith("{"):
        metadata, constraint_pos = top_metadata_and_constraints(text)
        if metadata.get("format") != "f2-sparse-ring-v1":
            raise ValueError("unexpected sparse checkpoint format")
        saved_parameters = int(metadata["nparameters"])
        if saved_parameters != nparameters:
            raise ValueError(
                f"checkpoint has {saved_parameters} parameters, not {nparameters}"
            )
    else:
        metadata, constraint_pos = {}, 0
    size_counts: Counter[int] = Counter()
    supports_by_size: dict[int, Counter[tuple[int, ...]]] = defaultdict(Counter)
    occurrences = [0] * nparameters
    parent = list(range(nparameters))
    rank = [0] * nparameters
    constraint_count = 0
    term_count = 0
    auxiliary_generators: Counter[int] = Counter()

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        if rank[ra] == rank[rb]:
            rank[ra] += 1

    for tag, polynomial in iter_json_array(text, constraint_pos):
        constraint_count += 1
        support: set[int] = set()
        for monomial, _coefficient in polynomial:
            term_count += 1
            for generator, exponent in monomial:
                if not exponent:
                    continue
                # The detached arrays use generator 0 for t and the next
                # nparameters generators for p0,...,p_(nparameters-1).
                if 1 <= generator <= nparameters:
                    support.add(generator - 1)
                else:
                    auxiliary_generators[generator] += 1
        ordered = tuple(sorted(support))
        size_counts[len(ordered)] += 1
        supports_by_size[len(ordered)][ordered] += 1
        for variable in ordered:
            occurrences[variable] += 1
        if ordered:
            for variable in ordered[1:]:
                union(ordered[0], variable)

    components: dict[int, list[int]] = defaultdict(list)
    for variable in range(nparameters):
        components[find(variable)].append(variable)
    component_sizes = sorted((len(items) for items in components.values()), reverse=True)
    repeated_supports = {
        str(size): [
            {"variables": list(support), "constraints": count}
            for support, count in counts.most_common(8)
        ]
        for size, counts in sorted(supports_by_size.items())
        if size <= 20
    }
    return {
        "source": str(source),
        "source_sha256": sha256(source),
        "branch": "2+1+1",
        "scope": "exact variable-support incidence of sparse necessary equations; not exclusion",
        "checkpoint_metadata": {
            key: metadata.get(key)
            for key in ("format", "branch", "quadratic", "nparameters")
        } if metadata else None,
        "nparameters": nparameters,
        "constraints": constraint_count,
        "terms": term_count,
        "support_size_histogram": {str(k): v for k, v in sorted(size_counts.items())},
        "distinct_support_sets_by_size": {
            str(k): len(v) for k, v in sorted(supports_by_size.items())
        },
        "most_repeated_supports": repeated_supports,
        "variable_occurrence_counts": {f"p{i}": n for i, n in enumerate(occurrences)},
        "incidence_component_sizes": component_sizes,
        "auxiliary_generator_term_counts": dict(sorted(auxiliary_generators.items())),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("nparameters", type=int)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    result = incidence(args.source, args.nparameters)
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        key: result[key]
        for key in (
            "source_sha256", "nparameters", "constraints", "terms",
            "support_size_histogram", "distinct_support_sets_by_size",
            "most_repeated_supports", "incidence_component_sizes",
            "auxiliary_generator_term_counts",
        )
    }, indent=2))
    print(f"report={args.report}")


if __name__ == "__main__":
    main()
