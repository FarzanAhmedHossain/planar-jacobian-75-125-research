#!/usr/bin/env python3
"""Stream the full d17 necessary system on the exact U1=0 slice.

The weight-49 certificate gives p1=-2*p0 on the terminal slice and
U1=-(p0/3)*(t-1)^5*(t^2-3*t+3). Thus U1=0 is exactly p0=p1=0. The d17
checkpoint also has the independently certified specialization p79=0.
This scanner applies all three specializations to every saved equation and
searches for a constant contradiction or a univariate gcd-1 certificate.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

from f2_75_125_211_streaming_checkpoint_scan import (
    JsonStream,
    finish_gcd_state,
    update_gcd_state,
)


EXPECTED_SOURCE_SHA256 = (
    "c205df23dc70fb51df2142c482c3145bc7fc25c6be95f28131818c9c06ca3812"
)
KILLED_NAMES = {"p0", "p1", "p79"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    source = Path("research/f2_75_125_sparse_ring_lift_211_roots17_dd_20261024_d17_p79zero_snapshot.json")
    report_path = Path("research/f2_75_125_sparse_ring_lift_211_roots17_dd_20261024_d17_p0p1zero_scan.json")
    source_hash = sha256(source)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise ValueError(f"d17 checkpoint SHA-256 mismatch: {source_hash}")

    reader = JsonStream(source)
    try:
        metadata, _ = reader.top_level_constraints()
        if metadata.get("format") != "f2-sparse-ring-v1" or metadata.get("branch") != "2+1+1":
            raise ValueError("expected the rational d17 2+1+1 checkpoint")
        generators = metadata["generators"]
        nparameters = int(metadata["nparameters"])
        generator_to_parameter = {
            generators.index(f"p{i}"): i for i in range(nparameters)
        }
        killed_ids = {generators.index(name) for name in KILLED_NAMES}

        size_counts: Counter[int] = Counter()
        occurrences = [0] * nparameters
        parent = list(range(nparameters))
        ranks = [0] * nparameters
        gcd_states = {
            i: {"count": 0, "gcd": None, "bezout": [], "tags": [],
                "equations": []}
            for i in range(nparameters)
        }
        constant_constraints = []
        row_count = nonzero_rows = zero_rows = term_count = surviving_terms = 0
        sparse_entries = auxiliary_terms = 0

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a: int, b: int) -> None:
            ra, rb = find(a), find(b)
            if ra == rb:
                return
            if ranks[ra] < ranks[rb]:
                ra, rb = rb, ra
            parent[rb] = ra
            if ranks[ra] == ranks[rb]:
                ranks[ra] += 1

        for tag, sparse_terms in reader.array_values():
            row_count += 1
            polynomial: dict[tuple[tuple[int, int], ...], Fraction] = defaultdict(Fraction)
            has_auxiliary = False
            for sparse, values in sparse_terms:
                term_count += 1
                if len(values) != 1:
                    raise ValueError("expected rational checkpoint coefficients")
                sparse_entries += len(sparse)
                if any(generator_id in killed_ids and exponent
                       for generator_id, exponent in sparse):
                    continue
                monomial = []
                for generator_id, exponent in sparse:
                    exponent = int(exponent)
                    if not exponent or generator_id in killed_ids:
                        continue
                    parameter_index = generator_to_parameter.get(generator_id)
                    if parameter_index is None:
                        has_auxiliary = True
                        auxiliary_terms += 1
                    else:
                        monomial.append((parameter_index, exponent))
                polynomial[tuple(sorted(monomial))] += Fraction(values[0])

            polynomial = {m: c for m, c in polynomial.items() if c}
            if not polynomial:
                zero_rows += 1
                continue
            nonzero_rows += 1
            surviving_terms += len(polynomial)
            support = sorted({i for monomial in polynomial for i, _ in monomial})
            size_counts[len(support)] += 1
            for variable in support:
                occurrences[variable] += 1
            if support:
                for variable in support[1:]:
                    union(support[0], variable)

            if not has_auxiliary and not support:
                value = sum(polynomial.values(), Fraction(0))
                if value:
                    constant_constraints.append({"tag": tag, "value": str(value)})
                continue
            if has_auxiliary or len(support) != 1:
                continue
            variable = support[0]
            coefficients: dict[int, Fraction] = defaultdict(Fraction)
            for monomial, coefficient in polynomial.items():
                exponent = sum(power for index, power in monomial if index == variable)
                coefficients[exponent] += coefficient
            equation = [Fraction(0)] * (max(coefficients, default=0) + 1)
            for exponent, coefficient in coefficients.items():
                equation[exponent] += coefficient
            if any(equation):
                update_gcd_state(gcd_states[variable], tag, equation)

        components: dict[int, list[int]] = defaultdict(list)
        for variable in range(nparameters):
            components[find(variable)].append(variable)
        component_sizes = sorted((len(items) for items in components.values()), reverse=True)
        univariate = {
            f"p{i}": finish_gcd_state(state)
            for i, state in gcd_states.items() if state["count"]
        }
        inconsistent_gcds = {
            name: value for name, value in univariate.items()
            if value["status"] == "univariate_inconsistent"
        }
        empty = bool(constant_constraints or inconsistent_gcds)

        result = {
            "source": str(source),
            "source_sha256": source_hash,
            "branch": metadata["branch"],
            "nparameters": nparameters,
            "forced_specialization": {"p79": "0", "p0": "0", "p1": "0"},
            "slice_interpretation": (
                "On the terminal weight-49 line p1=-2*p0, so p0=p1=0 is exactly U1=0."
            ),
            "checkpoint_rows": row_count,
            "nonzero_rows_after_specialization": nonzero_rows,
            "zero_rows_after_specialization": zero_rows,
            "input_terms": term_count,
            "surviving_sparse_terms": surviving_terms,
            "sparse_entries_read": sparse_entries,
            "auxiliary_term_occurrences": auxiliary_terms,
            "support_size_histogram": {str(k): v for k, v in sorted(size_counts.items())},
            "incidence_component_sizes": component_sizes,
            "variable_occurrence_counts": {f"p{i}": n for i, n in enumerate(occurrences)},
            "constant_constraints": constant_constraints,
            "univariate_gcds": univariate,
            "inconsistent_univariate_gcds": inconsistent_gcds,
            "necessary_constraint_system_empty": empty,
            "scope": (
                "Exact specialization of every saved d17 equation. Empty proves the "
                "U1=0 slice of this necessary checkpoint system is inconsistent; "
                "nonempty/undecided does not produce a branch exclusion."
            ),
        }
        report_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({key: result[key] for key in (
            "source_sha256", "checkpoint_rows", "nonzero_rows_after_specialization",
            "zero_rows_after_specialization", "surviving_sparse_terms",
            "support_size_histogram", "incidence_component_sizes",
            "constant_constraints", "inconsistent_univariate_gcds",
            "necessary_constraint_system_empty", "scope"
        )}, indent=2), flush=True)
        print(f"report={report_path}", flush=True)
    finally:
        reader.close()


if __name__ == "__main__":
    main()
