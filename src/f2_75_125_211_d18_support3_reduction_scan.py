#!/usr/bin/env python3
"""Propagate exact p13,p14,p15 relations through the full d18 checkpoint.

The strict-support Groebner ideal yields p14=-14*p15 and
p13=(452/5)*p15 on its zero set. This script verifies those basis relations,
substitutes them into every complete checkpoint equation, and scans the result
for exact constant or univariate-gcd contradictions.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import sympy as sp

from f2_75_125_211_streaming_checkpoint_scan import (
    JsonStream,
    finish_gcd_state,
    update_gcd_state,
)


TARGET_RELATIONS = {"p13": Fraction(452, 5), "p14": Fraction(-14, 1),
                    "p15": Fraction(1, 1)}
TARGET_NAMES = tuple(TARGET_RELATIONS)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_support3_relations(path: Path, source_hash: str) -> tuple[str, str]:
    report = json.loads(path.read_text(encoding="utf-8"))
    report_hash = sha256(path)
    if report.get("source_sha256") != source_hash:
        raise ValueError("support-three Groebner report belongs to another checkpoint")
    if report.get("selected_variables") != list(TARGET_NAMES):
        raise ValueError("unexpected support-three variable set")
    if report.get("unit_ideal") is not False:
        raise ValueError("expected the nonunit support-three ideal")
    basis = report.get("groebner_basis", [])
    if len(basis) != 4:
        raise ValueError(f"expected four exact Groebner basis elements, found {len(basis)}")
    p13, p14, p15 = sp.symbols("p13 p14 p15")
    polys = [sp.sympify(poly, locals={"p13": p13, "p14": p14, "p15": p15})
             for poly in basis]
    if sp.expand(polys[0] - (p14 + 14*p15)**5) != 0:
        raise ValueError("first basis element does not give p14=-14*p15")
    if sp.expand(polys[2].subs(p14, -14*p15)
                 - (5*p13 - 452*p15)**3/125) != 0:
        raise ValueError("third basis element does not give p13=(452/5)*p15")
    return report_hash, (
        "the exact support-three basis contains (p14+14*p15)^5 and, after "
        "p14=-14*p15, (5*p13-452*p15)^3/125"
    )


def main() -> None:
    source = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_state.json")
    support3_report = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support3_gb.json")
    report_path = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support3_reduction_scan.json")
    source_hash = sha256(source)
    expected_hash = "e2e35061e883437806d6486d8f61fdc4cc6b963c49655366a6795db210d874e7"
    if source_hash != expected_hash:
        raise ValueError(f"degree-18 checkpoint SHA-256 mismatch: {source_hash}")
    support3_hash, support3_claim = verify_support3_relations(
        support3_report, source_hash
    )

    reader = JsonStream(source)
    try:
        metadata, _ = reader.top_level_constraints()
        if (metadata.get("format") != "f2-sparse-ring-v1"
                or metadata.get("branch") != "2+1+1"
                or metadata.get("quadratic")):
            raise ValueError("expected the rational d18 2+1+1 checkpoint")
        generators = metadata["generators"]
        nparameters = int(metadata["nparameters"])
        target_ids = {generators.index(name): name for name in TARGET_NAMES}
        parameter_ids = {
            generators.index(f"p{i}"): f"p{i}"
            for i in range(nparameters)
            if f"p{i}" not in TARGET_RELATIONS
        }

        support_counts: Counter[int] = Counter()
        supports_by_size: dict[int, Counter[tuple[str, ...]]] = defaultdict(Counter)
        occurrences: Counter[str] = Counter()
        parents: dict[str, str] = {}
        ranks: Counter[str] = Counter()
        gcd_states: dict[str, dict] = defaultdict(lambda: {
            "count": 0, "gcd": None, "bezout": [], "tags": [], "equations": []
        })
        constant_constraints = []
        rows_scanned = nonzero_rows = zero_rows = input_terms = output_terms = 0

        def find(name: str) -> str:
            parents.setdefault(name, name)
            if parents[name] != name:
                parents[name] = find(parents[name])
            return parents[name]

        def union(a: str, b: str) -> None:
            ra, rb = find(a), find(b)
            if ra == rb:
                return
            if ranks[ra] < ranks[rb]:
                ra, rb = rb, ra
            parents[rb] = ra
            if ranks[ra] == ranks[rb]:
                ranks[ra] += 1

        for tag, sparse_terms in reader.array_values():
            rows_scanned += 1
            polynomial: dict[tuple[str, ...], Fraction] = defaultdict(Fraction)
            has_auxiliary = False
            for sparse, values in sparse_terms:
                input_terms += 1
                if len(values) != 1:
                    raise ValueError("expected rational checkpoint coefficients")
                coefficient = Fraction(values[0])
                powers: dict[str, int] = defaultdict(int)
                for generator, exponent in sparse:
                    exponent = int(exponent)
                    if not exponent:
                        continue
                    if generator in target_ids:
                        name = target_ids[generator]
                        powers["r"] += exponent
                        coefficient *= TARGET_RELATIONS[name] ** exponent
                    elif generator in parameter_ids:
                        powers[parameter_ids[generator]] += exponent
                    else:
                        has_auxiliary = True
                monomial = tuple(sorted(
                    name for name, exponent in powers.items()
                    for _ in range(exponent)
                ))
                polynomial[monomial] += coefficient

            polynomial = {monomial: coefficient
                          for monomial, coefficient in polynomial.items()
                          if coefficient}
            if not polynomial:
                zero_rows += 1
                continue
            nonzero_rows += 1
            output_terms += len(polynomial)
            support = sorted({name for monomial in polynomial for name in monomial})
            if has_auxiliary:
                raise ValueError(f"unexpected auxiliary generator remains in {tag}")
            support_counts[len(support)] += 1
            supports_by_size[len(support)][tuple(support)] += 1
            for variable in support:
                occurrences[variable] += 1
                find(variable)
            if support:
                for variable in support[1:]:
                    union(support[0], variable)

            if not support:
                value = sum(polynomial.values(), Fraction(0))
                if value:
                    constant_constraints.append({"tag": tag, "value": str(value)})
                continue
            if len(support) != 1:
                continue
            variable = support[0]
            degree_coefficients: dict[int, Fraction] = defaultdict(Fraction)
            for monomial, coefficient in polynomial.items():
                degree_coefficients[len(monomial)] += coefficient
            equation = [Fraction(0)] * (max(degree_coefficients, default=0) + 1)
            for degree, coefficient in degree_coefficients.items():
                equation[degree] += coefficient
            if any(equation):
                update_gcd_state(gcd_states[variable], tag, equation)

        components: dict[str, int] = Counter()
        for variable in list(parents):
            components[find(variable)] += 1
        component_sizes = sorted(components.values(), reverse=True)
        repeated_supports = {
            str(size): [
                {"variables": list(support), "constraints": count}
                for support, count in counts.most_common(8)
            ]
            for size, counts in sorted(supports_by_size.items()) if size <= 20
        }
        univariate = {name: finish_gcd_state(state)
                      for name, state in gcd_states.items() if state["count"]}
        inconsistent = {name: result for name, result in univariate.items()
                        if result["status"] == "univariate_inconsistent"}
        empty = bool(constant_constraints or inconsistent)
        result = {
            "source": str(source),
            "source_sha256": source_hash,
            "support3_basis_report": str(support3_report),
            "support3_basis_report_sha256": support3_hash,
            "branch": metadata["branch"],
            "checkpoint_parameters": nparameters,
            "forced_specialization": {
                "p14": "-14*r",
                "p13": "(452/5)*r",
                "p15": "r",
            },
            "basis_relation_evidence": support3_claim,
            "checkpoint_rows_scanned": rows_scanned,
            "nonzero_rows_after_reduction": nonzero_rows,
            "zero_rows_after_reduction": zero_rows,
            "input_sparse_terms": input_terms,
            "output_sparse_terms": output_terms,
            "support_size_histogram": {
                str(size): count for size, count in sorted(support_counts.items())
            },
            "most_repeated_supports": repeated_supports,
            "incidence_component_sizes": component_sizes,
            "constant_constraints": constant_constraints,
            "univariate_gcds": univariate,
            "inconsistent_univariate_gcds": inconsistent,
            "necessary_constraint_system_empty": empty,
            "scope": (
                "Exact propagation of radical consequences from the d18 strict-support "
                "Groebner ideal through every complete saved necessary equation. Empty "
                "would exclude this U1=0 chart slice; nonempty/undecided is not a branch "
                "exclusion."
            ),
        }
        report_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({key: result[key] for key in (
            "source_sha256", "support3_basis_report_sha256",
            "checkpoint_rows_scanned", "nonzero_rows_after_reduction",
            "zero_rows_after_reduction", "output_sparse_terms",
            "support_size_histogram", "incidence_component_sizes",
            "most_repeated_supports",
            "constant_constraints", "inconsistent_univariate_gcds",
            "necessary_constraint_system_empty", "scope",
        )}, indent=2), flush=True)
        print(f"report={report_path}", flush=True)
    finally:
        reader.close()


if __name__ == "__main__":
    main()
