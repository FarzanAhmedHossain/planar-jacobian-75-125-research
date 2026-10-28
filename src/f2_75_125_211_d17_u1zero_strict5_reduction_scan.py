#!/usr/bin/env python3
"""Push the exact p1=0 strict-five-variable consequences through d17.

On the terminal weight-49 line, p1=0 is the U1=0 slice and p0=0. The exact
five-variable F5B basis on that slice forces p16=-14*p17 and
p15=(452/5)*p17. This script verifies those consequences against the saved
basis, applies the resulting polynomial substitution to every full d17
necessary equation, and scans for a constant or univariate gcd-1 conflict.
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


EXPECTED_SOURCE_SHA256 = (
    "c205df23dc70fb51df2142c482c3145bc7fc25c6be95f28131818c9c06ca3812"
)
KILLED_NAMES = {"p0", "p1", "p79"}
REPLACED_NAMES = {"p15", "p16", "p17"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def check_basis_claim(report_path: Path) -> tuple[str, str]:
    report_hash = sha256(report_path)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("unit_ideal") is not False:
        raise ValueError("expected the exact nonunit p1=0 strict-five-variable basis")
    basis = report.get("groebner_basis", [])
    if len(basis) != 4:
        raise ValueError(f"expected four F5B basis elements, found {len(basis)}")
    p7, p15, p16, p17 = sp.symbols("p7 p15 p16 p17")
    polys = [sp.sympify(poly, locals={"p7": p7, "p15": p15, "p16": p16, "p17": p17})
             for poly in basis]
    if sp.expand(polys[0] - (p16 + 14*p17)**5) != 0:
        raise ValueError("first saved basis element does not give (p16+14*p17)^5")
    if sp.expand(polys[1].subs(p16, -14*p17)) != 0:
        raise ValueError("second saved basis element does not vanish on p16=-14*p17")
    if sp.expand(polys[2].subs(p16, -14*p17) - (5*p15 - 452*p17)**3/125) != 0:
        raise ValueError("third saved basis element gives a different p15 relation")
    if sp.expand(polys[3].subs(p16, -14*p17)) != 0:
        raise ValueError("fourth saved basis element does not vanish on p16=-14*p17")
    return report_hash, "F5B basis forces p16=-14*p17 and p15=(452/5)*p17 on its zero set"


def main() -> None:
    source = Path("research/f2_75_125_sparse_ring_lift_211_roots17_dd_20261024_d17_p79zero_snapshot.json")
    basis_path = Path("research/f2_75_125_sparse_ring_lift_211_roots17_dd_20261024_d17_strict5_p1zero_gb.json")
    report_path = Path("research/f2_75_125_sparse_ring_lift_211_roots17_dd_20261024_d17_u1zero_strict5_reduction_scan.json")
    source_hash = sha256(source)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise ValueError(f"d17 checkpoint SHA-256 mismatch: {source_hash}")
    basis_hash, basis_claim = check_basis_claim(basis_path)

    reader = JsonStream(source)
    try:
        metadata, _ = reader.top_level_constraints()
        if metadata.get("format") != "f2-sparse-ring-v1" or metadata.get("branch") != "2+1+1":
            raise ValueError("expected the rational d17 2+1+1 checkpoint")
        generators = metadata["generators"]
        nparameters = int(metadata["nparameters"])
        parameter_ids = {generators.index(f"p{i}"): f"p{i}" for i in range(nparameters)}
        killed_ids = {generators.index(name) for name in KILLED_NAMES}
        replaced_ids = {generators.index(name): name for name in REPLACED_NAMES}

        support_counts: Counter[int] = Counter()
        occurrences: Counter[str] = Counter()
        parents: dict[str, str] = {}
        ranks: Counter[str] = Counter()
        gcd_states: dict[str, dict] = defaultdict(lambda: {
            "count": 0, "gcd": None, "bezout": [], "tags": [], "equations": []
        })
        constant_constraints = []
        row_count = nonzero_rows = zero_rows = input_terms = output_terms = auxiliary_terms = 0

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
            row_count += 1
            polynomial: dict[tuple[tuple[str, int], ...], Fraction] = defaultdict(Fraction)
            has_auxiliary = False
            for sparse, values in sparse_terms:
                input_terms += 1
                if len(values) != 1:
                    raise ValueError("expected rational checkpoint coefficients")
                if any(index in killed_ids and exponent for index, exponent in sparse):
                    continue
                powers: dict[str, int] = defaultdict(int)
                coefficient = Fraction(values[0])
                for generator_id, exponent in sparse:
                    exponent = int(exponent)
                    if not exponent or generator_id in killed_ids:
                        continue
                    if generator_id in replaced_ids:
                        name = replaced_ids[generator_id]
                        powers["r"] += exponent
                        if name == "p15":
                            coefficient *= Fraction(452, 5) ** exponent
                        elif name == "p16":
                            coefficient *= Fraction(-14, 1) ** exponent
                    elif generator_id in parameter_ids:
                        powers[parameter_ids[generator_id]] += exponent
                    else:
                        has_auxiliary = True
                        auxiliary_terms += 1
                monomial = tuple(sorted((name, exponent) for name, exponent in powers.items() if exponent))
                polynomial[monomial] += coefficient

            polynomial = {m: c for m, c in polynomial.items() if c}
            if not polynomial:
                zero_rows += 1
                continue
            nonzero_rows += 1
            output_terms += len(polynomial)
            support = sorted({name for monomial in polynomial for name, _exp in monomial})
            support_counts[len(support)] += 1
            for variable in support:
                occurrences[variable] += 1
                find(variable)
            if support:
                for variable in support[1:]:
                    union(support[0], variable)

            if has_auxiliary or support:
                if has_auxiliary or len(support) != 1:
                    continue
                variable = support[0]
                coefficients: dict[int, Fraction] = defaultdict(Fraction)
                for monomial, coefficient in polynomial.items():
                    exponent = sum(power for name, power in monomial if name == variable)
                    coefficients[exponent] += coefficient
                equation = [Fraction(0)] * (max(coefficients, default=0) + 1)
                for exponent, coefficient in coefficients.items():
                    equation[exponent] += coefficient
                if any(equation):
                    update_gcd_state(gcd_states[variable], tag, equation)
                continue

            value = sum(polynomial.values(), Fraction(0))
            if value:
                constant_constraints.append({"tag": tag, "value": str(value)})

        components: dict[str, int] = Counter()
        for variable in list(parents):
            components[find(variable)] += 1
        component_sizes = sorted(components.values(), reverse=True)
        univariate = {
            name: finish_gcd_state(state)
            for name, state in gcd_states.items() if state["count"]
        }
        inconsistent_gcds = {name: value for name, value in univariate.items()
                             if value["status"] == "univariate_inconsistent"}
        empty = bool(constant_constraints or inconsistent_gcds)
        result = {
            "source": str(source),
            "source_sha256": source_hash,
            "strict5_basis_report": str(basis_path),
            "strict5_basis_report_sha256": basis_hash,
            "branch": metadata["branch"],
            "base_specialization": {"p79": "0", "p0": "0", "p1": "0"},
            "additional_radical_reduction": {
                "source": basis_claim,
                "p15": "(452/5)*r",
                "p16": "-14*r",
                "p17": "r",
            },
            "checkpoint_rows": row_count,
            "nonzero_rows_after_reduction": nonzero_rows,
            "zero_rows_after_reduction": zero_rows,
            "input_sparse_terms": input_terms,
            "output_sparse_terms": output_terms,
            "auxiliary_term_occurrences": auxiliary_terms,
            "support_size_histogram": {str(k): v for k, v in sorted(support_counts.items())},
            "incidence_component_sizes": component_sizes,
            "constant_constraints": constant_constraints,
            "univariate_gcds": univariate,
            "inconsistent_univariate_gcds": inconsistent_gcds,
            "necessary_constraint_system_empty": empty,
            "scope": (
                "This is the full d17 necessary system restricted to the exact U1=0 "
                "slice and the radical relations forced by its strict-five-variable "
                "subsystem. Empty is an exclusion of this U1=0 chart slice only; "
                "nonempty/undecided is not a branch exclusion."
            ),
        }
        report_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({key: result[key] for key in (
            "source_sha256", "strict5_basis_report_sha256", "checkpoint_rows",
            "nonzero_rows_after_reduction", "zero_rows_after_reduction",
            "output_sparse_terms", "support_size_histogram",
            "incidence_component_sizes", "constant_constraints",
            "inconsistent_univariate_gcds", "necessary_constraint_system_empty", "scope"
        )}, indent=2), flush=True)
        print(f"report={report_path}", flush=True)
    finally:
        reader.close()


if __name__ == "__main__":
    main()
