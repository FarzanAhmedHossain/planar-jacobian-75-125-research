#!/usr/bin/env python3
"""Propagate exact support-eight radical relations through the d18 system.

The d18 support-three ideal first gives p14=-14*r, p13=452*r/5, p15=r.
After that substitution, the complete 143-equation support-eight subsystem
has a Groebner basis containing L^3, where L is a linear form in p34..p39.
Thus L=0 on every solution. This script substitutes both consequences into
every complete checkpoint equation; it never projects away mixed monomials.
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
from f2_75_125_211_d18_support3_reduction_scan import (
    TARGET_NAMES as SUPPORT3_NAMES,
    verify_support3_relations,
)


SOURCE = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_state.json")
SUPPORT3_REPORT = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support3_gb.json")
SUPPORT8_REPORT = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support8_gb.json")
EXPECTED_SOURCE_SHA256 = "e2e35061e883437806d6486d8f61fdc4cc6b963c49655366a6795db210d874e7"
SUPPORT8_NAMES = ("p5", "p34", "p35", "p36", "p37", "p38", "p39", "r")
L_NAMES = ("p34", "p35", "p36", "p37", "p38", "p39")
L_COEFFICIENTS = {
    "p34": Fraction(224725),
    "p35": Fraction(1639725),
    "p36": Fraction(6284625),
    "p37": Fraction(15845845),
    "p38": Fraction(24712395),
    "p39": Fraction(-163942748),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_support8_relation(path: Path, source_hash: str) -> tuple[str, str]:
    report = json.loads(path.read_text(encoding="utf-8"))
    report_hash = sha256(path)
    if report.get("source_sha256") != source_hash:
        raise ValueError("support-eight report belongs to another checkpoint")
    if report.get("selected_variables") != list(SUPPORT8_NAMES):
        raise ValueError("unexpected support-eight variable set")
    if report.get("input_equations") != 143 or report.get("unit_ideal") is not False:
        raise ValueError("unexpected support-eight subsystem scope/status")
    basis = report.get("groebner_basis", [])
    if len(basis) != 3:
        raise ValueError(f"expected three support-eight basis elements, found {len(basis)}")
    symbols = sp.symbols("p5 p34 p35 p36 p37 p38 p39 r")
    local = dict(zip(SUPPORT8_NAMES, symbols))
    parsed = [sp.sympify(poly, locals=local) for poly in basis]
    linear = sum(Fraction(coefficient) * local[name]
                 for name, coefficient in L_COEFFICIENTS.items())
    cubic_scale = Fraction(1, 11348910401078125)
    if sp.expand(parsed[1] - cubic_scale * linear**3) != 0:
        raise ValueError("support-eight basis does not contain the claimed cube")
    return report_hash, f"basis element 2 is L^3/{cubic_scale.denominator}; L=0 on the solution set"


def linear_power(exponent: int) -> dict[tuple[tuple[str, int], ...], Fraction]:
    """Return the exact expansion of (p34 solved from L=0)^exponent."""
    if exponent < 0 or exponent > 3:
        raise ValueError(f"unexpected p34 exponent {exponent}; audited maximum is 3")
    replacement = {
        name: -coefficient / L_COEFFICIENTS["p34"]
        for name, coefficient in L_COEFFICIENTS.items()
        if name != "p34"
    }
    result: dict[tuple[tuple[str, int], ...], Fraction] = {(): Fraction(1)}
    for _ in range(exponent):
        next_result: dict[tuple[tuple[str, int], ...], Fraction] = defaultdict(Fraction)
        for monomial, coefficient in result.items():
            powers = dict(monomial)
            for name, scalar in replacement.items():
                expanded = dict(powers)
                expanded[name] = expanded.get(name, 0) + 1
                key = tuple(sorted((key, power) for key, power in expanded.items() if power))
                next_result[key] += coefficient * scalar
        result = {monomial: coefficient for monomial, coefficient in next_result.items()
                  if coefficient}
    return result


def main() -> None:
    source_hash = sha256(SOURCE)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise ValueError(f"degree-18 checkpoint SHA-256 mismatch: {source_hash}")
    support3_hash, support3_claim = verify_support3_relations(SUPPORT3_REPORT, source_hash)
    support8_hash, support8_claim = verify_support8_relation(SUPPORT8_REPORT, source_hash)

    reader = JsonStream(SOURCE)
    try:
        metadata, _ = reader.top_level_constraints()
        if (metadata.get("format") != "f2-sparse-ring-v1"
                or metadata.get("branch") != "2+1+1"
                or metadata.get("quadratic")):
            raise ValueError("expected the rational d18 2+1+1 checkpoint")
        generators = metadata["generators"]
        nparameters = int(metadata["nparameters"])
        support3_ids = {generators.index(name): name for name in SUPPORT3_NAMES}
        parameter_ids = {
            generators.index(f"p{i}"): f"p{i}"
            for i in range(nparameters)
            if f"p{i}" not in SUPPORT3_NAMES
        }
        powers_by_p34 = {power: linear_power(power) for power in range(4)}

        support_counts: Counter[int] = Counter()
        support_sets: Counter[tuple[str, ...]] = Counter()
        occurrences: Counter[str] = Counter()
        parent: dict[str, str] = {}
        rank: Counter[str] = Counter()
        constant_constraints = []
        gcd_states: dict[str, dict] = defaultdict(lambda: {
            "count": 0, "gcd": None, "bezout": [], "tags": [], "equations": []
        })
        rows_scanned = rows_nonzero = rows_zero = input_terms = output_terms = 0
        output_sparse_terms = 0

        def find(name: str) -> str:
            parent.setdefault(name, name)
            if parent[name] != name:
                parent[name] = find(parent[name])
            return parent[name]

        def union(first: str, second: str) -> None:
            left, right = find(first), find(second)
            if left == right:
                return
            if rank[left] < rank[right]:
                left, right = right, left
            parent[right] = left
            if rank[left] == rank[right]:
                rank[left] += 1

        for tag, sparse_terms in reader.array_values():
            rows_scanned += 1
            polynomial: dict[tuple[tuple[str, int], ...], Fraction] = defaultdict(Fraction)
            for sparse, values in sparse_terms:
                input_terms += 1
                if len(values) != 1:
                    raise ValueError("expected rational checkpoint coefficients")
                coefficient = Fraction(values[0])
                powers: dict[str, int] = defaultdict(int)
                for generator, raw_exponent in sparse:
                    exponent = int(raw_exponent)
                    if not exponent:
                        continue
                    if generator in support3_ids:
                        name = support3_ids[generator]
                        coefficient *= {
                            "p13": Fraction(452, 5),
                            "p14": Fraction(-14),
                            "p15": Fraction(1),
                        }[name] ** exponent
                        powers["r"] += exponent
                    elif generator in parameter_ids:
                        name = parameter_ids[generator]
                        if name == "p34":
                            powers["p34"] += exponent
                        else:
                            powers[name] += exponent
                    else:
                        raise ValueError(f"unexpected nonparameter generator in {tag}")

                p34_exponent = powers.pop("p34", 0)
                base = tuple(sorted((name, exponent) for name, exponent in powers.items()
                                    if exponent))
                for expanded, scalar in powers_by_p34[p34_exponent].items():
                    merged = defaultdict(int)
                    for name, exponent in base:
                        merged[name] += exponent
                    for name, exponent in expanded:
                        merged[name] += exponent
                    monomial = tuple(sorted((name, exponent) for name, exponent in merged.items()
                                            if exponent))
                    polynomial[monomial] += coefficient * scalar

            polynomial = {monomial: coefficient for monomial, coefficient in polynomial.items()
                          if coefficient}
            if not polynomial:
                rows_zero += 1
                continue
            rows_nonzero += 1
            output_terms += len(polynomial)
            output_sparse_terms += sum(len(monomial) for monomial in polynomial)
            support = tuple(sorted({name for monomial in polynomial for name, _ in monomial}))
            support_counts[len(support)] += 1
            support_sets[support] += 1
            for name in support:
                occurrences[name] += 1
            if support:
                for name in support[1:]:
                    union(support[0], name)
            if not support:
                value = polynomial.get((), Fraction(0))
                constant_constraints.append({"tag": tag, "constant": str(value)})
                continue
            if len(support) == 1:
                name = support[0]
                state = gcd_states[name]
                degree = max(dict(monomial).get(name, 0) for monomial in polynomial)
                equation = [Fraction(0)] * (degree + 1)
                for monomial, coefficient in polynomial.items():
                    equation[dict(monomial).get(name, 0)] += coefficient
                update_gcd_state(state, tag, equation)

        components: Counter[str] = Counter()
        for name in parent:
            components[find(name)] += 1
        univariate_gcds = {
            name: finish_gcd_state(state)
            for name, state in gcd_states.items()
        }
    finally:
        reader.close()

    report = {
        "source": str(SOURCE),
        "source_sha256": source_hash,
        "support3_report_sha256": support3_hash,
        "support3_relation_evidence": support3_claim,
        "support8_report_sha256": support8_hash,
        "support8_relation_evidence": support8_claim,
        "branch": metadata["branch"],
        "checkpoint_parameters": nparameters,
        "radical_specializations": {
            "p14": "-14*r",
            "p13": "(452/5)*r",
            "p15": "r",
            "p34": "(-1639725*p35-6284625*p36-15845845*p37-24712395*p38+163942748*p39)/224725",
        },
        "checkpoint_rows_scanned": rows_scanned,
        "nonzero_rows_after_reduction": rows_nonzero,
        "zero_rows_after_reduction": rows_zero,
        "input_sparse_terms": input_terms,
        "output_terms": output_terms,
        "output_sparse_entries": output_sparse_terms,
        "support_size_histogram": {str(size): count for size, count in sorted(support_counts.items())},
        "most_repeated_supports": [
            {"variables": list(support), "constraints": count}
            for support, count in support_sets.most_common(20)
        ],
        "variable_occurrence_counts": dict(sorted(occurrences.items())),
        "incidence_component_sizes": sorted(components.values(), reverse=True),
        "constant_constraints": constant_constraints,
        "univariate_gcds": univariate_gcds,
        "necessary_constraint_system_empty": bool(constant_constraints) or any(
            state["status"] == "univariate_inconsistent"
            for state in univariate_gcds.values()
        ),
        "scope": (
            "Exact scan after radical consequences from two complete support-only Groebner subsystems. "
            "All checkpoint equations are substituted in full; mixed monomials are retained. "
            "An empty reduced necessary system excludes only this saved U1=0 chart slice."
        ),
    }
    output = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support8_reduction_scan.json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        key: report[key] for key in (
            "source_sha256", "radical_specializations", "checkpoint_rows_scanned",
            "nonzero_rows_after_reduction", "zero_rows_after_reduction",
            "input_sparse_terms", "output_terms", "support_size_histogram",
            "incidence_component_sizes", "constant_constraints",
            "necessary_constraint_system_empty", "scope",
        )
    }, indent=2), flush=True)
    print(f"report={output}", flush=True)


if __name__ == "__main__":
    main()
