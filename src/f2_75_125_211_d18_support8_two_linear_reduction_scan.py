#!/usr/bin/env python3
"""Propagate both exact support-eight linear consequences through d18."""

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
LINE_REPORT = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support8_line_section.json")
EXPECTED_SOURCE_SHA256 = "e2e35061e883437806d6486d8f61fdc4cc6b963c49655366a6795db210d874e7"
SUPPORT8_NAMES = ("p5", "p34", "p35", "p36", "p37", "p38", "p39", "r")
L_COEFFICIENTS = {
    "p34": Fraction(224725), "p35": Fraction(1639725),
    "p36": Fraction(6284625), "p37": Fraction(15845845),
    "p38": Fraction(24712395), "p39": Fraction(-163942748),
}
M_COEFFICIENTS = {
    "p35": Fraction(375), "p36": Fraction(3750),
    "p37": Fraction(20625), "p38": Fraction(82650),
    "p39": Fraction(-3016650),
}
SUPPORT12_NAMES = (
    "p24", "p25", "p26", "p27", "p28", "p44", "p45", "p46",
    "p47", "p48", "p5", "r",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def multiply_monomials(*monomials: tuple[tuple[str, int], ...]) -> tuple[tuple[str, int], ...]:
    powers: dict[str, int] = defaultdict(int)
    for monomial in monomials:
        for name, exponent in monomial:
            powers[name] += exponent
    return tuple(sorted((name, exponent) for name, exponent in powers.items() if exponent))


def polynomial_power(
    polynomial: dict[tuple[tuple[str, int], ...], Fraction], exponent: int
) -> dict[tuple[tuple[str, int], ...], Fraction]:
    result = {(): Fraction(1)}
    for _ in range(exponent):
        next_result: dict[tuple[tuple[str, int], ...], Fraction] = defaultdict(Fraction)
        for left_monomial, left_coefficient in result.items():
            for right_monomial, right_coefficient in polynomial.items():
                next_result[multiply_monomials(left_monomial, right_monomial)] += (
                    left_coefficient * right_coefficient
                )
        result = {monomial: coefficient for monomial, coefficient in next_result.items()
                  if coefficient}
    return result


def verify_radical_relations(source_hash: str) -> tuple[str, str, str]:
    support8 = json.loads(SUPPORT8_REPORT.read_text(encoding="utf-8"))
    line = json.loads(LINE_REPORT.read_text(encoding="utf-8"))
    support8_hash = sha256(SUPPORT8_REPORT)
    line_hash = sha256(LINE_REPORT)
    if support8.get("source_sha256") != source_hash or line.get("source_sha256") != source_hash:
        raise ValueError("support-eight relation reports belong to another checkpoint")
    if line.get("support8_report_sha256") != support8_hash:
        raise ValueError("line-section report does not bind the current support-eight report")
    if support8.get("selected_variables") != list(SUPPORT8_NAMES):
        raise ValueError("unexpected support-eight variable set")
    basis_strings = support8.get("groebner_basis", [])
    if len(basis_strings) != 3 or support8.get("input_equations") != 143:
        raise ValueError("unexpected support-eight exact subsystem")

    p5, p34, p35, p36, p37, p38, p39, r = sp.symbols("p5 p34 p35 p36 p37 p38 p39 r")
    local = dict(zip(SUPPORT8_NAMES, (p5, p34, p35, p36, p37, p38, p39, r)))
    basis = [sp.sympify(expression, locals=local) for expression in basis_strings]
    linear = 224725*p34 + 1639725*p35 + 6284625*p36 + 15845845*p37 + 24712395*p38 - 163942748*p39
    if sp.expand(basis[1] - linear**3/sp.Integer(11348910401078125)) != 0:
        raise ValueError("support-eight basis does not contain the exact L^3 relation")
    p34_section = (
        -1639725*p35 - 6284625*p36 - 15845845*p37
        - 24712395*p38 + 163942748*p39
    ) / sp.Integer(224725)
    quadratic_section = sp.cancel(basis[0].subs(p34, p34_section))
    m_form = 375*p35 + 3750*p36 + 20625*p37 + 82650*p38 - 3016650*p39 + 251692*p5*r
    if sp.expand(quadratic_section - m_form**2/sp.Integer(63348862864)) != 0:
        raise ValueError("support-eight basis does not restrict to the exact M^2 relation")
    if sp.expand(basis[2].subs(p34, p34_section)) != 0:
        raise ValueError("third support-eight basis element did not vanish on L=0")
    return support8_hash, line_hash, (
        "L^3/11348910401078125=0 implies L=0; on L=0, "
        "M^2/63348862864=0 implies M=0"
    )


def replacement_polynomial(expression: sp.Expr, variables: tuple[sp.Symbol, ...]):
    polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
    result = {}
    names = tuple(str(variable) for variable in variables)
    for exponents, coefficient in polynomial.terms():
        monomial = tuple((name, exponent) for name, exponent in zip(names, exponents) if exponent)
        result[monomial] = Fraction(int(coefficient.p), int(coefficient.q))
    return result


def main() -> None:
    source_hash = sha256(SOURCE)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise ValueError(f"degree-18 checkpoint SHA-256 mismatch: {source_hash}")
    support3_hash, support3_claim = verify_support3_relations(SUPPORT3_REPORT, source_hash)
    support8_hash, line_hash, support8_claim = verify_radical_relations(source_hash)

    p5, p36, p37, p38, p39, r = sp.symbols("p5 p36 p37 p38 p39 r")
    p35_replacement = -(
        3750*p36 + 20625*p37 + 82650*p38 - 3016650*p39 + 251692*p5*r
    ) / sp.Integer(375)
    p34_replacement = sp.cancel((
        -1639725*p35_replacement - 6284625*p36 - 15845845*p37
        - 24712395*p38 + 163942748*p39
    ) / sp.Integer(224725))
    downstream_variables = (p5, p36, p37, p38, p39, r)
    p34_polynomial = replacement_polynomial(p34_replacement, downstream_variables)
    p35_polynomial = replacement_polynomial(p35_replacement, downstream_variables)
    p34_powers = {power: polynomial_power(p34_polynomial, power) for power in range(4)}
    p35_powers = {power: polynomial_power(p35_polynomial, power) for power in range(4)}

    reader = JsonStream(SOURCE)
    try:
        metadata, _ = reader.top_level_constraints()
        if metadata.get("format") != "f2-sparse-ring-v1" or metadata.get("branch") != "2+1+1" or metadata.get("quadratic"):
            raise ValueError("expected the rational d18 2+1+1 checkpoint")
        generators = metadata["generators"]
        nparameters = int(metadata["nparameters"])
        support3_ids = {generators.index(name): name for name in SUPPORT3_NAMES}
        parameter_ids = {
            generators.index(f"p{i}"): f"p{i}"
            for i in range(nparameters)
            if f"p{i}" not in SUPPORT3_NAMES
        }

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
        output_sparse_entries = 0
        max_joint_exponent = 0
        support12_equations = []

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
                        powers[parameter_ids[generator]] += exponent
                    else:
                        raise ValueError(f"unexpected nonparameter generator in {tag}")

                exponent34 = powers.pop("p34", 0)
                exponent35 = powers.pop("p35", 0)
                joint_exponent = exponent34 + exponent35
                max_joint_exponent = max(max_joint_exponent, joint_exponent)
                if joint_exponent > 3:
                    raise ValueError(f"uncertified joint exponent {joint_exponent} in {tag}")
                base = tuple(sorted((name, exponent) for name, exponent in powers.items() if exponent))
                for mon34, coeff34 in p34_powers[exponent34].items():
                    for mon35, coeff35 in p35_powers[exponent35].items():
                        monomial = multiply_monomials(base, mon34, mon35)
                        polynomial[monomial] += coefficient * coeff34 * coeff35

            polynomial = {monomial: coefficient for monomial, coefficient in polynomial.items()
                          if coefficient}
            if not polynomial:
                rows_zero += 1
                continue
            rows_nonzero += 1
            output_terms += len(polynomial)
            output_sparse_entries += sum(len(monomial) for monomial in polynomial)
            support = tuple(sorted({name for monomial in polynomial for name, _ in monomial}))
            support_counts[len(support)] += 1
            support_sets[support] += 1
            if support == tuple(sorted(SUPPORT12_NAMES)):
                support12_equations.append({
                    "tag": tag,
                    "terms": [
                        [list(monomial), str(coefficient)]
                        for monomial, coefficient in sorted(polynomial.items())
                    ],
                })
            for name in support:
                occurrences[name] += 1
            if support:
                for name in support[1:]:
                    union(support[0], name)
            if not support:
                constant_constraints.append({"tag": tag, "constant": str(polynomial.get((), 0))})
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
        univariate_gcds = {name: finish_gcd_state(state) for name, state in gcd_states.items()}
    finally:
        reader.close()

    report = {
        "source": str(SOURCE),
        "source_sha256": source_hash,
        "support3_report_sha256": support3_hash,
        "support3_relation_evidence": support3_claim,
        "support8_report_sha256": support8_hash,
        "support8_line_section_report_sha256": line_hash,
        "support8_relation_evidence": support8_claim,
        "branch": metadata["branch"],
        "checkpoint_parameters": nparameters,
        "radical_specializations": {
            "p14": "-14*r",
            "p13": "(452/5)*r",
            "p15": "r",
            "M": "375*p35+3750*p36+20625*p37+82650*p38-3016650*p39+251692*p5*r=0",
            "L": "224725*p34+1639725*p35+6284625*p36+15845845*p37+24712395*p38-163942748*p39=0",
            "p35": "-(3750*p36+20625*p37+82650*p38-3016650*p39+251692*p5*r)/375",
            "p34": str(p34_replacement),
        },
        "checkpoint_rows_scanned": rows_scanned,
        "nonzero_rows_after_reduction": rows_nonzero,
        "zero_rows_after_reduction": rows_zero,
        "input_sparse_terms": input_terms,
        "maximum_original_joint_exponent_p34_p35": max_joint_exponent,
        "output_terms": output_terms,
        "output_sparse_entries": output_sparse_entries,
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
            state["status"] == "univariate_inconsistent" for state in univariate_gcds.values()
        ),
        "scope": (
            "Exact scan after radical consequences from two complete support-only Groebner subsystems. "
            "All checkpoint equations are substituted in full; mixed monomials are retained. "
            "An empty reduced necessary system excludes only this saved U1=0 chart slice."
        ),
    }
    support12_input = {
        "source": str(SOURCE),
        "source_sha256": source_hash,
        "support3_report_sha256": support3_hash,
        "support8_report_sha256": support8_hash,
        "support8_line_section_report_sha256": line_hash,
        "branch": metadata["branch"],
        "selected_variables": list(SUPPORT12_NAMES),
        "selection_semantics": "Keep complete transformed equations whose entire support is exactly these 12 variables; no mixed terms are discarded.",
        "equation_count": len(support12_equations),
        "equations": support12_equations,
        "scope": "Exact necessary subsystem extracted after the two verified radical substitutions.",
    }
    support12_path = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_input.json")
    support12_path.write_text(json.dumps(support12_input, separators=(",", ":")) + "\n", encoding="utf-8")
    output = Path("research/f2_75_125_sparse_ring_lift_211_roots18_u1zero_support8_two_linear_reduction_scan.json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        key: report[key] for key in (
            "source_sha256", "radical_specializations", "checkpoint_rows_scanned",
            "nonzero_rows_after_reduction", "zero_rows_after_reduction",
            "input_sparse_terms", "maximum_original_joint_exponent_p34_p35",
            "output_terms", "support_size_histogram", "incidence_component_sizes",
            "constant_constraints", "necessary_constraint_system_empty", "scope",
        )
    }, indent=2), flush=True)
    print(f"report={output}", flush=True)
    print(f"support12_input={support12_path} equations={len(support12_equations)}", flush=True)


if __name__ == "__main__":
    main()
