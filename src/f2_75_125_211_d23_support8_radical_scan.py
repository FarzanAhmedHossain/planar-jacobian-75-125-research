#!/usr/bin/env python3
"""Propagate the exact d23 support-eight radical relations through the checkpoint."""

from __future__ import annotations

import argparse
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


NAMES = ("p5", "p15", "p34", "p35", "p36", "p37", "p38", "p39")
EXPECTED_SOURCE_SHA256 = "270f0eb737463c75d3bdfbcdfb3a303a2e7fa0af2da13f63a9e2d01ec5e04b59"


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


def polynomial_multiply(left, right):
    result: dict[tuple[tuple[str, int], ...], Fraction] = defaultdict(Fraction)
    for monomial_left, coefficient_left in left.items():
        for monomial_right, coefficient_right in right.items():
            result[multiply_monomials(monomial_left, monomial_right)] += (
                coefficient_left * coefficient_right
            )
    return {monomial: coefficient for monomial, coefficient in result.items()
            if coefficient}


def polynomial_powers(polynomial, maximum: int):
    result = {0: {(): Fraction(1)}}
    for exponent in range(1, maximum + 1):
        result[exponent] = polynomial_multiply(result[exponent - 1], polynomial)
    return result


def replacement_polynomial(expression: sp.Expr, variables: tuple[sp.Symbol, ...]):
    polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
    names = tuple(str(variable) for variable in variables)
    result = {}
    for exponents, coefficient in polynomial.terms():
        monomial = tuple((name, exponent) for name, exponent in zip(names, exponents)
                         if exponent)
        result[monomial] = Fraction(int(coefficient.p), int(coefficient.q))
    return result


def verify_radical_relations(report_path: Path, source_hash: str):
    report = json.loads(report_path.read_text(encoding="utf-8"))
    report_hash = sha256(report_path)
    if report.get("source_sha256") != source_hash:
        raise ValueError("support-eight Groebner report is bound to another checkpoint")
    if report.get("selected_variables") != list(NAMES):
        raise ValueError("unexpected support-eight variables")
    if report.get("input_equations") != 143 or report.get("unit_ideal") is not False:
        raise ValueError("unexpected support-eight subsystem scope/status")
    basis = report.get("groebner_basis", [])
    if len(basis) != 3:
        raise ValueError(f"expected a three-element basis, found {len(basis)}")

    p5, p15, p34, p35, p36, p37, p38, p39 = sp.symbols("p5 p15 p34 p35 p36 p37 p38 p39")
    local = dict(zip(NAMES, (p5, p15, p34, p35, p36, p37, p38, p39)))
    parsed = [sp.sympify(poly, locals=local) for poly in basis]
    linear = (
        224725*p34 + 1639725*p35 + 6284625*p36 + 15845845*p37
        + 24712395*p38 - 163942748*p39
    )
    if sp.expand(parsed[2] - linear**3/sp.Integer(11348910401078125)) != 0:
        raise ValueError("basis does not contain the exact L^3 relation")
    p34_section = -(
        1639725*p35 + 6284625*p36 + 15845845*p37
        + 24712395*p38 - 163942748*p39
    ) / sp.Integer(224725)
    m_form = (
        375*p35 + 3750*p36 + 20625*p37 + 82650*p38
        - 3016650*p39 + 251692*p5*p15
    )
    if sp.expand(parsed[0].subs(p34, p34_section)
                 - m_form**2/sp.Integer(63348862864)) != 0:
        raise ValueError("basis does not give the exact M^2 relation on L=0")
    return report_hash, (
        "The support-eight ideal contains L^3/11348910401078125 and, on L=0, "
        "M^2/63348862864; over C this forces L=M=0."
    )


def scan(source: Path, support8_report: Path, output: Path) -> dict:
    source_hash = EXPECTED_SOURCE_SHA256
    support8_hash, radical_evidence = verify_radical_relations(
        support8_report, source_hash
    )

    p5, p15, p34, p35, p36, p37, p38, p39 = sp.symbols("p5 p15 p34 p35 p36 p37 p38 p39")
    p35_expr = -(
        3750*p36 + 20625*p37 + 82650*p38 - 3016650*p39
        + 251692*p5*p15
    ) / sp.Integer(375)
    p34_expr = sp.cancel((
        -1639725*p35_expr - 6284625*p36 - 15845845*p37
        - 24712395*p38 + 163942748*p39
    ) / sp.Integer(224725))
    downstream = (p5, p15, p36, p37, p38, p39)
    p34_base = replacement_polynomial(p34_expr, downstream)
    p35_base = replacement_polynomial(p35_expr, downstream)
    p34_powers = {0: {(): Fraction(1)}}
    p35_powers = {0: {(): Fraction(1)}}

    reader = JsonStream(source)
    try:
        metadata, _ = reader.top_level_constraints()
        if (metadata.get("format") != "f2-sparse-ring-v1"
                or metadata.get("branch") != "2+1+1"
                or metadata.get("quadratic")):
            raise ValueError("expected rational d23 2+1+1 checkpoint")
        if int(metadata.get("nparameters", -1)) != 211:
            raise ValueError("unexpected parameter count")
        generators = metadata["generators"]
        parameter_ids = {index: name for index, name in enumerate(generators)
                         if name.startswith("p") and name[1:].isdigit()}
        target_ids = {generators.index("p34"): "p34",
                      generators.index("p35"): "p35"}

        support_counts: Counter[int] = Counter()
        support_sets: Counter[tuple[str, ...]] = Counter()
        occurrences: Counter[str] = Counter()
        parents: dict[str, str] = {}
        ranks: Counter[str] = Counter()
        gcd_states: dict[str, dict] = defaultdict(lambda: {
            "count": 0, "gcd": None, "bezout": [], "tags": [], "equations": []
        })
        constant_constraints = []
        target_exact_support = (
            "p15", "p24", "p25", "p26", "p27", "p28",
            "p44", "p45", "p46", "p47", "p48", "p5",
        )
        target_exact_support_rows = []
        rows_scanned = rows_nonzero = rows_zero = 0
        input_terms = output_terms = output_sparse_entries = 0
        maximum_joint_exponent = maximum_output_row_terms = 0

        def find(name: str) -> str:
            parents.setdefault(name, name)
            if parents[name] != name:
                parents[name] = find(parents[name])
            return parents[name]

        def union(first: str, second: str) -> None:
            left, right = find(first), find(second)
            if left == right:
                return
            if ranks[left] < ranks[right]:
                left, right = right, left
            parents[right] = left
            if ranks[left] == ranks[right]:
                ranks[left] += 1

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
                    if generator in target_ids:
                        powers[target_ids[generator]] += exponent
                    elif generator in parameter_ids:
                        powers[parameter_ids[generator]] += exponent
                    else:
                        raise ValueError(f"unexpected nonparameter generator in {tag}")

                exponent34 = powers.pop("p34", 0)
                exponent35 = powers.pop("p35", 0)
                joint_exponent = exponent34 + exponent35
                maximum_joint_exponent = max(maximum_joint_exponent, joint_exponent)
                while len(p34_powers) <= exponent34:
                    p34_powers[len(p34_powers)] = polynomial_multiply(
                        p34_powers[len(p34_powers)-1], p34_base
                    )
                while len(p35_powers) <= exponent35:
                    p35_powers[len(p35_powers)] = polynomial_multiply(
                        p35_powers[len(p35_powers)-1], p35_base
                    )
                base = tuple(sorted((name, exponent) for name, exponent in powers.items()
                                    if exponent))
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
            maximum_output_row_terms = max(maximum_output_row_terms, len(polynomial))
            output_terms += len(polynomial)
            output_sparse_entries += sum(len(monomial) for monomial in polynomial)
            support = tuple(sorted({name for monomial in polynomial for name, _ in monomial}))
            support_counts[len(support)] += 1
            support_sets[support] += 1
            if support == target_exact_support:
                target_exact_support_rows.append({
                    "tag": tag,
                    "terms": [
                        [list(monomial), str(coefficient)]
                        for monomial, coefficient in sorted(polynomial.items())
                    ],
                })
            for name in support:
                occurrences[name] += 1
                find(name)
            if support:
                for name in support[1:]:
                    union(support[0], name)
            if not support:
                constant = polynomial.get((), Fraction(0))
                if constant:
                    constant_constraints.append({"tag": tag, "value": str(constant)})
                continue
            if len(support) == 1:
                name = support[0]
                degree = max(dict(monomial).get(name, 0) for monomial in polynomial)
                equation = [Fraction(0)] * (degree + 1)
                for monomial, coefficient in polynomial.items():
                    equation[dict(monomial).get(name, 0)] += coefficient
                update_gcd_state(gcd_states[name], tag, equation)
            if rows_scanned % 100 == 0:
                print(f"rows={rows_scanned} input_terms={input_terms} "
                      f"output_terms={output_terms} max_joint_exponent={maximum_joint_exponent} "
                      f"max_row_terms={maximum_output_row_terms}", flush=True)
    finally:
        reader.close()

    components: Counter[str] = Counter()
    for name in parents:
        components[find(name)] += 1
    univariate_gcds = {
        name: finish_gcd_state(state)
        for name, state in gcd_states.items() if state["count"]
    }
    empty = bool(constant_constraints) or any(
        row["status"] == "univariate_inconsistent" for row in univariate_gcds.values()
    )
    report = {
        "source": str(source),
        "source_sha256": source_hash,
        "support8_report": str(support8_report),
        "support8_report_sha256": support8_hash,
        "support8_radical_evidence": radical_evidence,
        "branch": metadata["branch"],
        "parameters": int(metadata["nparameters"]),
        "radical_specializations": {
            "L": "224725*p34+1639725*p35+6284625*p36+15845845*p37+24712395*p38-163942748*p39=0",
            "M": "375*p35+3750*p36+20625*p37+82650*p38-3016650*p39+251692*p5*p15=0",
            "p35": str(p35_expr),
            "p34": str(p34_expr),
        },
        "checkpoint_rows_scanned": rows_scanned,
        "nonzero_rows_after_reduction": rows_nonzero,
        "zero_rows_after_reduction": rows_zero,
        "input_sparse_terms": input_terms,
        "maximum_original_joint_exponent_p34_p35": maximum_joint_exponent,
        "maximum_output_terms_in_one_row": maximum_output_row_terms,
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
        "necessary_constraint_system_empty": empty,
        "scope": (
            "Exact scan of every saved d23 equation after the radical consequences "
            "L=M=0 from a complete support-only Groebner basis. Mixed monomials are "
            "retained. A contradiction excludes only this saved U1=0 support-three-line chart."
        ),
    }
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    support_rows_output = output.with_name(output.stem + "_support12_equations.json")
    support_rows_output.write_text(json.dumps({
        "source_sha256": source_hash,
        "support8_report_sha256": support8_hash,
        "variables": list(target_exact_support),
        "equation_count": len(target_exact_support_rows),
        "equations": target_exact_support_rows,
    }, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in (
        "source_sha256", "support8_report_sha256", "checkpoint_rows_scanned",
        "nonzero_rows_after_reduction", "zero_rows_after_reduction",
        "input_sparse_terms", "maximum_original_joint_exponent_p34_p35",
        "maximum_output_terms_in_one_row", "output_terms", "support_size_histogram",
        "incidence_component_sizes", "constant_constraints", "univariate_gcds",
        "necessary_constraint_system_empty", "scope",
    )}, indent=2), flush=True)
    print(f"report={output}", flush=True)
    print(f"support12_equations={support_rows_output}", flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("support8_report", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    source_hash = sha256(args.source)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise ValueError(f"degree-23 checkpoint SHA-256 mismatch: {source_hash}")
    scan(args.source, args.support8_report, args.output)


if __name__ == "__main__":
    main()
