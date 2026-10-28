#!/usr/bin/env python3
"""Low-memory exact support/affine audit of a 2+2 root-10 state.

The root state is one necessary F2 subsystem over Q(sqrt(5)); this audit
streams its embedded constraint array, reports small-support equations, and
checks whether the affine-linear constraints alone are inconsistent. It
does not treat a consistent profile as an exclusion.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def constraint_array(path: Path):
    decoder = json.JSONDecoder()
    chunk_size = 1 << 20
    with path.open("r", encoding="utf-8") as stream:
        buffer = ""
        pos = 0
        eof = False

        def refill() -> bool:
            nonlocal buffer, pos, eof
            if eof:
                return False
            chunk = stream.read(chunk_size)
            if not chunk:
                eof = True
                return False
            buffer = buffer[pos:] + chunk
            pos = 0
            return True

        marker = '"constraints"'
        while marker not in buffer:
            if not refill():
                raise ValueError("missing top-level constraints array")
        pos = buffer.index(marker) + len(marker)
        while True:
            while pos < len(buffer) and (buffer[pos].isspace() or buffer[pos] == ":"):
                pos += 1
            if pos < len(buffer):
                break
            if not refill():
                raise ValueError("truncated constraints key")
        if buffer[pos] != "[":
            raise ValueError("constraints value is not an array")
        pos += 1

        while True:
            while True:
                while pos < len(buffer) and (buffer[pos].isspace() or buffer[pos] == ","):
                    pos += 1
                if pos < len(buffer):
                    break
                if not refill():
                    raise ValueError("unexpected EOF in constraints array")
            if buffer[pos] == "]":
                return
            while True:
                try:
                    item, end = decoder.raw_decode(buffer, pos)
                    break
                except json.JSONDecodeError:
                    if not refill():
                        raise
            yield item
            pos = end
            if pos > chunk_size:
                buffer = buffer[pos:]
                pos = 0


def coeff_pair(values) -> tuple[Fraction, Fraction]:
    if len(values) == 1:
        return Fraction(values[0]), Fraction(0)
    if len(values) == 2:
        return Fraction(values[0]), Fraction(values[1])
    raise ValueError(f"unexpected Q(sqrt(5)) coefficient encoding {values!r}")


def add_pairs(a: tuple[Fraction, Fraction], b: tuple[Fraction, Fraction]):
    return a[0] + b[0], a[1] + b[1]


def is_zero_pair(a: tuple[Fraction, Fraction]) -> bool:
    return a[0] == 0 and a[1] == 0


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("state", type=Path)
    parser.add_argument("branch_tag", type=str)
    parser.add_argument("nparameters", type=int)
    parser.add_argument("report", type=Path)
    parser.add_argument("--force-zero-variable", action="append", default=[])
    parser.add_argument("--force-zero-certificate", type=Path)
    parser.add_argument("--max-small-support", type=int, default=3)
    args = parser.parse_args()

    source_hash = sha256(args.state)
    field = QQ.algebraic_field(sp.sqrt(5))
    variable_names = [f"p{i}" for i in range(args.nparameters)]
    forced_zero = set()
    certificate = None
    certificate_vars = set()
    if args.force_zero_certificate:
        certificate = json.loads(
            args.force_zero_certificate.read_text(encoding="utf-8")
        )
        if certificate.get("source_sha256") != source_hash:
            raise ValueError("zero-substitution certificate does not match source state")
        certificate_names = certificate.get(
            "forced_zero_parameters_after_outer_linear_block",
            certificate.get(
                "forced_zero_parameters_after_core_block",
                certificate.get(
                    "forced_zero_parameters_after_core_degree2",
                    certificate.get(
                        "forced_zero_parameters_after_third_block",
                        certificate.get(
                            "forced_zero_parameters_after_second_block",
                            certificate.get("forced_zero_parameters", []),
                        ),
                    ),
                ),
            ),
        )
        certificate_vars = {
            int(name[1:]) for name in certificate_names
            if name.startswith("p") and name[1:].isdigit()
        }
        forced_zero.update(certificate_vars)
    for name in args.force_zero_variable:
        if not name.startswith("p") or not name[1:].isdigit():
            raise ValueError(f"invalid variable name {name!r}")
        index = int(name[1:])
        if not 0 <= index < args.nparameters:
            raise ValueError(f"forced variable {name} is outside the parameter range")
        forced_zero.add(index)

    zero_witnesses = {}
    if forced_zero:
        for saved in constraint_array(args.state):
            tag, encoded = saved
            monomials = {}
            for sparse, values in encoded:
                exp_map = {int(i): int(e) for i, e in sparse if e}
                parameter_exp = {
                    i - 1: e for i, e in exp_map.items()
                    if 1 <= i <= args.nparameters
                }
                exponent = tuple(parameter_exp.get(i, 0)
                                 for i in range(args.nparameters))
                monomials[exponent] = add_pairs(
                    monomials.get(exponent, (Fraction(0), Fraction(0))),
                    coeff_pair(values),
                )
            monomials = {e: c for e, c in monomials.items() if not is_zero_pair(c)}
            if len(monomials) != 1:
                continue
            exponent, coefficient = next(iter(monomials.items()))
            support = [i for i, power in enumerate(exponent) if power]
            if len(support) == 1 and exponent[support[0]] == 1:
                zero_witnesses.setdefault(support[0], {
                    "tag": tag,
                    "coefficient": [fraction_text(coefficient[0]),
                                   fraction_text(coefficient[1])],
                })
        missing = forced_zero - set(zero_witnesses) - certificate_vars
        if missing:
            raise ValueError(
                "requested zero substitution lacks a source one-term linear witness: "
                + ", ".join(variable_names[i] for i in sorted(missing))
            )

    support_hist: Counter[int] = Counter()
    degree_span_hist: Counter[str] = Counter()
    term_hist: Counter[int] = Counter()
    support_sets: Counter[tuple[int, ...]] = Counter()
    active_variables: set[int] = set()
    affine_rows = []
    affine_tags = []
    univariate = []
    small_support = []
    constant_only = []
    zero_count = 0
    term_count = 0
    constraint_count = 0

    sqrt5 = sp.sqrt(5)
    for saved in constraint_array(args.state):
        tag, encoded = saved
        constraint_count += 1
        poly: dict[tuple[int, ...], tuple[Fraction, Fraction]] = {}
        for sparse, values in encoded:
            exp_map = {int(i): int(e) for i, e in sparse if e}
            if any(i == 0 and e for i, e in exp_map.items()):
                raise AssertionError(f"unexpected t-dependence in constraint {tag}")
            parameter_exp = {
                i - 1: e for i, e in exp_map.items()
                if 1 <= i <= args.nparameters
            }
            other = [i for i in exp_map if i > args.nparameters]
            if other:
                raise AssertionError(f"non-parameter generators in {tag}: {other[:4]}")
            if any(parameter_exp.get(i, 0) for i in forced_zero):
                term_count += 1
                continue
            active_variables.update(parameter_exp)
            exponent = tuple(parameter_exp.get(i, 0) for i in range(args.nparameters))
            coeff = coeff_pair(values)
            poly[exponent] = add_pairs(poly.get(exponent, (Fraction(0), Fraction(0))), coeff)
            term_count += 1
        poly = {e: c for e, c in poly.items() if not is_zero_pair(c)}
        if not poly:
            zero_count += 1
            continue

        support = tuple(i for i in range(args.nparameters)
                        if any(exponent[i] for exponent in poly))
        degrees = {sum(exponent) for exponent in poly}
        support_hist[len(support)] += 1
        support_sets[support] += 1
        degree_span_hist[f"{min(degrees)}:{max(degrees)}"] += 1
        term_hist[len(poly)] += 1
        if not support:
            constant_only.append({"tag": tag, "constant": [fraction_text(x) for x in next(iter(poly.values()))]})
        if len(support) == 1:
            var = support[0]
            univariate.append({
                "tag": tag,
                "variable": variable_names[var],
                "terms": [
                    {"power": exponent[var], "coefficient": [fraction_text(c[0]), fraction_text(c[1])]}
                    for exponent, c in sorted(poly.items(), key=lambda item: item[0][var])
                ],
            })
        if len(support) <= args.max_small_support and len(small_support) < 100:
            small_support.append({
                "tag": tag,
                "support": [variable_names[i] for i in support],
                "term_count": len(poly),
                "degree_range": [min(degrees), max(degrees)],
                "terms": [
                    {"exponents": [[variable_names[i], exponent[i]] for i in support if exponent[i]],
                     "coefficient": [fraction_text(c[0]), fraction_text(c[1])]}
                    for exponent, c in sorted(poly.items())
                ],
            })

        if max(degrees) <= 1:
            row = [field.zero] * (args.nparameters + 1)
            for exponent, (rational, sqrt_coefficient) in poly.items():
                expression = sp.Rational(rational.numerator, rational.denominator)
                expression += sp.Rational(sqrt_coefficient.numerator, sqrt_coefficient.denominator) * sqrt5
                if sum(exponent) == 0:
                    column = 0
                else:
                    index = next(i for i, power in enumerate(exponent) if power)
                    if exponent[index] != 1 or sum(exponent) != 1:
                        raise AssertionError("degree-one monomial has unexpected exponent")
                    column = index + 1
                row[column] += field.convert(expression)
            affine_rows.append(row)
            affine_tags.append(tag)

    if affine_rows:
        matrix = DomainMatrix(affine_rows, (len(affine_rows), args.nparameters + 1), field)
        coefficient_matrix = DomainMatrix(
            [row[1:] for row in affine_rows],
            (len(affine_rows), args.nparameters),
            field,
        )
        _, augmented_pivots = matrix.rref()
        _, coefficient_pivots = coefficient_matrix.rref()
        affine_augmented_rank = len(augmented_pivots)
        affine_coefficient_rank = len(coefficient_pivots)
    else:
        affine_augmented_rank = affine_coefficient_rank = 0

    result = {
        "source": str(args.state),
        "source_sha256": source_hash,
        "branch": args.branch_tag,
        "scope": "exact streamed support and affine-constraint profile of one 2+2 root-10 necessary system, optionally after witnessed zero substitutions; not a whole-branch exclusion",
        "coefficient_field": "Q(sqrt(5))",
        "nparameters": args.nparameters,
        "small_support_cutoff": args.max_small_support,
        "constraint_count": constraint_count,
        "zero_constraint_count": zero_count,
        "term_count": term_count,
        "parameter_support_size_histogram": dict(sorted(support_hist.items())),
        "constraint_total_degree_span_histogram": dict(sorted(degree_span_hist.items())),
        "polynomial_term_count_histogram": dict(sorted(term_hist.items())),
        "most_repeated_supports": [
            {"variables": [variable_names[i] for i in support], "equation_count": count}
            for support, count in support_sets.most_common(30)
        ],
        "active_variables": [variable_names[i] for i in sorted(active_variables)],
        "witnessed_zero_substitutions": {
            variable_names[i]: (
                zero_witnesses[i] if i in zero_witnesses else {
                    "certificate": str(args.force_zero_certificate),
                    "certificate_sha256": sha256(args.force_zero_certificate),
                }
            )
            for i in sorted(forced_zero)
        },
        "affine_linear_system": {
            "equation_count": len(affine_rows),
            "coefficient_rank": affine_coefficient_rank,
            "augmented_rank": affine_augmented_rank,
            "inconsistent": affine_augmented_rank > affine_coefficient_rank,
            "equation_tags": affine_tags,
        },
        "constant_only_constraints": constant_only,
        "univariate_constraints": univariate,
        "small_support_constraints": small_support,
        "warning": "The equations are necessary and only extend through root distance 10. Consistency or a partial reduction does not exclude the F2 branch.",
    }
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
