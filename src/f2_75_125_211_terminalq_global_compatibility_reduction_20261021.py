#!/usr/bin/env python3
"""Exact row-space reduction after adding all distance-14 terminal-Q rows.

The 141-row depth-14 system is taken from its already pinned source. Ten
additional raw/postsolve terminal-Q equations are streamed from the completed
distance-14 checkpoint. All row operations are over Q and are applied to the
seven new-parameter columns and the forcing polynomial together.
"""

from __future__ import annotations

import hashlib
import json
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import ijson


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "f2_75_125_211_p43zero_depth14_work_flint_20261021.json"
TERMINAL_SOURCE = ROOT / "f2_75_125_sparse_ring_lift_211_roots14_p43zero_d14_terminal_state.json"
OUTPUT = ROOT / "f2_75_125_211_terminalq_global_compatibility_reduction_20261021.json"
SOURCE_SHA256 = "04e1297cd63fd845e958d930aff0e433f4a14cedb37bc18dc3d62549cc284555"
TERMINAL_SOURCE_SHA256 = "0412d7637384badd19ff22cc7ce4813636cdabd51897ec13508e0668b7e44272"
OLD_NAMES = [*(f"p{i}" for i in range(36)), *(f"p{i}" for i in range(37, 43))]
NEW_NAMES = [f"p{i}" for i in range(44, 51)]
WEIGHT_BY_PARAMETER = {
    **{0: 1, 1: 2, 2: 3, 3: 3, 4: 3},
    **{i: 4 for i in range(5, 13)},
    **{i: 5 for i in range(13, 20)},
    **{i: 6 for i in range(20, 25)},
    **{i: 7 for i in range(25, 30)},
    **{i: 8 for i in range(30, 36)},
    **{i: 9 for i in range(37, 43)},
}
OLD_WEIGHTS = [WEIGHT_BY_PARAMETER[int(name[1:])] for name in OLD_NAMES]

sys.path.insert(0, str(ROOT / "f2_75_125_flint_deps"))
from flint import fmpq, fmpq_mpoly_ctx, __version__ as FLINT_VERSION  # noqa: E402


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def as_mpoly(ctx, terms: dict[tuple[int, ...], Fraction]):
    return ctx.from_dict({
        exponent: fmpq(coefficient.numerator, coefficient.denominator)
        for exponent, coefficient in terms.items() if coefficient
    })


def to_fraction(value) -> Fraction:
    return Fraction(str(value))


def encode_poly(poly) -> list[list[object]]:
    result = []
    for exponents, coefficient in sorted(poly.to_dict().items()):
        sparse = [[i, int(power)] for i, power in enumerate(exponents) if power]
        result.append([sparse, str(coefficient)])
    return result


def weighted_degrees(poly) -> set[int]:
    return {
        sum(weight * power for weight, power in zip(OLD_WEIGHTS, exponents))
        for exponents in poly.to_dict()
    }


def read_generators(path: Path) -> list[str]:
    with path.open("rb") as stream:
        return list(ijson.items(stream, "generators.item"))


def read_rows(path: Path, generators: list[str], *, terminal_q_only: bool = False):
    old_index = {name: i for i, name in enumerate(OLD_NAMES)}
    new_index = {name: i for i, name in enumerate(NEW_NAMES)}
    rows = []
    max_new_degree = 0
    with path.open("rb") as stream:
        for tag, terms in ijson.items(stream, "constraints.item"):
            if terminal_q_only and not tag.startswith("terminal-Q-distance14:"):
                continue
            columns: list[dict[tuple[int, ...], Fraction]] = [
                {} for _ in range(len(NEW_NAMES) + 1)
            ]
            for sparse, values in terms:
                coefficient = Fraction(values[0])
                new_factors = [
                    (generators[int(index)], int(power))
                    for index, power in sparse
                    if generators[int(index)] in new_index
                ]
                degree = sum(power for _, power in new_factors)
                max_new_degree = max(max_new_degree, degree)
                if degree > 1 or (new_factors and (len(new_factors) != 1 or new_factors[0][1] != 1)):
                    raise AssertionError(f"{tag} is nonlinear in p44,...,p50")
                column = new_index[new_factors[0][0]] if new_factors else len(NEW_NAMES)
                exponents = [0] * len(OLD_NAMES)
                for index, power in sparse:
                    name = generators[int(index)]
                    if name in old_index:
                        exponents[old_index[name]] += int(power)
                    elif name in new_index:
                        continue
                    elif name == "t":
                        raise AssertionError(f"unexpected active t in {tag}")
                    else:
                        raise AssertionError(f"unclassified generator {name} in {tag}")
                monomial = tuple(exponents)
                columns[column][monomial] = columns[column].get(monomial, Fraction(0)) + coefficient
            rows.append((tag, columns))
    return rows, max_new_degree


def pivot_of(coefficients):
    for column, polynomial in enumerate(coefficients):
        terms = polynomial.to_dict()
        if terms:
            exponent = min(terms)
            return (column, exponent), terms[exponent]
    return None, None


def main() -> None:
    sha_base = sha256_file(SOURCE)
    sha_terminal = sha256_file(TERMINAL_SOURCE)
    if sha_base != SOURCE_SHA256 or sha_terminal != TERMINAL_SOURCE_SHA256:
        raise SystemExit(f"source hash mismatch: base={sha_base}; terminal={sha_terminal}")
    base_generators = read_generators(SOURCE)
    terminal_generators = read_generators(TERMINAL_SOURCE)
    if base_generators != terminal_generators:
        raise SystemExit("source generator arrays differ")

    base_rows, degree_base = read_rows(SOURCE, base_generators)
    terminal_rows, degree_terminal = read_rows(
        TERMINAL_SOURCE, terminal_generators, terminal_q_only=True
    )
    if len(base_rows) != 141 or len(terminal_rows) != 10:
        raise SystemExit(f"expected 141 base + 10 terminal-Q rows; got {len(base_rows)} + {len(terminal_rows)}")
    base_tags = {tag for tag, _ in base_rows}
    if base_tags.intersection(tag for tag, _ in terminal_rows):
        raise SystemExit("terminal-Q tags overlap the 141-row base system")

    ctx = fmpq_mpoly_ctx.get(tuple(OLD_NAMES), "lex")
    all_rows = base_rows + terminal_rows
    matrix = [[as_mpoly(ctx, column) for column in columns] for _, columns in all_rows]
    if max(degree_base, degree_terminal) != 1:
        raise SystemExit("unexpected maximum degree in the new parameters")

    basis = []
    compatibility = []
    redundant_rows = []
    for row_index, ((tag, _columns), polynomials) in enumerate(zip(all_rows, matrix)):
        coefficients = polynomials[:7]
        forcing = polynomials[7]
        combination = {row_index: Fraction(1)}
        for pivot, basis_coefficients, basis_forcing, basis_combination in basis:
            column, exponent = pivot
            factor = coefficients[column].to_dict().get(exponent, fmpq(0))
            if not factor:
                continue
            factor_q = to_fraction(factor)
            for j in range(7):
                coefficients[j] -= factor * basis_coefficients[j]
            forcing -= factor * basis_forcing
            for source_row, scalar in basis_combination.items():
                value = combination.get(source_row, Fraction(0)) - factor_q * scalar
                if value:
                    combination[source_row] = value
                else:
                    combination.pop(source_row, None)

        pivot, value = pivot_of(coefficients)
        if pivot is not None:
            inverse = fmpq(1) / value
            inverse_q = to_fraction(inverse)
            coefficients = [poly * inverse for poly in coefficients]
            forcing *= inverse
            combination = {i: c * inverse_q for i, c in combination.items() if c * inverse_q}
            basis.append((pivot, coefficients, forcing, combination))
            continue

        if not forcing:
            redundant_rows.append({"row": row_index, "tag": tag})
            continue
        degrees = weighted_degrees(forcing)
        if degrees != {14}:
            raise AssertionError(f"compatibility row {tag} has weighted degrees {degrees}")
        support = sorted({
            i for exponent, _coefficient in forcing.to_dict().items()
            for i, power in enumerate(exponent) if power
        })
        compatibility.append({
            "row": row_index,
            "tag": tag,
            "term_count": len(forcing.to_dict()),
            "support_variables": [OLD_NAMES[i] for i in support],
            "weighted_degree": 14,
            "forcing_polynomial": encode_poly(forcing),
            "row_combination": [[i, str(c)] for i, c in sorted(combination.items())],
        })

    content_histogram = Counter(len(item["support_variables"]) for item in compatibility)
    residual_term_histogram = Counter(item["term_count"] for item in compatibility)
    constant_residuals = [item["tag"] for item in compatibility if item["term_count"] == 1
                          and item["forcing_polynomial"][0][0] == []]
    result = {
        "source": SOURCE.name,
        "source_sha256": sha_base,
        "terminal_source": TERMINAL_SOURCE.name,
        "terminal_source_sha256": sha_terminal,
        "branch": "2+1+1",
        "coefficient_field": "Q",
        "flint_version": FLINT_VERSION,
        "base_constraint_count": len(base_rows),
        "terminal_q_constraint_count": len(terminal_rows),
        "total_constraint_count": len(all_rows),
        "old_parameters": OLD_NAMES,
        "new_parameters": NEW_NAMES,
        "max_new_parameter_degree": max(degree_base, degree_terminal),
        "constant_rational_coefficient_row_rank": len(basis),
        "basis_pivots": [
            {"row": min(combo), "pivot_column": pivot[0],
             "pivot_monomial": [[i, int(power)] for i, power in enumerate(pivot[1]) if power]}
            for pivot, _coefficients, _forcing, combo in basis
        ],
        "nonzero_compatibility_equation_count": len(compatibility),
        "redundant_zero_rows": redundant_rows,
        "compatibility_support_size_histogram": dict(sorted(content_histogram.items())),
        "compatibility_term_count_histogram": dict(sorted(residual_term_histogram.items())),
        "constant_nonzero_compatibility_equations": constant_residuals,
        "compatibility_equations": compatibility,
        "interpretation": "Every listed compatibility equation is an exact Q-linear combination of source constraints and is necessary on every specialization of the old parameters. The remaining coefficient-basis rows still form an affine system in p44,...,p50; this reduction is not by itself a branch exclusion.",
    }
    OUTPUT.write_text(json.dumps(result, separators=(",", ":")) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "source_sha256", "terminal_source_sha256", "base_constraint_count",
        "terminal_q_constraint_count", "constant_rational_coefficient_row_rank",
        "nonzero_compatibility_equation_count", "redundant_zero_rows",
        "compatibility_support_size_histogram", "compatibility_term_count_histogram",
        "constant_nonzero_compatibility_equations", "interpretation",
    )}, indent=2), flush=True)
    print(f"output={OUTPUT.name} bytes={OUTPUT.stat().st_size} sha256={sha256_file(OUTPUT)}", flush=True)


if __name__ == "__main__":
    main()
