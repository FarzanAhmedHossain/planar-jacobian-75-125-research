#!/usr/bin/env python3
"""Modular Macaulay certificates for the weighted support-12 subideal.

The input is the detached depth-13 necessary-equation array. A full column
rank over a prime field, with all rational denominators invertible there,
certifies the same rank over Q. The result is local to this support-12
subideal and is never promoted to exclusion of the full F2 branch.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

from flint import nmod_mat

from f2_75_125_211_degree_structure import streamed_array


WEIGHTS = (1, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4)
PRIME = 2_305_843_009_213_693_951  # 2^61 - 1


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def monomials_of_weight(weights: tuple[int, ...], target: int) -> list[tuple[int, ...]]:
    result: list[tuple[int, ...]] = []
    current = [0] * len(weights)

    def visit(index: int, remaining: int) -> None:
        if index == len(weights):
            if remaining == 0:
                result.append(tuple(current))
            return
        weight = weights[index]
        for exponent in range(remaining // weight + 1):
            current[index] = exponent
            visit(index + 1, remaining - exponent * weight)
        current[index] = 0

    visit(0, target)
    return result


def mod_fraction(value: Fraction, prime: int) -> int:
    denominator = value.denominator % prime
    if denominator == 0:
        raise ValueError(f"prime {prime} divides a source coefficient denominator")
    return (value.numerator % prime) * pow(denominator, -1, prime) % prime


def modular_rank_and_rows(rows: list[list[int]], ncols: int, prime: int) -> tuple[int, list[int], int | None]:
    matrix = nmod_mat(rows, prime)
    rank = matrix.rank()
    if rank != ncols:
        return rank, [], None
    reduced, reduced_rank = matrix.transpose().rref()
    if reduced_rank != ncols:
        raise AssertionError("transpose rank disagrees with matrix rank")
    independent_rows: list[int] = []
    for i in range(ncols):
        pivot = next((j for j in range(len(rows)) if int(reduced[i, j]) != 0), None)
        if pivot is None:
            raise AssertionError("missing pivot in full-row-rank transpose")
        independent_rows.append(pivot)
    minor = nmod_mat(
        [[rows[row][col] for col in range(ncols)] for row in independent_rows],
        prime,
    )
    determinant = int(minor.det())
    if determinant == 0:
        raise AssertionError("reported independent rows have zero determinant")
    return rank, independent_rows, determinant


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("incidence_report", type=Path)
    parser.add_argument("nparameters", type=int)
    parser.add_argument("support_size", type=int)
    parser.add_argument("receipt", type=Path)
    args = parser.parse_args()

    source_hash = sha256(args.source)
    incidence = json.loads(args.incidence_report.read_text(encoding="utf-8"))
    if source_hash != incidence["source_sha256"]:
        raise ValueError("incidence report does not match the exact source hash")
    supports = incidence["most_repeated_supports"][str(args.support_size)]
    if len(supports) != 1:
        raise ValueError("requested support size is not a unique support set")
    support = tuple(supports[0]["variables"])
    if len(support) != len(WEIGHTS):
        raise ValueError("the verified weights apply to exactly 12 variables")

    polynomials: dict[int, list[tuple[str, dict[tuple[int, ...], Fraction]]]] = {
        8: [], 12: []
    }
    term_count = 0
    for tag, saved in streamed_array(args.source):
        term_map: defaultdict[tuple[int, ...], Fraction] = defaultdict(Fraction)
        parameter_support: set[int] = set()
        for sparse, values in saved:
            exponent_map = {
                i - 1: e for i, e in sparse
                if 1 <= i <= args.nparameters and e
            }
            parameter_support.update(exponent_map)
            exponent = tuple(exponent_map.get(i, 0) for i in support)
            if len(values) != 1:
                raise ValueError(f"non-rational coefficient encoding in {tag}")
            term_map[exponent] += Fraction(values[0])
            term_count += 1
        if tuple(sorted(parameter_support)) != support:
            continue
        term_map = {e: c for e, c in term_map.items() if c}
        degrees = {
            sum(weight * exponent for weight, exponent in zip(WEIGHTS, monomial))
            for monomial in term_map
        }
        if len(degrees) != 1:
            raise AssertionError(f"weights do not grade {tag}: {sorted(degrees)}")
        degree = next(iter(degrees))
        if degree not in polynomials:
            raise AssertionError(f"unexpected weighted degree {degree} in {tag}")
        polynomials[degree].append((tag, term_map))

    if (len(polynomials[8]), len(polynomials[12])) != (112, 35):
        raise AssertionError(
            f"unexpected selected generator counts: "
            f"degree 8={len(polynomials[8])}, degree 12={len(polynomials[12])}"
        )

    mon4 = monomials_of_weight(WEIGHTS, 4)
    mon8 = monomials_of_weight(WEIGHTS, 8)
    mon12 = monomials_of_weight(WEIGHTS, 12)
    index8 = {e: i for i, e in enumerate(mon8)}
    index12 = {e: i for i, e in enumerate(mon12)}

    def dense_row(poly: dict[tuple[int, ...], Fraction], index: dict[tuple[int, ...], int]) -> list[int]:
        row = [0] * len(index)
        for exponent, coefficient in poly.items():
            row[index[exponent]] = mod_fraction(coefficient, PRIME)
        return row

    rows8 = [dense_row(poly, index8) for _tag, poly in polynomials[8]]
    rank8, independent8, det8 = modular_rank_and_rows(rows8, len(mon8), PRIME)

    rows12: list[list[int]] = []
    for _tag, poly in polynomials[12]:
        rows12.append(dense_row(poly, index12))

    for _tag, poly in polynomials[8]:
        for multiplier in mon4:
            product: defaultdict[tuple[int, ...], int] = defaultdict(int)
            for exponent, coefficient in poly.items():
                shifted = tuple(a + b for a, b in zip(exponent, multiplier))
                product[shifted] = (product[shifted] + mod_fraction(coefficient, PRIME)) % PRIME
            row = [0] * len(mon12)
            for exponent, coefficient in product.items():
                row[index12[exponent]] = coefficient
            rows12.append(row)

    if len(rows12) != 1491:
        raise AssertionError(f"unexpected degree-12 Macaulay row count {len(rows12)}")
    rank12, independent12, det12 = modular_rank_and_rows(rows12, len(mon12), PRIME)

    power_claims = {}
    if rank12 == len(mon12):
        for variable, weight in zip((f"p{i}" for i in support), WEIGHTS):
            power_claims[variable] = {
                "power": 12 // weight,
                "weighted_degree": 12,
                "membership_consequence": "the pure power lies in the support-12 subideal over Q",
            }

    result = {
        "source": str(args.source),
        "source_sha256": source_hash,
        "incidence_report_sha256": sha256(args.incidence_report),
        "branch": "2+1+1",
        "scope": "exact modular Macaulay rank certificate for the 12-variable subideal only; not a whole-branch exclusion",
        "coefficient_field": "Q, certified by reduction modulo a prime not dividing any source denominator",
        "prime": PRIME,
        "support_variables": [f"p{i}" for i in support],
        "positive_integer_weights": list(WEIGHTS),
        "selected_generators": {str(d): len(polynomials[d]) for d in (8, 12)},
        "source_terms_scanned": term_count,
        "weighted_monomial_counts": {"4": len(mon4), "8": len(mon8), "12": len(mon12)},
        "degree_8_generator_matrix": {
            "shape": [len(rows8), len(mon8)],
            "rank_mod_prime": rank8,
            "independent_generator_rows": independent8,
            "independent_minor_determinant_mod_prime": det8,
        },
        "degree_12_macaulay_matrix": {
            "shape": [len(rows12), len(mon12)],
            "construction": "all 112 degree-8 generators times the 13 weight-4 monomials, followed by the 35 degree-12 generators",
            "rank_mod_prime": rank12,
            "independent_generator_rows": independent12,
            "independent_minor_determinant_mod_prime": det12,
        },
        "pure_power_consequences": power_claims,
        "interpretation": "If the degree-12 matrix has full column rank, the support-12 ideal contains every weight-12 form and hence the displayed pure powers; with positive weights, its only common zero is the origin.",
        "warning": "This is a necessary-equation subideal in 12 of 75 carried parameters. Even forcing these 12 parameters to zero would not exclude the full 2+1+1 branch.",
    }
    args.receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
