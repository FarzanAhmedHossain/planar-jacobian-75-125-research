#!/usr/bin/env python3
"""Exact modular rank witness for the seven new depth-14 parameters.

The 141 retained equations are affine-linear in p44,...,p50.  This checker
evaluates their coefficient matrix and forcing column at one rational point
of the 42 old parameters.  A nonzero 8x8 minor of [A|b] proves that the
generic old-parameter point is incompatible with all seven new parameters.
The vanishing locus of that minor remains open, so this is a branch split,
not an exclusion.

Requires ijson for streaming the 204 MB exact checkpoint.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

import ijson


SOURCE = Path("research/f2_75_125_211_p43zero_depth14_work_flint_20261021.json")
OUT = Path("research/f2_75_125_211_depth14_newparam_rank_witness_20261021.json")
PRIME = 32749
NEW_NAMES = tuple(f"p{i}" for i in range(44, 51))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def coefficient_mod(text: str) -> int:
    value = Fraction(text)
    denominator = value.denominator % PRIME
    if denominator == 0:
        raise ValueError("chosen prime divides a checkpoint denominator")
    return (value.numerator % PRIME) * pow(denominator, -1, PRIME) % PRIME


def rank_profile(matrix: list[list[int]]) -> tuple[int, list[int], list[int]]:
    a = [row[:] for row in matrix]
    row_ids = list(range(len(a)))
    ncols = len(a[0]) if a else 0
    rank = 0
    pivot_rows: list[int] = []
    pivot_cols: list[int] = []
    for col in range(ncols):
        pivot = next((r for r in range(rank, len(a)) if a[r][col] % PRIME), None)
        if pivot is None:
            continue
        if pivot != rank:
            a[rank], a[pivot] = a[pivot], a[rank]
            row_ids[rank], row_ids[pivot] = row_ids[pivot], row_ids[rank]
        pivot_value = a[rank][col] % PRIME
        inverse = pow(pivot_value, -1, PRIME)
        for c in range(col, ncols):
            a[rank][c] = a[rank][c] * inverse % PRIME
        for r in range(rank + 1, len(a)):
            factor = a[r][col] % PRIME
            if factor:
                for c in range(col, ncols):
                    a[r][c] = (a[r][c] - factor * a[rank][c]) % PRIME
        pivot_rows.append(row_ids[rank])
        pivot_cols.append(col)
        rank += 1
        if rank == len(a):
            break
    return rank, pivot_rows, pivot_cols


def determinant_mod(matrix: list[list[int]]) -> int:
    a = [[value % PRIME for value in row] for row in matrix]
    n = len(a)
    if any(len(row) != n for row in a):
        raise ValueError("selected witness minor is not square")
    determinant = 1
    for col in range(n):
        pivot = next((r for r in range(col, n) if a[r][col]), None)
        if pivot is None:
            return 0
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            determinant = -determinant
        value = a[col][col] % PRIME
        determinant = determinant * value % PRIME
        inverse = pow(value, -1, PRIME)
        for row in range(col + 1, n):
            factor = a[row][col] * inverse % PRIME
            for j in range(col + 1, n):
                a[row][j] = (a[row][j] - factor * a[col][j]) % PRIME
            a[row][col] = 0
    return determinant % PRIME


def main() -> None:
    source_hash = sha256_file(SOURCE)
    with SOURCE.open("rb") as stream:
        names = list(ijson.items(stream, "generators.item"))
    name_to_index = {name: i for i, name in enumerate(names)}
    new_indices = {name_to_index[name]: j for j, name in enumerate(NEW_NAMES)}
    if len(new_indices) != 7:
        raise AssertionError("checkpoint does not contain p44 through p50")

    # Only the 42 old parameters that actually occur in the retained
    # constraints are specialized.  The checkpoint also contains many
    # unrelated ring generators, which must not be counted as old parameters.
    old_names = [*(f"p{i}" for i in range(36)), *(f"p{i}" for i in range(37, 43))]
    if any(name not in name_to_index for name in old_names):
        raise AssertionError("checkpoint is missing an expected active old parameter")
    # A deterministic hash-based point avoids the arithmetic progression used
    # in the first exploratory run, which landed on a low-rank specialization.
    old_values = {
        name: 1 + int.from_bytes(hashlib.sha256(f"depth14-witness-{name}".encode()).digest()[:8], "big") % (PRIME - 1)
        for name in old_names
    }
    if any(value == 0 for value in old_values.values()):
        raise AssertionError("specialization unexpectedly set an old parameter to zero")

    a_rows: list[list[int]] = []
    b_rows: list[int] = []
    row_tags: list[str] = []
    new_degree_max = 0
    with SOURCE.open("rb") as stream:
        for tag, terms in ijson.items(stream, "constraints.item"):
            row = [0] * 7
            forcing = 0
            for sparse, values in terms:
                coefficient = coefficient_mod(values[0])
                new_factors = [(i, exponent) for i, exponent in sparse if i in new_indices]
                new_degree = sum(exponent for _, exponent in new_factors)
                new_degree_max = max(new_degree_max, new_degree)
                if new_degree > 1:
                    raise AssertionError(f"constraint {tag} is nonlinear in p44,...,p50")
                monomial_value = 1
                for i, exponent in sparse:
                    if i in new_indices:
                        continue
                    name = names[i]
                    if name == "t":
                        continue
                    if name not in old_values:
                        raise ValueError(f"unclassified checkpoint variable {name}")
                    monomial_value = monomial_value * pow(old_values[name], exponent, PRIME) % PRIME
                value = coefficient * monomial_value % PRIME
                if new_factors:
                    index, exponent = new_factors[0]
                    if exponent != 1 or len(new_factors) != 1:
                        raise AssertionError(f"constraint {tag} has an invalid new-variable monomial")
                    row[new_indices[index]] = (row[new_indices[index]] + value) % PRIME
                else:
                    forcing = (forcing + value) % PRIME
            a_rows.append(row)
            b_rows.append(forcing)
            row_tags.append(tag)

    augmented = [row + [forcing] for row, forcing in zip(a_rows, b_rows)]
    rank_a, rows_a, cols_a = rank_profile(a_rows)
    rank_aug, rows_aug, cols_aug = rank_profile(augmented)
    a_det = determinant_mod([[a_rows[r][c] for c in cols_a] for r in rows_a]) if rank_a else 0
    aug_det = determinant_mod([[augmented[r][c] for c in cols_aug] for r in rows_aug]) if rank_aug else 0
    result = {
        "source": str(SOURCE),
        "source_sha256": source_hash,
        "branch": "2+1+1",
        "constraints": len(a_rows),
        "old_parameter_count": len(old_names),
        "new_parameters": list(NEW_NAMES),
        "max_new_parameter_degree": new_degree_max,
        "prime": PRIME,
        "specialization": old_values,
        "coefficient_matrix_rank_mod_prime": rank_a,
        "coefficient_minor_rows": rows_a,
        "coefficient_minor_columns": cols_a,
        "coefficient_minor_determinant_mod_prime": a_det,
        "augmented_matrix_rank_mod_prime": rank_aug,
        "augmented_minor_rows": rows_aug,
        "augmented_minor_columns": cols_aug,
        "augmented_minor_determinant_mod_prime": aug_det,
        "generic_incompatibility_witness": rank_a == 7 and rank_aug == 8 and a_det != 0 and aug_det != 0,
        "interpretation": "Only when rank(A)=7 and rank([A|b])=8 at this specialization do the two nonzero minors prove generic incompatibility: the coefficient and augmented minors remain nonzero on a nonempty Zariski-open set, so any solutions lie in their proper closed complement. This remains a branch split, not a branch exclusion.",
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
