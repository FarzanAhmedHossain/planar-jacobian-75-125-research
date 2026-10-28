#!/usr/bin/env python3
"""Test whether old-parameter coordinate axes lift through the full d14 block.

The 141 saved equations are affine-linear in p44,...,p50. For each axis in
the 42 old parameters, this restricts the original affine system to an exact
141-by-7 rational linear system. Weighted homogeneity makes its ranks constant
along the nonzero part of that axis. The check is an exact necessary-system
screen only, not an exclusion of the 2+1+1 branch.
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import f2_75_125_211_depth14_coordinate_axis_screen as stream  # noqa: E402

SOURCE = HERE / "f2_75_125_211_p43zero_depth14_work_flint_20261021.json"
OUTPUT = HERE / "f2_75_125_211_depth14_affine_axis_liftability.json"
OLD_NAMES = [*(f"p{i}" for i in range(36)), *(f"p{i}" for i in range(37, 43))]
NEW_NAMES = [f"p{i}" for i in range(44, 51)]
OLD_WEIGHTS_BY_INDEX = {
    **{0: 1, 1: 2, 2: 3, 3: 3, 4: 3},
    **{i: 4 for i in range(5, 13)},
    **{i: 5 for i in range(13, 20)},
    **{i: 6 for i in range(20, 25)},
    **{i: 7 for i in range(25, 30)},
    **{i: 8 for i in range(30, 36)},
    **{i: 9 for i in range(37, 43)},
}
OLD_WEIGHTS = {name: OLD_WEIGHTS_BY_INDEX[int(name[1:])] for name in OLD_NAMES}


def exact_rank(matrix: list[list[Fraction]]) -> int:
    if not matrix:
        return 0
    rows = [row[:] for row in matrix]
    nrows, ncols = len(rows), len(rows[0])
    rank = 0
    for col in range(ncols):
        pivot = next((r for r in range(rank, nrows) if rows[r][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_value = rows[rank][col]
        for r in range(rank + 1, nrows):
            if rows[r][col]:
                factor = rows[r][col] / pivot_value
                for c in range(col, ncols):
                    rows[r][c] -= factor * rows[rank][c]
        rank += 1
        if rank == nrows:
            break
    return rank


def main() -> None:
    source = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else SOURCE
    output = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else OUTPUT
    source_hash, offset = stream.source_hash_and_constraints_offset(source)
    meta = stream.metadata(source, offset)
    generator_index = {name: i for i, name in enumerate(meta["generators"])}
    old_index = {generator_index[name]: name for name in OLD_NAMES}
    new_index = {generator_index[name]: j for j, name in enumerate(NEW_NAMES)}

    rows: list[tuple[str, dict[str, list[Fraction]], dict[str, Fraction]]] = []
    matrices = {
        axis: [[Fraction(0) for _ in NEW_NAMES] for _ in range(141)]
        for axis in OLD_NAMES
    }
    forcings = {axis: [Fraction(0) for _ in range(141)] for axis in OLD_NAMES}
    coefficient_exponents = {
        axis: [[set() for _ in range(141)] for _ in NEW_NAMES]
        for axis in OLD_NAMES
    }
    forcing_exponents = {
        axis: [set() for _ in range(141)] for axis in OLD_NAMES
    }
    new_weights: dict[str, set[int]] = {name: set() for name in NEW_NAMES}
    equation_tags = []
    term_count = 0

    for row_num, item in enumerate(stream.constraint_items(source, offset)):
        if row_num >= 141 or not isinstance(item, list) or len(item) != 2:
            raise ValueError("unexpected affine-system row encoding")
        tag, encoded_terms = item
        equation_tags.append(tag)
        row_coeffs = {axis: [Fraction(0) for _ in NEW_NAMES] for axis in OLD_NAMES}
        row_forcing = {axis: Fraction(0) for axis in OLD_NAMES}
        row_coef_exps = {axis: [set() for _ in NEW_NAMES] for axis in OLD_NAMES}
        row_force_exps = {axis: set() for axis in OLD_NAMES}

        for sparse, values in encoded_terms:
            term_count += 1
            if len(values) != 1:
                raise ValueError("expected rational coefficients")
            coefficient = Fraction(values[0])
            old_support: dict[str, int] = defaultdict(int)
            new_variable = None
            new_power = 0
            for index, power in sparse:
                name = meta["generators"][index]
                if index in old_index:
                    old_support[old_index[index]] += power
                elif index in new_index:
                    if new_variable is not None:
                        raise ValueError(f"nonlinear new-parameter term at {tag}")
                    new_variable = NEW_NAMES[new_index[index]]
                    new_power = power
                else:
                    raise ValueError(f"unclassified generator {name} in {tag}")
            if new_variable is not None and new_power != 1:
                raise ValueError(f"new parameter is not affine at {tag}")

            old_weight = sum(OLD_WEIGHTS[name] * power for name, power in old_support.items())
            implied_new_weight = 14 - old_weight
            if new_variable is not None:
                if implied_new_weight <= 0:
                    raise ValueError(f"invalid weighted degree for {new_variable} in {tag}")
                new_weights[new_variable].add(implied_new_weight)
            elif old_weight != 14:
                raise ValueError(f"pure-old term in {tag} has weighted degree {old_weight}")

            if len(old_support) > 1:
                continue
            axes = OLD_NAMES if not old_support else [next(iter(old_support))]
            for axis in axes:
                old_power = old_support.get(axis, 0)
                if new_variable is None:
                    row_forcing[axis] += coefficient
                    row_force_exps[axis].add(old_power)
                else:
                    j = NEW_NAMES.index(new_variable)
                    row_coeffs[axis][j] += coefficient
                    row_coef_exps[axis][j].add(old_power)

        for axis in OLD_NAMES:
            matrices[axis][row_num] = row_coeffs[axis]
            forcings[axis][row_num] = row_forcing[axis]
            for j in range(len(NEW_NAMES)):
                coefficient_exponents[axis][j][row_num] = row_coef_exps[axis][j]
            forcing_exponents[axis][row_num] = row_force_exps[axis]
        if (row_num + 1) % 25 == 0:
            print(f"affine_axis_scan equations={row_num + 1} terms={term_count}", flush=True)

    if len(equation_tags) != 141:
        raise AssertionError(f"expected 141 affine rows, found {len(equation_tags)}")
    if any(len(weights) != 1 for weights in new_weights.values()):
        raise AssertionError(f"new parameter weights not uniquely determined: {new_weights}")

    axes = {}
    for axis in OLD_NAMES:
        for j, exponent_sets in enumerate(coefficient_exponents[axis]):
            nonempty = [powers for powers in exponent_sets if powers]
            if any(len(powers) != 1 for powers in nonempty):
                raise AssertionError(f"nonuniform coefficient scaling on {axis}, {NEW_NAMES[j]}")
            global_powers = set().union(*exponent_sets)
            if len(global_powers) > 1:
                raise AssertionError(
                    f"coefficient scaling varies across equations on {axis}, {NEW_NAMES[j]}"
                )
        nonempty_forcing = [powers for powers in forcing_exponents[axis] if powers]
        if any(len(powers) != 1 for powers in nonempty_forcing):
            raise AssertionError(f"nonuniform forcing scaling on {axis}")
        global_forcing_powers = set().union(*forcing_exponents[axis])
        if len(global_forcing_powers) > 1:
            raise AssertionError(f"forcing scaling varies across equations on {axis}")

        coefficient_rank = exact_rank(matrices[axis])
        augmented = [row + [forcing] for row, forcing in zip(matrices[axis], forcings[axis])]
        augmented_rank = exact_rank(augmented)
        compatible = coefficient_rank == augmented_rank
        axes[axis] = {
            "status": "d14_axis_lifts" if compatible else "d14_axis_incompatible",
            "coefficient_rank_Q": coefficient_rank,
            "augmented_rank_Q": augmented_rank,
            "coefficient_power_exponents": [
                sorted(set().union(*coefficient_exponents[axis][j]))
                for j in range(len(NEW_NAMES))
            ],
            "forcing_power_exponents": sorted(set().union(*forcing_exponents[axis])),
        }

    result = {
        "source": source.name,
        "source_sha256": source_hash,
        "branch": "2+1+1",
        "scope": "exact rank consistency of the original 141-row affine d14 system on each nonzero old-parameter coordinate axis; not full-branch exclusion",
        "equations": len(equation_tags),
        "sparse_terms": term_count,
        "old_variables": OLD_NAMES,
        "old_weights": [OLD_WEIGHTS[name] for name in OLD_NAMES],
        "new_variables": NEW_NAMES,
        "new_weights": {name: next(iter(weights)) for name, weights in new_weights.items()},
        "axes": axes,
        "axis_status_counts": {
            status: sum(entry["status"] == status for entry in axes.values())
            for status in sorted({entry["status"] for entry in axes.values()})
        },
    }
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "source_sha256": source_hash,
        "equations": len(equation_tags),
        "sparse_terms": term_count,
        "new_weights": result["new_weights"],
        "axis_status_counts": result["axis_status_counts"],
        "axes": {name: (entry["coefficient_rank_Q"], entry["augmented_rank_Q"], entry["status"])
                 for name, entry in axes.items()},
        "receipt": output.name,
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
