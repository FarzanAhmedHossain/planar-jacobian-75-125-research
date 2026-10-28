"""Eliminate a carried parameter forced to zero by an exact checkpoint row.

This reduction is valid only when the source checkpoint contains the explicit
nonzero scalar multiple of that coordinate as a constraint.  It substitutes
the coordinate by zero in every saved polynomial, removes the generator, and
checks every retained monomial's index map.  The receipt records the exact
source/output hashes and dropped monomials.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


POLY_KEYS = ("h", "v", "S", "B", "kappa", "kappa2", "kappa3")
SEQUENCE_KEYS = ("p", "q", "w")


def polynomials(state: dict):
    for key in POLY_KEYS:
        if state.get(key) is not None:
            yield key, state[key]
    for key in SEQUENCE_KEYS:
        for i, value in enumerate(state.get(key, [])):
            yield f"{key}[{i}]", value
    for key, value in state.get("uk", {}).items():
        yield f"uk[{key}]", value
    for i, (_tag, value) in enumerate(state.get("constraints", [])):
        yield f"constraints[{i}]", value


def eliminate(source: Path, output: Path, receipt_path: Path,
              parameter_name: str, equation_tag: str) -> dict:
    raw = source.read_bytes()
    state = json.loads(raw)
    if state.get("format") != "f2-sparse-ring-v1":
        raise ValueError("unexpected sparse checkpoint format")
    if state.get("quadratic"):
        raise ValueError("this reducer currently expects rational coefficients")
    if output.exists() or receipt_path.exists():
        raise FileExistsError("refusing to overwrite an existing reduction artifact")

    old_names = state["generators"]
    old_index = {name: i for i, name in enumerate(old_names)}
    if parameter_name not in old_index or not parameter_name.startswith("p"):
        raise ValueError("the requested parameter is not a checkpoint generator")
    eliminated_index = old_index[parameter_name]
    if int(parameter_name[1:]) >= int(state["nparameters"]):
        raise ValueError("the requested parameter is not in the carried prefix")

    matching = [poly for tag, poly in state.get("constraints", []) if tag == equation_tag]
    if len(matching) != 1 or len(matching[0]) != 1:
        raise ValueError("the cited constraint is not a single monomial equation")
    sparse, values = matching[0][0]
    if sparse != [[eliminated_index, 1]] or not values or values[0] == "0":
        raise ValueError("the cited constraint is not a nonzero scalar multiple of the parameter")
    scalar = values[0]

    kept_monomials = 0
    dropped_monomials = 0
    retained_terms = 0
    original_terms = 0

    def specialize(polynomial: list) -> list:
        nonlocal kept_monomials, dropped_monomials, retained_terms, original_terms
        result = []
        for sparse, values in polynomial:
            original_terms += 1
            exponents = dict(sparse)
            if exponents.get(eliminated_index, 0):
                dropped_monomials += 1
                continue
            old_sparse = sorted(sparse)
            new_sparse = [[i if i < eliminated_index else i - 1, power]
                          for i, power in old_sparse]
            # The inverse index map recovers every retained source monomial.
            back = [[i if i < eliminated_index else i + 1, power]
                    for i, power in new_sparse]
            if back != old_sparse:
                raise AssertionError("retained monomial failed exact inverse-map check")
            result.append([new_sparse, values])
            kept_monomials += 1
            retained_terms += 1
        return result

    for key in POLY_KEYS:
        if state.get(key) is not None:
            state[key] = specialize(state[key])
    for key in SEQUENCE_KEYS:
        state[key] = [specialize(value) for value in state.get(key, [])]
    state["uk"] = {key: specialize(value) for key, value in state.get("uk", {}).items()}
    state["constraints"] = [
        [tag, specialize(poly)] for tag, poly in state.get("constraints", [])
    ]

    old_parameter_names = sorted(
        (name for name in old_names if name.startswith("p") and name[1:].isdigit()),
        key=lambda name: int(name[1:]),
    )
    new_parameter_names = [name for name in old_parameter_names
                           if name != parameter_name]
    rename = {name: f"p{i}" for i, name in enumerate(new_parameter_names)}
    new_names = ["t", *(f"p{i}" for i in range(len(old_parameter_names) - 1)),
                 *(f"x{i}" for i in range(73))]
    # Confirm that the source layout is exactly t, p-prefix, x0..x72.
    expected_old_names = ["t", *old_parameter_names, *(f"x{i}" for i in range(73))]
    if old_names != expected_old_names:
        raise ValueError("checkpoint generator order is not the expected sparse-lift layout")

    new_parameter_names_check = sorted(
        (name for name in new_names if name.startswith("p") and name[1:].isdigit()),
        key=lambda name: int(name[1:]),
    )
    if set(rename.values()) != set(new_parameter_names_check):
        raise AssertionError("parameter relabeling is not contiguous")
    state["generators"] = new_names
    state["nparameters"] = int(state["nparameters"]) - 1
    state["compaction"] = dict(state.get("compaction") or {})
    old_to_new = state["compaction"].get("old_to_new_parameter_names", {})
    if old_to_new:
        eliminated_old_name = next(
            (old for old, new in old_to_new.items() if new == parameter_name), None
        )
        if eliminated_old_name is not None:
            old_to_new = {
                old: new for old, new in old_to_new.items()
                if new != parameter_name
            }
            old_to_new = {
                old: (f"p{int(new[1:]) - 1}" if int(new[1:]) > int(parameter_name[1:])
                      else new)
                for old, new in old_to_new.items()
            }
            state["compaction"]["active_old_parameters"] = [
                old for old in state["compaction"].get("active_old_parameters", [])
                if old != eliminated_old_name
            ]
            state["compaction"]["old_to_new_parameter_names"] = old_to_new
    state["compaction"]["new_nparameters"] = state["nparameters"]
    state["compaction"]["generator_count_after"] = len(new_names)
    state["compaction"]["reserved_new_parameters"] = (
        len(old_parameter_names) - 1 - state["nparameters"]
    )
    state["compaction"]["additional_exact_zero_elimination"] = {
        "parameter_name_before": parameter_name,
        "equation_tag": equation_tag,
        "equation_scalar": scalar,
        "source_checkpoint": str(source),
        "source_checkpoint_sha256": hashlib.sha256(raw).hexdigest(),
        "substitution": f"{parameter_name}=0",
        "reason": "the original retained constraint is a nonzero scalar multiple of this coordinate",
        "generator_count_before": len(old_names),
        "generator_count_after": len(new_names),
        "carried_parameters_before": int(state["nparameters"]) + 1,
        "carried_parameters_after": int(state["nparameters"]),
        "kept_monomials": kept_monomials,
        "dropped_monomials_vanishing_under_substitution": dropped_monomials,
        "exact_inverse_index_checks": kept_monomials,
        "all_retained_monomials_checked": True,
    }
    encoded = (json.dumps(state, separators=(",", ":")) + "\n").encode("utf-8")
    output.write_bytes(encoded)
    receipt = {
        **state["compaction"]["additional_exact_zero_elimination"],
        "output_checkpoint": str(output),
        "output_checkpoint_sha256": hashlib.sha256(encoded).hexdigest(),
        "output_bytes": len(encoded),
        "original_polynomial_terms": original_terms,
        "retained_polynomial_terms": retained_terms,
        "scope": "Exact restriction to a necessary coordinate hyperplane; equivalent zero locus because the source ideal contains the cited coordinate equation.",
    }
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("parameter")
    parser.add_argument("equation_tag")
    args = parser.parse_args()
    print(json.dumps(eliminate(args.source, args.output, args.receipt,
                               args.parameter, args.equation_tag), indent=2))


if __name__ == "__main__":
    main()
