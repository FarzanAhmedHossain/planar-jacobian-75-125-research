"""Exactly relabel a saved sparse lift to remove unused parameter variables.

Every monomial's coefficient and exponent are preserved.  Used old parameters
are relabelled contiguously, while a reserve of fresh parameter names remains
available for the next affine block.  A per-monomial inverse-map check makes
this a coordinate relabelling, not a specialization of the solution set.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


POLY_KEYS = ("h", "v", "S", "B", "kappa", "kappa2", "kappa3")
SEQUENCE_KEYS = ("p", "q", "w")


def polynomial_entries(state: dict):
    for key in POLY_KEYS:
        if state.get(key) is not None:
            yield key, state[key]
    for key in SEQUENCE_KEYS:
        for i, polynomial in enumerate(state.get(key, [])):
            yield f"{key}[{i}]", polynomial
    for key, polynomial in state.get("uk", {}).items():
        yield f"uk[{key}]", polynomial
    for i, (_tag, polynomial) in enumerate(state.get("constraints", [])):
        yield f"constraints[{i}]", polynomial


def compact(source: Path, output: Path, receipt_path: Path,
            reserve_new_parameters: int) -> dict:
    raw = source.read_bytes()
    state = json.loads(raw)
    if state.get("format") != "f2-sparse-ring-v1":
        raise ValueError("unexpected sparse checkpoint format")
    if state.get("compaction") is not None:
        raise ValueError("refusing to compact an already compacted checkpoint")
    if reserve_new_parameters < 0:
        raise ValueError("reserve must be nonnegative")

    old_names = state["generators"]
    old_index = {name: i for i, name in enumerate(old_names)}
    used_indices = set()
    monomial_count = 0
    term_count = 0
    for _path, polynomial in polynomial_entries(state):
        term_count += len(polynomial)
        for sparse, _values in polynomial:
            monomial_count += 1
            used_indices.update(index for index, _power in sparse)

    used_names = {old_names[index] for index in used_indices}
    active_old_parameters = sorted(
        (name for name in used_names
         if name.startswith("p") and name[1:].isdigit()),
        key=lambda name: int(name[1:]),
    )
    if not active_old_parameters:
        raise ValueError("checkpoint has no active parameter variables")
    old_nparameters = int(state["nparameters"])
    if any(int(name[1:]) >= old_nparameters for name in active_old_parameters):
        raise ValueError("an active parameter lies beyond the allocated parameter prefix")
    if "t" not in old_names:
        raise ValueError("checkpoint has no t generator")

    parameter_map = {
        old_name: f"p{new_index}"
        for new_index, old_name in enumerate(active_old_parameters)
    }
    capacity = len(active_old_parameters) + reserve_new_parameters
    new_names = ["t", *(f"p{i}" for i in range(capacity)),
                 *(f"x{i}" for i in range(73))]
    new_index = {name: i for i, name in enumerate(new_names)}

    # Preserve unknown names and their order, but move them behind the compact
    # parameter prefix.  No unknown is allowed to disappear from the ring.
    for i in range(73):
        if f"x{i}" not in old_index:
            raise ValueError(f"source checkpoint lacks x{i}")

    remapped_monomials = 0

    def remap_polynomial(polynomial: list) -> list:
        nonlocal remapped_monomials
        result = []
        for sparse, values in polynomial:
            new_sparse = []
            for old_generator_index, exponent in sparse:
                old_name = old_names[old_generator_index]
                if old_name == "t":
                    target_name = "t"
                elif old_name in parameter_map:
                    target_name = parameter_map[old_name]
                elif old_name.startswith("x") and old_name[1:].isdigit():
                    target_name = old_name
                else:
                    raise ValueError(f"unsupported active generator {old_name}")
                new_sparse.append([new_index[target_name], exponent])
            new_sparse.sort()

            # Verify this monomial maps back exactly, including its exponent
            # data.  The coefficient list is retained byte-for-byte.
            back = []
            for compact_index, exponent in new_sparse:
                compact_name = new_names[compact_index]
                if compact_name == "t":
                    source_name = "t"
                elif compact_name.startswith("p") and compact_name[1:].isdigit():
                    source_name = next(
                        old for old, new in parameter_map.items()
                        if new == compact_name
                    )
                else:
                    source_name = compact_name
                back.append([old_index[source_name], exponent])
            back.sort()
            if back != sparse:
                raise AssertionError("monomial relabelling failed its exact inverse check")
            result.append([new_sparse, values])
            remapped_monomials += 1
        return result

    for key in POLY_KEYS:
        if state.get(key) is not None:
            state[key] = remap_polynomial(state[key])
    for key in SEQUENCE_KEYS:
        state[key] = [remap_polynomial(poly) for poly in state.get(key, [])]
    state["uk"] = {key: remap_polynomial(poly)
                    for key, poly in state.get("uk", {}).items()}
    state["constraints"] = [
        [tag, remap_polynomial(poly)] for tag, poly in state.get("constraints", [])
    ]

    source_hash = hashlib.sha256(raw).hexdigest()
    state["generators"] = new_names
    state["nparameters"] = len(active_old_parameters)
    state["compaction"] = {
        "format": "exact-parameter-relabel-v1",
        "source_checkpoint": str(source),
        "source_checkpoint_sha256": source_hash,
        "old_nparameters": old_nparameters,
        "new_nparameters": len(active_old_parameters),
        "active_old_parameters": active_old_parameters,
        "old_to_new_parameter_names": parameter_map,
        "reserved_new_parameters": reserve_new_parameters,
        "generator_count_before": len(old_names),
        "generator_count_after": len(new_names),
        "term_count_before_and_after": term_count,
        "monomial_inverse_checks": remapped_monomials,
        "exact_relabeling_verified": remapped_monomials == monomial_count,
    }
    if not state["compaction"]["exact_relabeling_verified"]:
        raise AssertionError("not all monomials passed the exact inverse check")
    if output.exists() or receipt_path.exists():
        raise FileExistsError("refusing to overwrite a prior compact checkpoint or receipt")

    encoded = (json.dumps(state, separators=(",", ":")) + "\n").encode("utf-8")
    output.write_bytes(encoded)
    receipt = {
        **state["compaction"],
        "output_checkpoint": str(output),
        "output_checkpoint_sha256": hashlib.sha256(encoded).hexdigest(),
        "output_bytes": len(encoded),
        "term_count": term_count,
        "active_old_generator_indices": sorted(used_indices),
        "scope": "Exact coordinate relabeling only; no equations were added or removed.",
    }
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--reserve", type=int, default=72,
                        help="new affine parameters reserved for the next block")
    args = parser.parse_args()
    print(json.dumps(compact(args.source, args.output, args.receipt, args.reserve),
                     indent=2))


if __name__ == "__main__":
    main()
