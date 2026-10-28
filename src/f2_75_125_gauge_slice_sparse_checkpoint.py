"""Choose the exact x-translation gauge [t^2]P10 = 0 in a sparse checkpoint.

This is a coordinate slice for the source automorphism x -> x+c, not a
mathematical exclusion.  The input checkpoint must already contain the exact
depth-six 2+2 lift (preferably after corrected-power replay).  All stored
polynomials and constraints are specialized and reindexed exactly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ

from f2_75_125_sparse_ring_lift import SparseLift


def encode(poly, field):
    """Use the sparse checkpoint's exact coefficient encoding."""
    terms = []
    for powers, coefficient in sorted(poly.items()):
        sparse = [[i, power] for i, power in enumerate(powers) if power]
        values = [coefficient] if field == QQ else coefficient.to_list()
        terms.append([sparse, [str(value) for value in values]])
    return terms


def field_string(field, value):
    return str(sp.simplify(field.to_sympy(value)))


def compose_everywhere(lift, replacements):
    """Apply a simultaneous polynomial substitution to all checkpoint data."""
    for attribute in ("S", "kappa", "kappa2", "kappa3"):
        setattr(lift, attribute, getattr(lift, attribute).compose(replacements))
    if lift.B is not None:
        lift.B = lift.B.compose(replacements)
    lift.uk = {key: value.compose(replacements)
               for key, value in lift.uk.items()}
    lift.p = [value.compose(replacements) for value in lift.p]
    lift.q = [value.compose(replacements) for value in lift.q]
    lift.w = [value.compose(replacements) for value in lift.w]
    lift.constraints = [(tag, value.compose(replacements))
                        for tag, value in lift.constraints]
    for attribute in ("center_a", "center_b"):
        if hasattr(lift, attribute):
            setattr(lift, attribute,
                    getattr(lift, attribute).compose(replacements))


def active_parameter_indices(lift):
    polynomials = [lift.S, lift.kappa, lift.kappa2, lift.kappa3]
    if lift.B is not None:
        polynomials.append(lift.B)
    polynomials += list(lift.uk.values()) + lift.p + lift.q + lift.w
    polynomials += [value for _, value in lift.constraints]
    for attribute in ("center_a", "center_b"):
        if hasattr(lift, attribute):
            polynomials.append(getattr(lift, attribute))
    used = set()
    for polynomial in polynomials:
        for powers in polynomial:
            used.update(i - 1 for i, exponent in enumerate(powers)
                        if 1 <= i <= 160 and exponent)
    return sorted(used)


def verify_terminal_p(lift):
    """Recheck every saved P terminal position, not only the first 32."""
    positions = 0
    through_five = 0
    by_distance = []
    position_records = []
    for distance in range(len(lift.p)):
        bound = (71 - 5 * distance) // 12
        if bound < 0:
            continue
        chart = lift.centered_slice(distance, bound)
        checked_here = 0
        for degree in range(bound + 1):
            checked_here += 1
            positions += 1
            if distance <= 5:
                through_five += 1
            coefficient = lift.local_coefficient(chart, degree)
            if coefficient:
                raise AssertionError(
                    f"nonzero P terminal coefficient at ({distance},{degree})"
                )
            position_records.append({
                "distance": distance,
                "u_degree": degree,
                "coefficient_is_zero": True,
                "original_32_position": distance <= 5,
            })
        by_distance.append({"distance": distance, "positions": checked_here})
    expected_original = [
        [distance, degree]
        for distance in range(6)
        for degree in range((71 - 5 * distance) // 12 + 1)
    ]
    actual_original = [
        [record["distance"], record["u_degree"]]
        for record in position_records if record["original_32_position"]
    ]
    if through_five != 32 or actual_original != expected_original:
        raise AssertionError(f"expected 32 original P positions, got {through_five}")
    if positions != 36:
        raise AssertionError(f"expected 36 saved P positions at depth 6, got {positions}")
    return {"positions_checked": positions,
            "positions_through_distance5": through_five,
            "original_32_positions_explicit": expected_original,
            "all_positions_with_zero_status": position_records,
            "nonzero_positions": 0,
            "by_distance": by_distance}


def constraint_inventory(lift):
    inventory = []
    zero_tags = []
    for tag, polynomial in lift.constraints:
        if not polynomial:
            zero_tags.append(tag)
            continue
        encoded = encode(polynomial, lift.field)
        payload = json.dumps(encoded, separators=(",", ":"))
        degree = max((sum(powers[1:]) for powers in polynomial), default=0)
        inventory.append({
            "tag": tag,
            "terms": len(polynomial),
            "parameter_degree": degree,
            "sha256_of_exact_sparse_polynomial": hashlib.sha256(
                payload.encode("utf-8")).hexdigest(),
        })
    return inventory, zero_tags


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path,
                        help="exact depth-6 2+2 checkpoint; replayed state preferred")
    parser.add_argument("--output", type=Path, required=True,
                        help="new gauge-sliced checkpoint; must not already exist")
    parser.add_argument("--receipt", type=Path, required=True,
                        help="new gauge-slice receipt; must not already exist")
    args = parser.parse_args()

    source = args.checkpoint.resolve()
    output = args.output.resolve()
    receipt_path = args.receipt.resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    if source in (output, receipt_path) or output == receipt_path:
        raise ValueError("input, output checkpoint, and receipt paths must differ")
    if output.exists() or receipt_path.exists():
        raise FileExistsError("refusing to overwrite an existing output or receipt")

    print(f"loading exact checkpoint: {source}", flush=True)
    input_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    lift = SparseLift.from_checkpoint(source)
    print("checkpoint replay and original terminal-P identity checks passed", flush=True)
    if lift.tag != "2+2" or len(lift.p) != 7 or len(lift.q) != 7:
        raise ValueError("expected a depth-6 2+2 sparse-lift checkpoint")
    if len(lift.w) != 7:
        raise ValueError("expected seven saved W slices through distance 6")
    if lift.nparameters != 56:
        raise ValueError(f"expected 56 pre-slice parameters, got {lift.nparameters}")

    field = lift.field
    if field != QQ.algebraic_field(sp.sqrt(5)):
        raise ValueError("the 2+2 gauge slice requires QQ(sqrt(5))")

    # In the checkpoint ring, generator 0 is t and generator 46 is p45.
    # The serialized P10 beta confirms this zero-based parameter name.
    pivot_index = 46
    pivot = lift.parameters[45]
    beta = lift.local_coefficient(lift.p[5], 2)  # [t^2]P10
    zero_powers = (0,) * len(lift.ring.gens)
    pivot_powers = list(zero_powers)
    pivot_powers[pivot_index] = 1
    beta_terms = dict(beta.items())
    gamma = beta_terms.get(tuple(pivot_powers), field.zero)
    if not gamma:
        raise AssertionError("p45 has zero coefficient in beta")

    # The audited pivot must be a simple affine term, with no p45 powers or
    # mixed monomials.  This makes solving beta=0 exact and unambiguous.
    for powers in beta:
        exponent = powers[pivot_index]
        if exponent and (exponent != 1 or
                         any(powers[i] for i in range(len(powers))
                             if i != pivot_index)):
            raise AssertionError("beta is not affine-linear in p45")
    gamma_poly = lift.ring.ground_new(gamma)
    beta_remainder = beta - gamma_poly * pivot
    if any(powers[pivot_index] for powers in beta_remainder):
        raise AssertionError("beta remainder still depends on p45")

    expected_gamma = field.from_sympy(
        sp.Rational(260329251, 173750) * sp.sqrt(5)
        + sp.Rational(471620427, 173750)
    )
    if gamma != expected_gamma:
        raise AssertionError("p45 pivot coefficient differs from the audited value")
    print("verified nonzero affine p45 pivot in beta", flush=True)

    # The x-translation changes beta by 3*h(0)^3*c.  Verify this exact
    # nonzero source-orbit coefficient from the saved edge polynomial.
    h0 = dict(lift.h.items()).get(zero_powers, field.zero)
    source_shift_coefficient = field.convert(3) * h0**3
    expected_shift = field.from_sympy(483 + 216 * sp.sqrt(5))
    if not source_shift_coefficient or source_shift_coefficient != expected_shift:
        raise AssertionError("unexpected or zero x-translation coefficient")

    gamma_inverse = field.from_sympy(1 / field.to_sympy(gamma))
    p45_solution = -beta_remainder * lift.ring.ground_new(gamma_inverse)
    if any(powers[pivot_index] for powers in p45_solution):
        raise AssertionError("solved p45 expression is not independent of p45")

    constraints_before = len(lift.constraints)
    # First eliminate the old p45 exactly throughout all checkpoint fields.
    compose_everywhere(lift, {pivot: p45_solution})
    print("substituted the gauge solution through every checkpoint polynomial", flush=True)

    # Compact the remaining old parameter slots 46..55 to slots 45..54, so the
    # existing affine allocator can safely continue at parameter index 55.
    reindex = {lift.parameters[index]: lift.parameters[index - 1]
               for index in range(46, 56)}
    compose_everywhere(lift, reindex)
    lift.nparameters = 55

    beta_after = lift.local_coefficient(lift.p[5], 2)
    if beta_after:
        raise AssertionError("gauge-sliced checkpoint has nonzero [t^2]P10")
    if len(lift.constraints) != constraints_before:
        raise AssertionError("constraint list length changed during substitution")

    terminal_p = verify_terminal_p(lift)
    print("rechecked all 36 terminal-P positions after substitution", flush=True)
    nonzero_constraints, zero_constraint_tags = constraint_inventory(lift)
    active = active_parameter_indices(lift)
    if any(index >= lift.nparameters for index in active):
        raise AssertionError("active parameter lies beyond the allocator cursor")

    output.parent.mkdir(parents=True, exist_ok=True)
    lift.save_checkpoint(output)
    saved = json.loads(output.read_text(encoding="utf-8"))
    saved["source_translation_gauge"] = {
        "status": "coordinate gauge slice; not an exclusion",
        "normalization": "[t^2]P10 = 0",
        "source_automorphism": "x -> x+c",
        "source_shift_coefficient": field_string(field, source_shift_coefficient),
        "pivot_parameter_before_reindex": "p45",
        "pivot_coefficient_in_beta_before_reindex": field_string(field, gamma),
        "affine_coordinate_change_invertible": True,
        "pivot_eliminated_and_remaining_parameters_reindexed": True,
        "active_parameter_count": len(active),
    }
    output.write_text(json.dumps(saved, separators=(",", ":")) + "\n",
                      encoding="utf-8")
    output_hash = hashlib.sha256(output.read_bytes()).hexdigest()

    report = {
        "branch": lift.tag,
        "input_checkpoint": str(source),
        "input_sha256": input_hash,
        "output_checkpoint": str(output),
        "output_sha256": output_hash,
        "depth": len(lift.p) - 1,
        "gauge": {
            "kind": "source-translation coordinate slice; not an exclusion",
            "normalization": "beta=[t^2]P10=0",
            "source_automorphism": "x -> x+c",
            "source_shift_coefficient_3h0_cubed": field_string(
                field, source_shift_coefficient),
            "coefficient_of_old_p45_in_beta": field_string(field, gamma),
            "beta_terms_before_slice": len(beta),
            "beta_remainder_terms_before_slice": len(beta_remainder),
            "old_p45_solution_exact_sparse_terms": encode(p45_solution, field),
            "coordinate_change": (
                "beta=gamma*p45+beta_remainder, with all other allocated "
                "coordinates fixed; inverse p45=-beta_remainder/gamma "
                "on beta=0"
            ),
            "old_parameter_names": ["t", *(f"p{i}" for i in range(56))],
        },
        "parameters": {
            "allocated_before": 56,
            "allocator_cursor_after": lift.nparameters,
            "active_after": len(active),
            "active_indices_after": active,
        },
        "substitution": {
            "all_parameterized_checkpoint_fields": [
                "S", "B", "kappa", "kappa2", "kappa3", "uk", "p", "q", "w",
                "constraints", "center_a", "center_b"
            ],
            "checkpoint_affine_block_receipts_preserved": lift.receipt,
            "fixed_source_fields": ["h", "v", "u"],
        },
        "P_terminal_recheck": terminal_p,
        "constraints": {
            "count_before": constraints_before,
            "count_after": len(lift.constraints),
            "nonzero_count_after": len(nonzero_constraints),
            "input_constraint_list_empty": constraints_before == 0,
            "terminal_rows_enforced_by_checkpoint_identity_checks": True,
            "constraint_scope_note": (
                "No explicit retained constraint polynomials were saved in "
                "this checkpoint. SparseLift.from_checkpoint rejects the "
                "checkpoint unless every saved terminal-P slice vanishes; "
                "the helper then directly recomputes all 36 terminal-P "
                "positions after gauge substitution."
            ),
            "nonzero_constraints_with_tags_and_hashes": nonzero_constraints,
            "zero_after_slice_but_retained_with_tags": zero_constraint_tags,
            "exact_polynomials_retained_in_output_checkpoint": True,
        },
        "scope": (
            "exact specialization to the source-translation gauge slice; "
            "all stored p/q/w slices and retained constraints are substituted"
        ),
    }
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(report, indent=2) + "\n",
                            encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
