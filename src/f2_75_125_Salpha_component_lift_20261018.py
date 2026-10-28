"""Continue the dual-root F2 lift on one exact depth-seven component.

The source checkpoint's ideal was verified to be (F3*F4), up to a nonzero
scalar.  This driver adds F3=0 or F4=0 as a retained necessary condition,
then propagates it through the same affine and terminal reductions as the
dual-root lift.  A resulting checkpoint is still a necessary system; a
component is excluded only after an exact unit-ideal certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from f2_75_125_sparse_ring_lift import SparseLift
from f2_75_125_dual_root_terminal_lift_20261018 import (
    append_existing_alpha_rows,
    impose_next_alpha_rows,
    initialize_alpha_chart,
)


DEFAULT_CHECKPOINT = HERE / "f2_75_125_sparse_ring_lift_22_terminal_roots7_gauged_Salpha_slice_state.json"
FACTOR_LOG = HERE / "f2_75_125_sparse_ring_lift_22_Salpha_slice_factorcheck_verified.log"
FACTOR_RECEIPT = HERE / "f2_75_125_Salpha_slice_gb_receipt_20261018.json"


def read_factor(component: int) -> tuple[str, sp.Expr]:
    log = FACTOR_LOG.read_text(encoding="utf-8")
    factor_index = component - 1  # factorize(): _[1] is the scalar, _[2:3] are factors
    match = re.search(rf"^\s*_\[{factor_index}\]=(.*)$", log, re.MULTILINE)
    if not match:
        raise ValueError(f"factor _[{factor_index}] not found in verified Singular log")
    names = ("p1", "p7", "p8", "p9", "p16", "p17", "p18", "p19",
             "p20", "p27", "p28", "p29", "p30", "p31", "p32", "p33")
    namespace = {name: sp.Symbol(name) for name in names}
    namespace["a"] = sp.Symbol("a")
    source = match.group(1).strip()
    expression = sp.expand(sp.sympify(source, locals=namespace).subs(
        namespace["a"], sp.sqrt(5)
    ))
    return source, expression


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--component", type=int, choices=(3, 4), required=True)
    parser.add_argument("--resume", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--through", type=int, default=8)
    parser.add_argument(
        "--eliminate-component-parameter", action="store_true",
        help=("solve the verified factor equation for a parameter with "
              "nonzero constant coefficient, then substitute exactly"),
    )
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()

    base_receipt = json.loads(FACTOR_RECEIPT.read_text(encoding="utf-8"))
    expected_hash = base_receipt["input_state_sha256"]
    source_hash = hashlib.sha256(args.resume.read_bytes()).hexdigest()
    lift = SparseLift.from_checkpoint(
        args.resume, verify_saved_terminal_rows=False
    )
    compaction = getattr(lift, "compaction", None) or {}
    original_source_hash = compaction.get("source_checkpoint_sha256", source_hash)
    if original_source_hash != expected_hash:
        raise ValueError("source checkpoint hash does not match the verified factorization")

    factor_source, factor_expr = read_factor(args.component)
    old_to_new = compaction.get("old_to_new_parameter_names", {})
    if old_to_new:
        factor_expr = factor_expr.xreplace({
            sp.Symbol(old): sp.Symbol(new) for old, new in old_to_new.items()
        })
    if lift.tag != "2+2" or lift.field == QQ or args.component not in (3, 4):
        raise ValueError("component lift requires the quadratic-field F2 2+2 slice")
    component_poly = lift.ring.from_expr(factor_expr)
    if not component_poly:
        raise AssertionError("component factor decoded as zero")
    tag = f"Salpha-component-F{args.component}"
    matching_constraints = [value for existing, value in lift.constraints
                            if existing == tag]
    factor_log_hash = hashlib.sha256(FACTOR_LOG.read_bytes()).hexdigest()
    resume_receipt = None
    prior = None
    resume_name = args.resume.name
    if resume_name.endswith("_state.json"):
        resume_receipt = args.resume.with_name(
            resume_name[:-len("_state.json")] + ".json"
        )
        if resume_receipt.exists():
            prior = json.loads(resume_receipt.read_text(encoding="utf-8"))
    factor_already_eliminated = False
    if matching_constraints:
        if len(matching_constraints) != 1:
            raise ValueError("saved state has duplicate component-factor constraints")
        factor_was_present = True
        # A resumed checkpoint may already have passed through affine
        # parameter substitutions, so its tagged factor is the exact image
        # of the original factor, not necessarily its original expression.
        # Verify the producing receipt and use that pushed-forward equation.
        if resume_receipt is None:
            raise ValueError("cannot infer receipt for a component checkpoint")
        if prior is None:
            raise FileNotFoundError(f"component checkpoint receipt missing: {resume_receipt}")
        if (prior.get("component") != f"F{args.component}=0"
                or prior.get("component_tag") != tag
                or prior.get("factor_log_sha256") != factor_log_hash
                or prior.get("original_source_checkpoint_sha256") != expected_hash
                or prior.get("digest") != lift.digest()):
            raise ValueError("resumed component state does not match its exact source receipt")
        factor_equation = matching_constraints[0]
        factor_already_eliminated = not factor_equation
    elif (prior is not None
          and prior.get("component") == f"F{args.component}=0"
          and prior.get("component_tag") == tag
          and prior.get("factor_log_sha256") == factor_log_hash
          and prior.get("original_source_checkpoint_sha256") == expected_hash
          and prior.get("component_parameter_elimination") is not None
          and prior.get("digest") == lift.digest()):
        # Older component runs removed the now-zero factor row when saving.
        # The verified producing receipt preserves its provenance.
        factor_was_present = True
        factor_already_eliminated = True
        factor_equation = lift.ring.zero
    else:
        lift.constraints.append((tag, component_poly))
        factor_was_present = False
        factor_equation = component_poly

    elimination = None
    if args.eliminate_component_parameter:
        # In the compact parameterization, F3 is affine with a nonzero
        # coefficient in old p16 (new p4), while F4 is affine in old p27
        # (new p10).  Their coefficients lie in Q(sqrt(5)) and are units.
        # The quotient by this component equation is therefore exactly a
        # polynomial ring in the other parameters; no saturation is needed.
        if not factor_equation:
            # Earlier affine substitutions may already have made the tagged
            # factor identically zero on the retained affine chart.
            elimination = {
                "component_factor": "identically zero after verified resume",
                "zero_constraints_removed": 0,
            }
        else:
            preferred_old = {
            3: ("p16", "p17", "p18", "p19", "p20"),
            4: ("p27", "p28", "p29", "p30", "p31", "p32", "p33"),
            }[args.component]
            preferred_names = [old_to_new.get(name, name) for name in preferred_old]
            preferred = [value for name in preferred_names
                         for value in lift.parameters[:lift.nparameters]
                         if str(value) == name]
            candidates = list(dict.fromkeys(
                [*preferred, *lift.parameters[:lift.nparameters]]
            ))
            factor_expr_current = factor_equation.as_expr()
            elimination_data = None
            for parameter_generator in candidates:
                parameter_name = str(parameter_generator)
                parameter_symbol = sp.Symbol(parameter_name)
                factor_as_poly = sp.Poly(factor_expr_current, parameter_symbol)
                if factor_as_poly.degree() != 1:
                    continue
                coefficient = factor_as_poly.coeff_monomial(parameter_symbol)
                if coefficient.free_symbols or coefficient == 0:
                    continue
                coefficient_element = lift.field.from_sympy(coefficient)
                if not coefficient_element:
                    continue
                solution_expr = sp.expand(
                    -factor_as_poly.coeff_monomial(1) / coefficient
                )
                if parameter_symbol in solution_expr.free_symbols:
                    continue
                solution_poly = lift.ring.from_expr(solution_expr)
                replacement = {parameter_generator: solution_poly}
                if factor_equation.compose(replacement):
                    continue
                elimination_data = (parameter_name, coefficient, solution_expr,
                                    solution_poly, replacement, parameter_generator)
                break
            if elimination_data is None:
                raise ValueError(
                    f"no current parameter has a nonzero constant linear coefficient "
                    f"in the F{args.component} component equation"
                )
            (parameter_name, coefficient, solution_expr,
             solution_poly, replacement, parameter_generator) = elimination_data
            lift.replace_parameters(replacement)
            if factor_equation.compose(replacement):
                raise AssertionError("exact component substitution did not annihilate its factor")
            previous_constraints = len(lift.constraints)
            lift.constraints = [
                (name, value) for name, value in lift.constraints
                if value or name == tag
            ]
            elimination = {
                "preferred_source_parameters": list(preferred_old),
                "ring_parameter": parameter_name,
                "coefficient": str(coefficient),
                "solution_term_count": len(sp.Add.make_args(solution_expr)),
                "zero_constraints_removed": previous_constraints-len(lift.constraints),
            }
            print(json.dumps({"stage": "exact component parameter elimination",
                              **elimination}), flush=True)

    alpha, c0, center_a, center_b = initialize_alpha_chart(lift)
    if factor_was_present:
        # The matching checkpoint receipt and digest were verified above; its
        # alpha rows were already backfilled and propagated through this depth.
        # Re-expanding every historical root chart would only duplicate work.
        added = []
    else:
        added = append_existing_alpha_rows(lift, alpha, center_a, center_b)
    lift.alpha_center_a = center_a
    lift.alpha_center_b = center_b
    original_replace = lift.replace_parameters

    def replace_with_alpha(replacements):
        original_replace(replacements)
        lift.alpha_center_a = lift.alpha_center_a.compose(replacements)
        lift.alpha_center_b = lift.alpha_center_b.compose(replacements)

    lift.replace_parameters = replace_with_alpha
    started = time.monotonic()
    for distance in range(len(lift.p), args.through + 1):
        lift.block(distance)
        lift.impose_terminal_chart(distance)
        lift.impose_terminal_q(distance)
        lift.impose_root_values(distance)
        impose_next_alpha_rows(
            lift, distance, alpha, lift.alpha_center_a, lift.alpha_center_b
        )
        lift.save_checkpoint(args.checkpoint)
    lift.save_checkpoint(args.checkpoint)

    receipt = {
        "branch": lift.tag,
        "component": f"F{args.component}=0",
        "component_tag": tag,
        "factor_log_sha256": hashlib.sha256(FACTOR_LOG.read_bytes()).hexdigest(),
        "source_checkpoint": str(args.resume),
        "source_checkpoint_sha256": source_hash,
        "original_source_checkpoint_sha256": original_source_hash,
        "parameter_relabeling": old_to_new,
        "factor_expression": factor_source,
        "factor_added_before_continuation": True,
        "factor_present_in_resume_state": factor_was_present,
        "factor_already_eliminated_in_resume_state": factor_already_eliminated,
        "component_parameter_elimination": elimination,
        "backfilled_alpha_rows": len(added),
        "through_distance": args.through,
        "parameter_count": lift.nparameters,
        "constraint_count": len(lift.constraints),
        "constraint_term_count": sum(len(poly) for _, poly in lift.constraints),
        "blocks": lift.receipt,
        "digest": lift.digest(),
        "elapsed_seconds": time.monotonic() - started,
        "scope": "component-restricted necessary dual-root system; not an exclusion certificate",
    }
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt), flush=True)


if __name__ == "__main__":
    main()
