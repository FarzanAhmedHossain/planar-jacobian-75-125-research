"""Lift a sparse-graph point back to the original F2 component checkpoint."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import sympy as sp

import f2_75_125_quadratic_affine_profile as core
from f2_75_125_quadratic_affine_profile_stream import (
    iter_json_array,
    top_metadata_and_constraints,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_rational(value) -> Fraction:
    value = sp.Rational(value)
    return Fraction(int(value.p), int(value.q))


def parse_quad(text: str):
    root5 = sp.sqrt(5)
    expr = sp.expand(sp.sympify(text, locals={"sqrt": sp.sqrt}))
    poly = sp.Poly(expr, root5)
    return (
        parse_rational(poly.coeff_monomial(1)),
        parse_rational(poly.coeff_monomial(root5)),
    )


def load_rref(receipt):
    if "rref_field_elements" in receipt:
        return [[core.from_json(value) for value in row]
                for row in receipt["rref_field_elements"]]
    return [[parse_quad(value) for value in row] for row in receipt["rref"]]


def monomial_value(monomial, assignment):
    value = core.ONE
    for index, exponent in monomial:
        for _ in range(exponent):
            value = core.mul(value, assignment[index])
    return value


def monomials_for(component, stage_index, receipt, names):
    if "monomials" in receipt:
        return [tuple((int(i), int(e)) for i, e in monomial)
                for monomial in receipt["monomials"]]
    if component != "F3=0":
        raise ValueError("F4 stages must carry their monomial index lists")
    if stage_index == 0:  # support-six graph on p5,p6,p7,p8,U=p0^3,V=p0*p3
        return [
            ((names.index("p5"), 1),),
            ((names.index("p6"), 1),),
            ((names.index("p7"), 1),),
            ((names.index("p8"), 1),),
            ((names.index("p0"), 3),),
            tuple(sorted(((names.index("p0"), 1), (names.index("p3"), 1)))),
        ]
    if stage_index == 1:  # F3 distance-four graph
        return [
            ((names.index(f"p{i}"), 1),) for i in range(10, 17)
        ] + [
            ((names.index("p0"), 4),),
            tuple(sorted(((names.index("p0"), 2), (names.index("p3"), 1)))),
            tuple(sorted(((names.index("p0"), 1), (names.index("p8"), 1)))),
            ((names.index("p3"), 2),),
        ]
    raise ValueError(f"unknown legacy F3 receipt at stage {stage_index}")


def verify_constraints(path: Path, names, assignment):
    text = path.read_text(encoding="utf-8")
    metadata, pos = top_metadata_and_constraints(text)
    if metadata.get("generators") != names:
        raise ValueError(f"generator order changed in {path}")
    nonzero = []
    count = 0
    for tag, saved in iter_json_array(text, pos):
        count += 1
        value = core.ZERO
        for sparse, coefficient_values in saved:
            monomial = tuple((int(index), int(exponent))
                             for index, exponent in sparse if exponent)
            value = core.add(
                value,
                core.mul(core.from_json(coefficient_values),
                         monomial_value(monomial, assignment)),
            )
        if not core.iszero(value):
            nonzero.append({"tag": tag, "value": core.show(value)})
            if len(nonzero) >= 20:
                break
    return count, nonzero


def apply_stage(receipt_path: Path, checkpoint_path: Path, stage_index: int,
                component: str, names, assignment):
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("component") != component:
        raise ValueError(f"component mismatch in {receipt_path}")
    expected_hash = (receipt.get("source_checkpoint_sha256")
                     or receipt.get("checkpoint_sha256")
                     or receipt.get("source_sha256"))
    actual_hash = sha256(checkpoint_path)
    if expected_hash != actual_hash:
        raise ValueError(f"checkpoint hash mismatch for {receipt_path}")
    rref = load_rref(receipt)
    pivots = [int(value) for value in receipt["pivot_columns"]]
    monomials = monomials_for(component, stage_index, receipt, names)
    if len(rref) != len(pivots) or len(monomials) != len(rref[0]):
        raise ValueError(f"RREF dimensions mismatch in {receipt_path}")
    pivot_generators = []
    for row_index, pivot_column in enumerate(pivots):
        pivot_monomial = monomials[pivot_column]
        if len(pivot_monomial) != 1 or pivot_monomial[0][1] != 1:
            raise ValueError("a pivot column is not a single parameter")
        pivot_generator = pivot_monomial[0][0]
        value = core.ZERO
        for column, monomial in enumerate(monomials):
            if column in pivots:
                continue
            value = core.add(
                value,
                core.mul(rref[row_index][column],
                         monomial_value(monomial, assignment)),
            )
        assignment[pivot_generator] = core.sub(core.ZERO, value)
        pivot_generators.append(names[pivot_generator])
    return {
        "linear_receipt": str(receipt_path),
        "source_checkpoint": str(checkpoint_path),
        "source_checkpoint_sha256": actual_hash,
        "pivot_generators_solved": pivot_generators,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--component", choices=("F3=0", "F4=0"), required=True)
    parser.add_argument("--root-checkpoint", type=Path, required=True)
    parser.add_argument("--root-receipt", type=Path, required=True)
    parser.add_argument("--pregraph-checkpoint", type=Path, required=True)
    parser.add_argument("--linear-receipt", type=Path, action="append", required=True)
    parser.add_argument("--stage-checkpoint", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    if len(args.linear_receipt) != len(args.stage_checkpoint):
        raise ValueError("supply one source checkpoint for every linear receipt")

    root_receipt = json.loads(args.root_receipt.read_text(encoding="utf-8"))
    if root_receipt.get("component") != args.component or root_receipt.get("through_distance") != 9:
        raise ValueError("root sidecar does not certify the requested depth-nine component")
    root_hash = sha256(args.root_checkpoint)
    pregraph_text = args.pregraph_checkpoint.read_text(encoding="utf-8")
    pregraph = json.loads(pregraph_text)
    if pregraph.get("compaction", {}).get("source_checkpoint_sha256") != root_hash:
        raise ValueError("pregraph necessary-system checkpoint is not derived from root checkpoint")
    names = pregraph.get("generators")
    if not names:
        raise ValueError("generator list missing from pregraph checkpoint")
    root_text = args.root_checkpoint.read_text(encoding="utf-8")
    root_data, _root_pos = top_metadata_and_constraints(root_text)
    if root_data.get("generators") != names:
        raise ValueError("root and pregraph generator orders differ")
    if root_data.get("branch") != "2+2" or root_data.get("quadratic") is not True:
        raise ValueError("root checkpoint is not the quadratic F2 2+2 system")
    del root_text, pregraph_text

    assignment = [(index + 2, 0) for index in range(len(names))]
    stage_results = []
    for stage_index, (receipt, checkpoint) in enumerate(
            zip(args.linear_receipt, args.stage_checkpoint)):
        stage_results.append(apply_stage(
            receipt, checkpoint, stage_index, args.component, names, assignment
        ))

    # The common support-four graph for p1,p2, replayed on both components.
    p0 = assignment[names.index("p0")]
    p3 = assignment[names.index("p3")]
    A = (Fraction(672, 5), Fraction(288, 5))
    B = (Fraction(1728), Fraction(19008, 25))
    C = (Fraction(8), Fraction(8, 5))
    D = (Fraction(459, 10), Fraction(153, 10))
    p2 = core.mul(core.inv(A), core.sub(core.mul(p0, p0), core.mul(B, p3)))
    p1 = core.sub(core.ZERO, core.add(core.mul(C, p2), core.mul(D, p3)))
    assignment[names.index("p1")] = p1
    assignment[names.index("p2")] = p2

    p1p2_checks = [
        core.sub(core.add(core.mul(A, p2), core.mul(B, p3)), core.mul(p0, p0)),
        core.add(core.add(p1, core.mul(C, p2)), core.mul(D, p3)),
    ]
    if any(not core.iszero(value) for value in p1p2_checks):
        raise AssertionError("support-four p1,p2 graph failed exact replay")

    checkpoint_verifications = []
    for checkpoint in (args.root_checkpoint, args.pregraph_checkpoint,
                       *args.stage_checkpoint):
        count, nonzero = verify_constraints(checkpoint, names, assignment)
        if nonzero:
            raise AssertionError(
                f"lifted assignment fails {len(nonzero)} constraints in {checkpoint}: {nonzero[:5]}"
            )
        checkpoint_verifications.append({
            "checkpoint": str(checkpoint),
            "checkpoint_sha256": sha256(checkpoint),
            "constraint_count": count,
            "nonzero_constraint_count": 0,
            "all_constraints_exactly_satisfied": True,
        })
    root_verification = checkpoint_verifications[0]
    pregraph_verification = checkpoint_verifications[1]

    result = {
        "component": args.component,
        "branch": "2+2",
        "root_checkpoint": str(args.root_checkpoint),
        "root_checkpoint_sha256": root_hash,
        "root_sidecar": str(args.root_receipt),
        "pregraph_checkpoint": str(args.pregraph_checkpoint),
        "pregraph_checkpoint_sha256": sha256(args.pregraph_checkpoint),
        "field": "Q(sqrt(5))",
        "free_value_rule": "each uneliminated generator starts at its position index plus 2 in Q; exact graph pivots and then p1,p2 are back-substituted",
        "stage_back_substitutions": stage_results,
        "support_four_p1_p2_graph_exactly_satisfied": True,
        "original_root_constraint_count": root_verification["constraint_count"],
        "original_root_nonzero_constraint_count": 0,
        "original_root_constraints_exactly_satisfied": True,
        "pregraph_constraint_count": pregraph_verification["constraint_count"],
        "pregraph_nonzero_constraint_count": 0,
        "pregraph_constraints_exactly_satisfied": True,
        "checkpoint_verifications": checkpoint_verifications,
        "generator_values": {
            name: core.show(value) for name, value in zip(names, assignment)
        },
        "zero_generator_values": [
            name for name, value in zip(names, assignment) if core.iszero(value)
        ],
        "scope": (
            "Exact common zero of the saved polynomial equalities in the original "
            "depth-nine component checkpoint and every supplied reduced checkpoint. "
            "It is not a Keller pair and does not certify omitted inequalities, "
            "chart-open conditions, or source-to-carrier attachment."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "component", "original_root_constraint_count",
        "original_root_nonzero_constraint_count",
        "original_root_constraints_exactly_satisfied",
        "pregraph_constraint_count", "pregraph_nonzero_constraint_count",
        "stage_back_substitutions", "zero_generator_values", "scope",
    )}, indent=2))


if __name__ == "__main__":
    main()
