"""Restrict a saved 2+1+1 checkpoint to its exact support-three radical line.

The d19 support-three Groebner report proves that every zero of that complete
necessary subsystem lies on p13=452*p15/5, p14=-14*p15.  Substituting this
line in the full checkpoint is therefore an exact necessary-slice reduction;
it does not claim to exclude the branch.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp

import f2_75_125_sparse_ring_lift as lift_module


HERE = Path(__file__).resolve().parent
DEFAULT_SOURCE = HERE / "f2_75_125_sparse_ring_lift_211_roots19_u1zero_streamed_state.json"
DEFAULT_LINE_REPORT = HERE / "f2_75_125_211_support3_radical_line_20261025.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--line-report", type=Path, default=DEFAULT_LINE_REPORT)
    parser.add_argument(
        "--checkpoint", type=Path,
        default=HERE / "f2_75_125_sparse_ring_lift_211_roots19_u1zero_support3line_state.json",
    )
    parser.add_argument(
        "--receipt", type=Path,
        default=HERE / "f2_75_125_211_support3line_specialization_receipt.json",
    )
    parser.add_argument("--flint-threads", type=int, default=8)
    args = parser.parse_args()
    if args.flint_threads < 1:
        raise ValueError("--flint-threads must be positive")

    line = json.loads(args.line_report.read_text(encoding="utf-8"))
    if (line.get("branch") != "2+1+1, U1=0, d19 necessary system"
            or line.get("geometric_zero_set") !=
            "(p13,p14,p15) = (452*s/5,-14*s,s), s in algebraic closure"
            or not line.get("line_containment_verified")
            or not line.get("basis_restrictions_verified")):
        raise ValueError("the support-three report does not verify the expected line")

    source_hash = sha256(args.source)
    if source_hash != line.get("source_sha256"):
        raise ValueError("support-three line report is not bound to this checkpoint")

    lift_module.POLY_BACKEND = "flint-mpoly"
    import flint
    flint.ctx.threads = args.flint_threads
    print(f"flint_threads={flint.ctx.threads}", flush=True)

    lift = lift_module.SparseLift.from_checkpoint_streaming(args.source)
    if lift.tag != "2+1+1" or len(lift.p) - 1 != 19:
        raise ValueError("expected the d19 2+1+1 checkpoint")

    p13, p14, p15 = lift.parameters[13:16]
    replacements = {
        p13: lift.const(sp.Rational(452, 5)) * p15,
        p14: -14 * p15,
    }
    if (5 * p13 - 452 * p15).compose(replacements):
        raise AssertionError("p13 line equation did not vanish after substitution")
    if (p14 + 14 * p15).compose(replacements):
        raise AssertionError("p14 line equation did not vanish after substitution")

    lift.replace_parameters(replacements)
    lift.constraints = [(tag, polynomial) for tag, polynomial in lift.constraints if polynomial]
    if hasattr(lift, "last_free"):
        lift.last_free = tuple(
            parameter for parameter in lift.last_free
            if parameter not in (p13, p14)
        )
    lift.receipt.append({
        "step": "exact support-three line specialization",
        "report_sha256": sha256(args.line_report),
        "source_sha256": source_hash,
        "substitutions": ["p13=452*p15/5", "p14=-14*p15"],
        "scope": "necessary d19 2+1+1 U1=0 slice only; no branch exclusion",
    })

    args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    lift.save_checkpoint(args.checkpoint)
    report = {
        "status": "specialized_necessary_slice_checkpoint_written",
        "branch": lift.tag,
        "depth": len(lift.p) - 1,
        "source": str(args.source),
        "source_sha256": source_hash,
        "line_report": str(args.line_report),
        "line_report_sha256": sha256(args.line_report),
        "substitutions": ["p13=452*p15/5", "p14=-14*p15"],
        "retained_constraints": len(lift.constraints),
        "checkpoint": str(args.checkpoint),
        "checkpoint_sha256": lift.last_checkpoint_sha256,
        "scope": "This is a necessary-condition slice, not an exclusion of 2+1+1 or (75,125).",
    }
    args.receipt.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
