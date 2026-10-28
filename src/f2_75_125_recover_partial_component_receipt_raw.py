"""Recover exact component sidecars directly from sparse JSON checkpoints.

This avoids rebuilding a very large SymPy polynomial ring solely to compute
the checkpoint digest and read zero/size metadata.  The digest formatting is
the one used by SymPy's sparse-ring coefficients with python-flint enabled;
the implementation is checked against existing F3 and F4 receipts.

This writes continuation metadata only.  It is not a branch-exclusion
certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path


HERE = Path(__file__).resolve().parent
FACTOR_LOG = HERE / "f2_75_125_sparse_ring_lift_22_Salpha_slice_factorcheck_verified.log"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def flint_rational_repr(value: str) -> str:
    """Match repr(fmpq) for the canonical rational strings saved in JSON."""
    if "/" in value:
        numerator, denominator = value.split("/", 1)
        numerator, denominator = int(numerator), int(denominator)
    else:
        numerator, denominator = int(value), 1
    if denominator == 1:
        return str(numerator)
    return f"{numerator}/{denominator}"


def coefficient_repr(values: list[str], quadratic: bool) -> str:
    if not quadratic:
        return flint_rational_repr(values[0])
    coefficients = ", ".join(flint_rational_repr(value) for value in values)
    return f"ANP([{coefficients}], [1, 0, -5], QQ)"


def checkpoint_digest(data: dict) -> str:
    """Reproduce SparseLift.digest() without constructing RingElements."""
    if not data.get("quadratic"):
        # This task's component checkpoints are all over Q(sqrt(5)); refusing
        # other fields avoids silently applying an unverified repr convention.
        raise ValueError("raw digest is audited only for the quadratic field Q(sqrt(5))")
    generator_count = len(data["generators"])
    digest = hashlib.sha256()

    def add_polynomial(encoded: list) -> None:
        for sparse, values in encoded:
            powers = [0] * generator_count
            for old_index, exponent in sparse:
                powers[old_index] = exponent
            term_repr = (
                "(" + repr(tuple(powers)) + ", "
                + coefficient_repr(values, True) + ")"
            )
            digest.update(term_repr.encode())

    for sequence_name in ("p", "q", "w"):
        for polynomial in data[sequence_name]:
            add_polynomial(polynomial)
    for tag, polynomial in data["constraints"]:
        digest.update(tag.encode())
        add_polynomial(polynomial)
    return digest.hexdigest()


def recover(state_path: Path, template_path: Path, out_path: Path,
            component: int) -> dict:
    if out_path.exists():
        raise FileExistsError(f"refusing to overwrite {out_path}")
    if os.environ.get("SYMPY_GROUND_TYPES") != "flint":
        raise RuntimeError("run with SYMPY_GROUND_TYPES=flint to match the audited digest")

    template = json.loads(template_path.read_text(encoding="utf-8"))
    data = json.loads(state_path.read_text(encoding="utf-8"))
    expected_component = f"F{component}=0"
    expected_tag = f"Salpha-component-F{component}"
    factor_hash = sha256(FACTOR_LOG)
    if (template.get("component") != expected_component
            or template.get("component_tag") != expected_tag
            or template.get("factor_log_sha256") != factor_hash):
        raise ValueError("template is not the verified receipt for this component")
    if data.get("format") != "f2-sparse-ring-v1" or data.get("branch") != "2+2":
        raise ValueError("expected an F2 2+2 sparse-ring checkpoint")
    compaction = data.get("compaction") or {}
    original_hash = compaction.get("source_checkpoint_sha256")
    if original_hash != template.get("original_source_checkpoint_sha256"):
        raise ValueError("partial state and template have different source provenance")

    tagged = [poly for tag, poly in data["constraints"] if tag == expected_tag]
    if len(tagged) != 1:
        raise ValueError("partial state must contain exactly one component tag")
    if (template.get("component_parameter_elimination") is not None
            and tagged[0]):
        raise ValueError("component parameter was reportedly eliminated but its row is nonzero")
    numeric_steps = [entry["step"] for entry in data["receipt"]
                     if isinstance(entry.get("step"), int)]
    if not numeric_steps:
        raise ValueError("partial state receipt has no numeric distance")

    term_count = sum(len(poly) for _tag, poly in data["constraints"])
    receipt = dict(template)
    receipt.update({
        "source_checkpoint": str(state_path),
        "source_checkpoint_sha256": sha256(state_path),
        "factor_log_sha256": factor_hash,
        "original_source_checkpoint_sha256": original_hash,
        "through_distance": max(numeric_steps),
        "parameter_count": data["nparameters"],
        "constraint_count": len(data["constraints"]),
        "constraint_term_count": term_count,
        "blocks": data["receipt"],
        "digest": checkpoint_digest(data),
        "recovered_from_partial_checkpoint": True,
        "scope": ("recovered exact partial checkpoint sidecar for continuation; "
                  "not an exclusion certificate"),
    })
    out_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return {
        "state": str(state_path),
        "sidecar": str(out_path),
        "component": expected_component,
        "through_distance": receipt["through_distance"],
        "parameters": receipt["parameter_count"],
        "constraints": receipt["constraint_count"],
        "constraint_terms": term_count,
        "state_sha256": receipt["source_checkpoint_sha256"],
        "digest": receipt["digest"],
        "component_row_is_zero": not tagged[0],
        "scope": receipt["scope"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--component", type=int, choices=(3, 4), required=True)
    args = parser.parse_args()
    print(json.dumps(recover(args.state, args.template, args.out, args.component),
                     indent=2), flush=True)


if __name__ == "__main__":
    main()
