"""Verify that the newly extracted d19 low-support rows are redundant.

Every coefficient component is divided exactly by the fourth generator of
the previously computed exact support-12 Groebner basis.  The surrounding
receipts bind this calculation to the d19 checkpoint and to the exact reuse
of its support-12 equations from d18.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import sympy as sp
from flint import fmpq
from sympy.polys.domains import QQ


ROOT = Path("research")
SOURCE = ROOT / "f2_75_125_sparse_ring_lift_211_roots19_u1zero_streamed_state.json"
SMALL = ROOT / "f2_75_125_211_new_small_support_equations_20261025.json"
SUPPORT12_GB = ROOT / "f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_gb_4gen_weighted_flint.json"
SUPPORT12_INPUT = ROOT / "f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_input.json"
SUPPORT12_REUSE = ROOT / "f2_75_125_211_verify_d19_support12_reuse_20261025.json"
OUTPUT = ROOT / "f2_75_125_211_new_small_support_redundancy_20261025.json"
EXPECTED_SOURCE_SHA256 = "0f293c437b1cc18f379e3a8d0bb2e98647616e062fa321e91cfdc45ebe67926d"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if sha256(SOURCE) != EXPECTED_SOURCE_SHA256:
        raise ValueError("d19 checkpoint hash mismatch")

    small = json.loads(SMALL.read_text(encoding="utf-8"))
    if small.get("source_sha256") != EXPECTED_SOURCE_SHA256:
        raise ValueError("small-support extraction is not bound to d19")
    if small.get("source_sha256_rechecked") != EXPECTED_SOURCE_SHA256:
        raise ValueError("small-support source hash was not rechecked")
    if small.get("reduced_support_histogram_through_size17") != {
        "12": 142, "13": 1, "17": 2
    }:
        raise ValueError("unexpected low-support histogram")

    gb_report = json.loads(SUPPORT12_GB.read_text(encoding="utf-8"))
    support12_input = json.loads(SUPPORT12_INPUT.read_text(encoding="utf-8"))
    reuse = json.loads(SUPPORT12_REUSE.read_text(encoding="utf-8"))
    if gb_report.get("branch") != "2+1+1" or gb_report.get("unit_ideal"):
        raise ValueError("expected the saved nonunit 2+1+1 support-12 basis")
    if (gb_report.get("equation_count") != 4
            or gb_report.get("source_equation_count") != 142
            or not gb_report.get("exact_linear_reduction_verified")):
        raise ValueError("support-12 basis lacks its exact four-generator receipt")
    if support12_input.get("equation_count") != 142:
        raise ValueError("unexpected d18 support-12 input")
    if not reuse.get("same_142_equations_and_full_support12_group"):
        raise ValueError("d19 support-12 equations are not verified as the d18 group")

    research = Path(__file__).resolve().parent
    sys.path.insert(0, str(research))
    from f2_75_125_flint_mpoly_backend import FlintRing

    equation_records = small["equations"]
    old_variables = [
        name for name in equation_records[1]["support"]
        if name not in {"p29", "p36", "p37", "p38", "p39"}
    ]
    if set(old_variables) != set(support12_input["selected_variables"]):
        raise ValueError("the low-support equations use an unexpected old-variable ring")

    ring = FlintRing(old_variables)
    symbols = sp.symbols(" ".join(old_variables))
    symbol_map = dict(zip(old_variables, symbols))
    basis_text = gb_report["basis"][3]["polynomial"]
    basis_poly = sp.Poly(sp.sympify(basis_text, locals=symbol_map),
                         *symbols, domain=QQ)
    divisor = ring.from_dict(dict(basis_poly.terms()))
    divisor_text_hash = hashlib.sha256(basis_text.encode("utf-8")).hexdigest()

    verified_equations = []
    for equation in equation_records:
        fresh = (["p16"] if equation["tag"] == "16:1:t^6" else
                 ["p29", "p36", "p37", "p38", "p39"])
        components = [{} for _ in range(len(fresh) + 1)]
        for monomial, coefficient in equation["terms"]:
            powers = dict(monomial)
            present = [name for name in fresh if powers.get(name, 0)]
            if len(present) > 1 or any(powers.get(name, 0) > 1 for name in fresh):
                raise ValueError(f"fresh parameters are not affine in {equation['tag']}")
            component = fresh.index(present[0]) + 1 if present else 0
            exponents = tuple(powers.get(name, 0) for name in old_variables)
            components[component][exponents] = (
                components[component].get(exponents, fmpq(0))
                + fmpq(coefficient)
            )

        quotients = []
        for component_index, terms in enumerate(components):
            polynomial = ring.from_prevalidated_dict({
                exponents: coefficient
                for exponents, coefficient in terms.items()
                if coefficient
            })
            quotient, remainder = divmod(polynomial.data, divisor.data)
            if not remainder.is_zero():
                raise ArithmeticError(
                    f"component {component_index} of {equation['tag']} "
                    "is not divisible by the support-12 basis generator"
                )
            quotients.append({
                "parameter": "1" if component_index == 0 else fresh[component_index - 1],
                "input_term_count": len(polynomial),
                "quotient": str(quotient),
                "quotient_term_count": len(quotient),
                "remainder_zero": True,
            })
        verified_equations.append({
            "tag": equation["tag"],
            "full_support": equation["support"],
            "input_term_count": equation["term_count"],
            "fresh_parameters": fresh,
            "all_coefficient_components_divisible": True,
            "components": quotients,
        })

    result = {
        "scope": (
            "Exact componentwise divisibility proof that each of the three newly "
            "extracted support-13/support-17 equations lies in the extension of "
            "the previously computed support-12 ideal. This adds no restriction "
            "to that necessary-system slice and excludes no larger F2 branch."
        ),
        "d19_checkpoint_sha256": EXPECTED_SOURCE_SHA256,
        "small_support_equations_sha256": sha256(SMALL),
        "support12_groebner_report_sha256": sha256(SUPPORT12_GB),
        "support12_input_sha256": sha256(SUPPORT12_INPUT),
        "support12_reuse_receipt_sha256": sha256(SUPPORT12_REUSE),
        "support12_reuse_verified": True,
        "support12_equations": 142,
        "support12_ideal_unit": False,
        "divisor_basis_index_zero_based": 3,
        "divisor_polynomial_sha256": divisor_text_hash,
        "divisor_term_count": len(divisor),
        "divisor_total_degree": int(divisor.data.total_degree()),
        "low_support_equations_verified": len(verified_equations),
        "equations": verified_equations,
        "all_new_low_support_equations_redundant": len(verified_equations) == 3,
        "whole_2_plus_1_plus_1_excluded_percent": 0.0,
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "divisor_term_count": result["divisor_term_count"],
        "divisor_total_degree": result["divisor_total_degree"],
        "low_support_equations_verified": result["low_support_equations_verified"],
        "all_new_low_support_equations_redundant": result[
            "all_new_low_support_equations_redundant"
        ],
        "output": str(OUTPUT),
    }, indent=2))


if __name__ == "__main__":
    main()
