"""Exact alpha-root tuned-jet audit for the depth-six F2 (75,125) state.

This deliberately checks only conditions that follow from the global K6--K8
polynomiality identities at alpha=(3+sqrt(5))/2.  It does not impose an
alpha-specific terminal-chain chart, nor does it claim full branch exclusion.
The source profile was computed from the exact Q-replayed depth-six checkpoint.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parent
PROFILE = ROOT / "f2_75_125_depth6_k8_alpha_branch_profile_agent_20261018.json"
STATE = ROOT / "f2_75_125_sparse_ring_lift_22_terminal_roots6_q_replayed_state.json"
REPLAY = ROOT / "f2_75_125_sparse_ring_lift_22_terminal_roots6_q_replayed.json"
OUTPUT = ROOT / "f2_75_125_alpha_tuned_rootwise_receipt_20261018.json"


def profile_polynomial(profile: dict, name: str) -> sp.Expr:
    """Decode sparse ring terms; coefficients are [sqrt(5), rational]."""
    out = sp.Integer(0)
    sqrt5 = sp.sqrt(5)
    for monomial, coefficient_pair in profile["equations"][name]:
        if len(coefficient_pair) == 1:
            coefficient = sp.Rational(coefficient_pair[0])
        elif len(coefficient_pair) == 2:
            coefficient = (sp.Rational(coefficient_pair[0]) * sqrt5
                           + sp.Rational(coefficient_pair[1]))
        else:
            raise ValueError(f"unexpected coefficient encoding in {name}")
        term = coefficient
        # The checkpoint ring begins t,p0,p1,...; these local reports omit t
        # after evaluation and therefore sparse index 2 denotes parameter p1.
        for generator_index, exponent in monomial:
            if generator_index < 2:
                raise ValueError(f"unexpected non-parameter index in {name}")
            term *= sp.Symbol(f"p{generator_index - 1}") ** exponent
        out += term
    return sp.expand(out)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    profile = json.loads(PROFILE.read_text(encoding="utf-8"))
    replay = json.loads(REPLAY.read_text(encoding="utf-8"))
    if profile["input_state"].replace("\\", "/").endswith(STATE.name) is False:
        raise AssertionError("profile does not identify the Q-replayed checkpoint")
    if replay["output_sha256"] != sha256(STATE):
        raise AssertionError("Q-replayed checkpoint hash differs from its receipt")

    S = profile_polynomial(profile, "S_alpha")
    D1 = profile_polynomial(profile, "Dprime_alpha")
    E1 = profile_polynomial(profile, "Eprime_alpha")
    K8 = profile_polynomial(profile, "K8_tuned_contact_alpha")

    # This specialization keeps S(alpha) nonzero while satisfying the extra
    # double-root K7 derivative conditions and the tuned K8 contact equation.
    p = {f"p{i}": sp.Symbol(f"p{i}") for i in
         (1, 2, 7, 8, 9, 16, 17, 18, 19, 20, 27, 28, 29, 30, 31, 32, 33)}
    values = {
        p["p1"]: sp.Integer(1), p["p2"]: sp.Integer(0),
        p["p8"]: sp.Integer(0), p["p9"]: sp.Integer(0),
        p["p17"]: sp.Integer(0), p["p18"]: sp.Integer(0),
        p["p19"]: sp.Integer(0), p["p20"]: sp.Integer(0),
        p["p28"]: sp.Integer(0), p["p29"]: sp.Integer(0),
        p["p30"]: sp.Integer(0), p["p31"]: sp.Integer(0),
        p["p32"]: sp.Integer(0), p["p33"]: sp.Integer(0),
    }
    values[p["p7"]] = sp.solve(D1.subs(values), p["p7"])[0]
    values[p["p16"]] = sp.solve(E1.subs(values), p["p16"])[0]
    values[p["p27"]] = sp.solve(K8.subs(values), p["p27"])[0]

    residuals = {
        "S_alpha": sp.simplify(S.subs(values)),
        "K7_Dprime_alpha": sp.simplify(D1.subs(values)),
        "K7_Eprime_alpha": sp.simplify(E1.subs(values)),
        "K8_contact_alpha": sp.simplify(K8.subs(values)),
    }
    assert residuals["S_alpha"] != 0
    assert all(residuals[key] == 0 for key in
               ("K7_Dprime_alpha", "K7_Eprime_alpha", "K8_contact_alpha"))

    receipt = {
        "root": "alpha=(3+sqrt(5))/2",
        "source_profile": PROFILE.name,
        "source_profile_sha256": sha256(PROFILE),
        "source_state": STATE.name,
        "source_state_sha256": sha256(STATE),
        "q_replay_receipt": {
            "saved_corrected_power_relations_recovered": all(
                row["saved_Q_slice_recovered"]
                for row in replay["corrected_power_replay"]
            ),
            "terminal_Q_positions_checked": replay["terminal_Q_positions_checked"],
            "retained_constraint_count": replay["retained_constraint_count"],
        },
        "rootwise_conditions": {
            "K7_derivative_conditions": ["D'(alpha)=0", "E'(alpha)=0"],
            "K8_contact": "(S*Delta/h - 81*U4)(alpha)=0, after the K7 tuned double-root orders",
            "K6_value_conditions_D_alpha_E_alpha": "identically zero in the source profile",
        },
        "exact_specialization": {
            "p1": "1", "p2": "0", "p8": "0", "p9": "0",
            "p17": "0", "p18": "0", "p19": "0", "p20": "0",
            "p28": "0", "p29": "0", "p30": "0", "p31": "0",
            "p32": "0", "p33": "0",
            "p7": str(sp.simplify(values[p["p7"]])),
            "p16": str(sp.simplify(values[p["p16"]])),
            "p27": str(sp.simplify(values[p["p27"]])),
            "all_other_free_parameters": "0 (they do not occur in these equations)",
        },
        "exact_residuals": {key: str(value) for key, value in residuals.items()},
        "conclusion": (
            "S(alpha) != 0 is compatible with the alpha-local K7 derivative and "
            "K8 contact equations in the depth-six Q-replayed family. Those "
            "equations add nonzero rootwise conditions, but they do not force "
            "S(alpha)=0. This is not a global K6-K8 ideal proof or a branch exclusion."
        ),
        "scope": "alpha-root necessary conditions only; no alpha-specific terminal chart",
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
