"""Substitute the exact second-root weight-49 relation in a sparse checkpoint.

The t=1 chart has already reduced S(t) to a two-parameter family.  The
second eligible-root chart requires S(alpha)=0, which gives
p2=((sqrt(5)-5)/5)*p1.  This script applies that exact specialization to
every encoded polynomial in a sparse checkpoint using only Python's standard
library, preserving all other equations and bands.

The output describes a necessary slice of the 2+2 system; it is not a
whole-branch exclusion.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path


Quad = tuple[Fraction, Fraction]  # a + b*sqrt(5)
ZERO: Quad = (Fraction(0), Fraction(0))
RATIO: Quad = (Fraction(-1), Fraction(1, 5))
ALPHA: Quad = (Fraction(3, 2), Fraction(1, 2))


def add(x: Quad, y: Quad) -> Quad:
    return (x[0] + y[0], x[1] + y[1])


def mul(x: Quad, y: Quad) -> Quad:
    return (x[0] * y[0] + 5 * x[1] * y[1],
            x[0] * y[1] + x[1] * y[0])


def power(x: Quad, exponent: int) -> Quad:
    if exponent < 0:
        raise ValueError("negative field exponent")
    result = (Fraction(1), Fraction(0))
    for _ in range(exponent):
        result = mul(result, x)
    return result


def parse_coefficient(values: list[str]) -> Quad:
    if len(values) == 1:
        return (Fraction(values[0]), Fraction(0))
    if len(values) == 2:
        # AlgebraicField.to_list stores descending powers: sqrt coefficient,
        # then rational coefficient.
        return (Fraction(values[1]), Fraction(values[0]))
    raise AssertionError(f"unexpected quadratic coefficient encoding: {values}")


def encode_coefficient(value: Quad) -> list[str]:
    rational, sqrt_coefficient = value
    if not sqrt_coefficient:
        return [str(rational)]
    return [str(sqrt_coefficient), str(rational)]


def transform_polynomial(terms: list, p1_index: int, p2_index: int) -> tuple[list, int, int]:
    collected: dict[tuple[tuple[int, int], ...], Quad] = {}
    for sparse, values in terms:
        powers = {index: exponent for index, exponent in sparse}
        exponent = powers.pop(p2_index, 0)
        if exponent:
            powers[p1_index] = powers.get(p1_index, 0) + exponent
        monomial = tuple(sorted(powers.items()))
        coefficient = mul(parse_coefficient(values), power(RATIO, exponent))
        collected[monomial] = add(collected.get(monomial, ZERO), coefficient)

    result = []
    for monomial, coefficient in sorted(collected.items()):
        if coefficient == ZERO:
            continue
        sparse = [[index, exponent] for index, exponent in monomial if exponent]
        result.append([sparse, encode_coefficient(coefficient)])
    return result, len(terms), len(result)


def transform_state(data: dict, p1_index: int, p2_index: int) -> tuple[int, int]:
    before = after = 0

    def apply(terms: list | None) -> list | None:
        nonlocal before, after
        if terms is None:
            return None
        transformed, old_count, new_count = transform_polynomial(terms, p1_index, p2_index)
        before += old_count
        after += new_count
        return transformed

    for key in ("h", "v", "S", "B", "kappa", "kappa2", "kappa3"):
        data[key] = apply(data.get(key))
    data["uk"] = {key: apply(value) for key, value in data["uk"].items()}
    for key in ("p", "q", "w"):
        data[key] = [apply(value) for value in data[key]]
    data["constraints"] = [[tag, apply(poly)] for tag, poly in data["constraints"]]
    return before, after


def verify_no_p2(data: dict, p2_index: int) -> None:
    polynomials = [data.get(key) for key in ("h", "v", "S", "B", "kappa", "kappa2", "kappa3")]
    polynomials += list(data["uk"].values()) + data["p"] + data["q"] + data["w"]
    polynomials += [poly for _, poly in data["constraints"]]
    for polynomial in polynomials:
        if polynomial is None:
            continue
        if any(any(index == p2_index for index, _ in sparse)
               for sparse, _ in polynomial):
            raise AssertionError("eliminated p2 remains in a polynomial")


def evaluate_s_alpha(data: dict, t_index: int, p1_index: int) -> Quad:
    result = ZERO
    for sparse, values in data["S"]:
        powers = {index: exponent for index, exponent in sparse}
        if powers.get(p1_index, 0) != 1 or any(
            index not in (t_index, p1_index) for index in powers
        ):
            raise AssertionError(f"unexpected S monomial after slicing: {powers}")
        result = add(
            result,
            mul(parse_coefficient(values), power(ALPHA, powers.get(t_index, 0))),
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=Path(
        "research/f2_75_125_sparse_ring_lift_22_terminal_roots7_gauged_state.json"))
    parser.add_argument("--output", type=Path, default=Path(
        "research/f2_75_125_sparse_ring_lift_22_terminal_roots7_gauged_Salpha_slice_state.json"))
    parser.add_argument("--receipt", type=Path, default=Path(
        "research/f2_75_125_dual_root_checkpoint_slice_receipt_20261018.json"))
    args = parser.parse_args()

    source_bytes = args.source.read_bytes()
    data = json.loads(source_bytes)
    if data.get("format") != "f2-sparse-ring-v1" or data.get("branch") != "2+2":
        raise AssertionError("expected an F2 2+2 sparse-ring checkpoint")
    if not data.get("quadratic"):
        raise AssertionError("expected coefficient field Q(sqrt(5))")
    p1_index = data["generators"].index("p1")
    p2_index = data["generators"].index("p2")

    before, after = transform_state(data, p1_index, p2_index)
    verify_no_p2(data, p2_index)
    s_alpha = evaluate_s_alpha(data, data["generators"].index("t"), p1_index)
    if s_alpha != ZERO:
        raise AssertionError(f"S(alpha) did not vanish after substitution: {s_alpha}")
    data["second_root_slice"] = {
        "relation": "p2=((sqrt(5)-5)/5)*p1",
        "equation": "S(alpha)=0",
        "alpha": "(3+sqrt(5))/2",
        "scope": "necessary parameter slice from the second eligible root chart",
        "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
    }

    output_bytes = (json.dumps(data, separators=(",", ":")) + "\n").encode()
    args.output.write_bytes(output_bytes)
    receipt = {
        "branch": "2+2",
        "source": str(args.source),
        "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "output": str(args.output),
        "output_sha256": hashlib.sha256(output_bytes).hexdigest(),
        "eliminated_parameter": "p2",
        "substitution": "p2=((sqrt(5)-5)/5)*p1",
        "stored_term_count_before": before,
        "stored_term_count_after": after,
        "p2_occurrences_after": 0,
        "S_alpha_residual": "0",
        "active_parameter_reserve_unchanged": data["nparameters"],
        "scope": "all saved necessary equations restricted exactly to S(alpha)=0; not a branch exclusion",
    }
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
