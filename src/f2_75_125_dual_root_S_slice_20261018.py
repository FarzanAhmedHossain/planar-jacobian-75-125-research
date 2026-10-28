"""Replay the linear S(alpha)=0 cut in the F2 2+2 root charts.

This uses only the standard library and the exact saved d7 checkpoint.  The
checkpoint already incorporates the distinguished-root chart S(1)=0.  The
second eligible-root terminal chart forces S(alpha)=0.  Since the saved S is
linear in p1,p2, the result is an exact one-parameter slice of the weight-49
edge kernel.  This is a necessary slice, not a full-branch exclusion.
"""

from __future__ import annotations

import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path


Quad = tuple[Fraction, Fraction]  # a + b*sqrt(5)
ZERO: Quad = (Fraction(0), Fraction(0))
ONE: Quad = (Fraction(1), Fraction(0))
SQRT5: Quad = (Fraction(0), Fraction(1))
ALPHA: Quad = (Fraction(3, 2), Fraction(1, 2))


def add(x: Quad, y: Quad) -> Quad:
    return (x[0] + y[0], x[1] + y[1])


def mul(x: Quad, y: Quad) -> Quad:
    return (x[0] * y[0] + 5 * x[1] * y[1],
            x[0] * y[1] + x[1] * y[0])


def inv(x: Quad) -> Quad:
    norm = x[0] * x[0] - 5 * x[1] * x[1]
    if not norm:
        raise ZeroDivisionError("zero norm in Q(sqrt(5))")
    return (x[0] / norm, -x[1] / norm)


def neg(x: Quad) -> Quad:
    return (-x[0], -x[1])


def parse_coefficient(values: list[str], quadratic: bool) -> Quad:
    if quadratic:
        if len(values) == 1:
            return (Fraction(values[0]), Fraction(0))
        if len(values) != 2:
            raise AssertionError(f"expected Q(sqrt(5)) coefficient: {values}")
        # SparseLift serializes AlgebraicField.to_list() in descending powers.
        return (Fraction(values[1]), Fraction(values[0]))
    return (Fraction(values[0]), Fraction(0))


def fmt(x: Quad) -> str:
    a, b = x
    return f"({a}) + ({b})*sqrt(5)"


def evaluate(coefficients: dict[int, Quad], point: Quad) -> Quad:
    if not coefficients:
        return ZERO
    result = ZERO
    for degree in range(max(coefficients), -1, -1):
        result = mul(result, point)
        result = add(result, coefficients.get(degree, ZERO))
    return result


def main() -> None:
    source = Path(sys.argv[1] if len(sys.argv) > 1 else
                  "research/f2_75_125_sparse_ring_lift_22_terminal_roots7_gauged_state.json")
    data_bytes = source.read_bytes()
    data = json.loads(data_bytes)
    if data["branch"] != "2+2" or not data["quadratic"]:
        raise AssertionError("expected the exact Q(sqrt(5)) 2+2 checkpoint")
    names = data["generators"]

    # Read S(t) as a coefficient dictionary in the active parameters.
    by_parameter: dict[str, dict[int, Quad]] = {}
    for sparse, values in data["S"]:
        powers = {names[index]: exponent for index, exponent in sparse}
        parameter_names = [name for name in powers if name != "t"]
        if len(parameter_names) != 1 or powers[parameter_names[0]] != 1:
            raise AssertionError(f"S is not linear in one active parameter: {powers}")
        parameter = parameter_names[0]
        degree = powers.get("t", 0)
        coefficients = by_parameter.setdefault(parameter, {})
        if degree in coefficients:
            raise AssertionError("duplicate sparse S coefficient")
        coefficients[degree] = parse_coefficient(values, True)

    if set(by_parameter) != {"p1", "p2"}:
        raise AssertionError(f"unexpected active S parameters: {set(by_parameter)}")
    s1 = {name: evaluate(poly, ONE) for name, poly in by_parameter.items()}
    if any(value != ZERO for value in s1.values()):
        raise AssertionError(f"saved distinguished chart does not give S(1)=0: {s1}")
    sa = {name: evaluate(poly, ALPHA) for name, poly in by_parameter.items()}
    if any(value == ZERO for value in sa.values()):
        raise AssertionError("second-root evaluation unexpectedly loses a parameter")

    ratio = neg(mul(sa["p1"], inv(sa["p2"])))
    expected = (Fraction(-1), Fraction(1, 5))
    if ratio != expected:
        raise AssertionError(f"unexpected p2/p1 ratio: {fmt(ratio)}")

    receipt = {
        "source": str(source),
        "source_sha256": hashlib.sha256(data_bytes).hexdigest(),
        "branch": "2+2",
        "checkpoint_depth": len(data["p"]) - 1,
        "active_parameter_count_before_second_root_cut": data["nparameters"],
        "saved_S_parameters": sorted(by_parameter),
        "S1_coefficients": {name: fmt(value) for name, value in s1.items()},
        "Salpha_coefficients": {name: fmt(value) for name, value in sa.items()},
        "second_root_constraint": "S(alpha)=0, alpha=(3+sqrt(5))/2",
        "exact_slice": "p2=((sqrt(5)-5)/5)*p1",
        "scope": "one exact linear cut in the weight-49 kernel; not a full-branch exclusion",
    }
    output = Path("research/f2_75_125_dual_root_S_slice_receipt_20261018.json")
    output.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
