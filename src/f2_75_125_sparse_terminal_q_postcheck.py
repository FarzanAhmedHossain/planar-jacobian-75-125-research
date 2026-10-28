"""Exact Q-side terminal-halfspace audit from a sparse-lift checkpoint.

The checkpoint already satisfies the P terminal rows through its last band.
This reconstructs the same collapsed source shear u -> u+a*w+b*w^2 and
checks every Q terminal coefficient available in that checkpoint.
"""

from __future__ import annotations

import argparse
import json
from math import comb, factorial
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.rings import ring


def load_sparse(path: Path):
    data = json.loads(path.read_text(encoding="utf-8"))
    field = QQ.algebraic_field(sp.sqrt(5)) if data["quadratic"] else QQ
    R, *gens = ring(data["generators"], field)

    def decode(terms):
        result = {}
        for sparse, values in terms:
            powers = [0] * len(gens)
            for index, exponent in sparse:
                powers[index] = exponent
            if data["quadratic"]:
                value_expr = sum(
                    (sp.Rational(value) * sp.sqrt(5) ** (len(values) - 1 - i)
                     for i, value in enumerate(values)), sp.S.Zero)
                value = field.from_sympy(value_expr)
            else:
                value = field.convert(sp.Rational(values[0]))
            result[tuple(powers)] = value
        return R.from_dict(result)

    state = {
        key: decode(data[key]) for key in ("h", "v", "S", "B", "kappa")
    }
    state["kappa2"] = decode(data["kappa2"]) if "kappa2" in data else R.zero
    state["kappa3"] = decode(data["kappa3"]) if "kappa3" in data else R.zero
    state["constraints"] = [(tag, decode(poly))
                             for tag, poly in data.get("constraints", [])]
    state["uk"] = {int(k): decode(v) for k, v in data["uk"].items()}
    state["p"] = [decode(x) for x in data["p"]]
    state["q"] = [decode(x) for x in data["q"]]
    state["w"] = [decode(x) for x in data["w"]]
    state.update(R=R, gens=gens, field=field, data=data)
    return state


def coefficient(poly, degree, R):
    return R.from_dict({(0, *powers[1:]): value
                        for powers, value in poly.items()
                        if powers[0] == degree})


def local_slice(poly, lead, distance, bound, R):
    out = R.zero
    for powers, value in poly.items():
        exponent = 5 * powers[0] - lead + distance
        if exponent < 0:
            raise AssertionError("checkpoint monomial lies outside source support")
        for degree in range(min(exponent, bound) + 1):
            out += R.term_new((degree, *powers[1:]), value * comb(exponent, degree))
    return out


def centered_slice(series, lead, distance, bound, a, b, R, field):
    result = R.zero
    modulus = R.gens[0] ** (bound + 1)
    for i in range(distance + 1):
        deficit = distance - i
        for n in range(deficit // 2 + 1):
            m = deficit - 2 * n
            derivative = local_slice(series[i], lead, i, bound + m + n, R)
            for _ in range(m + n):
                derivative = derivative.diff(R.gens[0])
            a_power = R.one if m == 0 else a ** m
            b_power = R.one if n == 0 else b ** n
            result += derivative * a_power * b_power * field.convert(
                sp.Rational(1, factorial(m) * factorial(n)))
    return result.rem(modulus)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    s = load_sparse(args.checkpoint)
    R, field, p, q = s["R"], s["field"], s["p"], s["q"]
    if len(q) != len(p):
        raise AssertionError("P/Q checkpoint depths differ")
    depth = len(q) - 1

    Cpoly = local_slice(p[0], 15, 0, 6, R)
    zero = (0,) * len(R.gens)
    C = Cpoly[(6, *zero[1:])]
    if not C:
        raise AssertionError("terminal P leading coefficient vanished")
    A1 = local_slice(p[1], 15, 1, 5, R)
    a = -coefficient(A1, 5, R) / (6 * C)
    b = R.zero
    R2 = centered_slice(p, 15, 2, 5, a, b, R, field)
    b = -coefficient(R2, 5, R) / (6 * C)

    # Verify the full P-side terminal halfspace before using this checkpoint
    # as the locus on which the Q-side conditions are audited.  In particular,
    # this checks all 32 forbidden P positions through distance five, plus any
    # later terminal slices already carried by a deeper checkpoint.
    p_terminal_positions = 0
    p_terminal_positions_through5 = 0
    for distance in range(depth + 1):
        p_bound = (71 - 5 * distance) // 12
        if p_bound < 0:
            continue
        p_chart = centered_slice(p, 15, distance, p_bound, a, b, R, field)
        for degree in range(p_bound + 1):
            p_terminal_positions += 1
            if distance <= 5:
                p_terminal_positions_through5 += 1
            if coefficient(p_chart, degree, R):
                raise AssertionError(
                    f"saved P terminal coefficient is nonzero at ({distance},{degree})")

    checks = []
    for distance in range(depth + 1):
        bound = (119 - 5 * distance) // 12
        if bound < 0:
            continue
        chart = centered_slice(q, 25, distance, bound, a, b, R, field)
        positions = []
        nonzero = []
        for degree in range(bound + 1):
            value = coefficient(chart, degree, R)
            positions.append([distance, degree])
            if value:
                nonzero.append({"position": [distance, degree],
                                "terms": len(value),
                                "expression": str(value.as_expr())})
        checks.append({"distance": distance, "bound": bound,
                       "positions": len(positions), "nonzero": nonzero})

    report = {
        "branch": s["data"]["branch"],
        "depth": depth,
        "P_terminal_positions_checked": p_terminal_positions,
        "P_terminal_positions_checked_through5": p_terminal_positions_through5,
        "candidate_Q_terminal_positions": sum(x["positions"] for x in checks),
        "nonzero_positions": sum(len(x["nonzero"]) for x in checks),
        "checks": checks,
        "scope": "exact Q terminal coefficients on the checkpoint's P-side necessary locus",
    }
    output = json.dumps(report, indent=2)
    if args.receipt:
        args.receipt.write_text(output + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
