"""Extend the exact F2 2+2 lift with the second eligible double-root chart.

The base sparse lift imposes the F2 terminal chart at t=1.  In the 2+2
configuration h=(t-1)^2(t-alpha)^2, alpha=(3+sqrt(5))/2 is also a double
root, and every corresponding root of h(z^5)^3 meets the GGV first-shear
multiplicity cutoff.  This driver retains the original chart and imposes the
same necessary terminal support conditions in the alpha-root chart.

For the latter, write z=theta*(1+u), theta^5=alpha.  Every monomial in a
fixed band has a common nonzero theta factor times alpha^i*(1+u)^e, so the
vanishing equations are defined over Q(sqrt(5)); no degree-10 field is
needed.  This is still a necessary-condition lift, not an exclusion proof.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from math import comb, factorial
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ

from f2_75_125_sparse_ring_lift import SparseLift


def coeff_at(lift: SparseLift, polynomial, degree: int):
    return lift.ring.from_dict({
        (0, *powers[1:]): coefficient
        for powers, coefficient in polynomial.items()
        if powers[0] == degree
    })


def alpha_slice(lift: SparseLift, band, distance: int, bound: int,
                base_degree: int, alpha, *, q_support: bool = False):
    """Expand one band at z=theta*(1+u), omitting a common nonzero factor."""
    out = lift.ring.zero
    if q_support:
        low = -(-(25 - distance) // 5)
        high = min(25, (125 - distance) // 5)
    alpha_powers = [lift.field.one]
    for _ in range(max((powers[0] for powers in band), default=0)):
        alpha_powers.append(alpha_powers[-1] * alpha)
    for powers, coefficient in band.items():
        index = powers[0]
        if q_support and not (low <= index <= high):
            continue
        exponent = 5 * index - base_degree + distance
        if exponent < 0:
            raise AssertionError("source monomial fell outside its F2 support")
        scale = coefficient * alpha_powers[index]
        for degree in range(min(exponent, bound) + 1):
            monomial = (degree, *powers[1:])
            out += lift.ring.term_new(
                monomial, scale * comb(exponent, degree)
            )
    return out


def centered_alpha_slice(lift: SparseLift, bands, distance: int, bound: int,
                         center_a, center_b, base_degree: int, alpha,
                         *, q_support: bool = False):
    result = lift.ring.zero
    modulus = lift.t ** (bound + 1)
    for i in range(distance + 1):
        deficit = distance - i
        for n in range(deficit // 2 + 1):
            m = deficit - 2 * n
            local = alpha_slice(
                lift, bands[i], i, bound + m + n, base_degree, alpha,
                q_support=q_support,
            )
            for _ in range(m + n):
                local = local.diff(lift.t)
            # PolyElement deliberately rejects 0**0; the coefficient of an
            # absent center term uses the multiplicative identity instead.
            scale_a = lift.ring.one if m == 0 else center_a**m
            scale_b = lift.ring.one if n == 0 else center_b**n
            scale = scale_a * scale_b
            scale *= lift.field.convert(sp.Rational(1, factorial(m) * factorial(n)))
            result += local * scale
    return result.rem(modulus)


def initialize_alpha_chart(lift: SparseLift):
    sqrt5 = sp.sqrt(5)
    alpha = lift.field.from_sympy((sp.Integer(3) + sqrt5) / 2)
    p0 = alpha_slice(lift, lift.p[0], 0, 6, 15, alpha)
    leading = coeff_at(lift, p0, 6)
    zero = (0,) * len(lift.ring.gens)
    if set(leading) != {zero}:
        raise AssertionError("alpha-chart leading u^6 coefficient is not scalar")
    c0 = leading[zero]
    if not c0:
        raise AssertionError("alpha-chart leading coefficient vanished")

    first = alpha_slice(lift, lift.p[1], 1, 5, 15, alpha)
    center_a = -coeff_at(lift, first, 5) / (6 * c0)
    second_at_a = centered_alpha_slice(
        lift, lift.p, 2, 5, center_a, lift.ring.zero, 15, alpha
    )
    center_b = -coeff_at(lift, second_at_a, 5) / (6 * c0)
    return alpha, c0, center_a, center_b


def append_existing_alpha_rows(lift: SparseLift, alpha, center_a, center_b):
    """Retain alpha-chart P/Q rows already present in the resumed checkpoint."""
    tags = {tag for tag, _ in lift.constraints}
    added = []
    depth = min(len(lift.p), len(lift.q)) - 1
    for distance in range(depth + 1):
        p_bound = (71 - 5 * distance) // 12
        if p_bound >= 0:
            row = centered_alpha_slice(
                lift, lift.p, distance, p_bound, center_a, center_b, 15, alpha
            )
            for degree in range(p_bound + 1):
                tag = f"terminal-alpha-P-distance{distance}:u^{degree}"
                if tag in tags:
                    continue
                equation = coeff_at(lift, row, degree)
                if equation:
                    lift.constraints.append((tag, equation))
                    added.append((tag, equation))
        q_bound = (119 - 5 * distance) // 12
        if q_bound >= 0:
            row = centered_alpha_slice(
                lift, lift.q, distance, q_bound, center_a, center_b, 25,
                alpha, q_support=True,
            )
            for degree in range(q_bound + 1):
                tag = f"terminal-alpha-Q-distance{distance}:u^{degree}"
                if tag in tags:
                    continue
                equation = coeff_at(lift, row, degree)
                if equation:
                    lift.constraints.append((tag, equation))
                    added.append((tag, equation))
    return added


def impose_next_alpha_rows(lift: SparseLift, distance: int, alpha,
                           center_a, center_b):
    p_bound = (71 - 5 * distance) // 12
    if p_bound >= 0:
        expression = centered_alpha_slice(
            lift, lift.p, distance, p_bound, center_a, center_b, 15, alpha
        )
        lift.impose_congruence(
            expression, lift.t ** (p_bound + 1), lift.last_free,
            f"terminal-alpha-P-distance{distance}",
        )
        center_a, center_b = lift.alpha_center_a, lift.alpha_center_b
    q_bound = (119 - 5 * distance) // 12
    if q_bound >= 0:
        expression = centered_alpha_slice(
            lift, lift.q, distance, q_bound, center_a, center_b, 25, alpha,
            q_support=True,
        )
        lift.impose_congruence(
            expression, lift.t ** (q_bound + 1), lift.last_free,
            f"terminal-alpha-Q-distance{distance}",
        )
    return center_a, center_b


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resume", type=Path, required=True)
    parser.add_argument("--through", type=int, default=14)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()

    source_hash = hashlib.sha256(args.resume.read_bytes()).hexdigest()
    # This source checkpoint's original t=1 terminal rows were checked in its
    # producing run.  Skip only the base loader's duplicate 42-minute replay.
    lift = SparseLift.from_checkpoint(
        args.resume, verify_saved_terminal_rows=False
    )
    if lift.tag != "2+2" or lift.field == QQ:
        raise ValueError("the alpha-root terminal chart is for the 2+2 branch")
    alpha, c0, center_a, center_b = initialize_alpha_chart(lift)
    added = append_existing_alpha_rows(lift, alpha, center_a, center_b)
    initial_added = len(added)
    print(json.dumps({
        "stage": "alpha-root backfill",
        "depth": len(lift.p) - 1,
        "new_rows": initial_added,
        "new_terms": sum(len(poly) for _, poly in added),
        "parameters": lift.nparameters,
    }), flush=True)

    # The root-specific center is determined by bands 0, 1, and 2.  Later
    # affine substitutions must update it just like the distinguished chart.
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
        "component": ((lift.compaction or {}).get("exact_graph_reduction") or {}).get(
            "component"
        ),
        "component_graph_reduction": (lift.compaction or {}).get(
            "exact_graph_reduction"
        ),
        "root": "alpha=(3+sqrt(5))/2; z=theta*(1+u), theta^5=alpha",
        "root_multiplicity": 6,
        "first_shear_cutoff": "15/4",
        "initial_checkpoint": str(args.resume),
        "initial_checkpoint_sha256": source_hash,
        "saved_t1_terminal_rows_rechecked_by_loader": False,
        "through_distance": args.through,
        "alpha_leading_u6": str(lift.field.to_sympy(c0)),
        "backfilled_alpha_rows": initial_added,
        "parameter_count": lift.nparameters,
        "constraint_count": len(lift.constraints),
        "constraint_term_count": sum(len(poly) for _, poly in lift.constraints),
        "constraint_tags": [tag for tag, _ in lift.constraints],
        "blocks": lift.receipt,
        "digest": lift.digest(),
        "elapsed_seconds": time.monotonic() - started,
        "scope": "necessary dual-root support system; not an exclusion certificate",
    }
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt), flush=True)


if __name__ == "__main__":
    main()
