"""Extend the exact total-degree recurrence beyond the first 15 bands.

The support envelopes are the standardized F2 rectangles
    P: 0<=i<=15, 0<=j<=60, 5i-j<=15,
    Q: 0<=i<=25, 0<=j<=100, 5i-j<=25.

For a degree deficit d, the new P/Q coefficients enter the bracket linearly.
Nonresonant rows are solved exactly over Q. Resonant and unsupported rows are
retained as polynomial compatibility conditions in the free coefficients.
This is a necessary-system generator, not by itself an exclusion.
"""

from __future__ import annotations

import argparse
import json
import pickle
from collections import defaultdict
from pathlib import Path

import sympy as sp


def layer_keys(d: int, x_cap: int, y_cap: int) -> range:
    """Coefficient indices k for x^(x_cap-d+k)y^(y_cap-k)."""
    low = max(0, d - x_cap)
    high = min(y_cap, (5 * d) // 6)
    return range(low, high + 1) if low <= high else range(0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--through", type=int, default=20)
    parser.add_argument("--resume", type=Path)
    parser.add_argument("--checkpoint-dir", type=Path)
    parser.add_argument("--constraints-out", type=Path)
    args = parser.parse_args()
    if not 1 <= args.through <= 198:
        raise SystemExit("--through must be in [1,198]")

    if args.resume:
        with args.resume.open("rb") as checkpoint_file:
            checkpoint = pickle.load(checkpoint_file)
        if checkpoint.get("format_version") != 1:
            raise SystemExit("unsupported checkpoint format")
        start_d = int(checkpoint["through"])
        if start_d >= args.through:
            raise SystemExit("resume checkpoint already reaches --through")
        p_layers = checkpoint["p_layers"]
        q_layers = checkpoint["q_layers"]
        p_free = checkpoint["p_free"]
        q_free = checkpoint["q_free"]
        constraints = checkpoint["constraints"]
    else:
        start_d = 0
        p_layers = {0: {0: sp.Integer(1)}}
        q_layers = {0: {0: sp.Integer(1)}}
        p_free = []
        q_free = []
        constraints = []

    for d in range(start_d + 1, args.through + 1):
        p_now: dict[int, sp.Expr] = {}
        q_now: dict[int, sp.Expr] = {}
        if d <= 75:
            p_now = {k: sp.Integer(0) for k in layer_keys(d, 15, 60)}
        if d <= 125:
            q_now = {k: sp.Integer(0) for k in layer_keys(d, 25, 100)}

        # Coefficients from already chosen layer pairs with a+b=d.
        residual_terms: dict[int, list[sp.Expr]] = defaultdict(list)
        for a in range(1, min(75, d - 1) + 1):
            b = d - a
            if b < 1 or b > 125 or a not in p_layers or b not in q_layers:
                continue
            for pk, pc in p_layers[a].items():
                px, py = 15 - a + pk, 60 - pk
                for qk, qc in q_layers[b].items():
                    qx, qy = 25 - b + qk, 100 - qk
                    jac = px * qy - py * qx
                    if not jac:
                        continue
                    out_x = px + qx - 1
                    if 0 <= out_x <= 198 - d:
                        residual_terms[out_x].append(jac * pc * qc)
        residuals = {
            r: sp.expand(sum(terms, sp.Integer(0)))
            for r, terms in residual_terms.items()
        }

        # Map every new coefficient to its output x exponent.
        new_rows: dict[int, tuple[int, bool, bool]] = {}
        for k in set(p_now) | set(q_now):
            out_x = 39 - d + k
            if 0 <= out_x <= 198 - d:
                if out_x in new_rows:
                    raise AssertionError("two layer indices reached one output row")
                new_rows[out_x] = (k, k in p_now, k in q_now)

        if d == 198:
            # The constant Jacobian is allowed to be any nonzero scalar; it
            # is not an equation with a pre-normalized right-hand side.
            jacobian_constant = residuals.get(0, sp.Integer(0))
            free = sorted(jacobian_constant.free_symbols, key=str)
            poly = sp.Poly(jacobian_constant, *free) if free else sp.Poly(jacobian_constant)
            print(
                f"constant_jacobian_terms={len(poly.terms())} "
                f"constant_jacobian_degree={poly.total_degree()} "
                f"constant_jacobian_identically_zero={jacobian_constant == 0}"
            )
            break

        # Any old residual outside the new support is itself a necessary row.
        for out_x, residual in residuals.items():
            if out_x not in new_rows and residual != 0:
                constraints.append((d, out_x, residual, "unsupported"))

        for out_x, (k, has_p, has_q) in new_rows.items():
            residual = residuals.get(out_x, sp.Integer(0))
            factor = 4 * d - 5 * k
            if factor == 0:
                if residual != 0:
                    constraints.append((d, out_x, residual, "resonance"))
                if has_p:
                    sym = sp.Symbol(f"p_{d}_{k}")
                    p_now[k] = sym
                    p_free.append(sym)
                if has_q:
                    sym = sp.Symbol(f"q_{d}_{k}")
                    q_now[k] = sym
                    q_free.append(sym)
            elif has_p and has_q:
                psym = sp.Symbol(f"p_{d}_{k}")
                p_now[k] = psym
                q_now[k] = sp.expand(
                    (5 * psym - residual / (5 * factor)) / 3
                )
                p_free.append(psym)
            elif has_p:
                p_now[k] = sp.expand(residual / (25 * factor))
            elif has_q:
                q_now[k] = sp.expand(-residual / (15 * factor))

        p_layers[d] = p_now
        q_layers[d] = q_now

        rows_here = [c for c in constraints if c[0] == d]
        max_terms = 0
        max_degree = 0
        for _, _, polynomial, _ in rows_here:
            free = sorted(polynomial.free_symbols, key=str)
            if polynomial == 0:
                continue
            poly = sp.Poly(polynomial, *free) if free else sp.Poly(polynomial)
            max_terms = max(max_terms, len(poly.terms()))
            max_degree = max(max_degree, int(poly.total_degree()))
        print(
            f"d={d} P_terms={len(p_now)} Q_terms={len(q_now)} "
            f"new_free={len(p_free)+len(q_free)} "
            f"compat_rows={len(rows_here)} max_compat_terms={max_terms} "
            f"max_compat_degree={max_degree}"
        )

        if args.checkpoint_dir:
            args.checkpoint_dir.mkdir(parents=True, exist_ok=True)
            checkpoint_path = args.checkpoint_dir / f"frontier_state_d{d:03d}.pkl"
            temporary_path = checkpoint_path.with_suffix(".pkl.tmp")
            with temporary_path.open("wb") as checkpoint_file:
                pickle.dump(
                    {
                        "format_version": 1,
                        "through": d,
                        "p_layers": p_layers,
                        "q_layers": q_layers,
                        "p_free": p_free,
                        "q_free": q_free,
                        "constraints": constraints,
                    },
                    checkpoint_file,
                    protocol=pickle.HIGHEST_PROTOCOL,
                )
            temporary_path.replace(checkpoint_path)

    print(f"through={args.through}")
    print(f"free_P={len(p_free)} free_Q={len(q_free)}")
    print(f"compatibility_rows={len(constraints)}")
    for d, out_x, polynomial, kind in constraints:
        poly = sp.Poly(polynomial, *sorted(polynomial.free_symbols, key=str))
        print(
            f"constraint d={d} out_x={out_x} kind={kind} "
            f"terms={len(poly.terms())} total_degree={poly.total_degree()}"
        )

    if args.constraints_out:
        args.constraints_out.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "through": args.through,
            "free_P": len(p_free),
            "free_Q": len(q_free),
            "compatibility_rows": len(constraints),
            "constraints": [
                {
                    "deficit": d,
                    "output_x_exponent": out_x,
                    "kind": kind,
                    "exact_srepr": sp.srepr(polynomial),
                }
                for d, out_x, polynomial, kind in constraints
            ],
        }
        args.constraints_out.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
