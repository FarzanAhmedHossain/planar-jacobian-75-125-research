"""FLINT-backed exact continuation of the F2 total-degree recurrence.

This mirrors f2_75_125_total_degree_frontier_extend.py, but performs sparse
multivariate arithmetic with python-flint.  Checkpoints contain plain Python
term dictionaries because FLINT multivariate polynomials are not pickleable.
The generated equations are necessary conditions, not a branch exclusion.
"""

from __future__ import annotations

import argparse
import json
import pickle
from collections import defaultdict
from pathlib import Path
from typing import Any

import sympy as sp
from flint import fmpq, fmpq_mpoly_ctx


def layer_keys(d: int, x_cap: int, y_cap: int) -> range:
    low = max(0, d - x_cap)
    high = min(y_cap, (5 * d) // 6)
    return range(low, high + 1) if low <= high else range(0)


def mpoly_data(poly: Any) -> list[tuple[tuple[int, ...], int, int]]:
    return [
        (tuple(map(int, exponents)), int(coefficient.numerator), int(coefficient.denominator))
        for exponents, coefficient in poly.to_dict().items()
    ]


def mpoly_from_data(ctx: Any, data: list[tuple[tuple[int, ...], int, int]]) -> Any:
    return ctx.from_dict(
        {tuple(exponents): fmpq(numerator, denominator) for exponents, numerator, denominator in data}
    )


def make_variable_names(old_checkpoint: dict[str, Any], through: int) -> list[str]:
    names = {str(symbol) for symbol in old_checkpoint["p_free"] + old_checkpoint["q_free"]}
    for d in range(int(old_checkpoint["through"]) + 1, through + 1):
        if d <= 75:
            names.update(f"p_{d}_{k}" for k in layer_keys(d, 15, 60))
        if d <= 125:
            names.update(f"q_{d}_{k}" for k in layer_keys(d, 25, 100))
    return sorted(names)


def sympy_to_mpoly(expr: sp.Expr, ctx: Any, variable_index: dict[str, int]) -> Any:
    if expr == 0:
        return ctx.from_dict({})
    symbols = sorted(expr.free_symbols, key=str)
    if not symbols:
        rational = sp.Rational(expr)
        return ctx.from_dict({(0,) * len(variable_index): fmpq(int(rational.p), int(rational.q))})
    polynomial = sp.Poly(expr, *symbols, domain=sp.QQ)
    result: dict[tuple[int, ...], fmpq] = {}
    nvars = len(variable_index)
    positions = [variable_index[str(symbol)] for symbol in symbols]
    for powers, coefficient in polynomial.terms():
        exponents = [0] * nvars
        for position, power in zip(positions, powers):
            exponents[position] = int(power)
        result[tuple(exponents)] = fmpq(int(coefficient.p), int(coefficient.q))
    return ctx.from_dict(result)


def is_zero(poly: Any) -> bool:
    return not poly.to_dict()


def term_count(poly: Any) -> int:
    return len(poly.to_dict())


def total_degree(poly: Any) -> int:
    terms = poly.to_dict()
    return max((sum(map(int, powers)) for powers in terms), default=0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--through", type=int, default=26)
    parser.add_argument("--resume", type=Path, help="Original SymPy checkpoint")
    parser.add_argument("--resume-flint", type=Path, help="FLINT term-dictionary checkpoint")
    parser.add_argument("--checkpoint-dir", type=Path)
    parser.add_argument("--constraints-out", type=Path)
    args = parser.parse_args()
    if not 1 <= args.through <= 198:
        raise SystemExit("--through must be in [1,198]")
    if bool(args.resume) == bool(args.resume_flint):
        raise SystemExit("provide exactly one of --resume or --resume-flint")

    if args.resume_flint:
        print(f"loading FLINT checkpoint {args.resume_flint}", flush=True)
        with args.resume_flint.open("rb") as stream:
            checkpoint = pickle.load(stream)
        if checkpoint.get("format_version") != 1 or checkpoint.get("backend") != "python-flint":
            raise SystemExit("unsupported FLINT checkpoint format")
        start_d = int(checkpoint["through"])
        if start_d >= args.through:
            raise SystemExit("resume checkpoint already reaches --through")
        variable_names = checkpoint["variable_names"]
        if any(
            name not in set(variable_names)
            for d in range(start_d + 1, args.through + 1)
            for name in (
                [f"p_{d}_{k}" for k in layer_keys(d, 15, 60)] if d <= 75 else []
            )
            + ([f"q_{d}_{k}" for k in layer_keys(d, 25, 100)] if d <= 125 else [])
        ):
            raise SystemExit("checkpoint lacks variables needed for the requested --through")
        p_data = checkpoint["p_layers"]
        q_data = checkpoint["q_layers"]
        p_free = list(checkpoint["p_free"])
        q_free = list(checkpoint["q_free"])
        constraint_data = checkpoint["constraints"]
        ctx = fmpq_mpoly_ctx.get(tuple(variable_names))
        p_layers = {
            int(d): {int(k): mpoly_from_data(ctx, data) for k, data in layer.items()}
            for d, layer in p_data.items()
        }
        q_layers = {
            int(d): {int(k): mpoly_from_data(ctx, data) for k, data in layer.items()}
            for d, layer in q_data.items()
        }
        constraints = [
            (int(d), int(out_x), mpoly_from_data(ctx, data), str(kind))
            for d, out_x, data, kind in constraint_data
        ]
    else:
        print(f"loading SymPy checkpoint {args.resume}", flush=True)
        with args.resume.open("rb") as stream:
            source = pickle.load(stream)
        if source.get("format_version") != 1:
            raise SystemExit("unsupported source checkpoint format")
        start_d = int(source["through"])
        if start_d >= args.through:
            raise SystemExit("resume checkpoint already reaches --through")
        variable_names = make_variable_names(source, args.through)
        ctx = fmpq_mpoly_ctx.get(tuple(variable_names))
        variable_index = {name: i for i, name in enumerate(variable_names)}
        print(
            f"loaded through={start_d}; converting {len(source['p_free']) + len(source['q_free'])} existing free variables "
            f"to FLINT in a {len(variable_names)}-variable context",
            flush=True,
        )
        convert_cache: dict[int, tuple[sp.Expr, Any]] = {}

        def convert(expr: sp.Expr) -> Any:
            key = id(expr)
            cached = convert_cache.get(key)
            if cached is not None and cached[0] is expr:
                return cached[1]
            result = sympy_to_mpoly(expr, ctx, variable_index)
            convert_cache[key] = (expr, result)
            return result

        p_layers = {
            int(d): {int(k): convert(expr) for k, expr in layer.items()}
            for d, layer in source["p_layers"].items()
        }
        q_layers = {
            int(d): {int(k): convert(expr) for k, expr in layer.items()}
            for d, layer in source["q_layers"].items()
        }
        p_free = [str(symbol) for symbol in source["p_free"]]
        q_free = [str(symbol) for symbol in source["q_free"]]
        constraints = [
            (int(d), int(out_x), convert(expr), str(kind))
            for d, out_x, expr, kind in source["constraints"]
        ]
        print(
            f"converted layers: P={sum(map(len, p_layers.values()))} Q={sum(map(len, q_layers.values()))}; "
            f"constraints={len(constraints)}",
            flush=True,
        )

    variable_index = {name: i for i, name in enumerate(variable_names)}
    nvars = len(variable_names)
    zero = ctx.from_dict({})
    variable_cache: dict[str, Any] = {}

    def variable(name: str) -> Any:
        if name not in variable_cache:
            exponents = [0] * nvars
            exponents[variable_index[name]] = 1
            variable_cache[name] = ctx.from_dict({tuple(exponents): fmpq(1)})
        return variable_cache[name]

    print(f"continuing deficits {start_d + 1}..{args.through}", flush=True)
    for d in range(start_d + 1, args.through + 1):
        p_now: dict[int, Any] = {}
        q_now: dict[int, Any] = {}
        if d <= 75:
            p_now = {k: zero for k in layer_keys(d, 15, 60)}
        if d <= 125:
            q_now = {k: zero for k in layer_keys(d, 25, 100)}

        residuals: dict[int, Any] = {}
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
                        contribution = jac * pc * qc
                        residuals[out_x] = residuals.get(out_x, zero) + contribution
        residuals = {out_x: poly for out_x, poly in residuals.items() if not is_zero(poly)}

        new_rows: dict[int, tuple[int, bool, bool]] = {}
        for k in set(p_now) | set(q_now):
            out_x = 39 - d + k
            if 0 <= out_x <= 198 - d:
                if out_x in new_rows:
                    raise AssertionError("two layer indices reached one output row")
                new_rows[out_x] = (k, k in p_now, k in q_now)

        if d == 198:
            constant = residuals.get(0, zero)
            print(
                f"constant_jacobian_terms={term_count(constant)} "
                f"constant_jacobian_degree={total_degree(constant)} "
                f"constant_jacobian_identically_zero={is_zero(constant)}",
                flush=True,
            )
            break

        for out_x, residual in residuals.items():
            if out_x not in new_rows:
                constraints.append((d, out_x, residual, "unsupported"))

        for out_x, (k, has_p, has_q) in new_rows.items():
            residual = residuals.get(out_x, zero)
            factor = 4 * d - 5 * k
            if factor == 0:
                if not is_zero(residual):
                    constraints.append((d, out_x, residual, "resonance"))
                if has_p:
                    name = f"p_{d}_{k}"
                    p_now[k] = variable(name)
                    p_free.append(name)
                if has_q:
                    name = f"q_{d}_{k}"
                    q_now[k] = variable(name)
                    q_free.append(name)
            elif has_p and has_q:
                name = f"p_{d}_{k}"
                p_symbol = variable(name)
                p_now[k] = p_symbol
                q_now[k] = (5 * p_symbol - residual / (5 * factor)) / 3
                p_free.append(name)
            elif has_p:
                p_now[k] = residual / (25 * factor)
            elif has_q:
                q_now[k] = -residual / (15 * factor)

        p_layers[d] = p_now
        q_layers[d] = q_now
        rows_here = [constraint for constraint in constraints if constraint[0] == d]
        print(
            f"d={d} P_terms={len(p_now)} Q_terms={len(q_now)} "
            f"new_free={len(p_free) + len(q_free)} compat_rows={len(rows_here)} "
            f"max_compat_terms={max((term_count(row[2]) for row in rows_here), default=0)} "
            f"max_compat_degree={max((total_degree(row[2]) for row in rows_here), default=0)}",
            flush=True,
        )

        if args.checkpoint_dir:
            args.checkpoint_dir.mkdir(parents=True, exist_ok=True)
            checkpoint_path = args.checkpoint_dir / f"frontier_flint_state_d{d:03d}.pkl"
            temporary_path = checkpoint_path.with_suffix(".pkl.tmp")
            encoded = {
                "format_version": 1,
                "backend": "python-flint",
                "through": d,
                "variable_names": variable_names,
                "p_layers": {
                    int(degree): {int(k): mpoly_data(poly) for k, poly in layer.items()}
                    for degree, layer in p_layers.items()
                },
                "q_layers": {
                    int(degree): {int(k): mpoly_data(poly) for k, poly in layer.items()}
                    for degree, layer in q_layers.items()
                },
                "p_free": p_free,
                "q_free": q_free,
                "constraints": [
                    (degree, out_x, mpoly_data(poly), kind)
                    for degree, out_x, poly, kind in constraints
                ],
            }
            with temporary_path.open("wb") as stream:
                pickle.dump(encoded, stream, protocol=pickle.HIGHEST_PROTOCOL)
            temporary_path.replace(checkpoint_path)
            print(f"checkpoint={checkpoint_path} bytes={checkpoint_path.stat().st_size}", flush=True)

    print(f"through={args.through} free_P={len(p_free)} free_Q={len(q_free)} constraints={len(constraints)}")
    for d, out_x, polynomial, kind in constraints:
        print(
            f"constraint d={d} out_x={out_x} kind={kind} "
            f"terms={term_count(polynomial)} total_degree={total_degree(polynomial)}"
        )

    if args.constraints_out:
        args.constraints_out.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "through": args.through,
            "free_P": len(p_free),
            "free_Q": len(q_free),
            "compatibility_rows": len(constraints),
            "variable_names": variable_names,
            "constraints": [
                {
                    "deficit": d,
                    "output_x_exponent": out_x,
                    "kind": kind,
                    "terms": [
                        {"exponents": list(exponents), "numerator": numerator, "denominator": denominator}
                        for exponents, numerator, denominator in mpoly_data(polynomial)
                    ],
                }
                for d, out_x, polynomial, kind in constraints
            ],
        }
        args.constraints_out.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
