"""Exact positive control for the first 60 y-adic coefficient rows.

The script constructs rational P- and Q-layers satisfying [P,Q]=1 through
output y-degree 59, with the entire (3,5) top-edge data fixed by
f=(1+t^4)^3 and g=(1+t^4)^5. It is a finite-truncation witness, not a Keller
pair or an exclusion proof.
"""

from __future__ import annotations

import json

import sympy as sp


def coefficient(poly: sp.Expr, var: sp.Symbol, degree: int) -> sp.Expr:
    return sp.Poly(sp.expand(poly), var).nth(degree)


def main() -> None:
    x, t = sp.symbols("x t")
    p0 = x**3 + x
    q0 = x**5 + x
    p0x = sp.diff(p0, x)
    q0x = sp.diff(q0, x)
    assert sp.degree(p0, x) == 3
    assert sp.degree(q0, x) == 5
    assert sp.gcd(p0x, q0x) == 1

    h = 1 + t**4
    f = sp.expand(h**3)
    g = sp.expand(h**5)
    assert sp.degree(f, t) == 12
    assert sp.degree(g, t) == 20
    assert sp.expand(3 * f * sp.diff(g, t) - 5 * sp.diff(f, t) * g) == 0
    assert coefficient(f, t, 0) == sp.LC(p0, x)
    assert coefficient(g, t, 0) == sp.LC(q0, x)

    p: dict[int, sp.Expr] = {0: p0}
    q: dict[int, sp.Expr] = {0: q0}
    edge_layers = 0
    free_dimensions_before_pinning = 0
    free_dimensions_after_pinning = 0

    for n in range(60):
        layer = n + 1
        d = (layer + 15) // 5
        e = d + 2
        assert d == min(15, (layer + 15) // 5)

        tail = sp.Integer(0)
        for i in range(1, n + 1):
            tail += (n - i + 1) * sp.diff(p[i], x) * q[n - i + 1]
        for i in range(n):
            tail -= (i + 1) * p[i + 1] * sp.diff(q[n - i], x)

        rhs = ((1 if n == 0 else 0) - tail) / sp.Integer(n + 1)
        p_fixed = sp.Integer(0)
        q_fixed = sp.Integer(0)
        fixed_coefficients: dict[tuple[str, int], sp.Expr] = {}
        if layer % 5 == 0:
            edge_layers += 1
            k = layer // 5
            p_lead = coefficient(f, t, k)
            q_lead = coefficient(g, t, k)
            p_fixed = p_lead * x**d
            q_fixed = q_lead * x**e
            fixed_coefficients[("p", d)] = p_lead
            fixed_coefficients[("q", e)] = q_lead

        residual = sp.expand(rhs - p0x * q_fixed + p_fixed * q0x)

        unknowns: list[tuple[str, int]] = []
        columns: list[sp.Expr] = []
        for power in range(d + 1):
            if ("p", power) not in fixed_coefficients:
                unknowns.append(("p", power))
                columns.append(-x**power * q0x)
        for power in range(e + 1):
            if ("q", power) not in fixed_coefficients:
                unknowns.append(("q", power))
                columns.append(p0x * x**power)

        row_degrees = range(d + 5)
        matrix = sp.Matrix(
            [
                [coefficient(col, x, power) for col in columns]
                for power in row_degrees
            ]
        )
        target = sp.Matrix([coefficient(residual, x, power) for power in row_degrees])
        solution, parameters = matrix.gauss_jordan_solve(target)
        zero_parameters = {parameter: 0 for parameter in parameters}
        solution = solution.subs(zero_parameters)
        assert matrix * solution == target, f"failed layer {layer}"

        p_layer = p_fixed
        q_layer = q_fixed
        for value, (kind, power) in zip(solution, unknowns, strict=True):
            if kind == "p":
                p_layer += value * x**power
            else:
                q_layer += value * x**power

        p[layer] = sp.Poly(sp.expand(p_layer), x).as_expr()
        q[layer] = sp.Poly(sp.expand(q_layer), x).as_expr()
        assert sp.degree(p[layer], x) <= d
        assert sp.degree(q[layer], x) <= e

        # The linear map has kernel dimension d-1. Prescribing the two edge
        # leaders removes one further degree of freedom on each edge layer.
        free_dimensions_before_pinning += d - 1
        free_dimensions_after_pinning += d - 1 - (1 if layer % 5 == 0 else 0)

        row_check = sp.expand(
            (n + 1) * (p0x * q[layer] - p[layer] * q0x) + tail
            - (1 if n == 0 else 0)
        )
        assert row_check == 0, f"recurrence row {n} failed"

    assert edge_layers == 12
    assert free_dimensions_before_pinning == 462
    assert free_dimensions_after_pinning == 450

    print(
        json.dumps(
            {
                "top_edge_identity": "3*f*g' - 5*f'*g = 0",
                "p0_prime_q0_prime_gcd_degree": int(sp.degree(sp.gcd(p0x, q0x), x)),
                "y_output_rows_solved": 60,
                "last_output_y_degree": 59,
                "prescribed_edge_layers": edge_layers,
                "syzygy_dimension_before_edge_pinning": free_dimensions_before_pinning,
                "remaining_affine_freedom_after_edge_pinning": free_dimensions_after_pinning,
                "all_exact_rational": True,
                "scope": "positive-control truncation only; not a full Keller pair",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
