"""Compare the q=4,s=2 K12/K14/K15 residues with the depth-6 2+2 state.

This is a diagnostic only. It decodes the exact checkpoint into Q(sqrt(5))[p]
Taylor jets, specializes the global U1=0 stratum by p1=p2=0, and records which
rootwise congruences are identities in the saved family. No primary state is
modified and no exclusion is claimed.
"""

from __future__ import annotations

import json
from math import comb
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.rings import ring


HERE = Path(__file__).resolve().parent
STATE = HERE / "f2_75_125_sparse_ring_lift_22_terminal_roots6_q_replayed_state.json"
OUT = HERE / "f2_75_125_depth6_q4s2_k12k14k15_compare_agent_20261018.json"
SQRT5 = sp.sqrt(5)
ALPHA = (3 + SQRT5) / 2


def main() -> None:
    data = json.loads(STATE.read_text(encoding="utf-8"))
    field = QQ.algebraic_field(SQRT5)
    names = [f"p{i}" for i in range(56)]
    P, *ps = ring(names, field)
    old_names = data["generators"]
    p_index = {name: i for i, name in enumerate(names)}
    zero = P.zero

    def scalar(value):
        return P.ground_new(field.convert(value))

    def coefficient(values):
        expr = sum(
            (sp.Rational(value) * SQRT5 ** (len(values) - 1 - i)
             for i, value in enumerate(values)),
            sp.S.Zero,
        )
        return field.from_sympy(expr)

    def jet(terms, tau, order=0, kill_u1=False):
        result = {}
        tau_k = field.from_sympy(sp.sympify(tau))
        for sparse, values in terms:
            t_degree = 0
            monomial = [0] * len(ps)
            discard = False
            for old_index, exponent in sparse:
                name = old_names[old_index]
                if name == "t":
                    t_degree = exponent
                else:
                    pi = p_index[name]
                    if kill_u1 and pi in (1, 2) and exponent:
                        discard = True
                        break
                    monomial[pi] = exponent
            if discard or t_degree < order:
                continue
            scale = field.convert(comb(t_degree, order)) * tau_k ** (t_degree - order)
            key = tuple(monomial)
            result[key] = result.get(key, field.zero) + coefficient(values) * scale
        return P.from_dict({m: c for m, c in result.items() if c})

    def stats(poly):
        return {"terms": len(poly), "is_zero": not bool(poly)}

    def field_expr(value):
        return sp.sstr(field.to_sympy(value))

    # Root values of the saved state, then the exact U1=0 specialization.
    roots = (("1", sp.Integer(1)), ("alpha", ALPHA))
    summary = {
        "input_state": str(STATE.relative_to(HERE.parent)),
        "state_constraints": len(data.get("constraints", [])),
        "parameter_count": data["nparameters"],
        "terminal_P_rows": "resume verifier rechecks distances 0..5 (32 forbidden coefficients)",
        "global_U1_zero_specialization": "p1=p2=0; S depends only on p1,p2",
        "roots": {},
        "scope": "exact necessary-condition comparison only; no exclusion",
    }

    for root_name, tau in roots:
        print(f"evaluating root {root_name}", flush=True)
        raw_S = jet(data["S"], tau)
        raw_B = jet(data["B"], tau)
        raw_U3 = jet(data["uk"]["3"], tau)
        raw_D = 3 * raw_B - raw_S**2
        raw_E = 27 * raw_U3 - raw_S**3
        rec = {
            "raw_root_values": {
                "S": stats(raw_S),
                "D=3B-S^2": stats(raw_D),
                "E=27U3-S^3": stats(raw_E),
            }
        }

        h_jets = [jet(data["h"], tau, j, kill_u1=True) for j in range(6)]
        c = h_jets[2]
        c1 = h_jets[3]
        if any(x for j, x in enumerate(h_jets) if j < 2):
            raise AssertionError("h did not have a double root")
        c_value = next(iter(c.values()))
        c1_value = next(iter(c1.values())) if c1 else field.zero
        tau_value = field.from_sympy(sp.sympify(tau))
        r0prime = 2 * tau_value - (1 + field.from_sympy(ALPHA))

        u = {
            k: [jet(data["uk"][str(k)], tau, j, kill_u1=True)
                for j in range(0, 7)]
            for k in range(2, 7)
        }
        profile = {
            "U2_coefficients_orders_0_to_3": [stats(u[2][j]) for j in range(4)],
            "U3_coefficients_orders_0_to_3": [stats(u[3][j]) for j in range(4)],
            "U4_coefficients_orders_0_to_1": [stats(u[4][j]) for j in range(2)],
            "U5_coefficients_orders_0_to_1": [stats(u[5][j]) for j in range(2)],
        }

        # q=4,s=2 leading coefficients; the audit's q=4,s=2 equations apply
        # only when d and a are nonzero. No division by d,a or residual factors
        # is done here.
        d = u[2][4]
        a = u[4][2]
        f = u[6][0] / scalar(tau_value)
        A = scalar(c_value**3) * a
        F = scalar(c_value**6) * f
        e12 = (81*A**3 - 162*A**2*d**2 + 486*A*F*d + 63*A*d**4
               - 729*F**2 - 108*F*d**3 - 7*d**6)
        e14 = (324*A**3*d - 729*A**2*F - 378*A**2*d**3
               + 972*A*F*d**2 + 126*A*d**5 - 729*F**2*d
               - 189*F*d**4 - 13*d**7)
        sel = {
            "d_q4_leading": stats(d),
            "a_s2_leading": stats(a),
            "K12_cleared_selector_residual": stats(e12),
            "K14_cleared_selector_residual": stats(e14),
            "base_quotient_r0_divides_3C4_minus_B2sq_at_root": stats(
                3 * a / scalar(c_value) - (d / scalar(c_value**2))**2
            ),
            "base_quotient_r0_divides_27U6_minus_tB2cube_at_root": stats(
                27 * u[6][0] - scalar(tau_value) * (d / scalar(c_value**2))**3
            ),
        }
        print(f"root {root_name}: selector remainders constructed", flush=True)

        # Conditional post-cube quotients at this root. These are Taylor
        # quotient values assuming the K12/K14 base numerators vanish there.
        # All denominators below are nonzero field constants (root and h lead).
        b2 = d / scalar(c_value**2)
        b2_1 = (u[2][5] - scalar(2*c_value*c1_value)*b2) / scalar(c_value**2)
        b3 = u[3][4] / scalar(c_value**2)
        c4 = a / scalar(c_value)
        c4_1 = (u[4][3] - scalar(c1_value)*c4) / scalar(c_value)
        L = (3*c4_1 - 2*b2*b2_1) / scalar(r0prime)
        T = u[5][2] / scalar(c_value)
        N = (27*u[6][1] - b2**3 - scalar(3*tau_value)*b2**2*b2_1) / scalar(r0prime)
        k15 = scalar(tau_value)*b3*N**2 + scalar(18*tau_value)*L*T*N + scalar(243)*T**3
        print(f"root {root_name}: K15 polynomial constructed", flush=True)
        sel["conditional_quotient_values"] = {
            "B3": stats(b3), "L": stats(L), "T": stats(T), "N": stats(N),
            "K15_r0_divisibility_residual": stats(k15),
            "K15_zero_cases_are_undivided": True,
            "K14_next_equation_coefficient_of_U8": stats(scalar(54)*N),
            "K14_next_equation_without_U8": stats(
                -scalar(18)*b3*T*N - scalar(tau_value)*L**2*N - scalar(81)*L*T**2
            ),
            "K14_next_when_N_zero_reduces_to_L_T_squared": stats(-scalar(81)*L*T**2),
        }
        rec["U1_zero_q4s2_profile"] = profile
        rec["K12_K14_selectors"] = sel
        summary["roots"][root_name] = rec

    OUT.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
