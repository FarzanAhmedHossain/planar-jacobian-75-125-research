"""Exact local valuation split from the first five Q-bands in F2 2+1+1.

Combines the V3 and V4 polynomiality conditions with the GGV terminal order
at the unique double root, then checks the unique V5 tie coefficient exactly.
The output is a necessary local split, not a whole-branch exclusion.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "F2-75-125-band-model.md"
OUTPUT = HERE / "f2_75_125_211_first5_band_root_split_20261021_v2.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite {OUTPUT}")

    u, c, a, d, e, f = sp.symbols("u c a d e f")
    h = c * u**2
    U1 = a * u**5
    U2 = d * u
    U3 = e
    U4 = f
    N5 = (
        810 * h**9 * (U1 * U4 + U2 * U3)
        - 135 * h**6 * (U1**2 * U3 + U1 * U2**2)
        + 60 * h**3 * U1**3 * U2
        - 7 * U1**5
    )
    coefficient_u19 = sp.expand(N5).coeff(u, 19)
    expected = 135 * c**6 * d * (6 * c**3 * e - a * d)
    if sp.factor(coefficient_u19 - expected) != 0:
        raise AssertionError("the V5 order-19 tie coefficient changed")

    # Orders of the six summands of N5 at a root where ord(h)=e.
    # The entries correspond to U1*U4, U2*U3, U1^2*U3,
    # U1*U2^2, U1^3*U2, and U1^5.
    def n5_orders(e_order: int, k: int, q: int, r: int, s: int = 0) -> list[int]:
        return [
            9 * e_order + k + s,
            9 * e_order + q + r,
            6 * e_order + 2 * k + r,
            6 * e_order + k + 2 * q,
            3 * e_order + 3 * k + q,
            5 * k,
        ]

    simple_root_orders = n5_orders(e_order=1, k=2, q=1, r=0)
    if simple_root_orders != [11, 10, 10, 10, 10, 10]:
        raise AssertionError("the simple-root V5 divisibility check changed")

    double_root_census = {
        "q>=2": n5_orders(e_order=2, k=5, q=2, r=0),
        "q=1,k=5,r=0": n5_orders(e_order=2, k=5, q=1, r=0),
        "q=1,k=5,r>=1": n5_orders(e_order=2, k=5, q=1, r=1),
        "q=1,k>=6,r=0": n5_orders(e_order=2, k=6, q=1, r=0),
        "q=1,k>=6,r>=1": n5_orders(e_order=2, k=6, q=1, r=1),
    }
    minimum_indices = {
        label: [i for i, order in enumerate(orders) if order == min(orders)]
        for label, orders in double_root_census.items()
    }
    if minimum_indices != {
        "q>=2": [1],  # all orders are at least 20; index 1 attains 20
        "q=1,k=5,r=0": [1, 3],
        "q=1,k=5,r>=1": [3],
        "q=1,k>=6,r=0": [1],
        "q=1,k>=6,r>=1": [1, 3],
    }:
        raise AssertionError("the double-root V5 case split changed")
    if min(double_root_census["q>=2"]) != 20:
        raise AssertionError("the q>=2 case no longer clears the h^10 threshold")
    if min(double_root_census["q=1,k>=6,r>=1"]) != 20:
        raise AssertionError("the k>=6,r>=1 case no longer clears the h^10 threshold")

    B1, C1, E0 = sp.symbols("B1 C1 E0")
    A = c**2 * B1
    D = c * C1
    if sp.factor(A * D - 6 * c**3 * E0 - c**3 * (B1 * C1 - 6 * E0)) != 0:
        raise AssertionError("the parameterized V5 tie relation changed")

    simple_root_minima = {
        "h^9 U1 U4": 11,
        "h^9 U2 U3": 10,
        "h^6 U1^2 U3": 10,
        "h^6 U1 U2^2": 10,
        "h^3 U1^3 U2": 10,
        "U1^5": 10,
    }
    if min(simple_root_minima.values()) < 10:
        raise AssertionError("a simple root has a forbidden V5 pole")

    receipt = {
        "source": str(SOURCE),
        "source_sha256": sha256(SOURCE),
        "branch": "2+1+1",
        "field": "characteristic zero",
        "root_pattern": "one double root and two simple roots of h",
        "V3_condition": "h^4 divides 18 h^3 U1 U2 - U1^3",
        "V3_root_bound": "ord_alpha(U1) >= ceil(4*ord_alpha(h)/3)",
        "terminal_root_condition": "ord_double_root(U1) >= 5",
        "deduced_U1": "U1=(t-1)^5 R(t)^2 B(t), degree(B) <= 2",
        "V4_condition": "h^7 divides 54 h^6 U1 U3 + 27 h^6 U2^2 - 9 h^3 U1^2 U2 + U1^4",
        "V4_deduction": "ord_alpha(U2) >= 1 at every root of h",
        "deduced_U2": "U2=(t-1) R(t) C(t), degree(C) <= 8",
        "V5_condition": "h^10 divides N5",
        "simple_root_V5_minimum_orders": simple_root_minima,
        "simple_root_V5_order_census": simple_root_orders,
        "double_root_V5_order_census": double_root_census,
        "double_root_V5_minimum_term_indices_zero_based": minimum_indices,
        "double_root_V5_cases": [
            "ord(U2) >= 2: all V5 terms have order >= 20",
            "ord(U2)=1, ord(U1)>=6, ord(U3)>=1: all terms have order >= 20",
            "ord(U2)=1, ord(U1)=5, ord(U3)=0: order-19 terms tie and require A*D=6*c^3*E",
        ],
        "order19_coefficient": "135*c^6*D*(6*c^3*E - A*D)",
        "tie_after_substitution": "B(1)*C(1)=6*U3(1)",
        "exhaustive_double_root_parameter_strata": [
            "C(1)=0 (ord(U2)>=2)",
            "C(1)!=0, B(1)=0, U3(1)=0",
            "C(1)!=0, B(1)!=0, B(1)*C(1)=6*U3(1)",
        ],
        "U1_free_coefficients_for_fixed_R": 3,
        "U2_free_coefficients_for_fixed_R": 9,
        "scope_warning": (
            "These are necessary first-five-band local conditions. They do not "
            "exclude the 2+1+1 branch or the degree pair."
        ),
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "EXACT_V5_LOCAL_SPLIT_REPLAY_PASSED",
        "branch": receipt["branch"],
        "deduced_U1": receipt["deduced_U1"],
        "deduced_U2": receipt["deduced_U2"],
        "order19_coefficient": str(coefficient_u19),
        "output": str(OUTPUT),
        "output_sha256": sha256(OUTPUT),
    }, indent=2))


if __name__ == "__main__":
    main()
