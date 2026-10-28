#!/usr/bin/env python3
"""Exact K9 tied-pole check for the explicit 2+1+1 local survivor.

For the explicit K3--K8 survivor from
f2_75_125_211_nonzero_all_deep_k8_witness_20261021.py, this checker:

* reconstructs the K9 first-integral coefficient
    [w^9]F^(5/3) + kappa_1/t [w^4]F^(4/3);
* proves by weighted valuations that arbitrary polynomial slices U6..U9
  contribute no pole at any root of h; and
* solves exactly for a single kappa_1 that removes all poles at h=0.

This tests one explicit local point, not the whole branch.  Its only purpose
is to determine whether this point lifts through K9 polynomiality.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from itertools import product
from math import factorial
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_nonzero_all_deep_k8_witness_20261021.json"
OUT = HERE / "f2_75_125_211_exact_k9_tied_residue_check_20261021.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def partitions(total: int, max_part: int):
    """Weighted partitions as multiplicities indexed by positive weight."""
    if max_part == 0:
        return [()] if total == 0 else []
    out = []

    def rec(weight: int, left: int, current: list[int]):
        if weight > max_part:
            if left == 0:
                out.append(tuple(current))
            return
        for multiplicity in range(left // weight + 1):
            current.append(multiplicity)
            rec(weight + 1, left - multiplicity * weight, current)
            current.pop()

    rec(1, total, [])
    return out


def generalized_binomial(alpha: Fraction, n: int) -> Fraction:
    value = Fraction(1)
    for i in range(n):
        value *= alpha - i
    return value / factorial(n)


def coefficient(alpha: Fraction, multiplicities: tuple[int, ...]) -> Fraction:
    count = sum(multiplicities)
    value = generalized_binomial(alpha, count) * factorial(count)
    for multiplicity in multiplicities:
        value /= factorial(multiplicity)
    return value


def reduced_parts(expr: sp.Expr, t: sp.Symbol, alpha: sp.Expr):
    expr = sp.cancel(sp.together(expr))
    numerator, denominator = sp.fraction(expr)
    pnum = sp.Poly(sp.expand(numerator), t, extension=alpha)
    pden = sp.Poly(sp.expand(denominator), t, extension=alpha)
    return pnum, pden


def coefficient_numerator(
    total: int,
    alpha: Fraction,
    slices: list[sp.Expr],
    h: sp.Expr,
    t: sp.Symbol,
    clear_power: int,
    clear_t: bool,
) -> sp.Expr:
    """Return h^clear_power (and optionally t) times [w^total]F^alpha."""
    base = 5 if alpha == Fraction(5, 3) else 4
    result = sp.Integer(0)
    for mult in partitions(total, len(slices)):
        if sum((i + 1) * m for i, m in enumerate(mult)) != total:
            continue
        count = sum(mult)
        exponent = clear_power + base - 3 * count
        if exponent < 0:
            raise AssertionError("clearing power is too small")
        term = sp.Rational(coefficient(alpha, mult).numerator,
                           coefficient(alpha, mult).denominator) * h**exponent
        for index, multiplicity in enumerate(mult):
            if multiplicity:
                term *= slices[index] ** multiplicity
        result += term
    if clear_t:
        result *= t
    return sp.expand(result)


def poly_remainder_coefficients(expr: sp.Expr, t: sp.Symbol, modulus: sp.Poly,
                                count: int, extension: sp.Expr) -> list[sp.Expr]:
    poly = sp.Poly(sp.expand(expr), t, extension=extension)
    _, remainder = sp.div(poly, modulus)
    return [sp.simplify(remainder.nth(i)) for i in range(count)]


def solve_k10_local_lift(t: sp.Symbol, h: sp.Expr, R: sp.Expr, U1: sp.Expr,
                         U2: sp.Expr, U3: sp.Expr, U4: sp.Expr, U5_over_t: sp.Expr,
                         alpha_field: sp.Expr):
    """Solve the exact K10 pole equations allowing the full supported U6."""
    unknown_u6 = sp.symbols("u6_0:12")
    kappa1 = sp.Symbol("kappa_1")
    known_orders = {
        1: (5, 3), 2: (3, 2), 3: (2, 2), 4: (0, 1), 5: (0, 1)
    }
    later_slice_terms = []
    for mult in partitions(10, 10):
        if sum((i + 1) * m for i, m in enumerate(mult)) != 10:
            continue
        if not any(i >= 6 and m for i, m in enumerate(mult)):
            continue
        count = sum(mult)
        for root_type, h_order in (("double", 2), ("simple", 1)):
            val = h_order * (5 - 3 * count)
            for i, multiplicity in enumerate(mult):
                if multiplicity and i < 5:
                    val += multiplicity * known_orders[i + 1][
                        0 if root_type == "double" else 1
                    ]
            later_slice_terms.append({
                "partition": "+".join(
                    str(i + 1) for i, m in enumerate(mult) for _ in range(m)
                ),
                "root_type": root_type,
                "minimum_order_with_unknown_slices_nonnegative": val,
            })
    if any(row["minimum_order_with_unknown_slices_nonnegative"] < 0
           for row in later_slice_terms):
        raise AssertionError("some U7..U10 terms can affect the K10 principal part")
    u6_poly = sum(value * t**i for i, value in enumerate(unknown_u6))
    # At distance six the support gives P_9=t^2 U6 with deg(U6)<=11;
    # after the t^3 near-edge normalization F6=U6/t.
    slices6 = [U1, U2, U3, U4, U5_over_t, u6_poly / t]
    main = coefficient_numerator(10, Fraction(5, 3), slices6, h, t,
                                 clear_power=25, clear_t=True)
    correction = coefficient_numerator(5, Fraction(4, 3),
                                       [U1, U2, U3, U4, U5_over_t], h, t,
                                       clear_power=25, clear_t=False)
    # The full regularity condition is h^25 | t*h^25*G10.  The kappa_2
    # term is proportional to h^28/t after this clearing, hence is zero
    # modulo h^25 in the local ring at the roots of h.
    base = sp.expand(main.subs({**{u: 0 for u in unknown_u6}, kappa1: 0}))
    # Main coefficient is linear in the U6 slice; correction is multiplied
    # by the one shared global resonance parameter kappa_1.
    pieces = [sp.expand(sp.diff(main, u)) for u in unknown_u6]
    pieces.append(correction)
    modulus = sp.Poly(sp.expand(h**25), t, extension=alpha_field)
    remainder_count = modulus.degree()
    base_rem = poly_remainder_coefficients(base, t, modulus, remainder_count, alpha_field)
    columns = [poly_remainder_coefficients(piece, t, modulus,
                                           remainder_count, alpha_field)
               for piece in pieces]
    matrix = sp.Matrix(remainder_count, len(columns),
                       lambda row, col: columns[col][row])
    rhs = sp.Matrix([-value for value in base_rem])
    augmented = matrix.row_join(rhs)
    rank = matrix.rank()
    augmented_rank = augmented.rank()
    if augmented_rank > rank:
        status = "NO_SUPPORTED_U6_OR_KAPPA1_LIFTS_THIS_PROFILE_THROUGH_K10"
        witness = None
        solution = None
        # A left-null vector with nonzero pairing against rhs is a compact
        # exact inconsistency witness in the chosen remainder basis.
        for vector in matrix.T.nullspace():
            pairing = sp.simplify((vector.T * rhs)[0])
            if pairing != 0:
                column_pairings = [sp.simplify(value) for value in vector.T * matrix]
                if any(value != 0 for value in column_pairings):
                    raise AssertionError("proposed K10 left-null vector does not annihilate the coefficient matrix")
                if sp.simplify((vector.T * rhs)[0]) == 0:
                    raise AssertionError("proposed K10 left-null vector has zero forcing pairing")
                witness = {
                    "functional_basis": "coefficient of t^i in remainder modulo h^25, i=0,...,99",
                    "functional_vector_length": remainder_count,
                    "left_null_vector_nonzero_entries": [
                        [i, str(value)] for i, value in enumerate(vector) if value != 0
                    ],
                    "pairings_with_all_12_U6_columns_and_kappa1_column": [
                        str(value) for value in column_pairings
                    ],
                    "pairing_with_constant_forcing": str(
                        sp.simplify((vector.T * sp.Matrix(base_rem))[0])
                    ),
                    "pairing_with_rhs": str(pairing),
                }
                break
        if witness is None:
            # Nullspace of matrix.T lives in the equation-row space; this
            # fallback solves the transposed homogeneous system directly.
            left_null = matrix.T.nullspace()
            if left_null:
                vector = left_null[0]
                witness = {
                    "left_null_vector_nonzero_entries": [
                        [i, str(value)] for i, value in enumerate(vector) if value != 0
                    ],
                    "pairing_with_rhs": str(sp.simplify((vector.T * rhs)[0])),
                }
    else:
        status = "SUPPORTED_U6_AND_KAPPA1_LIFT_EXISTS_THROUGH_K10"
        witness = None
        solutions = sp.linsolve((matrix, rhs), [*unknown_u6, kappa1])
        solution = [str(sp.simplify(value)) for value in next(iter(solutions))]
        # Verify every exact coefficient remainder after substitution.
        substitutions = dict(zip([*unknown_u6, kappa1], next(iter(solutions))))
        verified = all(sp.simplify(value.subs(substitutions)) == 0
                       for value in base_rem)
        if not verified:
            raise AssertionError("the reported K10 lift failed exact remainder replay")
    return {
        "profile_lift_status": status,
        "supported_U6_degree_cap": 11,
        "unknown_count_including_kappa1": len(pieces),
        "remainder_equation_count": remainder_count,
        "coefficient_matrix_rank_over_Qsqrt_minus3": rank,
        "augmented_matrix_rank_over_Qsqrt_minus3": augmented_rank,
        "later_U7_to_U10_terms_checked": len(later_slice_terms),
        "minimum_later_slice_term_order": min(
            row["minimum_order_with_unknown_slices_nonnegative"]
            for row in later_slice_terms
        ),
        "all_U7_to_U10_terms_regular_at_roots_of_h": True,
        "exact_solution_if_consistent": solution,
        "inconsistency_witness_if_inconsistent": witness,
    }


def main() -> None:
    receipt = json.loads(SOURCE.read_text(encoding="utf-8"))
    # Use the exact values recorded by the witness checker, in K=Q(sqrt(-3)).
    t = sp.Symbol("t")
    s = sp.sqrt(-3)
    h = (t - 1) ** 2 * (t**2 - 3 * t + 3)
    R = t**2 - 3 * t + 3
    U1 = (t - 1) ** 5 * R**3
    U2 = (t - 1) ** 3 * R**2
    U3 = (t - 1) ** 2 * R**2
    d0 = sp.Rational(1, 6) - s / 18
    d1 = sp.Rational(1, 9) + 8 * s / 27
    U4 = R * (d0 + d1 * (t - 1))
    U5_band = t * R
    U5_over_t = R
    slices = [U1, U2, U3, U4, U5_over_t]

    # At the double root t=1, h has order two. At either simple root of R,
    # h has order one. These orders follow exactly from the factorization.
    known_orders = {
        1: (5, 3),
        2: (3, 2),
        3: (2, 2),
        4: (0, 1),
        5: (0, 1),
    }

    # For every K9 partition involving U6..U9, test its worst-case valuation
    # with ord(Uj)>=0. Such terms must be regular to be irrelevant to the pole
    # calculation for the explicit U1..U5 assignment.
    high_slice_terms = []
    for mult in partitions(9, 9):
        if sum((i + 1) * m for i, m in enumerate(mult)) != 9:
            continue
        if not any(i >= 5 and m for i, m in enumerate(mult)):
            continue
        count = sum(mult)
        exp_h = 5 - 3 * count
        for root_type, h_order in (("double", 2), ("simple", 1)):
            val = h_order * exp_h
            for i, m in enumerate(mult):
                if m and i < 5:
                    val += m * known_orders[i + 1][0 if root_type == "double" else 1]
            high_slice_terms.append({
                "partition": "+".join(str(i + 1) for i, m in enumerate(mult) for _ in range(m)),
                "root_type": root_type,
                "minimum_order_with_unknown_slices_nonnegative": val,
            })
    if any(row["minimum_order_with_unknown_slices_nonnegative"] < 0 for row in high_slice_terms):
        raise AssertionError("some unknown U6..U9 terms can affect the K9 pole part")

    # Construct the known-slice main term. Any U6..U9 terms were proved
    # regular, so they cannot change its principal part at h=0.
    h_symbol = sp.Symbol("h")
    main_part = sp.Integer(0)
    for mult in partitions(9, 5):
        if sum((i + 1) * m for i, m in enumerate(mult)) != 9:
            continue
        count = sum(mult)
        exponent = 5 - 3 * count
        term = sp.Rational(coefficient(Fraction(5, 3), mult).numerator,
                           coefficient(Fraction(5, 3), mult).denominator) * h_symbol**exponent
        for index, multiplicity in enumerate(mult):
            if multiplicity:
                term *= slices[index] ** multiplicity
        main_part += term

    correction_part = sp.Integer(0)
    for mult in partitions(4, 4):
        if sum((i + 1) * m for i, m in enumerate(mult)) != 4:
            continue
        count = sum(mult)
        exponent = 4 - 3 * count
        term = sp.Rational(coefficient(Fraction(4, 3), mult).numerator,
                           coefficient(Fraction(4, 3), mult).denominator) * h_symbol**exponent
        for index, multiplicity in enumerate(mult):
            if multiplicity:
                term *= slices[index] ** multiplicity
        correction_part += term

    main_part = sp.cancel(main_part.subs(h_symbol, h))
    correction_part = sp.cancel(correction_part.subs(h_symbol, h) / t)
    kappa = sp.Symbol("kappa_1")
    G9 = sp.cancel(sp.together(main_part + kappa * correction_part))

    ph = sp.Poly(sp.expand(h), t, extension=s)
    main_num, main_den = reduced_parts(main_part, t, s)
    corr_num, corr_den = reduced_parts(correction_part, t, s)
    common_den = main_den * corr_den
    common_h_power = 0
    while True:
        quotient, remainder = sp.div(common_den, ph)
        if remainder:
            break
        common_den = quotient
        common_h_power += 1
    if common_h_power == 0:
        status = "REGULAR_FOR_ALL_KAPPA1"
        possible_kappas = []
        remainder_coefficients = []
        remainder_relation_rank = 0
    else:
        # On a common denominator, the pole-cancellation condition is
        # rem_h(N0 + kappa_1*N1)=0.  It is a linear system over Q(sqrt(-3)).
        n0 = main_num * corr_den
        n1 = corr_num * main_den
        _, rem0 = sp.div(n0, ph**common_h_power)
        _, rem1 = sp.div(n1, ph**common_h_power)
        equations = [sp.simplify(rem0.nth(i) + kappa * rem1.nth(i))
                     for i in range(max(rem0.degree(), rem1.degree()) + 1)]
        nonzero_pairs = [(sp.simplify(rem0.nth(i)), sp.simplify(rem1.nth(i)))
                         for i in range(len(equations))
                         if rem0.nth(i) != 0 or rem1.nth(i) != 0]
        if not nonzero_pairs:
            status = "REGULAR_FOR_ALL_KAPPA1"
            possible_kappas = []
        else:
            a0, b0 = next(((a, b) for a, b in nonzero_pairs if b != 0), (None, None))
            if b0 is None:
                status = "NO_KAPPA1_REMOVES_ALL_POLES"
                possible_kappas = []
            else:
                candidate = sp.simplify(-a0 / b0)
                compatible = all(sp.simplify(a + candidate * b) == 0
                                 for a, b in nonzero_pairs)
                if compatible:
                    status = "REGULAR_FOR_SOME_KAPPA1"
                    possible_kappas = [str(candidate)]
                else:
                    status = "NO_KAPPA1_REMOVES_ALL_POLES"
                    possible_kappas = []
        remainder_coefficients = [
            f"{sp.simplify(a)} + kappa_1*({sp.simplify(b)})"
            for a, b in nonzero_pairs
        ]
        remainder_relation_rank = 1 if status == "REGULAR_FOR_SOME_KAPPA1" else 2

    # The unnormalized support band is U5=t*R; the formal near-edge series
    # coefficient used in F is U5/t=R.  Keep both names explicit in receipt.
    result = {
        "source": SOURCE.name,
        "source_sha256": sha256(SOURCE),
        "field": "Q(sqrt(-3))",
        "carrier": "h=(t-1)^2(t^2-3t+3)",
        "fixed_slices": {
            "U1": "(t-1)^5 R^3",
            "U2": "(t-1)^3 R^2",
            "U3": "(t-1)^2 R^2",
            "U4": "R*(d0+d1*(t-1))",
            "U5_band": "t*R",
            "U5_over_t": "R",
            "R": "t^2-3t+3",
            "d0": str(d0),
            "d1": str(d1),
        },
        "unknown_U6_to_U9_terms": {
            "partition_count_checked": len(high_slice_terms),
            "minimum_valuation": min(row["minimum_order_with_unknown_slices_nonnegative"] for row in high_slice_terms),
            "all_regular_at_every_root_of_h": True,
            "detail": high_slice_terms,
        },
        "K9_expression": "[w^9]F^(5/3) + kappa_1/t * [w^4]F^(4/3)",
        "common_denominator_h_power": common_h_power,
        "K9_polynomiality_status_for_this_profile": status,
        "kappa_1_values_if_any": possible_kappas,
        "nonzero_exact_remainder_coefficients": remainder_coefficients[:80],
        "remainder_coefficient_count": len(remainder_coefficients),
        "linear_remainder_relation_rank": remainder_relation_rank,
        "scope": "one explicit K3-K8 local survivor; no whole-branch exclusion follows",
    }
    result["K10_local_lift"] = solve_k10_local_lift(
        t, h, R, U1, U2, U3, U4, U5_over_t, s
    )
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": status,
        "kappa_1_values": possible_kappas,
        "nonzero_remainder_coefficient_count": len(remainder_coefficients),
        "unknown_slice_terms_checked": len(high_slice_terms),
        "minimum_unknown_slice_term_order": min(row["minimum_order_with_unknown_slices_nonnegative"] for row in high_slice_terms),
        "K10_local_lift": result["K10_local_lift"],
        "receipt": str(OUT),
        "script_sha256": sha256(Path(__file__)),
        "receipt_sha256": sha256(OUT),
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
