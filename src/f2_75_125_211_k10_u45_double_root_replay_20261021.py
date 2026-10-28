#!/usr/bin/env python3
"""Fast exact replay of the decisive u^45 coefficient in one K10 lift test.

This isolates the t=1 double-root obstruction for the explicit K3--K8
2+1+1 survivor.  Arithmetic is performed in Q[s], s^2=-3, using pairs of
fractions, so no symbolic matrix rank is needed.  The check is pointwise in
the earlier slice coefficients; it is not a branch-wide exclusion.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from math import factorial
from pathlib import Path


HERE = Path(__file__).resolve().parent
K9_RECEIPT = HERE / "f2_75_125_211_exact_k9_tied_residue_check_20261021.json"
WITNESS_RECEIPT = HERE / "f2_75_125_211_nonzero_all_deep_k8_witness_20261021.json"
OUT = HERE / "f2_75_125_211_k10_u45_double_root_replay_20261021.json"
LENGTH = 46  # coefficients u^0,...,u^45
Alg = tuple[Fraction, Fraction]  # a + b*s, where s^2=-3
ZERO: Alg = (Fraction(0), Fraction(0))
ONE: Alg = (Fraction(1), Fraction(0))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def add(a: Alg, b: Alg) -> Alg:
    return (a[0] + b[0], a[1] + b[1])


def mul(a: Alg, b: Alg) -> Alg:
    return (a[0] * b[0] - 3 * a[1] * b[1],
            a[0] * b[1] + a[1] * b[0])


def scale(a: Alg, q: Fraction) -> Alg:
    return (a[0] * q, a[1] * q)


def p_add(a: list[Alg], b: list[Alg]) -> list[Alg]:
    return [add(x, y) for x, y in zip(a, b)]


def p_mul(a: list[Alg], b: list[Alg]) -> list[Alg]:
    out = [ZERO for _ in range(LENGTH)]
    for i, x in enumerate(a):
        if x == ZERO:
            continue
        for j, y in enumerate(b[:LENGTH - i]):
            if y != ZERO:
                out[i + j] = add(out[i + j], mul(x, y))
    return out


def p_scale(a: list[Alg], q: Fraction) -> list[Alg]:
    return [scale(x, q) for x in a]


def p_pow(a: list[Alg], exponent: int) -> list[Alg]:
    out = [ONE] + [ZERO for _ in range(LENGTH - 1)]
    for _ in range(exponent):
        out = p_mul(out, a)
    return out


def p_shift(a: list[Alg], amount: int) -> list[Alg]:
    return [ZERO for _ in range(amount)] + a[:LENGTH - amount]


def partitions(total: int, max_part: int) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []

    def rec(weight: int, left: int, current: list[int]) -> None:
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


def partition_coefficient(alpha: Fraction, multiplicities: tuple[int, ...]) -> Fraction:
    count = sum(multiplicities)
    value = generalized_binomial(alpha, count) * factorial(count)
    for multiplicity in multiplicities:
        value /= factorial(multiplicity)
    return value


def label(multiplicities: tuple[int, ...]) -> str:
    return "+".join(str(weight)
                     for weight, multiplicity in enumerate(multiplicities, 1)
                     for _ in range(multiplicity))


def format_alg(value: Alg) -> str:
    a, b = value
    if b == 0:
        return str(a)
    if a == 0:
        return f"({b})*sqrt(-3)"
    sign = "+" if b > 0 else "-"
    return f"{a} {sign} ({abs(b)})*sqrt(-3)"


def main() -> None:
    k9_receipt = json.loads(K9_RECEIPT.read_text(encoding="utf-8"))
    witness = json.loads(WITNESS_RECEIPT.read_text(encoding="utf-8"))
    if k9_receipt["source_sha256"] != sha256(WITNESS_RECEIPT):
        raise AssertionError("K9 receipt does not match the K3--K8 witness receipt")
    lift = k9_receipt["K10_local_lift"]
    if not (lift["later_U7_to_U10_terms_checked"] == 14
            and lift["minimum_later_slice_term_order"] >= 0
            and lift["all_U7_to_U10_terms_regular_at_roots_of_h"]):
        raise AssertionError("the omitted U7--U10 term valuations are not certified")

    # Local chart u=t-1.  Here R=1-u+u^2 and h=u^2 R.
    R = [ONE, (Fraction(-1), Fraction(0)), ONE] + [ZERO] * (LENGTH - 3)
    h = p_shift(R, 2)
    t = [ONE, ONE] + [ZERO] * (LENGTH - 2)
    U1 = p_shift(p_pow(R, 3), 5)
    U2 = p_shift(p_pow(R, 2), 3)
    U3 = p_shift(p_pow(R, 2), 2)
    # The two K8 coefficient equations for this ansatz are explicit.  The
    # first is quadratic in d0; after using it, the second is linear in d1.
    s = (Fraction(0), Fraction(1))
    d0_minus = add((Fraction(1, 6), Fraction(0)), scale(s, Fraction(-1, 18)))
    d0_plus = add((Fraction(1, 6), Fraction(0)), scale(s, Fraction(1, 18)))

    def k8_equations(d0: Alg, d1: Alg) -> tuple[Alg, Alg]:
        d0sq = mul(d0, d0)
        e4 = scale(add(add(scale(d0sq, Fraction(27)),
                           scale(d0, Fraction(-9))), ONE), Fraction(5, 243))
        numerator = add(add(add(add(scale(d0sq, Fraction(324)),
                                     scale(mul(d0, d1), Fraction(-162))),
                                 scale(d0, Fraction(-90))),
                             scale(d1, Fraction(27))),
                        (Fraction(17), Fraction(0)))
        e5 = scale(numerator, Fraction(-5, 729))
        return e4, e5

    def d1_from_k8(d0: Alg) -> Alg:
        # E5=0 and E4=0 reduce to
        # (27-162*d0)*d1=-(18*d0+5); this denominator is nonzero at both roots.
        denominator = add((Fraction(27), Fraction(0)), scale(d0, Fraction(-162)))
        rhs = add(scale(d0, Fraction(-18)), (Fraction(-5), Fraction(0)))
        if denominator == ZERO:
            raise AssertionError("the reduced K8 equation is not uniquely solvable for d1")
        # Divide in Q[s], s^2=-3.  For a+bs, 1/(a+bs)=(a-bs)/(a^2+3b^2).
        norm = denominator[0] ** 2 + 3 * denominator[1] ** 2
        inverse = (denominator[0] / norm, -denominator[1] / norm)
        return mul(rhs, inverse)

    k8_candidates = [(d0_minus, d1_from_k8(d0_minus)),
                     (d0_plus, d1_from_k8(d0_plus))]
    expected_k8_candidates = [
        ((Fraction(1, 6), Fraction(-1, 18)),
         (Fraction(1, 9), Fraction(8, 27))),
        ((Fraction(1, 6), Fraction(1, 18)),
         (Fraction(1, 9), Fraction(-8, 27))),
    ]
    if k8_candidates != expected_k8_candidates:
        raise AssertionError(f"unexpected K8 candidates: {k8_candidates}")
    for d0, d1 in k8_candidates:
        if k8_equations(d0, d1) != (ZERO, ZERO):
            raise AssertionError("an asserted K8 candidate fails its coefficient equations")

    # Compute [u^45] t*h^25*[w^10]F^(5/3) for both K8 candidates directly
    # by truncated polynomial arithmetic in Q[s], s^2=-3.
    forcing_by_candidate: list[dict[str, object]] = []
    main_partition_count = 0
    for d0, d1 in k8_candidates:
        linear = [d0, d1] + [ZERO] * (LENGTH - 2)
        U4 = p_mul(R, linear)
        F5 = R[:]  # U5/t, the normalized fifth formal coefficient
        slices = [U1, U2, U3, U4, F5]
        known_main = ZERO
        candidate_partition_count = 0
        for mult in partitions(10, 5):
            if sum(weight * m for weight, m in enumerate(mult, 1)) != 10:
                continue
            count = sum(mult)
            exponent = 25 + 5 - 3 * count
            if exponent < 0:
                raise AssertionError("unexpected negative cleared h exponent")
            term = p_mul(t, p_pow(h, exponent))
            for index, multiplicity in enumerate(mult):
                if multiplicity:
                    term = p_mul(term, p_pow(slices[index], multiplicity))
            known_main = add(known_main,
                             scale(term[45], partition_coefficient(Fraction(5, 3), mult)))
            candidate_partition_count += 1
        main_partition_count = candidate_partition_count
        forcing_by_candidate.append({
            "d0": format_alg(d0),
            "d1": format_alg(d1),
            "u45_forcing_coefficient": format_alg(known_main),
            "exact_pair_Q_plus_Qsqrt_minus3": [str(known_main[0]), str(known_main[1])],
        })

    # A valuation-only bound suffices to show the other supported terms
    # cannot affect this coefficient.  ord_u(U1,...,U5)=(5,3,2,0,0).
    known_orders = {1: 5, 2: 3, 3: 2, 4: 0, 5: 0}
    u6_orders: list[dict[str, object]] = []
    for mult in partitions(10, 10):
        if not (mult[5] and not any(mult[6:])):
            continue
        count = sum(mult)
        cleared_order = 2 * (25 + 5 - 3 * count)
        cleared_order += sum(m * known_orders[w]
                             for w, m in enumerate(mult, 1) if w <= 5)
        u6_orders.append({"partition": label(mult),
                          "cleared_numerator_order_lower_bound": cleared_order})
    if min(row["cleared_numerator_order_lower_bound"] for row in u6_orders) != 48:
        raise AssertionError("unexpected U6 contribution order")

    kappa1_orders: list[dict[str, object]] = []
    for mult in partitions(5, 5):
        if sum(weight * m for weight, m in enumerate(mult, 1)) != 5:
            continue
        count = sum(mult)
        cleared_order = 2 * (25 + 4 - 3 * count)
        cleared_order += sum(m * known_orders[w]
                             for w, m in enumerate(mult, 1))
        kappa1_orders.append({"partition": label(mult),
                              "cleared_numerator_order_lower_bound": cleared_order})
    if min(row["cleared_numerator_order_lower_bound"] for row in kappa1_orders) != 51:
        raise AssertionError("unexpected kappa_1 correction order")

    expected_forcings = [
        (Fraction(1, 1458), Fraction(-5, 4374)),
        (Fraction(1, 1458), Fraction(5, 4374)),
    ]
    actual_forcings = [tuple(row["exact_pair_Q_plus_Qsqrt_minus3"])
                       for row in forcing_by_candidate]
    if actual_forcings != [tuple(map(str, pair)) for pair in expected_forcings]:
        raise AssertionError(f"unexpected K10 u^45 forcing coefficients: {actual_forcings}")

    higher_profile = witness["higher_order_profile_candidate"]
    k8_solution = higher_profile["solutions_with_exact_global_h3_divisibility"][0]
    if (k8_solution["d0"] != "1/6 - sqrt(3)*I/18"
            or k8_solution["d1"] != "1/9 + 8*sqrt(3)*I/27"):
        raise AssertionError("the selected K8 branch does not match the local field assignment")

    result = {
        "status": "EXACT_U45_K10_LOCAL_NO_LIFT",
        "field": "Q[s], s^2=-3",
        "root_chart": "u=t-1; R=1-u+u^2; h=u^2 R",
        "source_K9_receipt_sha256": sha256(K9_RECEIPT),
        "source_K3_K8_witness_receipt_sha256": sha256(WITNESS_RECEIPT),
        "K10_polynomiality_requires": "[u^j]N10=0 for j=0,...,49 because ord_u(h^25)=50",
        "coefficient_checked": "[u^45]N10",
        "cleared_numerator": "N10=t*h^25*[w^10]F^(5/3)+kappa_1*h^25*[w^5]F^(4/3)+kappa_2*h^28/t plus later-slice terms",
        "known_slice_partition_count": main_partition_count,
        "K8_coefficient_equations": {
            "order_4": "5*(27*d0^2-9*d0+1)/243=0",
            "order_5": "-5*(324*d0^2-162*d0*d1-90*d0+27*d1+17)/729=0",
            "exact_solution_argument": "order_4 has exactly d0=1/6 +/- sqrt(-3)/18; modulo order_4, order_5 is linear in d1 with nonzero coefficient 27-162*d0, giving exactly two conjugate candidates",
        },
        "K8_admissible_candidates_and_K10_forcing": forcing_by_candidate,
        "U6_partitions": u6_orders,
        "minimum_cleared_U6_order": 48,
        "kappa1_partitions": kappa1_orders,
        "minimum_cleared_kappa1_order": 51,
        "kappa2_cleared_order": 56,
        "U7_to_U10_partition_count": lift["later_U7_to_U10_terms_checked"],
        "minimum_cleared_U7_to_U10_order": 50 + lift["minimum_later_slice_term_order"],
        "conclusion": "the K8 equations reduce this two-parameter U4 ansatz to exactly two conjugate candidates, and both have a nonzero u^45 coefficient, so neither has a supported K10 polynomial lift at t=1",
        "scope": "the subfamily with fixed h,U1,U2,U3,U5/t and U4=R*(d0+d1*u); not a component or the full 2+1+1 branch",
        "source_assignment_family": higher_profile["assignment_family"],
        "source_K8_solutions": [
            {"d0": row["d0"], "d1": row["d1"]}
            for row in higher_profile["solutions_with_exact_global_h3_divisibility"]
        ],
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({**result, "receipt": str(OUT),
                      "script_sha256": sha256(Path(__file__)),
                      "receipt_sha256": sha256(OUT)}, indent=2), flush=True)


if __name__ == "__main__":
    main()
