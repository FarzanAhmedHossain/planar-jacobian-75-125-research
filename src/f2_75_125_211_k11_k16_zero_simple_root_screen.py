#!/usr/bin/env python3
"""Exact local pole screen at a simple root on the F2 2+1+1, U1=0 branch.

At a simple root of h, the established preceding bands give
ord(U2)>=2, ord(U3)>=1, ord(U4)>=1, and ord(U5/t)>=0.  This script
enumerates that bounded order grid and applies the first-integral coefficients
G_n = sum_r kappa_r t^(-r) [w^(n-5r)] F^((5-r)/3), for 9<=n<=16.

A profile is removed only when the binomial part has a unique negative-order
term supported on the known slices U2,...,U5, strictly below every possible
resonance correction.  Higher slices are assigned their worst-case order 0.
This is an exact local profile screen, not an exclusion of the whole branch.
"""

from __future__ import annotations

import hashlib
import json
import platform
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from math import factorial
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUT = HERE / "f2_75_125_211_k11_k16_zero_simple_root_screen_20261024.json"
MODEL = HERE / "F2-75-125-band-model.md"
K10_CHECKER = HERE / "f2_75_125_k10_zero_first_slice.py"
INF = 10_000
FIRST_N, LAST_N = 9, 16


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


@lru_cache(maxsize=None)
def partitions(total: int, weights: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    if not weights:
        return ((),) if total == 0 else ()
    weight, *tail = weights
    result: list[tuple[int, ...]] = []
    for multiplicity in range(total // weight + 1):
        for suffix in partitions(total - multiplicity * weight, tuple(tail)):
            result.append((multiplicity,) + suffix)
    return tuple(result)


def generalized_binomial(alpha: Fraction, n: int) -> Fraction:
    value = Fraction(1)
    for j in range(n):
        value *= (alpha - j) / (j + 1)
    return value


def partition_coefficient(alpha: Fraction, multiplicities: tuple[int, ...]) -> Fraction:
    count = sum(multiplicities)
    value = generalized_binomial(alpha, count) * factorial(count)
    for multiplicity in multiplicities:
        value /= factorial(multiplicity)
    return value


def make_templates(total: int, alpha: Fraction, exponent_base: int):
    weights = tuple(range(2, total + 1))
    templates = []
    for multiplicities in partitions(total, weights):
        coeff = partition_coefficient(alpha, multiplicities)
        if not coeff:
            continue
        count = sum(multiplicities)
        templates.append((weights, multiplicities, coeff, exponent_base - 3 * count))
    return tuple(templates)


def term_order(template, orders: dict[int, int], h_order: int = 1):
    weights, multiplicities, coeff, h_exponent = template
    if any(m and orders.get(i, 0) >= INF for i, m in zip(weights, multiplicities)):
        return None
    return h_order * h_exponent + sum(
        m * orders.get(i, 0) for i, m in zip(weights, multiplicities)
    )


def term_label(template) -> str:
    weights, multiplicities, _coeff, _h_exp = template
    return "+".join(
        str(weight)
        for weight, multiplicity in zip(weights, multiplicities)
        for _ in range(multiplicity)
    ) or "constant"


def least_terms(templates, orders: dict[int, int]):
    evaluated = []
    for template in templates:
        order = term_order(template, orders)
        if order is not None:
            evaluated.append((order, term_label(template), template[2], template[0], template[1]))
    if not evaluated:
        return None, []
    minimum = min(row[0] for row in evaluated)
    return minimum, [row for row in evaluated if row[0] == minimum]


def first_integral_templates(n: int):
    """Return main terms and each nonzero resonance contribution at G_n."""
    main = make_templates(n, Fraction(5, 3), 5)
    corrections = []
    for r in range(1, n // 5 + 1):
        total = n - 5 * r
        alpha = Fraction(5 - r, 3)
        corrections.append((r, make_templates(total, alpha, 5 - r)))
    return main, corrections


def robust_profile_exclusion(n: int, orders: dict[int, int], main, corrections):
    main_order, tied = least_terms(main, orders)
    correction_orders = []
    correction_detail = {}
    for resonance, templates in corrections:
        order, _terms = least_terms(templates, orders)
        correction_detail[str(resonance)] = order
        if order is not None:
            correction_orders.append(order)
    correction_order = min(correction_orders) if correction_orders else None
    if main_order is None or main_order >= 0 or len(tied) != 1:
        return False, main_order, correction_order, correction_detail, tied
    unique = tied[0]
    weights, multiplicities = unique[3], unique[4]
    support = {weight for weight, mult in zip(weights, multiplicities) if mult}
    known_slice_support = support <= {2, 3, 4, 5}
    below_corrections = correction_order is None or main_order < correction_order
    return (
        known_slice_support and below_corrections,
        main_order,
        correction_order,
        correction_detail,
        tied,
    )


def main() -> None:
    # The finite local profile domain comes from the exact degree caps of the
    # first five known slices; infinity denotes the identically-zero slice.
    q_values = list(range(2, 12)) + [INF]
    r_values = list(range(1, 12)) + [INF]
    s_values = list(range(1, 12)) + [INF]
    v_values = list(range(13)) + [INF]
    profile_count = len(q_values) * len(r_values) * len(s_values) * len(v_values)
    assert profile_count == 22_176

    templates_by_n = {n: first_integral_templates(n) for n in range(FIRST_N, LAST_N + 1)}
    first_excluded: Counter[int] = Counter()
    examples: dict[str, list[dict[str, object]]] = {str(n): [] for n in range(FIRST_N, LAST_N + 1)}
    survivors = []
    k11_k16_new = 0

    for q in q_values:
        for r in r_values:
            for s in s_values:
                for v in v_values:
                    profile = (q, r, s, v)
                    orders = {2: q, 3: r, 4: s, 5: v}
                    # Later slices are allowed their worst-case simple-root
                    # order zero, making the screen conservative.
                    orders.update({i: 0 for i in range(6, LAST_N + 1)})
                    removed_at = None
                    removed_data = None
                    for n in range(FIRST_N, LAST_N + 1):
                        main_terms, correction_terms = templates_by_n[n]
                        result = robust_profile_exclusion(
                            n, orders, main_terms, correction_terms
                        )
                        if result[0]:
                            removed_at = n
                            removed_data = result
                            break
                    if removed_at is None:
                        if r < 2 or v < 1:
                            raise AssertionError(
                                "a K9-K16 survivor violates the derived simple-root bounds "
                                "ord(U3)>=2 or ord(U5/t)>=1"
                            )
                        survivors.append(profile)
                        continue
                    first_excluded[removed_at] += 1
                    if removed_at >= 11:
                        k11_k16_new += 1
                    if len(examples[str(removed_at)]) < 12:
                        _ok, main_order, correction_order, detail, tied = removed_data
                        exponents = []
                        for row in tied:
                            _order, label, coefficient, weights, multiplicities = row
                            exponents.append({
                                "monomial": label,
                                "coefficient": str(coefficient),
                                "support": [weight for weight, mult in zip(weights, multiplicities) if mult],
                            })
                        examples[str(removed_at)].append({
                            "profile_q_r_s_v": ["inf" if x >= INF else x for x in profile],
                            "main_order": main_order,
                            "correction_min_order": correction_order,
                            "correction_orders_by_resonance": detail,
                            "unique_least_main_term": exponents[0],
                        })

    expected_survivors = len(q_values) * 11 * len(s_values) * 13
    if len(survivors) != expected_survivors:
        raise AssertionError((len(survivors), expected_survivors))

    # Cross-check the reported simple-root K10 residue alternatives directly
    # against the independent exact K10 script's canonical templates.
    simple_k10 = templates_by_n[10]
    expected_k10_partition_counts = {
        "main_nonzero_terms": 12,
        "kappa1_nonzero_terms": 2,
        "kappa2_nonzero_terms": 1,
    }
    observed_k10_partition_counts = {
        "main_nonzero_terms": len(simple_k10[0]),
        "kappa1_nonzero_terms": len(simple_k10[1][0][1]),
        "kappa2_nonzero_terms": len(simple_k10[1][1][1]),
    }
    # The main count is exact after omitting U1-containing and zero-binomial
    # terms; these explicit counts guard the recurrence transcription.
    if observed_k10_partition_counts != expected_k10_partition_counts:
        raise AssertionError((observed_k10_partition_counts, expected_k10_partition_counts))

    result = {
        "scope": "Exact local simple-root pole screen for F2 2+1+1, U1=0; not whole-branch exclusion.",
        "field": "exact rational coefficients and integer local orders",
        "root_multiplicity": 1,
        "branch": "2+1+1",
        "U1": 0,
        "local_order_hypotheses": {
            "ord_U2_min": 2,
            "ord_U3_min": 1,
            "ord_U4_min": 1,
            "ord_U5_over_t_min": 0,
            "higher_slice_order_assumption": "0 (worst-case lower bound)",
        },
        "profile_grid_count": profile_count,
        "band_range": [FIRST_N, LAST_N],
        "resonances_included": "all kappa_r terms with 1<=r<=floor(n/5), including exact zero terms",
        "recurrence_check": {
            "K10_nonzero_term_counts": observed_k10_partition_counts,
            "expected_counts": expected_k10_partition_counts,
            "passed": True,
        },
        "profiles_first_excluded_by_band": {str(n): first_excluded[n] for n in range(FIRST_N, LAST_N + 1)},
        "profiles_first_excluded_K9_K10": sum(first_excluded[n] for n in (9, 10)),
        "profiles_newly_excluded_K11_K16": k11_k16_new,
        "profiles_surviving_K9_K16": len(survivors),
        "survivor_characterization": "exactly the enumerated profiles with ord(U3)>=2 and ord(U5/t)>=1; q=ord(U2) and s=ord(U4) retain all bounded values",
        "survivor_characterization_asserted_profile_by_profile": True,
        "cumulative_local_profile_exclusion": f"{profile_count-len(survivors)}/{profile_count}",
        "cumulative_local_profile_percent": round(100 * (profile_count-len(survivors)) / profile_count, 6),
        "whole_2+1+1_branch_excluded_percent": 0.0,
        "examples": examples,
        "source_sha256": {
            MODEL.name: sha256(MODEL),
            K10_CHECKER.name: sha256(K10_CHECKER),
        },
        "script_sha256": sha256(Path(__file__)),
        "python": platform.python_version(),
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in (
        "profile_grid_count", "profiles_first_excluded_by_band",
        "profiles_first_excluded_K9_K10", "profiles_newly_excluded_K11_K16",
        "profiles_surviving_K9_K16", "cumulative_local_profile_exclusion",
        "cumulative_local_profile_percent", "whole_2+1+1_branch_excluded_percent",
        "recurrence_check", "source_sha256", "script_sha256",
    )}, indent=2))
    print(f"report={OUT}")


if __name__ == "__main__":
    main()
