"""Exact audit for the shared, compact series-power cache in the F2 lift.

This checks that the cached coefficients used by the optimized block are
identical to independent direct coefficient extraction over QQ.  It verifies
an implementation optimization only; it is not a branch-exclusion proof.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp
from sympy.polys.domains import QQ


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--receipt",
        type=Path,
        default=Path(__file__).with_name(
            "f2_75_125_sparse_ring_power_cache_audit_20261024.json"
        ),
    )
    args = parser.parse_args()

    research = Path(__file__).resolve().parent
    sys.path.insert(0, str(research))
    import f2_75_125_sparse_ring_lift as lift

    lift.POLY_BACKEND = "flint-mpoly"
    ring, *gens = lift.ring(["t", "a", "b"], QQ)
    t, a, b = gens
    obj = lift.SparseLift.__new__(lift.SparseLift)
    obj.ring = ring
    checks = 0

    for distance in (2, 5, 10, 15, 18):
        p = [
            (i + 1) * a + (i + 2) * b * t + (i % 3 + 1)
            for i in range(distance + 1)
        ]
        q = [
            (2 * i + 1) * b + (i + 1) * a * t + (i % 2 + 1)
            for i in range(distance)
        ] + [ring.zero]

        p_requests = {5: (None, (distance,)), 4: (None, (distance,))}
        q_requests = {3: (None, (distance,)), 2: (None, (distance,))}
        if distance >= 5:
            p_requests[3] = (distance - 5, ())
            p_requests[2] = (distance - 5, ())
        if distance >= 10:
            q_requests[2] = (distance - 10, (distance,))
        if distance >= 15:
            p_requests[4] = (None, (distance, distance - 15))

        p_cache = lift.SparseLift.selected_series_power_coefficients(
            obj, p, p_requests, distance
        )
        q_cache = lift.SparseLift.selected_series_power_coefficients(
            obj, q, q_requests, distance
        )

        for exponent, (prefix, indices) in p_requests.items():
            needed = set(indices)
            if prefix is not None:
                needed.update(range(prefix + 1))
            direct = lift.SparseLift.series_power_coefficients(
                obj, p, exponent, distance
            )
            for degree in needed:
                assert p_cache[exponent][degree] == direct[degree]
                checks += 1

        for exponent, (prefix, indices) in q_requests.items():
            needed = set(indices)
            if prefix is not None:
                needed.update(range(prefix + 1))
            direct = lift.SparseLift.series_power_coefficients(
                obj, q, exponent, distance
            )
            for degree in needed:
                assert q_cache[exponent][degree] == direct[degree]
                checks += 1

        q_distance = a - 2 * b + (distance + 1) * t
        q_full = q[:distance] + [q_distance]
        q2_with_solved_term = q_cache[2][distance] + 2 * q[0] * q_distance
        assert q2_with_solved_term == lift.SparseLift.series_power_coefficient(
            obj, q_full, 2, distance
        )
        checks += 1

        if distance >= 5:
            target = distance - 5
            cached = sum(
                (p_cache[3][i] * q[target - i] for i in range(target + 1)),
                ring.zero,
            )
            direct = sum(
                (
                    lift.SparseLift.series_power_coefficient(obj, p, 3, i)
                    * q[target - i]
                    for i in range(target + 1)
                ),
                ring.zero,
            )
            assert cached == direct
            checks += 1

        if distance >= 10:
            target = distance - 10
            cached = sum(
                (p[i] * q_cache[2][target - i] for i in range(target + 1)),
                ring.zero,
            )
            direct = sum(
                (
                    p[i]
                    * lift.SparseLift.series_power_coefficient(
                        obj, q, 2, target - i
                    )
                    for i in range(target + 1)
                ),
                ring.zero,
            )
            assert cached == direct
            checks += 1

        if distance >= 15:
            assert p_cache[4][distance - 15] == (
                lift.SparseLift.series_power_coefficient(obj, p, 4, distance - 15)
            )
            checks += 1

    # Verify the degree-d unknown-slice split used by new_q: compute all
    # known slices in a parameter-only context, then restore the only possible
    # current-slice terms, 5*p0^4*p_d in P^5 and 4*p0^3*p_d in P^4.
    full_ring, *full_gens = lift.ring(["t", "a", "b", "x0"], QQ)
    t_full, a_full, b_full, x_full = full_gens
    affine_ring = full_ring.prefix(3)
    projected = lift.SparseLift.__new__(lift.SparseLift)
    projected.ring = full_ring
    projected.affine_ring = affine_ring
    projected_checks = 0

    for distance in (2, 5, 10, 15, 18):
        known_p = [
            (i + 1) * a_full + (i + 2) * b_full * t_full + (i % 3 + 1)
            for i in range(distance)
        ]
        p_current = (distance + 2) * x_full + (distance + 3) * a_full * t_full
        full_p = [*known_p, p_current]
        small_p = [value.project_to_ring(affine_ring) for value in known_p]
        small_p.append(affine_ring.zero)

        p_requests = {5: (None, (distance,)), 4: (None, (distance,))}
        if distance >= 5:
            p_requests[3] = (distance - 5, ())
            p_requests[2] = (distance - 5, ())
        if distance >= 15:
            p_requests[4] = (None, (distance, distance - 15))

        projected_cache = lift.SparseLift.selected_series_power_coefficients(
            projected, small_p, p_requests, distance,
            polynomial_ring=affine_ring,
        )
        for exponent, (prefix, indices) in p_requests.items():
            requested_degrees = set(indices)
            if prefix is not None:
                requested_degrees.update(range(prefix + 1))
            direct = lift.SparseLift.series_power_coefficients(
                projected, full_p, exponent, distance
            )
            for degree in requested_degrees:
                value = projected_cache[exponent][degree].embed_into(full_ring)
                if exponent == 5 and degree == distance:
                    value += 5 * full_p[0]**4 * p_current
                if exponent == 4 and degree == distance:
                    value += 4 * full_p[0]**3 * p_current
                assert value == direct[degree]
                projected_checks += 1

    if checks != 116:
        raise AssertionError(f"expected 116 exact equalities, got {checks}")

    source = Path(lift.__file__).resolve()
    report = {
        "result": "pass",
        "scope": "exact coefficient-cache parity only; not a branch exclusion",
        "coefficient_field": "QQ",
        "backend": "python-flint fmpq_mpoly",
        "distances": [2, 5, 10, 15, 18],
        "exact_equalities": checks,
        "projected_affine_equalities": projected_checks,
        "lift_source": str(source),
        "lift_source_sha256": digest(source),
        "audit_script_sha256": digest(Path(__file__).resolve()),
        "sympy_version": sp.__version__,
    }
    args.receipt.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"receipt={args.receipt}")


if __name__ == "__main__":
    main()
