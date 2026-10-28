"""Generate and audit the exact sparse Jacobian system for the F2 (75,125) envelope.

This script is deliberately a front-end model checker, not an exclusion proof.
It expands every bilinear coefficient in [P,Q]=1 for the support sets derived
from the standard-pair rectangles and the F2 regular corner of the already
standardized pair, as recorded in F2-75-125-band-model.md. It verifies
term-by-term agreement of the original (x,y) Jacobian and the (t,z) band
formula, and reports a digest of the resulting sparse equations.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from collections.abc import Iterator


Point = tuple[int, int]
Row = tuple[int, int]  # (z-band exponent, t exponent)


def support(i_max: int, j_max: int, edge_max: int) -> list[Point]:
    """The lattice points in 0<=i<=i_max, 0<=j<=j_max, 5i-j<=edge_max."""
    return [
        (i, j)
        for i in range(i_max + 1)
        for j in range(j_max + 1)
        if 5 * i - j <= edge_max
    ]


def band(point: Point) -> int:
    i, j = point
    return 5 * i - j


def iter_bilinear_terms(
    p: list[Point], q: list[Point]
) -> Iterator[tuple[int, int, Point, Point, int]]:
    """Yield each nonzero sparse term as (band, t_degree, p_point, q_point, coeff)."""
    for p_point in p:
        i, j = p_point
        k = band(p_point)
        for q_point in q:
            u, v = q_point
            ell = band(q_point)
            c_xy = i * v - j * u
            c_band = k * u - ell * i
            assert c_xy == c_band
            if c_xy == 0:
                continue

            t_exp = i + u - 1
            z_exp = k + ell - 4
            x_exp = i + u - 1
            y_exp = j + v - 1
            assert t_exp >= 0 and x_exp >= 0 and y_exp >= 0
            assert z_exp == 5 * x_exp - y_exp
            yield z_exp, t_exp, p_point, q_point, c_xy


def audit_system() -> dict[str, object]:
    p = support(15, 60, 15)
    q = support(25, 100, 25)
    assert len(p) == 586
    assert len(q) == 1576
    assert (3, 0) in p and (15, 60) in p
    assert (5, 0) in q and (25, 100) in q

    p_ids = {point: index for index, point in enumerate(p)}
    q_ids = {point: index for index, point in enumerate(q)}
    row_terms: Counter[Row] = Counter()
    by_band_terms: Counter[int] = Counter()
    singleton_candidate: dict[Row, tuple[Point, Point, int]] = {}
    top_edge_terms: list[tuple[int, int, int, int]] = []
    digest = hashlib.sha256()
    all_terms = 0

    # In x,y, [x^i y^j, x^u y^v] has coefficient i*v-j*u and exponent
    # (i+u-1,j+v-1).  In t=xy^5,z=y^-1 these become
    # t^(i+u-1) z^(k+l-4), with coefficient k*u-l*i.  The iterator checks
    # this identity term by term and is reusable by downstream eliminators.
    for z_exp, t_exp, p_point, q_point, c_xy in iter_bilinear_terms(p, q):
        i, j = p_point
        u, v = q_point
        k = band(p_point)
        ell = band(q_point)
        row = (z_exp, t_exp)
        if row_terms[row] == 0:
            singleton_candidate[row] = (p_point, q_point, c_xy)
        else:
            singleton_candidate.pop(row, None)
        row_terms[row] += 1
        by_band_terms[z_exp] += 1
        all_terms += 1
        digest.update(
            f"{z_exp},{t_exp},{p_ids[p_point]},{q_ids[q_point]},{c_xy};".encode()
        )

        if k == 15 and ell == 25:
            # Write i=3+a, u=5+b on the two forced edge supports.
            a, b = i - 3, u - 5
            assert 0 <= a <= 12 and 0 <= b <= 20
            assert t_exp == 7 + a + b
            assert c_xy == 5 * (3 * b - 5 * a)
            top_edge_terms.append((a, b, t_exp, c_xy))

    # The z^36 row is contributed only by the unique maximal band pair
    # (k,l)=(15,25), since k<=15 and l<=25 throughout the support envelope.
    top_rows = [row for row in row_terms if row[0] == 36]
    top_t_degrees = sorted(tp for z, tp in top_rows)
    assert top_t_degrees == list(range(8, 39))
    assert len(top_edge_terms) == sum(1 for a in range(13) for b in range(21) if 3 * b - 5 * a)
    assert all(z == 36 for z, _ in top_rows)

    # The constant right-hand side occurs in exactly one output coefficient.
    rhs_row = (0, 0)
    assert rhs_row in row_terms
    assert min(z for z, _ in row_terms) == -159
    assert max(z for z, _ in row_terms) == 36
    assert min(t for _, t in row_terms) == 0
    # The only possible t^39 term is the bracket of (15,60) and (25,100);
    # its coefficient is 15*100-60*25=0, so the largest nonzero row is t^38.
    assert max(t for _, t in row_terms) == 38

    largest_row, largest_count = max(row_terms.items(), key=lambda item: item[1])
    stream_digest = digest.hexdigest()
    assert all_terms == 912_713
    assert len(row_terms) == 3_790
    assert stream_digest == "2e36049f836f07e455fa076279a3fa41eef70d9f30ddc48ce2d03965ac19117e"
    assert not singleton_candidate

    terms_by_band = dict(sorted(by_band_terms.items()))
    rows_by_band = Counter(z for z, _ in row_terms)
    assert rows_by_band[36] == 31
    assert by_band_terms[36] == 268
    forced_p_endpoints = {(3, 0), (15, 60)}
    forced_q_endpoints = {(5, 0), (25, 100)}
    forced_endpoint_singletons = [
        {
            "row": [z, t],
            "P_monomial": list(p_point),
            "Q_monomial": list(q_point),
            "coefficient": coefficient,
            "rhs": 1 if (z, t) == rhs_row else 0,
        }
        for (z, t), (p_point, q_point, coefficient) in sorted(singleton_candidate.items())
        if p_point in forced_p_endpoints and q_point in forced_q_endpoints
    ]
    return {
        "P_variables": len(p),
        "Q_variables": len(q),
        "total_variables": len(p) + len(q),
        "P_band_range": [min(map(band, p)), max(map(band, p))],
        "Q_band_range": [min(map(band, q)), max(map(band, q))],
        "nonzero_equation_rows": len(row_terms),
        "nonzero_bilinear_terms": all_terms,
        "zero_bilinear_pairs": len(p) * len(q) - all_terms,
        "output_z_band_range": [min(z for z, _ in row_terms), max(z for z, _ in row_terms)],
        "output_t_degree_range": [min(t for _, t in row_terms), max(t for _, t in row_terms)],
        "largest_row": {"z_band": largest_row[0], "t_degree": largest_row[1], "terms": largest_count},
        "constant_jacobian_row": {"z_band": rhs_row[0], "t_degree": rhs_row[1], "rhs": 1},
        "top_edge_nonzero_terms": len(top_edge_terms),
        "top_edge_band": 36,
        "top_edge_equation_rows": rows_by_band[36],
        "top_edge_total_terms": terms_by_band[36],
        "top_edge_t_degree_range": [min(tp for _, _, tp, _ in top_edge_terms), max(tp for _, _, tp, _ in top_edge_terms)],
        "bilinear_term_stream_sha256": stream_digest,
        "singleton_equation_rows": len(singleton_candidate),
        "singleton_rows_with_both_forced_edge_endpoints": forced_endpoint_singletons,
        "coefficient_identity": "i*v-j*u = (5*i-j)*u-(5*u-v)*i",
        "equation_identity": "[P,Q]=1 iff each (z_band,t_degree) row equals its Kronecker-delta RHS",
    }


if __name__ == "__main__":
    print(json.dumps(audit_system(), indent=2, sort_keys=True))
