"""Exact support consequences of the Magnus-IV cube-root split for (75,125).

The checker uses only necessary supports:
  Supp(P) subset S_P, Supp(Z) subset the branch rectangle B, P = R^3 + Z,
  and the normalized corner coefficient r_(5,20)=1.

It verifies that the branchwise coefficient equations of R^3 outside
S_P union B force every R coefficient outside 5*i-j <= 5 to vanish.  It then
counts the resulting branchwise (R,Z,Q) variables.  This is a support
reduction, not an exclusion of the Keller pair.
"""

from collections import defaultdict
from fractions import Fraction
from itertools import combinations_with_replacement


def convex_hull(points):
    pts = sorted(set(points))

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def in_convex_polygon(poly, point):
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    return all(
        cross(poly[i], poly[(i + 1) % len(poly)], point) >= 0
        for i in range(len(poly))
    )


def magnus_vertices():
    a, b, m, n = 3, 5, 15, 60
    corner = Fraction(a, a + b)
    slope = (Fraction(n) - corner) / (Fraction(m) - corner)
    lower = lambda x: slope * (x - corner) + corner
    upper = lambda x: Fraction(n, m) * x
    y_limit = Fraction(a - 1, a) * n
    points = [
        (x, y)
        for x in range(m + 1)
        for y in range(n + 1)
        if 0 <= y < y_limit and lower(x) <= y <= upper(x)
    ]
    points.remove((0, 0))  # excluded by the noncomposition argument
    return points


def branch_hull_support(vertex, sp_vertices):
    x, y = vertex
    rectangle_vertices = [(0, 0), (x, 0), (x, y), (0, y)]
    hull = convex_hull(sp_vertices + rectangle_vertices)
    return {
        (i, j)
        for i in range(6)
        for j in range(21)
        if in_convex_polygon(hull, (3 * i, 3 * j))
    }


def forced_zero_cube_coefficients(r_support, allowed_rhs):
    """Propagate only unconditional pure-power consequences of R^3=rhs.

    For an exponent outside allowed_rhs, the coefficient in R^3 is zero.
    After the corner coefficient r_(5,20) is set to one and previously forced
    zero variables are substituted, a single pure-power monomial c*r^k=0
    forces r=0 in an integral domain.  Products of different variables are
    deliberately left branched and are not used to assert either variable is
    zero.
    """
    coefficients = defaultdict(list)
    for triple in combinations_with_replacement(sorted(r_support), 3):
        exponent = (
            sum(point[0] for point in triple),
            sum(point[1] for point in triple),
        )
        if exponent in allowed_rhs:
            continue
        distinct = len(set(triple))
        scalar = 6 if distinct == 3 else (3 if distinct == 2 else 1)
        coefficients[exponent].append((triple, scalar))

    corner = (5, 20)
    zeros = set()
    witnesses = []
    contradiction = False
    while True:
        newly_forced = {}
        for exponent, terms in sorted(coefficients.items()):
            surviving = defaultdict(int)
            for monomial, scalar in terms:
                if any(variable in zeros for variable in monomial):
                    continue
                reduced = tuple(variable for variable in monomial if variable != corner)
                surviving[reduced] += scalar
            if len(surviving) != 1:
                continue
            reduced = next(iter(surviving))
            if not reduced:
                contradiction = True
                break
            if len(set(reduced)) == 1:
                variable = reduced[0]
                if variable not in zeros:
                    newly_forced.setdefault(
                        variable,
                        (exponent, len(reduced), surviving[reduced]),
                    )
        if contradiction or not newly_forced:
            break
        for variable, (exponent, power, scalar) in sorted(newly_forced.items()):
            witnesses.append((variable, exponent, power, scalar, len(zeros)))
        zeros |= set(newly_forced)
    return coefficients, zeros, contradiction, witnesses


def main():
    p_support = {
        (i, j)
        for i in range(16)
        for j in range(61)
        if 5 * i - j <= 15
    }
    q_terms = sum(
        1
        for i in range(26)
        for j in range(101)
        if 5 * i - j <= 25
    )
    assert len(p_support) == 586
    assert q_terms == 1576

    sp_vertices = convex_hull(p_support)
    assert sp_vertices == [(0, 0), (3, 0), (15, 60), (0, 60)]

    expected_vertices = [
        (1, 3), (1, 4), (2, 7), (2, 8), (3, 12), (4, 16),
        (5, 20), (6, 24), (7, 28), (8, 32), (9, 36),
    ]
    vertices = magnus_vertices()
    assert vertices == expected_vertices

    r_edge_support = {
        (i, j)
        for i in range(6)
        for j in range(21)
        if 5 * i - j <= 5
    }
    assert len(r_edge_support) == 76

    print("vertex R_hull_terms forbidden_R3_coefficients forced_R_zeros Z_terms model_variables")
    for vertex in vertices:
        x, y = vertex
        r_support = branch_hull_support(vertex, sp_vertices)
        assert r_edge_support <= r_support

        z_rectangle = {
            (i, j)
            for i in range(x + 1)
            for j in range(y + 1)
        }
        allowed_rhs = p_support | z_rectangle
        coeffs, forced_zeros, contradiction, witnesses = forced_zero_cube_coefficients(
            r_support, allowed_rhs
        )
        assert not contradiction
        assert forced_zeros == r_support - r_edge_support

        # Once Supp(R) is in the scaled F2 polygon, every exponent in R^3 is
        # in S_P.  Therefore Z=P-R^3 also has support in S_P.
        z_support = z_rectangle & p_support
        assert vertex in z_support
        forced_r_free = len(r_edge_support) - 1  # normalize r_(5,20)=1
        model_variables = forced_r_free + len(z_support) + q_terms
        print(
            f"({x},{y}) {len(r_support)} {len(coeffs)} {len(forced_zeros)} "
            f"{len(z_support)} {model_variables}"
        )
        for variable, exponent, power, scalar, old_zero_count in witnesses:
            print(
                f"  force={variable[0]},{variable[1]} via R3 coefficient "
                f"at {exponent}: {scalar}*r_{variable[0]},{variable[1]}^{power}=0 "
                f"after {old_zero_count} prior zeros"
            )

    print("result=every Magnus branch forces Supp(R) into 5*i-j<=5")
    print("scope=exact support-cube consequences only; no [P,Q]=1 exclusion")


if __name__ == "__main__":
    main()
