"""Enumerate the Magnus-IV innermost-vertex region for the F2 (75,125) box."""

from fractions import Fraction


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


def main() -> None:
    a, b, m, n = 3, 5, 15, 60
    c = Fraction(a, a + b)
    slope = (Fraction(n) - c) / (Fraction(m) - c)
    lower = lambda x: slope * (x - c) + c
    upper = lambda x: Fraction(n, m) * x
    y_limit = Fraction(a - 1, a) * n

    points = [
        (x, y)
        for x in range(m + 1)
        for y in range(n + 1)
        if 0 <= y < y_limit and lower(x) <= y <= upper(x)
    ]
    points.remove((0, 0))

    # Branch-local support count for P=R^3+Z.  If S_P is the F2 outer support
    # polygon and B is the branch rectangle for Z, then
    #   3 N(R) = N(R^3) subset conv(Supp(P) union Supp(Z))
    #                          subset conv(S_P union B).
    # Enumerate lattice points of the scaled hull inside R's original box.
    # Q has the 1,576-term outer F2 support envelope.
    sp_vertices = [(0, 0), (3, 0), (15, 60), (0, 60)]
    q_terms = sum(
        1
        for i in range(26)
        for j in range(101)
        if 5 * i - j <= 25
    )
    assert q_terms == 1576

    print(f"a={a} b={b} m={m} n={n}")
    print(f"C={c} slope={slope} lower=(53*x-15)/13 upper=4*x y<{y_limit}")
    print(f"lattice_points={points}")
    print(f"count={len(points)}")
    print("vertex z_terms R_allowed_terms model_variables_after_hull_bound")
    for x, y in points:
        z_terms = (x + 1) * (y + 1)
        hull = convex_hull(sp_vertices + [(0, 0), (x, 0), (x, y), (0, y)])
        r_terms = [
            (i, j)
            for i in range(6)
            for j in range(21)
            if in_convex_polygon(hull, (3 * i, 3 * j))
        ]
        assert (5, 20) in r_terms
        print(f"({x},{y}) {z_terms} {len(r_terms)} {len(r_terms)-1+z_terms+q_terms}")


if __name__ == "__main__":
    main()
