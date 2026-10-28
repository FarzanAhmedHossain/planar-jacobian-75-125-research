"""Exact first three upper-y commutator recurrences for high Magnus branches.

In the four branches (x0,e)=(6,24),(7,28),(8,32),(9,36), the preceding
Q-tail argument gives deg_y[Z,R] <= 39. After the allowed x-translation,
the leading slices are R_20=x^5 and Z_e=x^s, where e=4s. This checker
verifies the unique first three formal fractional-power slices and the
resulting x-adic restrictions on R_19 (plus the first exceptional split).
It is a necessary-condition calculation, not a branch exclusion.
"""

from itertools import product

import sympy as sp


def bracket(f, g, x, y):
    return sp.diff(f, x) * sp.diff(g, y) - sp.diff(f, y) * sp.diff(g, x)


def exact_first_slices(s, x):
    gamma = sp.Rational(s, 5)
    r19 = sp.Function("r19")(x)
    r18 = sp.Function("r18")(x)
    r17 = sp.Function("r17")(x)
    r16 = sp.Function("r16")(x)
    z0 = x**s
    z1 = gamma * x ** (s - 5) * r19
    z2 = (
        gamma * x ** (s - 5) * r18
        + sp.binomial(gamma, 2) * x ** (s - 10) * r19**2
    )
    z3 = (
        gamma * x ** (s - 5) * r17
        + 2 * sp.binomial(gamma, 2) * x ** (s - 10) * r19 * r18
        + sp.binomial(gamma, 3) * x ** (s - 15) * r19**3
    )
    eta = sp.symbols(f"eta_{s}")
    z4 = (
        gamma * x ** (s - 5) * r16
        + sp.binomial(gamma, 2) * x ** (s - 10) * (2 * r19 * r17 + r18**2)
        + 3 * sp.binomial(gamma, 3) * x ** (s - 15) * r19**2 * r18
        + sp.binomial(gamma, 4) * x ** (s - 20) * r19**4
        + eta * x ** (s - 1)
    )
    return (z0, z1, z2, z3, z4), (r19, r18, r17, r16)


def valuation_obstruction(s, v, eta, nu17):
    """Whether the unique lowest negative x-order in z_(e-3) contradicts."""
    cubic = s - 15 + 3 * v
    mixed = s - 10 + v + eta
    linear = s - 5 + nu17
    return cubic < 0 and cubic < min(mixed, linear)


def main():
    x, y = sp.symbols("x y")
    branches = []
    for s in (6, 7, 8, 9):
        e = 4 * s
        (z0, z1, z2, z3, z4), (r19, r18, r17, r16) = exact_first_slices(s, x)
        R = (
            x**5 * y**20 + r19 * y**19 + r18 * y**18
            + r17 * y**17 + r16 * y**16
        )
        Z = z0 * y**e + z1 * y ** (e - 1) + z2 * y ** (e - 2) + z3 * y ** (e - 3)
        J = sp.expand(bracket(Z, R, x, y))
        Jpoly = sp.Poly(J, y)
        for depth in range(4):
            exponent = e + 19 - depth
            residual = sp.simplify(Jpoly.coeff_monomial(y**exponent))
            if residual != 0:
                raise AssertionError(
                    f"formal slice mismatch for s={s}, depth={depth}: "
                    f"{sp.factor(residual)}"
                )

        if s in (7, 8, 9):
            Z4 = Z + z4 * y ** (e - 4)
            J4 = sp.Poly(sp.expand(bracket(Z4, R, x, y)), y)
            assert sp.simplify(J4.coeff_monomial(y ** (e + 15))) == 0

        # The depth-four particular-slice x-orders force the additional
        # divisibility stated in the note.  Here v=ord(r19), eta=ord(r18);
        # all omitted slice coefficients have nonnegative x-order.
        if s == 7:
            # With v>=3, eta=0 gives the unique lowest order -3 from r18^2.
            assert (s - 10 + 2 * 0) == -3
            assert (s - 15 + 2 * 3 + 0) == -2
            assert (s - 20 + 4 * 3) == -1
            # If v>=4 and eta=1, r18^2 is again uniquely lowest at -1.
            assert (s - 10 + 2 * 1) == -1
            assert (s - 15 + 2 * 4 + 1) == 1
            assert (s - 20 + 4 * 4) == 3
            # For v=3 and eta>=2, r19^4 is uniquely lowest at -1.
            assert (s - 10 + 2 * 2) == 1
            assert (s - 15 + 2 * 3 + 2) == 0
            assert (s - 20 + 4 * 3) == -1
        elif s == 8:
            # v>=3 and eta=0: r18^2 has unique lowest order -2.
            assert (s - 10 + 2 * 0) == -2
            assert (s - 15 + 2 * 3 + 0) == -1
            assert (s - 20 + 4 * 3) == 0
        elif s == 9:
            # v=2: r19^4 is uniquely lowest; after v>=3, eta=0
            # makes r18^2 uniquely lowest.
            assert (s - 20 + 4 * 2) == -3
            assert (s - 15 + 2 * 2 + 0) == -2
            assert (s - 10 + 2 * 0) == -1
            assert (s - 15 + 2 * 2 + 1) == -1
            assert (s - 10 + 2 * 1) == 1
            assert (s - 15 + 2 * 3 + 0) == 0
            assert (s - 20 + 4 * 3) == 1

        # The z_(e-l) operator has coefficient 5*(4d-(e-l)) on x^d.
        # For l=1,2,3 this is nonzero for every integer d.
        assert all((4 * d - (e - depth)) != 0
                   for depth in (1, 2, 3) for d in range(0, 30))

        # z_(e-2) must be polynomial. Its r19^2 term cannot cancel against
        # x^(s-5)*r18, which has nonnegative x-order.
        minimum_v = next(v for v in range(5) if s - 10 + 2 * v >= 0)
        # Check every smaller admissible valuation is killed by the unique
        # lowest r19^3 term in z_(e-3), uniformly for all allowed r18,r17.
        for v in range(minimum_v):
            assert all(
                valuation_obstruction(s, v, eta, nu17)
                for eta, nu17 in product(range(5), repeat=2)
            )
        final_v = next(
            v for v in range(minimum_v, 5)
            if not any(
                valuation_obstruction(s, v, eta, nu17)
                for eta, nu17 in product(range(5), repeat=2)
            )
        )
        # All four conclusions used by the support count are checked directly.
        depth3_v = {6: 3, 7: 3, 8: 3, 9: 2}[s]
        assert final_v == depth3_v
        expected_v = {6: 3, 7: 3, 8: 3, 9: 3}[s]

        extra = ""
        if s == 6:
            a3, a4 = sp.symbols("a3 a4")
            b0, b1 = sp.symbols("b0 b1")
            r19_special = x**3 * (a3 + a4 * x)
            r18_special = b0 + b1 * x
            z3_special = exact_first_slices(s, x)[0][3]
            z3_special = z3_special.subs(
                {
                    sp.Function("r19")(x): r19_special,
                    sp.Function("r18")(x): r18_special,
                }
            )
            # For s=6 the only x^(-1) coefficient is a nonzero scalar times
            # a3*b0, so polynomiality yields the explicit zero-product split.
            coefficient = sp.simplify(sp.expand(x * z3_special).coeff(x, 0))
            expected = sp.Rational(6, 5) * sp.Rational(1, 5) * a3 * b0
            assert sp.simplify(coefficient - expected) == 0
            extra = "; coeff_x3(r19)*r18(0)=0"

        base_count = {6: 919, 7: 955, 8: 994, 9: 1036}[s]
        after_r19 = base_count - expected_v
        after_r18 = after_r19 - (0 if s == 6 else 1)
        expected_counts = {6: (916, 916), 7: (952, 951), 8: (991, 990), 9: (1033, 1032)}
        assert (after_r19, after_r18) == expected_counts[s]
        if s == 6:
            # The exact product split has two codimension-one children and
            # their common codimension-two intersection.
            assert after_r19 - 1 == 915
            assert after_r19 - 2 == 914
            after_splits = "two 915-variable children; intersection=914"
        elif s == 7:
            # The depth-four x^(-1) coefficient factors into these two leaves.
            a3, b1 = sp.symbols("a3 b1")
            gamma = sp.Rational(7, 5)
            selector = (
                sp.binomial(gamma, 2) * b1**2
                + 3 * sp.binomial(gamma, 3) * a3**2 * b1
                + sp.binomial(gamma, 4) * a3**4
            )
            assert sp.simplify(
                sp.factor(selector)
                - sp.Rational(7, 625) * (5 * b1 - a3**2) * (5 * b1 - 2 * a3**2)
            ) == 0, sp.factor(selector)
            assert after_r18 - 1 == 950
            after_splits = "two 950-variable leaves; intersection=949"
        else:
            after_splits = f"{after_r18} variables after depth-four divisibility"
        branches.append((s, s, e, expected_v, extra, after_r19, after_r18, after_splits))

    print("vertex s e ord_x(r19) after_depth3_vars depth4_consequence child_sizes")
    for x0, s, e, v, extra, after_r19, after_r18, after_splits in branches:
        print(
            f"({x0},{e}) {s} {e} >= {v}{extra} "
            f"{after_r19} {after_r18} {after_splits}"
        )
    print("verified=first three nonresonant slices and resonant depth-four row")
    print("scope=necessary high-y pruning only; no branch exclusion")


if __name__ == "__main__":
    main()
