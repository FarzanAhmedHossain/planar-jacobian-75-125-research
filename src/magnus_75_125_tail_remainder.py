"""Top-slice consequences of R^2 | 1-[Z,S] in the high Magnus branches."""

import sympy as sp


def bracket(f, g, x, y):
    return sp.diff(f, x) * sp.diff(g, y) - sp.diff(f, y) * sp.diff(g, x)


def main():
    x, y = sp.symbols("x y")
    ztop = sp.Function("z_e")(x)
    slow = sp.Function("s_60")(x)
    zlo = sp.Function("z_lower")(x)
    slo = sp.Function("s_lower")(x)

    for e in (24, 28, 32, 36):
        Z = ztop * y**e + zlo * y ** (e - 1)
        S = slow * y**60 + slo * y**59
        top = sp.Poly(sp.expand(bracket(Z, S, x, y)), y).coeff_monomial(y ** (e + 59))
        assert sp.simplify(top - (60 * sp.diff(ztop, x) * slow - e * ztop * sp.diff(slow, x))) == 0
        next_slice = sp.Poly(sp.expand(bracket(Z, S, x, y)), y).coeff_monomial(y ** (e + 58))
        expected_next = (
            59 * sp.diff(ztop, x) * slo - e * ztop * sp.diff(slo, x)
            + 60 * sp.diff(zlo, x) * slow - (e - 1) * zlo * sp.diff(slow, x)
        )
        assert sp.simplify(next_slice - expected_next) == 0

    r19fun, r18fun = sp.Function("r19")(x), sp.Function("r18")(x)
    Rhigh = x**5 * y**20 + r19fun * y**19 + r18fun * y**18
    R2 = sp.Poly(sp.expand(Rhigh**2), y)
    assert sp.simplify(R2.coeff_monomial(y**40) - x**10) == 0
    assert sp.simplify(R2.coeff_monomial(y**39) - 2 * x**5 * r19fun) == 0

    mu = sp.symbols("mu", nonzero=True)
    base_counts = {24: 919, 28: 955, 32: 994, 36: 1036}
    print("vertex forced_S60_power forced_S59_power local_reduction_count leaf_counts")
    for e in (24, 28, 32, 36):
        s = e // 4
        S_coeffs = sp.symbols(f"S60_0:{18}")
        s60 = sum(S_coeffs[i] * x**i for i in range(18))
        z_e = mu * x**s
        leading = sp.expand(60 * sp.diff(z_e, x) * s60 - e * z_e * sp.diff(s60, x))
        assert sp.expand(
            leading
            - 4 * s * mu * sum((15 - i) * S_coeffs[i] * x ** (s + i - 1) for i in range(18))
        ) == 0

        forced = tuple(i for i in range(18) if s + i - 1 < 10)
        assert forced == tuple(range(11 - s))
        assert all(15 - i != 0 for i in forced)

        # For the next y coefficient, write R^2*T with top quotient slice
        # t_(e+19).  Its only known correction is 2*x^5*r19*t_(e+19),
        # divisible by x^8 because x^3|r19.
        a = sp.symbols(f"a0:{2}")
        b = sp.symbols(f"b0:{s + 7}")
        c = sp.symbols("c0:17")
        d = sp.symbols("d0:8")
        r19 = x**3 * (a[0] + a[1] * x)
        s60 = x ** (11 - s) * sum(b[i] * x**i for i in range(s + 7))
        s59 = sum(c[i] * x**i for i in range(17))
        ttop = sum(d[i] * x**i for i in range(8))
        ztop = mu * x**s
        zprev = mu * sp.Rational(s, 5) * x ** (s - 5) * r19
        s59_operator = 59 * sp.diff(ztop, x) * s59 - e * ztop * sp.diff(s59, x)
        cross = 60 * sp.diff(zprev, x) * s60 - (e - 1) * zprev * sp.diff(s60, x)
        quotient_correction = 2 * x**5 * r19 * ttop
        cross_poly = sp.Poly(sp.expand(cross), x)
        correction_poly = sp.Poly(sp.expand(quotient_correction), x)
        assert all(cross_poly.nth(k) == 0 for k in range(8))
        assert all(correction_poly.nth(k) == 0 for k in range(8))

        forced59 = tuple(i for i in range(17) if s + i - 1 < 8)
        assert forced59 == tuple(range(max(0, 9 - s)))
        for i in forced59:
            assert 59 - 4 * i != 0
            assert sp.Poly(s59_operator, x).nth(s + i - 1) == mu * s * (59 - 4 * i) * c[i]

        reduction = 11 - s
        if e == 24:
            # x^3|r19, x^5|S60, x^3|S59, and a zero-product split.
            count = base_counts[e] - 3 - reduction - len(forced59)
            assert count == 908
            leaves = "907, 907; intersection 906"
        elif e == 28:
            # x^3|r19, x|r18, x^4|S60, x^2|S59, and selector branch.
            count = base_counts[e] - 3 - 1 - reduction - 1 - len(forced59)
            assert count == 944
            leaves = "944, 944; intersection 943"
        elif e == 32:
            count = base_counts[e] - 3 - 1 - reduction - len(forced59)
            assert count == 986
            leaves = "986"
        else:
            count = base_counts[e] - 3 - 1 - reduction - len(forced59)
            assert count == 1030
            leaves = "1030"

        s59_power = max(0, 9 - s)
        print(f"({s},{e}) x^{11-s}|S_60 x^{s59_power}|S_59 {count} {leaves}")

    print("identity=R^2 | 1-[Z,S] forces the y^(e+59) slice divisible by x^10")
    print("next_slice=the y^(e+58) remainder forces x^(9-s)|S_59 when s<9")
    print("scope=necessary support pruning only; no branch exclusion")


if __name__ == "__main__":
    main()
