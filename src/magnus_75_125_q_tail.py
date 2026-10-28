"""Exact high-y tail consequences of P=R^3+Z for the (75,125) model.

Inputs are the branchwise Magnus support model, deg_y(R)=20,
deg_y(Z)<=36, deg_y(Q)=100, and the top-slice relation
q_100=b*r_20^5 with b != 0.  The checker verifies the symbolic bracket
identities and the exact branch arithmetic behind the reduced Q-tail model.
The accompanying note gives the degree argument that turns these identities
into necessary support restrictions.  The final count also uses the exact
F2/Magnus edge bridge: four normalized R-edge coefficients are fixed by h,
the ell=15 edge of Z vanishes, and the ell=25 edge of the Q-tail S vanishes.
"""

import sympy as sp


def partitions(n, smallest=1):
    if n == 0:
        yield ()
        return
    for part in range(smallest, n + 1):
        for tail in partitions(n - part, part):
            yield (part,) + tail


def jacobian_bracket(f, g, x, y):
    return sp.diff(f, x) * sp.diff(g, y) - sp.diff(f, y) * sp.diff(g, x)


def main():
    x, y = sp.symbols("x y")
    r = sp.Function("r")(x, y)
    a = sp.Function("a")(x)
    k = sp.symbols("k", integer=True, nonnegative=True)

    # For a coefficient A(x), the exact bracket identity is
    # [R,A(x)R^k] = -A'(x)R^k R_y.
    bracket_identity = sp.simplify(
        jacobian_bracket(r, a * r**k, x, y)
        + sp.diff(a, x) * r**k * sp.diff(r, y)
    )
    assert bracket_identity == 0

    Zxy, Sxy = sp.Function("Z")(x, y), sp.Function("S")(x, y)
    b, c = sp.symbols("b c")
    q_decomposition_identity = sp.expand(
        jacobian_bracket(Zxy, b * r**5 + c * r**4 + Sxy, x, y)
        - jacobian_bracket(Zxy, Sxy, x, y)
        - (5 * b * r**4 + 4 * c * r**3)
        * jacobian_bracket(Zxy, r, x, y)
    )
    assert q_decomposition_identity == 0

    # Check the leading coefficient of [Z,R^5] at y^(e+99).
    leading_identities = {}
    for e in (24, 28, 32, 36):
        # Keep the leading coefficients as functions of x.  Terms one
        # y-degree lower cannot contribute to y^(e+99).
        rfun, zfun = sp.Function("r20")(x), sp.Function("z_e")(x)
        rlo, zlo = sp.symbols("r19 z_e_minus_1")
        Rf = rfun * y**20 + rlo * y**19
        Zf = zfun * y**e + zlo * y ** (e - 1)
        actual = sp.Poly(
            sp.expand(jacobian_bracket(Zf, Rf**5, x, y)), y
        ).coeff_monomial(y ** (e + 99))
        expected_fun = 5 * rfun**4 * (
            20 * rfun * sp.diff(zfun, x) - e * zfun * sp.diff(rfun, x)
        )
        assert sp.simplify(actual - expected_fun) == 0
        # The coefficient bound forces this expression to vanish, hence
        # (z_e^20 / r_20^e)' = 0 in characteristic zero.
        assert sp.simplify(
            sp.diff(zfun**20 / rfun**e, x)
            - (zfun**19 / rfun ** (e + 1))
            * (20 * rfun * sp.diff(zfun, x) - e * zfun * sp.diff(rfun, x))
        ) == 0
        leading_identities[e] = "20*r*z' - e*z*r' = 0"

    # The outer Q region has 1,576 monomials; the y<=60 support envelope has
    # 708.  Since Q=b*R^5+c*R^4+S and the exact edge bridge gives
    # Q_edge=b*(R_edge)^5, the c*R^4 term (weight at most 20) cannot meet the
    # weight-25 edge, so all 13 edge coefficients of S in this tail vanish.
    q_support = {
        (i, j)
        for i in range(26)
        for j in range(101)
        if 5 * i - j <= 25
    }
    q_low_tail = {(i, j) for i, j in q_support if j <= 60}
    q_tail_edge = {(i, j) for i, j in q_low_tail if 5 * i - j == 25}
    q_free_tail = q_low_tail - q_tail_edge
    assert len(q_support) == 1576
    assert len(q_low_tail) == 708
    assert q_tail_edge == {(i, 5 * i - 25) for i in range(5, 18)}
    assert len(q_tail_edge) == 13
    assert len(q_free_tail) == 695

    # The cube-root edge bridge fixes all five ell=5 coefficients of R;
    # r_(5,20)=1 is the corner normalization, leaving 71 other R variables.
    r_support = {
        (i, j)
        for i in range(6)
        for j in range(21)
        if 5 * i - j <= 5
    }
    r_edge = {(i, j) for i, j in r_support if 5 * i - j == 5}
    assert len(r_support) == 76
    assert r_edge == {(i, 5 * i - 5) for i in range(1, 6)}
    r_free_after_edge = len(r_support) - 1 - (len(r_edge) - 1)
    assert r_free_after_edge == 71

    p_support = {
        (i, j)
        for i in range(16)
        for j in range(61)
        if 5 * i - j <= 15
    }

    branches = [
        (1, 3, 8), (1, 4, 10), (2, 7, 24), (2, 8, 27),
        (3, 12, 52), (4, 16, 80), (5, 20, 111),
        (6, 24, 145), (7, 28, 182), (8, 32, 222), (9, 36, 265),
    ]
    # The raw support relaxation has 75 normalized R coefficients and 708 S
    # coefficients.  The exact edge bridge reduces these to 71 and 695;
    # it also fixes every Z coefficient on the edge ell=15 in the vertex box.
    print("vertex Z_terms Z_edge_pins raw_vars edge_refined_vars normalized_edge_refined_vars")
    for x0, e, z_terms in branches:
        z_support = {
            (i, j)
            for i in range(x0 + 1)
            for j in range(e + 1)
            if (i, j) in p_support
        }
        assert len(z_support) == z_terms
        z_edge = {(i, 5 * i - 15) for i in range(3, x0 + 1)
                  if 5 * i - 15 <= e}
        assert z_edge <= z_support
        if e > 20:
            allowed = [
                part for part in partitions(5)
                if all((e * m) % 20 == 0 for m in part)
            ]
        else:
            allowed = list(partitions(5))
        if e in (24, 28, 32, 36):
            assert allowed == [(5,)]
            assert x0 == e // 4
        variables = 75 + z_terms + 708 + 2
        edge_refined_variables = r_free_after_edge + (z_terms - len(z_edge)) + len(q_free_tail) + 2
        normalized_edge_refined_variables = edge_refined_variables
        if e in (24, 28, 32, 36):
            # The allowed x-translation x -> x + alpha preserves all support
            # envelopes and sends rho=(x-alpha)^5 to x^5.  This kills the five
            # nonleading coefficients of rho and the x0 lower coefficients
            # in the top y-slice z_e=mu*(x-alpha)^(e/4).
            normalized_edge_refined_variables -= 5 + x0
        print(
            f"({x0},{e}) {z_terms} {len(z_edge)} "
            f"{variables} {edge_refined_variables} "
            f"{normalized_edge_refined_variables}"
        )

    print("edge_R_free=71 (76 supported coefficients, one corner normalization, four h-pins)")
    print("edge_Z_pins=ell15 intersection with each branch rectangle")
    print("edge_S_pins=13 points (i,5*i-25), i=5..17; S_free_tail=695")
    print("high_slice_identity=Q=b*R^5+c*R^4+S, deg_y(S)<=60")
    print("high_branch_identity=deg_y([Z,R])<=39")
    print("low_y_congruence=R^2 divides 1-[Z,S]")
    print("high_vertices_force=r_20=scalar*(x-alpha)^5; z_e=scalar*(x-alpha)^(e/4)")
    print("after_x_translation=last four branches have r_20=x^5 and z_e=mu*x^(e/4)")
    print("scope=exact high-y support consequences; not a Keller-pair exclusion")


if __name__ == "__main__":
    main()
