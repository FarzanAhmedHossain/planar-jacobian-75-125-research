"""Audit a complete zero/nonzero split in the s=9, u0=u1=0 tail leaf.

The depth-4 and depth-5 fixed quotient coefficients give two homogeneous
linear equations in (u2,c2).  Their determinant is relevant only on the
nonzero (u2,c2) branch; the zero vector is retained as a separate leaf.
"""

from __future__ import annotations

import sympy as sp


def main() -> None:
    a3, rho1, u2, c2, u3, c3, t = sp.symbols(
        "a3 rho1 u2 c2 u3 c3 t"
    )

    # Depth-4 y^(e+56) fixed coefficient, after u0=u1=0.
    row4 = (
        (1828*a3**3 - 7620*a3*rho1)*u2
        + (-1740*a3**2 + 7075*rho1)*c2
    )

    # Depth-5 x^4 fixed coefficient (discarding its nonzero scalar 9*mu/625).
    row5 = (
        (9932*a3**4 - 42420*a3**2*rho1 + 3600*rho1**2)*u2
        + (-9420*a3**3 + 38975*a3*rho1)*c2
    )

    matrix = sp.Matrix([
        [sp.diff(row4, u2), sp.diff(row4, c2)],
        [sp.diff(row5, u2), sp.diff(row5, c2)],
    ])
    determinant = sp.factor(matrix.det())
    expected = -360 * (
        -172*a3**6 + 2925*a3**4*rho1
        - 26100*a3**2*rho1**2 + 70750*rho1**3
    )
    assert sp.expand(determinant - expected) == 0

    # Exact branch logic: for a homogeneous 2x2 system, either its vector is
    # zero or the coefficient matrix is singular.  Do not assert det=0 alone.
    assert sp.expand(row4.subs({u2: 0, c2: 0})) == 0
    assert sp.expand(row5.subs({u2: 0, c2: 0})) == 0

    # On the nonzero-vector leaf with a3=0, singularity forces rho1=0.
    assert sp.expand(expected.subs(a3, 0) - (-360*70750*rho1**3)) == 0

    # If the first vector vanishes, the next depth-4/depth-5 rows give a
    # second homogeneous system in (u3,c3).
    row4_next = (
        (1700*a3**3 - 7325*a3*rho1)*u3
        + (-1620*a3**2 + 6925*rho1)*c3
    )
    row5_next = (
        (2305*a3**4 - 9975*a3**2*rho1 + 500*rho1**2)*u3
        + (-2195*a3**3 + 9375*a3*rho1)*c3
    )
    matrix_next = sp.Matrix([
        [sp.diff(row4_next, u3), sp.diff(row4_next, c3)],
        [sp.diff(row5_next, u3), sp.diff(row5_next, c3)],
    ])
    determinant_next = sp.factor(matrix_next.det())
    expected_next = -50 * (
        -52*a3**6 + 2115*a3**4*rho1
        - 24300*a3**2*rho1**2 + 69250*rho1**3
    )
    assert sp.expand(determinant_next - expected_next) == 0
    assert sp.expand(row4_next.subs({u3: 0, c3: 0})) == 0
    assert sp.expand(row5_next.subs({u3: 0, c3: 0})) == 0

    cubic2 = 70750*t**3 - 26100*t**2 + 2925*t - 172
    cubic3 = 69250*t**3 - 24300*t**2 + 2115*t - 52
    discriminant2 = sp.factor(sp.discriminant(cubic2, t))
    discriminant3 = sp.factor(sp.discriminant(cubic3, t))
    assert discriminant2 == -762376168125000
    assert discriminant3 == 17335015875000
    assert sp.gcd(cubic2, sp.diff(cubic2, t)) == 1
    assert sp.gcd(cubic3, sp.diff(cubic3, t)) == 1
    assert sp.gcd(cubic2, cubic3) == 1

    print("branch_split=(u2,c2)=(0,0) OR ((u2,c2)!=(0,0) AND det=0)")
    print(f"determinant={determinant}")
    print("nonzero_leaf_a3_zero_implies=rho1=0")
    print(f"nonzero_first_pair_a3_nonzero_ratio_equation={cubic2}=0, t=rho1/a3^2")
    print(f"first_ratio_cubic_discriminant={discriminant2}")
    print("(u2,c2)=0 branch_split=(u3,c3)=0 OR ((u3,c3)!=(0,0) AND det_next=0)")
    print(f"determinant_next={determinant_next}")
    print(f"nonzero_second_pair_a3_nonzero_ratio_equation={cubic3}=0, t=rho1/a3^2")
    print(f"second_ratio_cubic_discriminant={discriminant3}")
    print("ratio_cubics_have_common_root=False")
    print("scope=necessary splits only; neither child is excluded")


if __name__ == "__main__":
    main()
