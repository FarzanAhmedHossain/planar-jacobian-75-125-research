"""Exact arithmetic of the F2 terminal type-I edge and its local selectors.

The attachment uses GGHV Theorem 2.20, its removal of type-III corners,
Theorem 3.1 and equation (3.16). This arithmetic checker does not replace
those theorem proofs and does not exclude the full degree pair.
"""

from __future__ import annotations

from fractions import Fraction as F
import json

import sympy as sp


def main():
    m, n, ell, aa, bb, kk = 3, 5, 5, 7, 2, 1
    assert (m+n)*bb*kk-n*(bb*ell-aa) == kk
    assert (m+n)*bb*kk-m*(bb*ell-aa) != kk
    p_end = (F(21,5), F(6))
    p_start = (F(4,5), F(1))
    q_end = (F(7), F(10))
    q_start = (F(1,5), F(0))
    weight = lambda point: 25*point[0]-17*point[1]
    assert weight(p_end) == weight(p_start) == 3
    assert weight(q_end) == weight(q_start) == 5
    assert 3+5 == 25-17
    assert ell % 25 != 0  # terminal type I.a is unavailable by Theorem 3.1
    # Type-III rho divides ell. Between initial slope 1/5 and terminal
    # slope 17/25, the only nontrivial possibilities are 2/5 and 3/5.
    collapsed = [(rho,sigma) for rho in range(1,ell+1)
                 for sigma in range(1,rho)
                 if ell%rho==0 and sp.gcd(rho,sigma)==1
                 and F(1,5)<F(sigma,rho)<F(17,25)]
    assert collapsed == [(5,2),(5,3)]
    # y -> y+c*x^(-sigma/rho), x=X^5 and u=X*y-1 give
    # u -> u+c*w^(5*sigma/rho-1), w=X^-1.
    shifts = [5*F(sigma,rho)-1 for rho,sigma in collapsed]
    assert shifts == [1,2]

    # x^i y^j = X^(5i-j)(1+u)^j. Factoring X^15 (or X^25)
    # yields the stated s. After the initial shear, x=X^5 and
    # Y=u/X, so a monomial X^k u^j has (25,-17)-weight 5k-12j.
    for i in range(26):
        for j in range(101):
            k=5*i-j
            assert 25*F(k+j,5)-17*j == 5*k-12*j
    p_forbidden=[(s,j) for s in range(15) for j in range(6)
                 if 5*s+12*j<72]
    q_forbidden=[(s,j) for s in range(24) for j in range(10)
                 if 5*s+12*j<120]
    assert len(p_forbidden)==53 and len(q_forbidden)==136
    assert [(s,j) for s in range(15) for j in range(7)
            if 5*s+12*j==72] == [(0,6),(12,1)]
    assert [(s,j) for s in range(25) for j in range(11)
            if 5*s+12*j==120] == [(0,10),(12,5),(24,0)]
    # This halfspace has no terms of ordinary total degree <6; the
    # only degree-six monomial is u^6.
    assert all(5*s+12*j<72 for s in range(6) for j in range(6-s))
    assert [(s,6-s) for s in range(7) if 5*s+12*(6-s)>=72] == [(0,6)]

    u,w,C,A,a,b,kappa = sp.symbols("u w C A a b kappa")
    # The generic ordinary tangent sextic must be a pure power. Higher
    # translations b*w^2 cannot change its ordinary homogeneous part.
    diagonal = [sp.binomial(6,s)*A**s/(6**s*C**(s-1)) for s in range(7)]
    tangent=sum(diagonal[s]*u**(6-s)*w**s for s in range(7))
    assert sp.expand(tangent-C*(u+A*w/(6*C))**6)==0
    centered=sp.expand(tangent.subs(u,u-A*w/(6*C)+b*w**2))
    total6=sum(coef*u**ij[0]*w**ij[1]
               for ij,coef in sp.Poly(centered,u,w).terms() if sum(ij)==6)
    assert sp.expand(total6-C*u**6)==0
    expected = {2:5*A**2/(12*C),3:5*A**3/(54*C**2),
                4:5*A**4/(432*C**3),5:A**5/(1296*C**4),
                6:A**6/(46656*C**5)}
    assert all(sp.simplify(diagonal[s]-v)==0 for s,v in expected.items())
    # The tuned quadratic-cube local face cannot be one of the omitted
    # multiplicity-six collapses. Integer-power shifts leave kappa^3
    # as its u^0*w^3 coefficient, forbidden by the terminal halfspace.
    tuned=sp.expand(((u+a*w+b*w**2)**2+kappa*w)**3)
    assert sp.Poly(tuned,u,w).coeff_monomial(w**3)==kappa**3

    # Exact Jacobian equations on the terminal boundary are consistent.
    z,AA,BB,CC,DD,EE=sp.symbols("z AA BB CC DD EE")
    f=AA*z**6+BB*z
    g=CC*z**10+DD*z**5+EE
    bracket=sp.expand((3*f*sp.diff(g,z)-5*sp.diff(f,z)*g)/25)
    expected_bracket=((25*BB*CC-15*AA*DD)*z**10
                      +(10*BB*DD-30*AA*EE)*z**5-5*BB*EE)/25
    assert sp.expand(bracket-expected_bracket)==0
    solution={CC:-9*AA**2/BB**3,DD:-15*AA/BB**2,EE:-5/BB}
    assert sp.simplify(bracket.subs(solution)-1)==0

    print(json.dumps({"terminal_normal":[25,-17],
                      "terminal_edge_weights":{"P":3,"Q":5},
                      "removable_typeIII_slopes":["2/5","3/5"],
                      "local_recentring":"u -> u+a*w+b*w^2",
                      "forbidden_coefficient_rows":{"P":53,"Q":136},
                      "necessary_orders_U1_to_U6":[5,4,3,2,1,0],
                      "tangent_cone":"C*(u+A*w/(6*C))^6",
                      "tuned_quadratic_cube":"excluded at the distinguished F2 root if kappa != 0",
                      "terminal_boundary_Jacobian":"consistent exact witness",
                      "status":"PASS",
                      "scope":"conditional on the cited GGHV terminal-chain attachment; no degree-pair exclusion"},indent=2))


if __name__=="__main__":
    main()
