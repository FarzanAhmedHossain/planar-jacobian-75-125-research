"""Universal small-jet checks for the F2 terminal a,b elimination.

Checks every forbidden coefficient through s=5 against direct substitution
of generic polynomials, and verifies the compatibility of source x-centering.
"""

from __future__ import annotations

from math import factorial
import json

import sympy as sp


def center_row(polys, s, u, a, b):
    row = sp.Integer(0)
    for i in range(s+1):
        for n in range((s-i)//2+1):
            m = s-i-2*n
            row += a**m*b**n*sp.diff(polys[i],u,m+n)/(factorial(m)*factorial(n))
    return sp.expand(row)


def through_w5(expr,w):
    pp=sp.Poly(sp.expand(expr),w)
    return sp.expand(sum(pp.nth(s)*w**s for s in range(6)))


def main():
    u,w,a,b,c=sp.symbols("u w a b c")
    polys=[]
    for i in range(6):
        # Keep A1 generic, including its forbidden coefficients. This
        # arithmetic audit must not assume the S(1)=0 selector itself.
        low=6 if i==0 else 0
        polys.append(sum(sp.Symbol(f"c{i}_{j}")*u**j for j in range(low,9)))
    C=sp.Symbol("c0_6")
    A=sp.Symbol("c1_5")
    a_value=-A/(6*C)
    r2=sp.expand(polys[2]+a_value*sp.diff(polys[1],u)
                 +a_value**2*sp.diff(polys[0],u,2)/2)
    b_value=-r2.coeff(u,5)/(6*C)

    # Direct polynomial substitution is independent of center_row's
    # differentiated-coefficient construction. Degree 8 is universal for
    # all forbidden jets through s=5: j+m+n <= j+s <= 8 there.
    original=sum(polys[i]*w**i for i in range(6))
    direct=sp.Poly(sp.expand(original.subs(u,u+a*w+b*w**2)),w,u)
    checked=0
    for s in range(6):
        formula=center_row(polys,s,u,a,b)
        max_j=(71-5*s)//12
        assert max_j+s<=8
        for j in range(max_j+1):
            assert sp.expand(direct.coeff_monomial(w**s*u**j)-formula.coeff(u,j))==0
            checked+=1

    row1=center_row(polys,1,u,a,b)
    assert sp.simplify(row1.coeff(u,5).subs(a,a_value))==0
    row2=center_row(polys,2,u,a,b)
    assert sp.simplify(row2.coeff(u,5).subs({a:a_value,b:b_value}))==0
    for j in range(5):
        assert sp.simplify(row2.coeff(u,j).subs(a,a_value)-r2.coeff(u,j))==0
        assert sp.diff(row2.coeff(u,j),b)==0
    expected_b=-(sp.Symbol("c2_5")+6*a_value*sp.Symbol("c1_6")
                  +21*a_value**2*sp.Symbol("c0_7"))/(6*C)
    assert sp.simplify(b_value-expected_b)==0

    # Verify the explicit row 3--5 coefficients, including mixed a,b terms.
    expected3=(polys[3]+a*sp.diff(polys[2],u)+b*sp.diff(polys[1],u)
               +a*a*sp.diff(polys[1],u,2)/2+a*b*sp.diff(polys[0],u,2)
               +a**3*sp.diff(polys[0],u,3)/6)
    expected4=(polys[4]+a*sp.diff(polys[3],u)+b*sp.diff(polys[2],u)
               +a*a*sp.diff(polys[2],u,2)/2+a*b*sp.diff(polys[1],u,2)
               +a**3*sp.diff(polys[1],u,3)/6+b*b*sp.diff(polys[0],u,2)/2
               +a*a*b*sp.diff(polys[0],u,3)/2+a**4*sp.diff(polys[0],u,4)/24)
    expected5=(polys[5]+a*sp.diff(polys[4],u)+b*sp.diff(polys[3],u)
               +a*a*sp.diff(polys[3],u,2)/2+a*b*sp.diff(polys[2],u,2)
               +a**3*sp.diff(polys[2],u,3)/6+b*b*sp.diff(polys[1],u,2)/2
               +a*a*b*sp.diff(polys[1],u,3)/2+a**4*sp.diff(polys[1],u,4)/24
               +a*b*b*sp.diff(polys[0],u,3)/2
               +a**3*b*sp.diff(polys[0],u,4)/6+a**5*sp.diff(polys[0],u,5)/120)
    assert sp.expand(center_row(polys,3,u,a,b)-expected3)==0
    assert sp.expand(center_row(polys,4,u,a,b)-expected4)==0
    assert sp.expand(center_row(polys,5,u,a,b)-expected5)==0

    # Source x->x+c has v=(1+c*w^5)^(1/5). Modulo w^6,
    # v=1+c*w^5/5, v^15=1+3*c*w^5, w/v=w (mod w^6).
    source_shifted=through_w5((1+3*c*w**5)*original.subs(u,u+c*(1+u)*w**5/5),w)
    delta5=3*polys[0]+(1+u)*sp.diff(polys[0],u)/5
    assert sp.expand(source_shifted-original-c*w**5*delta5)==0
    assert all(sp.expand(source_shifted).coeff(w,i)==polys[i] for i in range(5))
    assert all(sp.expand(delta5).coeff(u,j)==0 for j in range(5))
    shifted_polys=[*polys[:5],sp.expand(polys[5]+c*delta5)]
    for s in range(6):
        difference=sp.expand(center_row(shifted_polys,s,u,a,b)-center_row(polys,s,u,a,b))
        for j in range((71-5*s)//12+1):
            assert difference.coeff(u,j)==0

    print(json.dumps({"generic_forbidden_rows_through_s5_checked":checked,
                      "generic_u_degree":8,
                      "a":"-A1_u5/(6*A0_u6)",
                      "b":"-(A2_u5+6*a*A1_u6+21*a^2*A0_u7)/(6*A0_u6)",
                      "mixed_Taylor_coefficients_s3_to_s5":"PASS",
                      "A1_low_coefficients":"retained as generic variables; not assumed zero",
                      "x_centering":"a,b unchanged; forbidden rows unchanged through s5",
                      "status":"PASS",
                      "scope":"exact local centering arithmetic; no exclusion"},indent=2))


if __name__=="__main__":
    main()
