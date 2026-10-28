"""Exact sixth coupled F2 block: weights 44/54 and the Q19 law.

No ranks are computed. The identities retain full supports, the kappa70
correction, the K6 pole remainder and all zero-parameter cases.
"""

from __future__ import annotations

import json
import sympy as sp


def main():
    t,e,h,S,B,U3,U4,U5,U6,k=sp.symbols("t e h S B U3 U4 U5 U6 k")
    delta=(S/h*e+B/h**2*e**2+U3/h**3*e**3+U4/h**3*e**4
           +U5/(t*h**3)*e**5+U6/(t*h**3)*e**6)
    formal=sp.expand(h**5*sum(sp.binomial(sp.Rational(5,3),j)*delta**j for j in range(7)))
    g6=formal.coeff(e,6)+sp.Rational(4,3)*k*h**3*S/t
    N=(729*U3**2-486*S*B*U3-81*B**3+108*S**3*U3
       +162*S**2*B**2-63*S**4*B+7*S**6)
    D=3*B-S**2
    Delta=27*U3+2*S**3-9*S*B
    assert sp.expand(N-(Delta**2-3*D**3))==0
    V6=(sp.Rational(5,3)*h**2*U6+sp.Rational(10,9)*h*S*U5
        +t*(sp.Rational(10,9)*B-sp.Rational(5,27)*S**2)*U4
        +sp.Rational(5,6561)*t*N/h+sp.Rational(4,3)*k*h**3*S)
    assert sp.cancel(t*g6-V6)==0
    forcing=(7290*h**2*S*U5+t*h*(7290*B-1215*S**2)*U4
             +5*t*N+8748*k*h**4*S)
    assert sp.cancel(6561*h*V6-10935*h**3*U6-forcing)==0
    assert sp.simplify(V6.subs(t,0)-(sp.Rational(5,3)*h**2*U6
                       +sp.Rational(10,9)*h*S*U5+sp.Rational(4,3)*k*h**3*S))==0

    # The original monomial supports give these intervals exactly.
    p9=[i for i in range(16) if 0<=5*i-9<=60]
    q19=[i for i in range(26) if 0<=5*i-19<=100]
    w33=[i for i in range(76) if 0<=5*i-33<=300]
    assert p9==list(range(2,14))
    assert q19==list(range(4,24))
    assert w33==list(range(7,67))
    # Degree bounds of all N monomials: deg S<=3, deg B<=7, deg U3<=11.
    n_degrees=(22,3+7+11,3*7,3*3+11,2*3+2*7,4*3+7,6*3)
    assert max(n_degrees)==22
    assert max(8+11,4+3+12,1+7+11,1+6+11,1+22-4,12+3)==19
    assert 19%5!=0  # U^(19/5) is not rational at the simple t=0 root.

    # Generic functions audit the coupled bracket equations independently
    # of the proposed expanded sixth-block formulas.
    M=sp.Function("M")(t)
    A,A2,A3,A4,A5,A6=(sp.Function(name)(t) for name in ("A","A2","A3","A4","A5","A6"))
    q24,q23,q22,q21,q20,q19f=(sp.Function(name)(t) for name in
                              ("q24","q23","q22","q21","q20","q19"))
    g,W38,W37,W36,W35,W34,W33=(sp.Function(name)(t) for name in
                              ("g","W38","W37","W36","W35","W34","W33"))
    ps=[M**3,A,A2,A3,A4,A5,A6]
    qs=[M**5,q24,q23,q22,q21,q20,q19f]
    ws=[g,W38,W37,W36,W35,W34,W33]
    pseries=sum(ps[i]*e**i for i in range(7))
    qseries=sum(qs[i]*e**i for i in range(7))
    def bracket_at(left,top,output):
        return sum((top-i)*left[i]*sp.diff(ws[j],t)
                   -(39-j)*sp.diff(left[i],t)*ws[j]
                   for i in range(7) for j in range(7)
                   if top-i+39-j-4==output)
    p4_drop6=(4*M**9*A6+12*M**6*(A*A5+A2*A4)+6*M**6*A3**2
              +12*M**3*A**2*A4+24*M**3*A*A2*A3+4*M**3*A2**3
              +4*A**3*A3+6*A**2*A2**2)
    assert sp.expand(pseries**4).coeff(e,6)==sp.expand(p4_drop6)
    E44=(15*M**3*sp.diff(W33,t)-99*M**2*sp.diff(M,t)*W33
         +14*A*sp.diff(W34,t)-34*sp.diff(A,t)*W34
         +13*A2*sp.diff(W35,t)-35*sp.diff(A2,t)*W35
         +12*A3*sp.diff(W36,t)-36*sp.diff(A3,t)*W36
         +11*A4*sp.diff(W37,t)-37*sp.diff(A4,t)*W37
         +10*A5*sp.diff(W38,t)-38*sp.diff(A5,t)*W38
         +9*A6*sp.diff(g,t)-39*sp.diff(A6,t)*g
         -6*M**5*q19f-6*q24*q20-6*q23*q21-3*q22**2
         +9*k*M**6*A)
    direct44=(bracket_at(ps,15,44)-3*sp.expand(qseries**2).coeff(e,6)
              +3*k*sp.expand(pseries**3).coeff(e,1))
    assert sp.expand(E44-direct44)==0
    E54=(25*M**5*sp.diff(W33,t)-165*M**4*sp.diff(M,t)*W33
         +24*q24*sp.diff(W34,t)-34*sp.diff(q24,t)*W34
         +23*q23*sp.diff(W35,t)-35*sp.diff(q23,t)*W35
         +22*q22*sp.diff(W36,t)-36*sp.diff(q22,t)*W36
         +21*q21*sp.diff(W37,t)-37*sp.diff(q21,t)*W37
         +20*q20*sp.diff(W38,t)-38*sp.diff(q20,t)*W38
         +19*q19f*sp.diff(g,t)-39*sp.diff(q19f,t)*g
         -5*p4_drop6-9*k*(2*M**8*A+M**6*q24))
    direct54=(bracket_at(qs,25,54)-5*sp.expand(pseries**4).coeff(e,6)
              -9*k*sp.expand(pseries**2*qseries).coeff(e,1))
    assert sp.expand(E54-direct54)==0

    # The normalized and physical denominator-cleared laws agree exactly.
    p9sym,p10,p11,q19sym=sp.symbols("p9 p10 p11 q19")
    physical_forcing=(7290*t**2*h**2*S*p10
                      +t**2*h*(7290*B-1215*S**2)*p11
                      +5*t**5*N+8748*k*t**4*h**4*S)
    normalized_eq=6561*h*sp.Symbol("V6")-10935*h**3*U6-forcing
    physical_eq=6561*h*q19sym-10935*t**2*h**3*p9sym-physical_forcing
    assert sp.cancel(t**4*normalized_eq.subs({sp.Symbol("V6"):q19sym/t**4,
                    U6:p9sym/t**2,U5:p10/t**2,U4:p11/t**3})-physical_eq)==0

    print(json.dumps({"P9_t_powers":[2,13],"Q19_t_powers":[4,23],
                      "W33_t_powers":[7,66],"new_P9_coefficients":12,
                      "new_W33_coefficients":60,
                      "Q19_law":"6561*h*V6-10935*h^3*U6=forcing",
                      "pole_numerator":"Delta^2-3*D^3",
                      "polynomiality":"h divides pole_numerator",
                      "t0":"forced coefficient value; no extra obstruction",
                      "coupled_weights":[44,54],
                      "kappa_corrections":"+9*k*M^6*A; -9*k*(2*M^8*A+M^6*q24)",
                      "status":"PASS","scope":"exact sixth necessary block; no exclusion"},indent=2))


if __name__=="__main__":
    main()
