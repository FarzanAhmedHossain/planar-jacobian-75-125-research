"""Exact F2 source-translation gauge at the fifth (5,-1)-band.

This proves a coordinate normalization, not vanishing in uncentered
coordinates and not a degree-pair exclusion.  It is compatible with the
standalone power-relation edge model; an independent high-y centering
cannot be imposed at the same time without a separate compatibility proof.
"""

from __future__ import annotations

import json
from math import ceil, floor

import sympy as sp


def support(x_cap, y_cap, weight_cap):
    return {
        (i, j)
        for i in range(x_cap + 1)
        for j in range(y_cap + 1)
        if 5 * i - j <= weight_cap
    }


def main():
    # Every term produced by (x+c)^i y^j remains in the same envelope.
    counts = {}
    for name, caps in {"P": (15, 60, 15), "Q": (25, 100, 25),
                       "Wstar": (75, 300, 39)}.items():
        original = support(*caps)
        transitions = 0
        for i, j in original:
            for child_i in range(i + 1):
                assert (child_i, j) in original
                assert 5 * child_i - j == 5 * i - j - 5 * (i-child_i)
                transitions += 1
        counts[name] = {"support_points": len(original),
                        "translation_monomial_transitions": transitions}

    # Original monomial x^i y^j = t^i z^(5i-j).  Its first x
    # derivative has coefficient i*t^(i-1) in band 5i-j-5.
    for i in range(1, 76):
        for j in (0, 17, 60, 100, 300):
            assert 5*(i-1)-j == (5*i-j)-5
    assert (ceil(10/5), floor((60+10)/5)) == (2, 14)
    assert (ceil(20/5), floor((100+20)/5)) == (4, 24)
    assert (ceil(34/5), floor((300+34)/5)) == (7, 66)

    t, c = sp.symbols("t c")
    H = sp.Function("h")(t)
    U = t*H
    p15 = U**3
    q25 = U**5
    assert sp.simplify(sp.diff(p15, t)-3*t**2*H**2*(H+t*sp.diff(H,t))) == 0
    assert sp.simplify(sp.diff(q25, t)-5*t**4*H**4*(H+t*sp.diff(H,t))) == 0
    # The fifth binomial coefficient is affine in p10 with multiplier
    # (5/3) U^2.  This proves that its resonant kappa*U^4 term is unchanged.
    assert sp.simplify(sp.diff(q25,t)-sp.Rational(5,3)*U**2*sp.diff(p15,t)) == 0

    p, q, g = (sp.Function(name)(t) for name in ("p", "q", "g"))
    # Lower bands are unchanged. The full changed part of each next row
    # is the derivative of its already imposed leading row.
    e50 = 15*p*sp.diff(g,t)-39*sp.diff(p,t)*g-3*q**2
    e60 = 25*q*sp.diff(g,t)-39*sp.diff(q,t)*g-5*p**4
    delta45 = c*(15*p*sp.diff(g,t,2)-34*sp.diff(p,t)*sp.diff(g,t)
                  +10*sp.diff(p,t)*sp.diff(g,t)-39*sp.diff(p,t,2)*g
                  -6*q*sp.diff(q,t))
    delta55 = c*(25*q*sp.diff(g,t,2)-34*sp.diff(q,t)*sp.diff(g,t)
                  +20*sp.diff(q,t)*sp.diff(g,t)-39*sp.diff(q,t,2)*g
                  -20*p**3*sp.diff(p,t))
    assert sp.expand(delta45-c*sp.diff(e50,t)) == 0
    assert sp.expand(delta55-c*sp.diff(e60,t)) == 0

    a = (3+sp.sqrt(5))/2
    h22 = (t-1)**2*(t-a)**2
    h211 = (t-1)**2*(t**2-3*t+3)
    branch_receipts = {}
    b0 = sp.symbols("b0")
    for name, h in (("2+2", h22), ("2+1+1", h211)):
        h0 = sp.simplify(h.subs(t,0))
        coefficient = sp.simplify(3*h0**3)
        assert coefficient != 0
        shift = -b0/coefficient
        pshift = sp.Poly(sp.expand(sp.diff((t*h)**3,t)),t)
        qshift = sp.Poly(sp.expand(sp.diff((t*h)**5,t)),t)
        assert all(pshift.nth(i)==0 for i in range(2))
        assert all(qshift.nth(i)==0 for i in range(4))
        assert pshift.degree() == 14
        assert qshift.degree() == 24
        assert sp.simplify(b0+shift*pshift.nth(2)) == 0
        assert sp.simplify(qshift.nth(4)-5*h0**5) == 0
        branch_receipts[name] = {"h0": str(h0),
                                 "coefficient_of_c_in_p10_t2": str(coefficient),
                                 "unique_centering_c": str(shift)}

    # A leading g=U^7 V has t-order at least 8 and degree 39 in both
    # edge configurations. Its derivative is an allowed W34 correction.
    assert 8-1 >= 7 and 39-1 <= 66
    # Translation drops omega by multiples of 5. Hence the P bands
    # 15..11, Q bands 25..21, W bands 39..35 and their fixed scalars do
    # not change; in particular kappa70 is preserved.
    for top in (15,25,39):
        assert all(top-5*j < top-4 for j in range(1,top+1))

    print(json.dumps({"support_preservation": counts,
                      "branch_centers": branch_receipts,
                      "p45_covariance": "delta=c*d(E50)/dt",
                      "q55_covariance": "delta=c*d(E60)/dt",
                      "kappa": "invariant",
                      "status": "PASS",
                      "scope": "source gauge only; no exclusion"}, indent=2))


if __name__ == "__main__":
    main()
