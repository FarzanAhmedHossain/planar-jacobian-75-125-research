# Source and coordinate audit for the surviving F2 branch

Date: 2026-10-18. This is a new independent audit note. It does not exclude
the two-factor branch or the degree pair (75,125).

## A missing fifth-band coefficient has an available source gauge

The full F2 support permits
`P10 in span(t^2,...,t^14)` and `Q20 in span(t^4,...,t^24)`.
Consequently the existing fifth-band script's `P10=t^3 U5`, `deg U5<=11`,
does not describe the uncentered coefficient space. No inspected root
polynomiality condition proves its omitted `t^2` coefficient zero.

There is, however, an exact additional coordinate normalization. Apply the
source automorphism `T_c: x -> x+c, y -> y`. For each permitted monomial
`x^i y^j`, every translated child is `x^k y^j`, `0<=k<=i`, so the rectangle
and inequality `5i-j<=15` (respectively 25) remain satisfied. Degrees,
upper-right coefficients, the top edge and its quartic `h` remain unchanged.
The Jacobian constant is preserved, since `J(T_c)=1`.

In the band chart `x=t z^5`, `y=z^-1`, differentiation with respect to `x`
lowers the band weight by five and differentiates its coefficient in `t`.
With `U=t h`, `P15=U^3`, `Q25=U^5`, this gives exactly

```
P10 -> P10 + c (U^3)'
Q20 -> Q20 + c (U^5)'.
```

The coefficient of `t^2` in `(U^3)'` is `3 h(0)^3`, nonzero. Thus if `b0`
denotes the old `t^2` coefficient of `P10`, the unique choice
`c=-b0/(3h(0)^3)` gives `P10=t^3 U5`, `deg U5<=11`.
All bands `P15,...,P11` and `Q25,...,Q21` are fixed by this translation.
The fifth binomial coefficient changes by
`(5/3)U^2 c(U^3)'=c(U^5)'`; hence its additional `kappa U^4` term has the
same `kappa`. In particular the resonant `t^4` term of `Q20` remains in the
existing script after this centering.

The corrected power relation transforms by the same source automorphism,
with all scalar correction coefficients unchanged. If its leading band is
`W39=g`, then `W39,...,W35` are fixed and
`W34 -> W34+c g'`. Here `ord_0(g)>=8`, `deg(g)=39`, so `g'` is permitted by
the relaxed `W34` support `span(t^7,...,t^66)`.

The covariance of the two fifth-block equations is particularly direct.
Set

```
E50=15 P15 g' -39 P15' g -3 Q25^2,
E60=25 Q25 g' -39 Q25' g -5 P15^4.
```

Under the displayed shifts the full changes of the weight-45 and weight-55
rows are respectively `c E50'` and `c E60'`. Both vanish after the already
imposed leading equations. The scalar `kappa70` is invariant as well.
Thus the gauge-fixed fifth-band model can represent every source-translation
orbit of its standalone edge system. It must be explicitly called a
coordinate normalization, rather than a forced zero coefficient.

This gauge is independent of the other derivations' choice to move the
root of the upper-y carrier to zero. Those two x-centerings generally
consume the same source freedom and cannot be imposed together without
transporting the resulting coefficients. This note establishes the edge
gauge only. It does not identify the centered edge variables with the
separately centered high-y/Magnus variables.

The exact checker is
`f2_x_translation_fifth_band_audit_20261018.py`. On 2026-10-18 it passed:
support closure for all 586 P points, 1,576 Q points and 11,558 relaxed W*
points; all 291,423 translated monomial transitions; both functional
covariance identities; and the nonzero centering coefficients for both
quartic edge configurations. For `h=(t-1)^2(t-alpha)^2`,
`alpha=(3+sqrt(5))/2`, the coefficient is `483+216sqrt(5)`; its conjugate
is also nonzero. For `h=(t-1)^2(t^2-3t+3)` it is 81.

## What the classical results actually supply

[Appelgate–Onishi, 1985](https://doi.org/10.1016/0022-4049(85)90099-4)
prove invertibility when either component degree has at most two prime
factors. Both 75 and 125 have three factors counted with multiplicity;
this theorem does not apply. Dividing the numerical degrees by 25 does
not produce a polynomial Keller pair of degrees 3 and 5. No such
degree-dividing operation is justified by a common leading carrier.

[Lee–Li, Magnus IV](https://arxiv.org/html/2408.01279v1), Lemma 4.12 and
Theorem 5.8/Corollary 5.10, apply after the exact-rectangle normalization:
coprime `2<=a<b`, `a|m,n`, monic top-right coefficient, matching rectangle
ratio and constant Jacobian. They supply a pre-generator and the narrow
inner-vertex region. The workspace's `(a,b,m,n)=(3,5,15,60)` satisfies the
arithmetic hypotheses. Their special result in section 8 assumes `a=n`,
which fails here (`3!=60`). Conjectures C/D/E are not theorem inputs.
A genuine nonlinear composition `P=alpha(W)` would contradict the Keller
condition: `1=alpha'(W)[W,Q]` forces `alpha'` constant. Finite agreement
with a cube at an edge does not establish that composition.

Nagata's 1988/1989 titles are confirmed in primary bibliographies, but I
did not retrieve and inspect his original full theorem text in this audit.
Accordingly no new exclusion or source bridge is being credited to him.
The gcd is 25, outside the familiar prime or small-gcd shortcut in any case.

## A more targeted primary-source lead

[Makar-Limanov–Trakhtenberg, Properties of a Jacobian mate,
MPIM 2024-33](https://archive.mpim-bonn.mpg.de/5148/1/mpim-preprint_2024-33.pdf),
printed page 18, lists four necessary resolution possibilities for `D=75`.
Exactly one displayed row has `lambda0=5/3`:

```
v0=(15,60), v1'= (3,0), v1=3(7/5,2),
v2'=(4/5,1),
phi0=x(x y^5-r1)^2 q0(x y^5), deg q0=2.
```

The F2 situation therefore survives this independent primary enumeration.
The second, principal endpoint has now been transported into the ramified F2
chart. Put `t=xy^5`, write `P15=a t^3 h^3`, and use the row's
`phi0 proportional to xh`, so `h=(t-r1)^2 q0(t)`; scaling sets `r1=1`.
In the first-root chart `y=x^(-1/5)+Y`, `X=x^(1/5)`, `u=XY`, `w=X^(-1)`,
a monomial `x^mu Y^nu` maps to `X^(5mu-nu)u^nu`. Factoring `P=X^15 F`
and writing `5mu-nu=15-s`, `j=nu`, the endpoint line
`25mu-17nu=3` becomes `5s+12j=72`; its supporting side gives
`5s+12j>=72`.

With the permitted earlier type-III shears retained, these are conditions
on centered coefficients, not on the raw `U_s` separately:
`T_s=[w^s]F(u+aw+bw^2,w)`. The coefficients below the half-space must
vanish, while the boundary endpoint requires `[u]T_12 != 0`. This is the
same terminal half-space and endpoint already encoded by the GGV terminal
normal, so the Makar-Limanov--Trakhtenberg row supplies no new independent
equation or exclusion. It remains conditional on the `lambda0=5/3` row being
the relevant `D=75` resolution; the primary table lists other rows too.

## Proof requirements still outstanding

The finite outer envelope is theorem-facing only subject to its cited
GGV/Makar-Limanov inputs and the written minimality-free standardization
extraction. Its custody/arithmetic checkers do not prove those sources.
The surviving quartic edge cases are `2+2` and `2+1+1`, with their exact
algebraic root configurations recorded in the main working note.
Excluding either still requires an identity contradicting every permissible
point of its *necessary* coefficient system, including vanishing slices,
resonant scalars and any coordinate normalizations used. A result inside
the narrower normalized carrier system needs its separate source/support
attachment. None of these requirements is fulfilled by a tuned local cube
or by the gauge result proved here.
