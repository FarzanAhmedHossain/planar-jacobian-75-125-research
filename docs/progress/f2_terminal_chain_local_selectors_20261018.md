# Transporting the F2 terminal corner into local band coefficients

Date: 2026-10-18. This establishes additional necessary local equations
at the distinguished F2 double root, conditional on the same front-end
standardization and F2 complete-chain attachment used by the main model.
It excludes the tuned nonzero quadratic-cube local face there. It does
not exclude the complete two-factor branch or (75,125).

## Primary-source inputs, with scope

Work over an algebraically closed characteristic-zero field with a standard
`(3,5)`-pair in the F2 row. In
[GGHV 2017](https://arxiv.org/pdf/1708.07936), Theorem 2.20, the table on
printed page 25, Theorem 3.1, and equations (3.16)–(3.17), the row has
`A1=(7/5,2)`, `l1=5`, `k=1`, and terminal type I. Its proof on printed
page 14 removes intervening type-III corners while leaving `l1` fixed.
For type I.a, Theorem 3.1 requires `rho|l1`; for type I.b equation (3.16)
specifies the two starting points. These statements concern actual pairs
and have no minimality hypothesis in these sections.

[GGV 2016](https://arxiv.org/pdf/1401.1784), Propositions 5.17–5.18 on
printed pages 27–28, specify the actual automorphisms: `x` is fixed and
`y -> y+lambda*x^(sigma/rho)`. In type III the edge polynomial is a
single power of a linear factor and `rho|l`; type II allows multiple
factors and changes the complete chain. Thus the discarded corners are
literal source shears, not an unspecified replacement of the pair.

## The terminal normal is (25,-17)

For F2, equation (3.17) reads

```
(3+5)*2 - 5*(5*2-7) = 1.
```

The swapped ordered assignment fails this equality. Hence the terminal
type-I.b starts are

```
st(P)=(4/5,1), st(Q)=(1/5,0),
en(P)=3*A1=(21/5,6), en(Q)=5*A1=(7,10).
```

The line through the P endpoints has primitive normal `(25,-17)` and
P-weight 3; the Q-weight is 5. The sum 8 equals `25-17`, as required
for a nonzero constant edge Jacobian. Type I.a cannot occur at this
normal because `25` does not divide `l1=5`.

Before that terminal edge, the only possible removable type-III slopes
between `1/5` and `17/25` have `rho|5`; they are `2/5` and `3/5`.
Consequently the entire relevant source map is

```
x -> x,
y -> y + x^(-1/5) + a*x^(-2/5) + b*x^(-3/5),
```

after normalizing the first chosen ramified root to 1. Either `a` or `b`
may vanish when the corresponding corner was absent. All terms are in
`K[x^(+/-1/5),y]`. A slope `4/5` lies after the terminal slope and cannot
occur in this preceding chain. This retains all permitted collapsed
edges; simply assuming that the first shear alone suffices would be wrong.

## A finite local coefficient system

Set `x=X^5`, `w=X^-1`, and use the first ramified chart
`y=X^-1*(1+u)`. For every original monomial,

```
x^i y^j = X^(5i-j) (1+u)^j.
```

Write

```
P=X^15 F(u,w), F=sum_s A_s(u) w^s,
Q=X^25 G(u,w), G=sum_s B_s(u) w^s.
```

The F2 halfspace ensures `s>=0`; these are finite polynomials in `u,w`.
In terms of the old band coefficient `p_(15-s)(t)`,

```
A_s(u)=(1+u)^(-(15-s))*p_(15-s)((1+u)^5).
```

This is a polynomial despite its displayed denominator, since each
original monomial has nonnegative y exponent. In particular

```
A0=h((1+u)^5)^3;
A_s=(1+u)^s*U_s((1+u)^5), 1<=s<=4;
A5=U5((1+u)^5);
A6=(1+u)*U6((1+u)^5).
```

The remaining source shears are exactly `u -> u+a*w+b*w^2`.
In the resulting chart a monomial `X^k u^j` has `(25,-17)`-weight
`5k-12j`. Therefore the full terminal support conditions are

```
[w^s u^j] F(u+a*w+b*w^2,w)=0 whenever 5s+12j<72,
[w^s u^j] G(u+a*w+b*w^2,w)=0 whenever 5s+12j<120.
```

There are 53 candidate P coefficient positions and 136 Q positions.
Some equations vanish already from the top edge and prior band rows.
The two scalar variables `a,b` must be retained or eliminated exactly.
The P boundary positions are `(s,j)=(0,6),(12,1)`; Q has
`(0,10),(12,5),(24,0)`. The terminal starting coefficients at P `(12,1)`
and Q `(24,0)` are nonzero. These conditions impose the full terminal
halfspace, including information through original band distance 14 for P
and 23 for Q. They are stronger than the initial double-root statement.

## An immediate local exclusion and its zero branch

The P halfspace contains no monomial of ordinary `(u,w)` total degree
less than 6, and its only degree-six monomial is `u^6`. The recentering
has invertible linear part `(u,w)->(u+a*w,w)`, so the uncentered F must
have ordinary multiplicity 6 and tangent cone a sixth power.

As `t-1=(1+u)^5-1=5u+O(u^2)`, this gives

```
ord_1 U1>=5, ord_1 U2>=4, ord_1 U3>=3,
ord_1 U4>=2, ord_1 U5>=1.
```

In particular, writing `U1=h^2 S`, the previously surviving tuned case
`S(1)!=0` is impossible at the distinguished root. Its weighted local
face is proportional to `(u^2+k*w)^3` with `k!=0`. Such a polynomial has
two distinct roots of multiplicity three in the Puiseux chart; it is
type II, rather than a removable type-III sixth power. More directly,
after either permitted integer shift its coefficient of `w^3 u^0` is
still `k^3!=0`, contrary to the terminal halfspace. The `S(1)=0` branch
is retained.

There are also exact tangent-cone equations. Let `c=[(t-1)^2]h!=0`,
`C=c^3`, `A=[(t-1)^5]U1`. Once the lower-order vanishing above holds,
the leading diagonal coefficients of U2,...,U6 must be

```
[(t-1)^4]U2 = 5 A^2/(12 C),
[(t-1)^3]U3 = 5 A^3/(54 C^2),
[(t-1)^2]U4 = 5 A^4/(432 C^3),
[(t-1)^1]U5 = A^5/(1296 C^4),
U6(1)         = A^6/(46656 C^5).
```

Equivalently, the ordinary sextic is `C*(v+A*w/(6*C))^6` in the
`v=t-1` leading coordinate. Factors of five cancel in these displayed
relations when transporting to the ramified u coordinate. Division here
is only by the proved nonzero `C` and characteristic-zero constants.
No assumption that `A!=0` is made. On `U1=0`, all displayed diagonal
coefficients vanish, strengthening the local orders to

```
ord_1 U2>=5, ord_1 U3>=4, ord_1 U4>=3,
ord_1 U5>=2, ord_1 U6>=1.
```

These conditions are at the root selected by the F2 complete chain.
They have not been proved at every double root simultaneously. In `2+2`
the other double root still needs separate treatment; in `2+1+1` the
simple roots still need separate treatment.

## Boundary consistency and next calculation

The terminal boundary itself is not contradictory. With
`z=x^(17/25)*y`, its two edge polynomials can be written

```
Pedge=x^(3/25)*(A0*z^6+B0*z),
Qedge=x^(5/25)*(C0*z^10+D0*z^5+E0).
```

The edge Jacobian is `(3 f g'-5 f' g)/25`, so for Jacobian one it requires

```
25 B0 C0-15 A0 D0=0,
10 B0 D0-30 A0 E0=0,
-B0 E0/5=1.
```

For nonzero `A0,B0`, the assignments
`C0=-9 A0^2/B0^3`, `D0=-15 A0/B0^2`, `E0=-5/B0` satisfy them exactly.
For a general nonzero Jacobian scalar, replace the rightmost 1 by that
scalar. Thus the useful next step is imposing all 189 coefficient
positions jointly with the existing algebraic edge family, retaining
`a,b`, and testing whether the new root contact constraints collapse
its surviving U1/U2 parameter spaces. No terminal-edge contradiction
alone is claimed.

The accompanying
`f2_principal_edge_transport_audit_20261018.py` passed on 2026-10-18. It
checks the terminal endpoint arithmetic, all removable slopes, monomial
weight conversion, forbidden positions, tangent-cone identities including
the zero-A branch, the tuned-face obstruction, and the consistent terminal
boundary witness. Its source dependency remains the cited actual-pair
chain theorem and the current front-end F2 attachment.
