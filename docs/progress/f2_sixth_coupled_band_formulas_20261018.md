# Exact sixth coupled edge block and Q19 polynomiality

Date: 2026-10-18. These are necessary identities in the existing normalized
F2 edge model. They do not exclude a carrier or degree pair.

## Full support and notation

Put `M=t h`, `P15=M^3`, `Q25=M^5`. For `s=1,...,6`, denote the actual
band `P_(15-s)` by `A_s`, and write `q24,...,q19` for actual Q bands.
Thus

```
A1=t^3 U1=t^3 h^2 S,
A2=t^3 U2=t^3 h B,
A3=t^3 U3, A4=t^3 U4,
A5=t^2 U5, A6=P9=t^2 U6.
```

The exact supports are

```
P9 in span(t^2,...,t^13): deg U6<=11;
Q19=t^4 V6, deg V6<=19;
W33 in span(t^7,...,t^66).
```

The degree bound 12 on U6 would be a support relaxation: it introduces
`x^14 y^61` outside P's rectangle. If used, the Q19 law and its degree
cap force the extra leading coefficient zero. There are exactly 12 new
P9 coefficients and 60 new W33 coefficients.

## Denominator-cleared Q19 equation

Use the full fifth slice, including its possible constant U5 coefficient:

```
F=h^3+h^2 S w+h B w^2+U3 w^3+U4 w^4+(U5/t)w^5+(U6/t)w^6.
```

The first resonant scalar is kappa, normalized so that the actual Q20
adds `kappa M^4`. Through distance six the formal solution is the
binomial `P^(5/3)+kappa P^(4/3)`; no new resonance occurs. Define

```
N=729 U3^2-486 S B U3-81 B^3+108 S^3 U3
  +162 S^2 B^2-63 S^4 B+7 S^6,
D=3B-S^2, Delta=27U3+2S^3-9SB.
```

Then `N=Delta^2-3D^3`, and exact coefficient extraction gives

```
V6=(5/3)h^2 U6+(10/9)h S U5
   +t((10/9)B-(5/27)S^2)U4
   +(5/6561)t N/h+(4/3)kappa h^3 S.
```

Equivalently, the full polynomial equation linear in the two new bands is

```
6561 h V6-10935 h^3 U6
 =7290 h^2 S U5+t h(7290B-1215S^2)U4
  +5t N+8748 kappa h^4 S.
```

In actual physical-band variables this is

```
6561 h q19-10935 t^2 h^3 A6
 =7290 t^2 h^2 S A5+t^2 h(7290B-1215S^2)A4
  +5t^5 N+8748 kappa t^4 h^4 S.
```

These forms divide by no possibly zero coefficient, slice, S, B, kappa,
or parameter. All rational constants are valid in characteristic zero.

## Equivalence and efficient sparse implementation

Since `h(0)!=0`, t and h are coprime. Every term except `t N/h` in the
V6 expression is polynomial. Therefore polynomiality is exactly `h|N`,
equivalently `h|W6` in the main note's scaling `W6=5N/6561`.
With the exact older degree caps, `deg N<=22`, and when h divides N all
terms in V6 have degree at most 19. Thus that upper cap adds no further
restriction. The apparent t inverse in the normalized fifth and sixth
slices is cleared by the factor t defining V6. At t=0 the law forces

```
V6(0)=(5/3)h(0)^2 U6(0)+(10/9)h(0)S(0)U5(0)
      +(4/3)kappa h(0)^3 S(0),
```

which is a coefficient value, rather than another solvability obstruction.
In the x-centered gauge U5(0)=0, its middle term vanishes legitimately.

For uniqueness, the homogeneous physical sixth Q equation is
`5M q'-19M'q=0`; a nonzero solution would be `q=c M^(19/5)`. M has a
simple zero at t=0, so this has nonintegral valuation and cannot be a
rational function. Thus no homogeneous mode or zero branch was lost.

For the sparse solver, divide N by the fixed quartic h:

```
N=h J+R, deg_t R<4.
```

Retain every coefficient of R=0 as an equation in the OLD parameters.
Then define

```
q19=(5/3)M^2 A6+(10/9)t^2 h S A5
    +t^2((10/9)B-(5/27)S^2)A4
    +(5/6561)t^5 J+(4/3)kappa t^4 h^3 S.
```

Together with R=0 this is equivalent to the denominator-cleared equation.
It leaves only the 72 new P9/W33 variables and gives a constant matrix
for them. Dropping the forcing-only R equations would break equivalence.
Polynomial division by h uses its fixed nonzero leading coefficient in
the coefficient field, with no division by an old parameter.

## Coupled equations at weights 44 and 54

Let `g=W39`, and let W38,...,W34 be the previously solved bands. The
corrected power relation has top correction `-3 kappa P^3 Q`, giving

```
[P,W*]=3Q^2-3 kappa P^3+(lower weight terms),
[Q,W*]=5P^4+9 kappa P^2 Q+(lower weight terms).
```

Other correction products have brackets of weight below 44 and 54,
respectively, so they do not enter this block. With primes denoting
t derivatives, the new P equation is

```
0=15M^3 W33'-99M^2 M' W33
 +14A1 W34'-34A1' W34
 +13A2 W35'-35A2' W35
 +12A3 W36'-36A3' W36
 +11A4 W37'-37A4' W37
 +10A5 W38'-38A5' W38
 + 9A6 g'-39A6' g
 -6M^5 q19-6q24 q20-6q23 q21-3q22^2
 +9 kappa M^6 A1.
```

Set

```
C6=4M^9 A6+12M^6(A1 A5+A2 A4)+6M^6 A3^2
  +12M^3 A1^2 A4+24M^3 A1 A2 A3+4M^3 A2^3
  +4A1^3 A3+6A1^2 A2^2.
```

This is exactly the coefficient at distance six in P^4. The new Q
equation is

```
0=25M^5 W33'-165M^4 M' W33
 +24q24 W34'-34q24' W34
 +23q23 W35'-35q23' W35
 +22q22 W36'-36q22' W36
 +21q21 W37'-37q21' W37
 +20q20 W38'-38q20' W38
 +19q19 g'-39q19' g
 -5C6-9 kappa(2M^8 A1+M^6 q24).
```

The kappa terms are the distance-one coefficients of the top correction
products. Omitting them, or replacing them by their top coefficients
from the previous block, changes the equations.

`f2_sixth_coupled_band_formula_audit_20261018.py` passed on 2026-10-18.
It independently checks the Q19 binomial/resonance expansion, the
compressed pole numerator, normalized/physical linear laws, all support
intervals and degree caps, and both coupled equations against direct
bracket and power-series convolution. It performs no rank solve and
claims no branch exclusion.
