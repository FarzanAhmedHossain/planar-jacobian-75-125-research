# Exact elimination of the two terminal centering scalars

Date: 2026-10-18. This note checks the local algebra used to impose the
terminal F2 halfspace. The source attachment is proved separately in
`f2_terminal_chain_local_selectors_20261018.md`. No degree-pair exclusion
is claimed here.

Write `F(u,w)=sum_i A_i(u)w^i`, with `ord A0=6`, `C=[u^6]A0!=0`.
The first terminal selector `S(1)=0` gives `ord A1>=5`.

The coefficient of `u^5 w` after replacing `u` by `u+a w+b w^2` is
`[u^5]A1+6C*a`. Thus the unique choice is

```
a=-[u^5]A1/(6C).
```

Set `R2=A2+a A1'+a^2 A0''/2`. The full second row is `R2+b A0'`.
Because `A0'` has order five, its coefficients of `u^0,...,u^4` must
vanish independently of b. The coefficient of `u^5` then determines

```
b=-[u^5]R2/(6C)
 =-([u^5]A2+6a[u^6]A1+21a^2[u^7]A0)/(6C).
```

These formulas retain the `a=0` and `b=0` cases. The only denominator
requiring a geometric nonvanishing condition is the already nonzero C.

For every s, the exact centered row is

```
T_s=sum_(0<=i<=s) sum_(m,n>=0, m+2n=s-i)
      a^m b^n A_i^(m+n)/(m! n!).
```

This follows from the formal Taylor expansion, with derivative in u. The
equations imposed are `[u^j]T_s=0` when `5s+12j<72`. In particular the
last three rows through distance five are

```
T3=A3+a A2'+b A1'+a^2 A1''/2+ab A0''+a^3 A0'''/6;

T4=A4+a A3'+b A2'+a^2 A2''/2+ab A1''+a^3 A1'''/6
   +b^2 A0''/2+a^2 b A0'''/2+a^4 A0''''/24;

T5=A5+a A4'+b A3'+a^2 A3''/2+ab A2''+a^3 A2'''/6
   +b^2 A1''/2+a^2 b A1'''/2+a^4 A1''''/24
   +ab^2 A0'''/2+a^3 b A0''''/6+a^5 A0'''''/120.
```

The forbidden u powers are `0..5` in T1/T2, `0..4` in T3/T4, and
`0..3` in T5. Relations already used to determine a,b should be
substituted rather than counted as additional independent equations.

## Compatibility with the source x-centering gauge

Let `P_c(x,y)=P(x+c,y)`. In the original first ramified chart define
`v=(1+c w^5)^(1/5)` as its characteristic-zero formal binomial series,
with constant coefficient 1. Direct substitution gives the exact identity

```
F_c(u,w)=v^15 F(v(1+u)-1,w/v).
```

Modulo `w^6`, this is

```
F_c=F+c w^5*(3A0+(1+u)A0'/5).
```

Hence A0,...,A4 are unchanged, so a and b are unchanged. The change in
A5 has u order at least five; it cannot affect its forbidden powers
`u^0,...,u^3`. Consequently the source gauge setting `[t^2]P10=0`
commutes with all terminal equations through distance five in the precise
sense that those equations and their a,b values are invariant.

The compatibility extends to the whole terminal halfspace. Put
`phi(w)=a w+b w^2`, `Fbar(u,w)=F(u+phi(w),w)`. Then

```
Fbar_c(u,w)=v^15 Fbar(ut,wt),
wt=w/v,
ut=v(1+u+phi(w))-1-phi(w/v).
```

Here `wt=w+O(w^6)` and

```
ut=u+(c/5)(1+u)w^5+(2ac/5)w^6+(3bc/5)w^7+O(w^10).
```

For the positive local weights `weight(w)=5`, `weight(u)=12`, every
correction in wt has weight greater than 5, every correction in ut has
weight greater than 12, and `v^15` is a unit of constant coefficient 1.
Thus substitution preserves the ideals of series of weight at least 72
and at least 120. The same argument applies to Q, with multiplier `v^25`.
This is compatibility of the x-centering and terminal shear gauges. It
does not authorize simultaneously imposing a separate, independently
chosen upper-y carrier centering.

The accompanying checker `f2_terminal_centering_identities_20261018.py`
compares the Taylor row formula with direct substitution of generic
polynomials through every forbidden coefficient at s<=5, checks the
mixed factorial constants, a and b formulas, and x-gauge invariance
through that distance. Degree eight jets suffice universally for those
positions because `j+m+n<=j+s<=8`.

The checker was strengthened on 2026-10-18 to keep all A1 coefficients
of u degrees zero through eight as independent variables. It therefore
does not assume away the S(1)=0 selector while auditing the 32 positions.
The only input order restriction is A0 order six, already supplied by the
attached double root of h. The strengthened run passed all 32 positions
and reports explicitly that the low A1 coefficients were retained. This
checks coefficient transport and source-gauge invariance; a separate
solver receipt must still show that its resulting equations were imposed.
