# Leading-form factor audit for the degree pair `(75,125)`

## Exact common-power reduction

Work in characteristic zero and extend coefficients to an algebraic closure.
Let `A=P_75` and `B=Q_125` be the top ordinary homogeneous forms of a
hypothetical Keller pair `(P,Q)`. The degree-198 part of `[P,Q]` is `[A,B]`,
so `[A,B]=0`. Writing `A=y^75 a(x/y)` and `B=y^125 b(x/y)` gives

```text
[A,B] = y^198 (125 a' b - 75 a b')
       = 25 y^198 (5 a' b - 3 a b').
```

Thus `(a^5/b^3)'=0`, hence `A^5=c B^3` for some nonzero scalar `c`. Unique
factorization and `gcd(5,3)=1` give

```text
A = c_A H^3,     B = c_B H^5,     deg(H)=25.
```

The constants can be absorbed by independent nonzero target scalings. This is
a statement about top forms; lower-degree terms in `P` and `Q` remain
essential to the Keller equations.

## Exhaustive geometric factor strata

Over the algebraic closure, write

```text
H = c product_i L_i^(a_i),    sum_i a_i=25,
```

where the `L_i` are distinct linear forms. The number of distinct factors is
the number of points at infinity of the degree form.

- One factor: `P_75` has one point at infinity. Do not infer automorphy from
  Shastri's printed pointwise reading of Proposition 2. It cites Abhyankar's
  Theorem 19.4, whose hypothesis is universally quantified over every
  Jacobian pair; it does not establish automorphy for this one pair. A source
  shear `x -> x+y^k` can make any polynomial's ordinary top form a pure power,
  while preserving the Keller condition on the pair, so a pointwise theorem
  of that form would immediately settle the unresolved conjecture. The
  one-factor stratum remains open in this audit absent a separate argument.
- At least three factors: Abhyankar's two-point theorem says each component of
  a plane Jacobian pair has at most two points at infinity. This stratum is
  empty for a Keller pair.
- Exactly two factors: after an invertible linear source change,
  `H=x^a y^(25-a)` for some `1<=a<=24`. This is the only ordinary leading-form
  stratum left by these theorems.

In particular, `H=x^25+y^25` has 25 distinct linear factors over the
algebraic closure, so it is excluded by the two-point theorem. Over `Q` its
cyclotomic factorization has three irreducible factors; the geometric count
used here is over the algebraic closure.

## Relation to the GGHV/F2 branch

The factor count above is an ordinary top-form split. GGHV's 2017 enumeration
lists the degree-125 F2 member `(m,n)=(3,5)` with corner `A_0=(5,20)`; reaching
that normalized `5+20` carrier from the full two-factor stratum uses the
published Newton-chain reduction. It does not follow from the two-point
theorem alone. The F2 local `2+2` and `2+1+1` labels are finer residual
factorization cases; they should not be counted as the one/two/three-plus
factor strata of `H`.

Primary references:

- Abhyankar, “Conics as Axes and Jacobian Problem,” especially the statement
  that a Jacobian pair component has at most two points at infinity:
  <https://www.indianmathsoc.org/ms/mathstudent2007.pdf>
- Shastri, “A Discriminant Criterion for the Two-Dimensional Jacobian
  Problem,” Proposition 2 as printed; its pointwise scope is not accepted here
  because it cites Abhyankar's universally quantified Theorem 19.4:
  <https://ir.library.osaka-u.ac.jp/repo/ouka/all/6794/ojm27_01_11.pdf>
- Guccione–Guccione–Horruitiner–Valqui, “Some Algorithms Related to the
  Jacobian Conjecture,” §6 and the F2 `(3,5)` row:
  <https://arxiv.org/pdf/1708.07936>

## Scope

This closes only the three-or-more-factor ordinary leading-form stratum. The
one-factor and exactly-two-factor strata remain open here. It does not close
either F2 local branch. The current whole-branch exclusion percentages
therefore remain 0% for F2 `2+1+1` and 0% for paused F2 `2+2`.
