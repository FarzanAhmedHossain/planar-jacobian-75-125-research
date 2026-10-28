# Classical degree-criterion audit for the F2 pair `(75,125)`

This audit checks whether the classical Magnus, Applegate–Onishi, and Nagata
degree criteria directly remove the F2 candidate. They do not; the arithmetic
explains why the candidate is outside those criteria.

## Arithmetic of the pair

\[
75=3\cdot5\cdot5,\qquad 125=5\cdot5\cdot5,\qquad
\gcd(75,125)=25.
\]
The prime-factor counts here are with multiplicity, as in the component-degree
criterion.

| Criterion | Check at `(75,125)` | Outcome |
| --- | --- | --- |
| Magnus: coprime component degrees (and prime-degree corollaries) | gcd is 25; both degrees are composite | Not applicable |
| Applegate–Onishi degree criterion: one component degree has at most two prime factors | Both component degrees have three prime factors with multiplicity | Not applicable |
| Applegate–Onishi–Nagata gcd criterion: gcd at most 8 or prime | gcd is 25 | Not applicable |
| Abhyankar divisibility criterion | Neither 75 divides 125 nor 125 divides 75 | Not applicable |
| GGV exclusion of gcd `2p` and separate analysis of gcd 16 | gcd is 25 | Not applicable |

The original Applegate–Onishi article states the component-degree criterion;
the corrected/common-power structure remains valuable for organizing the
edge calculation, but it does not by itself settle this pair. The GGV 2017
algorithm places the work at the level of candidate Newton-polygon corners
and families. Thus the F2 edge-band equations and their exact elimination are
still required here.

## Relevance to the active local calculation

The common-power conclusion for the leading forms is already built into the
F2 setup: the relevant edge forms are proportional to powers with coprime
exponents 3 and 5. The quartic edge factor then splits into the `2+2` and
`2+1+1` root patterns. No classical degree test above removes either factor
pattern at gcd 25.

This is a scope audit, not a proof that the pair survives the full Jacobian
equations. The exclusion coverage for both F2 branches is unchanged.

## Ordinary leading-form factor split

For a hypothetical pair of degrees `(75,125)`, the top forms have the shape
`P_75 = a H^3`, `Q_125 = b H^5`, with `deg(H)=25`. Factor `H` over `C` and
let `r` be its number of distinct linear factors. The case `r=1` means that
`P` has one point at infinity. Shastri's 1990 Proposition 2 prints a
pointwise equivalence between this condition for a Jacobian pair and
automorphy, citing Abhyankar's Theorem 19.4. But Theorem 19.4 is a universally
quantified equivalence: the implication from *every* Jacobian pair having one
point at infinity to the Jacobian conjecture is equivalent to the
conjecture. It does not provide a pointwise implication for one pair. The
distinction matters because a high-degree source shear `x -> x+y^k` can make
the ordinary top form of any polynomial a pure power while preserving the
constant-Jacobian condition when applied to both components. The pointwise
reading would therefore convert the unresolved conjecture into a direct
exclusion; this audit does not rely on it. Jung--van der Kulk degree
divisibility applies only after automorphy has been proved.

Abhyankar's two-point theorem gives `r<=2`, so `r>=3` is excluded. Both the
one-factor and exactly-two-factor ordinary top-form strata remain unless
separately settled. The GGV F2 labels `2+2` and `2+1+1` refer to roots of a
separate quartic edge polynomial, so both F2 subbranches still require their
own lower-band exclusion. Current whole-branch coverage remains 0% for each.

Sources: [Shastri (1990), Proposition 2](https://ir.library.osaka-u.ac.jp/repo/ouka/all/6794/ojm27_01_11.pdf);
[Abhyankar, *Lectures on Expansion Techniques in Algebraic Geometry*,
Corollary 18.15 / Theorem 19.4](https://mathweb.tifr.res.in/Documents/Publications/Lectures/tifr57.pdf);
[Jung--van der Kulk degree divisibility](https://math.ac.vn/uploads/files/0302209.pdf).

## Sources

- Applegate and Onishi, [*The Jacobian conjecture in two variables* (1985)](https://doi.org/10.1016/0022-4049(85)90099-4).
- Guccione, Guccione, Horruitiner, and Valqui, [*Some algorithms related to the Jacobian Conjecture* (2017)](https://arxiv.org/abs/1708.07936).
- Guccione, Guccione, Horruitiner, and Valqui, [*On the shape of possible counterexamples to the Jacobian Conjecture*](https://arxiv.org/abs/1401.1784).
- The Magnus/Applegate–Onishi–Nagata gcd statement is summarized with its exact threshold in [this mathematical source](https://www.impan.pl/shop/publication/transaction/download/product/88176).

## Characteristic-pair sieve at gcd 25

There is one additional classical reduction that is sharper than the direct
degree criteria above. In Abhyankar's characteristic-pair setup, after the
usual generic linear-coordinate normalization to make the pair monic in one
variable, the first gcd entries are
`d_1 = 75` and `d_2 = gcd(75,125) = 25`. Each new characteristic exponent
strictly lowers this gcd, and the only divisor chain from 25 to 1 is
`25 -> 1` or `25 -> 5 -> 1`. Thus the characteristic-pair count is at most
three.

If `h = 2`, the two-characteristic-pair theorem implies the pair is
automorphic, which contradicts the nondivisibility of 75 and 125. If `h = 3`,
then necessarily `d_3 = 5` and `d_4 = 1`: `d_3 = 1` would have terminated the
sequence at `h = 2`. Abhyankar's reported “plus epsilon” result covers
`h = 3, d_3 <= 4`, so the only characteristic-pair signature not covered by
these results is precisely `h = 3, d_3 = 5`.

In the residual signature the characteristic sequence begins
`m_0 = -75, m_1 = -125, m_2, m_3`, with gcds
`d_1 = 75, d_2 = 25, d_3 = 5, d_4 = 1`. Consequently
`m_2 = 5u` with `5 ∤ u`, while `5 ∤ m_3`; the two later characteristic pairs
in Abhyankar's convention are `(u,5)` and `(m_3,5)`. This is an explicit
constraint to target in an approximate-root calculation. The `5` in `d_3=5`
is a gcd invariant; no equivalence with the factor multiplicity `5+20` of the
standardized F2 carrier has been established.

The same Oberwolfach report phrases the stronger claim for even `d_3` as a
result it proposes to prove, so this audit does **not** treat that statement as
an established theorem. The residual `h = 3, d_3 = 5` signature is a useful
target for extracting approximate-root/semigroup constraints and linking them
to the two-factor `5+20` F2 carrier. It does not by itself exclude that carrier
or change either quartic-edge branch's whole-branch coverage.

Sources: Abhyankar, [“Recent Progress on the Jacobian Conjecture” (Oberwolfach report)](https://ems.press/content/serial-article-files/45993), especially the introduction's definition of the gcd sequence and discussion of the two-characteristic-pair and plus-epsilon cases; Abhyankar, [“Some Thoughts on the Jacobian Conjecture, Part I”](https://doi.org/10.1016/j.jalgebra.2006.03.051).

## Appelgate–Onishi polygon lemmas in reduced ratio 3:5

For a hypothetical pair of degrees `(75,125)`, write `75=25·3` and `125=25·5`. The reduced exponents `m=3`, `n=5` are coprime, so the Appelgate–Onishi directional common-power lemmas apply: in every direction, the leading forms are proportional to `h^3` and `h^5`; their Newton polygons are homothetic, `W_P=3W` and `W_Q=5W`, for an integral convex polygon `W`. Nowicki–Nakai give explicit proofs of these lemmas and note that the original proof of the polygon statement needed clarification.

This is a genuine structural restriction and is the source of the common-power band equations used in the F2 model. It does not exclude either quartic root pattern: the `2+2` and `2+1+1` labels describe multiplicities of roots in one transformed edge factor, while the homothety statement constrains the two full Newton polygons. The ordinary top-degree two-point theorem is also not a substitute here, because this F2 quartic edge is not the ordinary total-degree leading form.

Source: Nowicki–Nakai, [“On Appelgate–Onishi’s Lemmas” (1988)](https://www.researchgate.net/publication/265368034_On_Appelgate-Onishi%27s_Lemmas), which states Lemmas A and B and proves the lattice-polygon homothety; Appelgate–Onishi, [“The Jacobian conjecture in two variables” (1985)](https://doi.org/10.1016/0022-4049(85)90099-4).

## Moskowicz variable-degree criteria at the standardized F2 corner

The strengthened variable-degree criteria do not close the surviving F2 carrier. Its support bounds give the actual partial degrees and leading-coefficient degrees

\[
n=\deg_y P=60,\quad u=\deg_x(a_n)=15,\qquad
r=\deg_y Q=100,\quad v=\deg_x(c_r)=25.
\]

Thus Moskowicz's notation has

\[
A=\gcd(n,u)=15,\quad C=\gcd(r,v)=25,\quad \gcd(A,C)=5,
\]

and normalized pairs

\[
(\tilde n,\tilde u)=(4,1)=(\tilde r,\tilde v).
\]

Theorem 2.4's unequal-normalized-pair hypothesis therefore fails. The exceptional-case Theorem 2.6 also does not apply: neither `A` nor `C` is in `{1,4} ∪ primes`, and `gcd(A,C)=5` is not 1 or 2. Theorem 2.7 fails for the same prime-factor reason. Proposition 4.1 only forces the ordinary two-factor multiplicity gcd to exceed 2; it permits the remaining gcd 5. These exact arithmetic checks explain why this Magnus/Appelgate–Onishi extension constrains, but does not exclude, either F2 quartic branch.

For the ordinary common form, write the exactly-two-factor case as
`H=L^a M^(25-a)` with `1 <= a <= 24` and nonproportional linear forms
`L,M`. The multiplicity condition `gcd(a,25-a)>2` then forces
`gcd(a,25)=5`, so, up to swapping the factors, the only remaining
multiplicity patterns are `(5,20)` and `(10,15)`. This is a finite restriction
on the ordinary top-form stratum only; it is not an exclusion of either pattern
and it does not identify them with the F2 quartic-edge `2+2`/`2+1+1` branches.

Source: Moskowicz, [*A Variation on Magnus' Theorem and its Generalizations*](https://arxiv.org/html/1810.08202v2), Theorems 2.4, 2.6, 2.7 and Proposition 4.1. The F2 support box and its nonzero upper-right corners supply the displayed actual partial degrees; no leading-form-only degree substitution is used here.

## McKay–Wang unequal-exponent condition on the two-factor top form

McKay–Wang's Newton-edge theorem requires the distinct factor exponents on an eligible negative-slope edge to be unequal. In the remaining ordinary two-factor forms, `H=L^a M^(25-a)` with `a=5` or `10` up to swapping, the exponents in `P_75` are `(3a, 3(25-a))`, namely `(15,60)` or `(30,45)`, and those in `Q_125` are `(5a, 5(25-a))`, namely `(25,100)` or `(50,75)`. Both are unequal. Thus this necessary edge condition is met in each of the two permitted multiplicity patterns and gives no exclusion; it does not replace the exact F2 lower-coefficient calculation.

Source: McKay and Wang, [“A Note on the Jacobian Condition and Two Points at Infinity”](https://doi.org/10.2307/2047856), Theorem 1.

### Applicability check for the F2 quartic edge

The tempting inference “`2+1+1` has three distinct roots, so McKay–Wang's two-factor bound excludes it” does not follow from Theorem 1. That theorem assumes an edge of **negative slope** in the Newton polygon of a polynomial pair in ordinary polynomial coordinates. The standardized F2 edge used here has endpoints `(3,0)` and `(15,60)`, hence slope `+5` in the `(x,y)` exponent diagram. Its quartic `h(t)` is written after the Laurent/monomial chart `t=xy^5`, `z=y^{-1}`; it is not the ordinary negative-slope face polynomial required by the theorem. Transporting the theorem through this chart would need a separate result. Therefore the root count of `h` alone gives no McKay–Wang exclusion of either F2 quartic branch.

## Moskowicz's stated limit at `(75,125)`

Moskowicz's Proposition 4.1 proves only the normalized multiplicity-gcd bound `d>2`. The paper's application section then explicitly says that this proposition does not discard the exceptional degree pairs `(75,125)` and `(64,224)`. This confirms that the Magnus-style multiplicity sieve is a reduction to the hard cases, not an exclusion of this one.

Source: Moskowicz, [*A Variation on Magnus' Theorem and its Generalizations*](https://arxiv.org/html/1810.08202v2), Proposition 4.1 and the discussion immediately following it.

## Direct Magnus, Nagata, and Applegate–Onishi degree tests

Nagata's revised two-variable account states Magnus's coprime-degree criterion, and gives the Applegate–Onishi cases as: prime `gcd(deg f,deg g)` forces one degree to divide the other; or one component degree is a product of two primes. Here `gcd(75,125)=25`, and both component degrees have three prime factors counted with multiplicity (`75=3·5·5`, `125=5·5·5`). These degree-only cases do not apply. The same paper notes the correction/completion of the Applegate–Onishi argument; invoking those results cannot discharge this pair without additional hypotheses beyond its degree arithmetic.

Source: M. Nagata, [*Two-dimensional Jacobian Conjecture* (revised account)](https://repository.kulib.kyoto-u.ac.jp/server/api/core/bitstreams/9ef8e868-5526-4830-b19f-543c0af09e7c/content), Magnus theorem and Applegate–Onishi theorem, PDF pp. 5–6.
