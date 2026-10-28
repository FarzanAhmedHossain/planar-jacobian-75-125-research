# Scope audit: one point at infinity and the `(75,125)` pair

Date: 2026-10-28.

## Finding

Do not count the ordinary one-factor case `H=L^25` as excluded by reading
Shastri's 1990 Proposition 2 pointwise. Shastri's text does print the
equivalence for a particular complex Jacobian pair, but it cites Abhyankar's
Theorem 19.4. In Abhyankar's source, Theorem 19.4 equates *universally
quantified assertions*: for every Jacobian pair, the pair is automorphic;
for every Jacobian pair, the first component has one point at infinity; and
related global conditions. The global equivalence does not establish the
pointwise implication for one selected pair.

## Shear check

This scope distinction is forced by a simple source automorphism. For any
nonconstant polynomial
`F(x,y)=sum c[i,j] x^i y^j`, choose `k>1` sufficiently large and let
`i0` be the largest `x` exponent and `j0` the largest `y` exponent among
terms with that `x` exponent. Then the unique highest-total-degree term of
`F(x+y^k,y)` is `c[i0,j0] y^(k*i0+j0)`. Thus the ordinary degree form of the
transformed polynomial is a pure power of one linear form. A single `k` can
be chosen for both components of any fixed pair. The source shear has
Jacobian one, so applying it to both entries preserves a constant Jacobian.

Consequently, if “one point at infinity” for a component of a Jacobian pair
were a pointwise automorphy criterion in arbitrary coordinates, every
hypothetical planar Keller pair could first be sheared into that condition,
which would prove the still-open planar Jacobian conjecture. The cited global
equivalence cannot be used in that way. The ordinary factor-count argument
therefore closes `r>=3` by Abhyankar's two-point theorem but does not close
`r=1` from this citation. This note does not assert a different resolution
of the one-factor case.

## Consequence for active `(75,125)` work

For a hypothetical pair of exact degrees `(75,125)`, the top forms are
`P_75=aH^3`, `Q_125=bH^5`, `deg(H)=25`. Abhyankar's two-point theorem
excludes every `H` with at least three distinct linear factors. The `r=1`
and `r=2` ordinary factor strata require valid arguments. This ordinary
factor split is distinct from the later F2 quartic root labels `2+1+1` and
`2+2`; the conditional GGV reduction to the F2 carrier must be audited on
its own hypotheses and remains the active coefficient-level route.

This source-scope correction changes no computational coverage:

- F2 `2+1+1, U1=0`: 0% branchwide exclusion.
- F2 `2+1+1, U1!=0`: 0% branchwide exclusion.
- F2 `2+2`: 0% branchwide exclusion; it remains paused.

The current finite-profile figures are separate screens, not branchwide
fractions: 9,334/12,128 (76.962401%) pruned in the recorded `U1=0`
double-root grid, and 22,830/51,870 (44.013881%) pruned in the recorded
`U1!=0` finite grid. For each F2 quartic branch, 1 of 2 coarse contact labels
is removed (50% of that two-label set only).

## Primary sources

- Shastri, “A Discriminant Criterion for the Two-Dimensional Jacobian
  Problem,” Osaka J. Math. 27 (1990), Proposition 2 and its remark:
  <https://ir.library.osaka-u.ac.jp/repo/ouka/all/6794/ojm27_01_11.pdf>
- Abhyankar, “Lectures on Expansion Techniques in Algebraic Geometry,”
  Theorem 19.4 (universal equivalence) and the two-point theorem:
  <https://mathweb.tifr.res.in/Documents/Publications/Lectures/tifr57.pdf>
- Abhyankar, “Conics as Axes and Jacobian Problem,” pp. 172–173, explaining
  the distinction between the global one-point assertion and the two-point
  theorem for individual Jacobian pairs:
  <https://indianmathsoc.org/ms/mathstudent2007.pdf>
