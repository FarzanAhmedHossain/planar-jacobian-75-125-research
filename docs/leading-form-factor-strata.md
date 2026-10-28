# Degree-25 leading-form factor strata for `(75,125)`

## Why a degree-25 form appears

Let `P,Q` be a hypothetical Keller pair over a characteristic-zero field, with total degrees 75 and 125. Their highest homogeneous components have zero bracket, because the nonzero constant `[P,Q]` has degree zero. After passing to an algebraic closure and choosing an affine coordinate on the projective line, this gives the relation

\[
3 f g' - 5 f' g=0,
\]

where `f` and `g` represent the degree forms. It implies `g^3 = c f^5` for a nonzero scalar `c`. Unique factorization then gives

\[
P_{75}=aH^3,\qquad Q_{125}=bH^5,
\qquad \deg H=25.
\]

This concerns the two leading homogeneous components. The literal pair `(H^3,H^5)` has Jacobian zero and is not itself a Keller pair; lower-degree terms of `P,Q` are essential.

## The exhaustive factor-count split

Over an algebraically closed field, the binary homogeneous form `H` splits into linear forms. Counting *distinct* linear factors gives exactly three mutually exclusive possibilities:

1. **One distinct factor:** `H=c L^25`.
2. **Exactly two distinct factors:** `H=c L^e M^(25-e)` for some `1 <= e <= 24` and nonproportional linear forms `L,M`.
3. **At least three distinct factors:** `H=c product_i L_i^{e_i}` with at least three pairwise nonproportional linear forms and positive exponents summing to 25.

The at-most-two-points-at-infinity theorem says that a plane Jacobian pair has degree forms with at most two distinct irreducible factors. Since the factors of `H` are exactly the distinct factors shared by `H^3` and `H^5`, possibility 3 is excluded for any hypothetical Keller pair. The theorem does not by itself decide possibilities 1 or 2.

## The specific choice `H=x^25+y^25`

In characteristic zero,

\[
x^{25}+y^{25}=\prod_{\zeta^{25}=-1}(x-\zeta y),
\]

and the 25 roots of `z^25=-1` are distinct. Thus this fixed `H` is a member of possibility 3, already ruled out by the two-points theorem. It cannot at the same time be a one-factor or two-factor example. Those are different choices of `H` and different strata of the search.

## When separate branch exclusions would be sufficient

Yes, the split can support an overall exclusion argument, with a crucial scope condition:

- the three cases above must be proved to cover every hypothetical pair after all allowed normalizations;
- the third case is already excluded by the cited theorem;
- the one-factor branch must be excluded for every admissible `H` and every compatible choice of lower-degree coefficients;
- the exactly-two-factor branch must likewise be excluded in full, including all multiplicity patterns `e=1,...,24`, coefficient parameters, and any distinct normalized subcases;
- any computer-assisted contradiction must be an exact proof for the entire branch (or a certified finite cover), not just a numerical sample or a bounded truncation;
- the reductions from an actual Keller pair to the equations tested in each branch must be proved and preserve all possible cases.

If those obligations are met, the exhaustive split implies that no pair of degrees `(75,125)` exists. Calling a calculation “local” does not change the quantifier: if it only excludes a neighborhood or a selected local model, then the rest of that branch remains open. In particular, showing that `x^25+y^25` is impossible does not eliminate all degree-25 forms.

## Do not conflate with the F2 quartic edge types

The degree-25 `H` factor count is a top-form classification. After a separate Newton-chain/F2 reduction, a quartic edge polynomial may have residual multiplicity patterns `2+2` or `2+1+1`. That quartic is not the original `H`; those are finer carrier cases inside the surviving reduction. The corrected-power selector gives necessary algebraic conditions for those cases but does not by itself exclude either entire branch or the degree pair. See [the F2 band model](F2-75-125-band-model.md) and [the edge-selector note](edge-selector.md).

## Sources

- Abhyankar, “Conics as Axes and Jacobian Problem,” *The Mathematics Student* 76 (2007), pp. 171–183. The exposition explains the correspondence between distinct degree-form factors and points at infinity, and the two-point bound. [PDF](https://www.indianmathsoc.org/ms/mathstudent2007.pdf)
- Abhyankar, *Expansion Techniques in Algebraic Geometry*, TIFR Lectures on Mathematics and Physics 57 (1977), cited for the two-points theorem and its broader consequences. [TIFR volume](https://mathweb.tifr.res.in/sites/default/files/publications/ln/tifr57.pdf)
