# Branch-coverage ledger for the `(75,125)` problem

This ledger separates global top-form strata from narrower normalized F2 coefficient models. A local certificate can close only the exact domain it covers.

| Layer | Cases | Current status in this archive | What remains for an overall exclusion |
|---|---|---|---|
| Degree-25 top form `H` | One distinct factor | Open in general | Exclude every admissible pure power and all compatible lower terms, or invoke a theorem whose hypotheses truly cover this pair |
| Degree-25 top form `H` | Exactly two distinct factors | Open in general | Cover all multiplicity splits `e+(25-e)`, all coefficients and lower terms, and prove the reduction to each tested model |
| Degree-25 top form `H` | Three or more distinct factors | Excluded by the at-most-two-points-at-infinity theorem | No additional search needed for this stratum |
| Fixed form `H=x^25+y^25` | 25 distinct factors | Excluded as a possible leading form | This removes this form only; the one- and two-factor strata are different forms |
| Normalized F2 quartic edge | `2+2` | Selector leaves an algebraic candidate family | Continue higher coefficient layers and prove any model-to-Keller-pair bridge |
| Normalized F2 quartic edge | `2+1+1` | Selector leaves an algebraic candidate family | Continue higher coefficient layers and prove any model-to-Keller-pair bridge |
| Remaining F2 coefficient space | Other classified branches and lower coefficients | Partially constrained by dated reports | Audit all branches, parameters, and lower terms; no whole surviving branch is declared eliminated here |

## Coverage standard

For each claimed branch exclusion, record (1) the exact branch definition, (2) the coefficient field and normalizations, (3) proof that every actual hypothetical pair in that branch maps into the model, (4) exact equations and parameter domains, (5) a checkable contradiction or a certified finite decomposition, and (6) proof that no exceptional or degenerate case was dropped. A Groebner-basis computation is useful only when the ideal, saturation, coefficient field, and certificate establish emptiness for the full branch.
