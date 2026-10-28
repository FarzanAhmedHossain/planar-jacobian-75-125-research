# Research status

**Snapshot date: 2026-10-29**

## Degree-pair claim

The planar Jacobian pair `(75,125)` remains open in this repository. No unconditional exclusion of the pair is claimed.

## Leading-form split

For a hypothetical pair in characteristic zero, write its degree forms as `P_75 = a H^3` and `Q_125 = b H^5`, where `deg H = 25`. Over the algebraic closure, count distinct linear factors of `H`:

- One factor: `H = c L^25`. Not excluded in general here.
- Two factors: `H = c L^e M^(25-e)`, with `1 <= e <= 24`. Not excluded in general here.
- At least three factors: excluded by the at-most-two-points-at-infinity theorem for a plane Jacobian pair.

Thus `H = x^25 + y^25` is not a surviving example: it has 25 distinct factors over the algebraic closure and is in the excluded third stratum. The literal pair `(H^3,H^5)` has zero Jacobian; it represents only the leading forms of the hypothetical map.

## Computational results

The compact top-form and corrected-power checkers have replayable exact outputs in `results/`. The normalized F2 band notes record necessary coefficient equations, selected local eliminations, and surviving configurations. In particular, the quartic edge selector leaves `2+2` and `2+1+1` parameter types; this computation does not exclude either entire type as a global Keller pair, nor bridge every local carrier result to the original pair without the required reduction hypotheses.

The historical reports preserve their dated state. A relaxed degree-25 lift stopped after a FLINT allocation failure; that is a resource failure, not an algebraic contradiction. No computation is described as currently running.

## Scope rule for future claims

A branchwise proof would exclude the degree pair only if the branch split covers every possible hypothetical pair and each remaining branch is shown to have no admissible solution, with every normalization and coefficient-system bridge justified. Excluding one sample `H`, one finite truncation, or one local parameter neighborhood does not establish that conclusion.
