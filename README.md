# Planar Jacobian Conjecture Research: Exact Reductions and Computational Records at the (75,125) Frontier

This repository collects source code, exact-computation records, and proof-scope notes for the planar Jacobian degree pair `(75,125)`.

## Current mathematical status

The pair `(75,125)` is **not excluded** here. For a hypothetical Keller pair with these degrees, the leading forms have the shape

\[
P_{75}=aH^3,\qquad Q_{125}=bH^5,\qquad \deg H=25,
\]

over an algebraic closure in characteristic zero. The choice `H = x^25 + y^25` has 25 distinct linear factors, so it belongs to the three-or-more-factor stratum. Abhyankar's at-most-two-points-at-infinity theorem rules out that whole stratum. The remaining one-factor and two-factor leading-form strata require separate arguments. Some normalized F2 computations impose additional exact necessary conditions and exclude selected local subcases, but no whole surviving stratum or the degree pair is excluded by this repository.

See [the factor-strata note](docs/leading-form-factor-strata.md), [the F2 band model](docs/F2-75-125-band-model.md), and [the dated status](STATUS.md).

## Repository map

- `src/`: research Python and Singular scripts copied from the direct project workspace.
- `src/singular/`: Singular input scripts.
- `docs/`: mathematical derivations, proof-scope notes, and dated progress reports.
- `results/`: small replayable certificates for the top-form audit and corrected-power edge selector.
- `data/json/`, `data/pkl/`: raw research records tracked with Git LFS; use `git lfs pull` after cloning. `data/MANIFEST.csv` records their source paths, sizes, and SHA-256 hashes.
- `references/`: bibliography and links to public source work.

The exact file inventory and exclusions are documented in [repository scope](docs/repository-scope.md). Raw records are research data, not proof certificates by themselves. Pickle files should only be loaded from trusted provenance.

## Replaying the compact exact checks

The two small audit scripts use Python 3 and SymPy 1.14.0:

```powershell
py -3 src/f2_75_125_top_form_branch_audit.py
py -3 src/f2_75_125_power_relation_edge_certificate.py
```

Their outputs are summarized in `results/`. They verify the stated finite algebraic identities and factor counts; they do not establish a global Jacobian-pair exclusion.

## Related work and attribution

This archive is a research notebook and computational supplement. It does not claim priority over or reproduce the independent public projects listed in [references](references/README.md), including Will Strinz's degree-`(75,125)` work and Roy van Rijn's Jacobian research repository. See each source for its own current claims and status.

## Updating the archive

Add new dated reports under `docs/progress/`, preserve the source and a replayable certificate for computational claims, and update `STATUS.md` only after checking the corresponding artifacts. Keep LFS data pointers and their manifest synchronized when adding raw checkpoints.
