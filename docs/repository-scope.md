# Repository scope and provenance

## Included

- 342 direct Python source files from the research workspace.
- 30 direct Singular input files.
- A SHA-256 manifest for 398 direct JSON research records and 29 pickle records. The raw payloads remain in the source workspace and are not included in this public snapshot because the initial multi-gigabyte LFS transfer repeatedly failed DNS and made only limited progress; the small curated checker JSON files in `results/` are included.
- 14 small Markdown progress notes, the F2 band model, compact checker outputs, and all 30 selected direct Singular scripts.

`data/MANIFEST.csv` records original relative paths, intended repository paths, sizes, and SHA-256 hashes for the raw JSON/pickle files. Some records may be partial, superseded, or exploratory; their inventory is for reproducibility and provenance, not an assertion that every artifact is a final proof certificate. Python pickle is executable serialization: load only records from a trusted source.

## Excluded

- Nested third-party Git repositories, including Roy van Rijn's and Will Strinz's projects. Their work is cited in `references/README.md`; no source content is claimed or copied from them.
- Vendored dependency trees, the profiler installation, environment caches, and generated Python bytecode.
- Transient logs, temporary/partial outputs, profiler traces, and machine-specific PowerShell launchers that contain local paths or process IDs.
- The 1.17 MB operational branch-coverage log with stale worker IDs, memory snapshots, and repeated process polling; its mathematical coverage distinctions are summarized in `docs/progress/branch-coverage-ledger.md`.
- Superseded temporary fragments, unless a later report explicitly cites them as a reproducibility input.

The selection avoids silently mirroring collaborators' or upstream authors' repositories and keeps this archive focused on the user's direct working tree. File counts were measured at the 2026-10-29 snapshot; update this note and the manifest when the archive changes.
