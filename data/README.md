# Raw research records

The raw JSON and pickle records are preserved byte-for-byte and tracked with Git LFS. Run `git lfs pull` after cloning to materialize non-empty records. Two source `.pkl` files are zero bytes and remain ordinary empty Git blobs. `MANIFEST.csv` records source-relative paths, repository paths, sizes, and SHA-256 digests for all 427 records.

These files include checkpoints and exploratory states. Their presence does not make them theorem certificates. Python pickle files can execute code when loaded, so only load them if you trust the provenance.
