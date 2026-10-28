# Raw research records

`MANIFEST.csv` records the source-relative paths, intended repository paths, sizes, and SHA-256 digests for 398 raw JSON files and 29 pickle records (including two zero-byte `.pkl` files). The exact source files remain in the original research workspace. They are not included in this public snapshot because the initial multi-gigabyte Git LFS upload repeatedly lost DNS and transferred only a small fraction of the data. The small, curated checker outputs are included separately in `results/`.

To add the raw payload later, restore the source files at the manifest paths and upload them with Git LFS after confirming that the account has enough storage and the transfer path is stable. The `.gitattributes` rules are ready for that update.

These files include checkpoints and exploratory states. Their presence does not make them theorem certificates. Python pickle files can execute code when loaded, so only load them if you trust the provenance.
