# C25/C26 checkpoint release and D27 status — 2026-11-06

## Public checkpoint archive

The C25/C26 checkpoint archive was published on GitHub under tag `f2-75-125-c25-c26-20261106`: [view the release](https://github.com/FarzanAhmedHossain/planar-jacobian-75-125-research/releases/tag/f2-75-125-c25-c26-20261106).

The release has 11 uploaded assets:

- C25: one checkpoint manifest and three ordered state parts, totaling 4,769,476,271 bytes. No standalone C25 receipt was located in the OCI handoff bucket or local archive.
- C26: one checkpoint manifest, five ordered state parts totaling 7,832,198,299 bytes, and the 71,834-byte D26 receipt.

The ordered C26 parts reassemble to the manifest-verified compressed checkpoint with SHA-256 `3d1f5e4aaafab2dd38d4b0d955f8730a9989c5f6cc299a39c007c245ce7034c5`; the gzip integrity check passed. The receipt SHA-256 is `ffbcdd20f142d11bb92757502be5bdc3b8c6c2aa73cc4a7f3e050b9f9e46a929`. GitHub also displays its automatically generated source-code archives; those are separate from the checkpoint assets.

## D27 snapshot

The latest read-only OCI snapshot recorded here is from 2026-11-06 07:55 UTC. The support-relaxed branch run is the $2+1+1$, depth-26, $U_1=0$ computation. It had entered distance 27, stage `unknown-series`, with $w_{\mathrm{count}}=60$ and $u_{\mathrm{count}}=12$. No C27 checkpoint or receipt existed at that snapshot.

At that sample, the worker used about 58.8 GB of resident memory, with about 0.6 GB available on the host. This is a resource-constrained computation and requires continued monitoring; the snapshot does not establish that the run has since completed or failed.

## Mathematical scope

These checkpoints and process snapshots are computational records for a necessary-condition system. They do not show that a Keller map of degrees $(75,125)$ exists, and they do not exclude the $2+1+1$ branch or the degree pair. The repository's [current status](../../STATUS.md) remains the authority for the proof-scope boundary.
