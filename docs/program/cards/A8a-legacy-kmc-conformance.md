# A8a-legacy-kmc-conformance — execute the curated 1999 C model against golden runs

- status: ready
- track: A (geochemistry / archaeology)
- priority: P2
- machine: any
- depends: A8-thesis-archive-intake
- claimed-by: —

## Objective

Turn A8's curated C source and five representative run fixtures into an executable conformance harness. Preserve the historical model semantics first; modernization belongs in the Rust port after behavior is mechanically captured.

1. Build `legacy/thesis-archive/c-model/source/` with a pinned compiler/container and document every compatibility flag or minimal source patch.
2. Replay the five `legacy/thesis-archive/golden-runs/` inputs in disposable directories.
3. Compare `results.dat`, `surfAl.out`, and `surfSi.out` structurally and numerically. Distinguish deterministic byte parity, PRNG/compiler drift, and true behavioral mismatch.
4. Emit a machine-readable conformance report suitable as a porting oracle for `petra/kmc-rs`.

## Constraints

- Never write to `/home/vsletten/Documents/personal/science/thesis`.
- Do not normalize a mismatch away. Preserve compiler, libc, architecture, seed, and floating-point provenance.
- If the historical source does not compile, freeze the unmodified failure first, then make the smallest reviewed compatibility patch in the curated copy.

## Acceptance

- Reproducible build command/toolchain lock.
- Five fixture outcomes with per-output hashes and semantic comparison.
- At least one sabotage test proves the comparator catches a perturbed trajectory.
- STATUS.md and this card updated in the delivering PR.

## Progress

- 2026-09-02 — card filed from A8 intake after five compact, cross-host/parameter-span golden runs were curated. (hermes-custom-build-001; profile=workstation)
