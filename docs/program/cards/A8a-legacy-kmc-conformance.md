# A8a-legacy-kmc-conformance — execute the curated 1999 C model against golden runs

- status: active
- track: A (geochemistry / archaeology)
- priority: P2
- machine: any
- depends: A8-thesis-archive-intake
- claimed-by: hermes-macbot-zero

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
- Diffusion IDs 24–27 are historically inactive: `isActive` returns `FALSE`, so `new_evtList` never schedules them even though `doReaction` can `diffuse()`. Do not implement diffusion activity/target selection or delete the unreachable entries in the curated byte copies; pin the as-implemented disable.

## Acceptance

- Reproducible build command/toolchain lock.
- Five fixture outcomes with per-output hashes and semantic comparison.
- At least one sabotage test proves the comparator catches a perturbed trajectory.
- STATUS.md and this card updated in the delivering PR.

## Progress

- 2026-09-02 — The live pinned `linux/amd64` Docker run exposed a genuine archived portability wart: all five unshimmed fixtures raised `SIGSEGV`. AddressSanitizer traced independent post-buffer and pre-buffer initialization reads to `rxnlist.c:93` and `lattice.c:271`. Added a content-locked external `malloc`/`free` linker wrapper with 64-byte prefix/suffix slack, preserving curated source bytes and documenting the defect as a Rust-port follow-up rather than silently correcting it. The unshimmed failure and ASAN evidence are frozen; the shimmed canonical five-fixture run is in progress. (hermes-macbot-zero; profile=laptop)
- 2026-09-02 — Second independent review found two remaining fail-open paths: drift candidates could still make the canonical gate green and missing inputs could escape before report initialization. The canonical gate now requires five exact byte-parity outcomes; the sabotage control exercises that gate; and all normal report-mode setup/input/control failures emit five ordered fail-closed outcomes. Nine tests plus an actual missing-input report probe pass. (hermes-macbot-zero; profile=laptop)
- 2026-09-02 — Independent review blocked the first checkpoint on drift over-classification, container/worktree coupling, non-gating controls, incomplete provenance, exception handling, and process cleanup. All six classes were repaired: per-output divergence never infers cause; non-canonical runs can emit only a narrowly evidenced drift *candidate*; source/input manifests, architecture, compiler target, libc, container digest, and seeds are recorded; diffusion/cross-host/sabotage controls now gate; setup/future failures still emit five outcomes; run roots are collision-safe; and compiler/model process groups have bounded timeouts. Eight focused tests and two consecutive five-fixture one-second fail-closed rehearsals pass. (hermes-macbot-zero; profile=laptop)
- 2026-09-02 — Claimed atomically at `f76a73d`; froze the modern-clang failure before compatibility work. Added a digest-pinned GCC/Python container recipe, a source-byte-preserving forced-include compatibility layer, structural/numeric comparator, diffusion-disable audit, historical cross-host control, and sabotage regression. Six focused tests pass; five full 5,000,000-step replays are next. (hermes-macbot-zero; profile=laptop)
- 2026-09-02 — card filed from A8 intake after five compact, cross-host/parameter-span golden runs were curated. (hermes-custom-build-001; profile=workstation)
