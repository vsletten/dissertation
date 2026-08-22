# A5p0-aging-observables — field-lab flagship Phase 0

- status: done
- track: A (geochemistry flagship)
- priority: P1
- machine: any
- depends: B4
- claimed-by: hermes-laptop

## Objective
The three observables the field-lab discrepancy study needs
(docs/scoping/field-lab-discrepancy.md Phase 0), built on B4's
framework: rate-spectra exporter (already partly in B4 — extend to
long-time series), BET-vs-geometric surface-area accounting, and
per-site exposure-age tracking. Then the smoke experiment: kaolinite
deck with a finite defect inventory, watch the bulk rate decline as the
inventory depletes — the qualitative H1 signature.

## Acceptance
- Observables tested; smoke experiment shows rate-vs-time decline with
  defect depletion, written up in docs/program/results/ with plots.

## Progress
- 2026-08-22T01:05:11+00:00 — hermes-laptop — Acceptance complete: 96 Petra workspace tests pass, including 20k-step kaolinite/Kossel bitwise parity; Clippy `-D warnings`, Rustfmt, diff checks, Python AST, SVG XML, and the seed-42 paranoid smoke are green. The smoke's 33 fast sites exhaust by step 40 while bulk propensity falls 330.0→1.15 (286.96×); results, summary CSV, and plot are committed under `docs/program/results/`.
- 2026-08-22T00:58:02+00:00 — hermes-laptop — Atomic `agents/A5p0-aging-observables` claim won. Implemented core-owned exposure timestamps/export, multi-state solid aliases, projected-geometric versus BET-like area accounting, and a reduced finite-defect kaolinite smoke deck. The seed-42 paranoid smoke exhausted 33 fast sites by step 40 and cut bulk propensity 286.96×; focused observable and production-kaolinite golden gates are green.
- 2026-08-22 — hermes-laptop — unblocked by B4's ensemble/observables framework; exposure-age core attachment contract is in `petra/docs/OBSERVABLES.md`.
- 2026-08-18 — card created (fable).

## Result

Delivered the complete Phase-0 observable slice:

- Rate spectra are cadence-sampled over long trajectories and exported in deterministic site order.
- Surface accounting now separates projected geometric footprint from the BET-like exposed-site census.
- Core-owned per-site exposure timestamps survive continuous exposure, clear on burial/removal, and restart on re-exposure without consuming simulation RNG.
- Multi-state aliases define mineral solids across heterogeneous site kinds; conflicting exposure definitions fail closed at compile time.
- Both production `kaolinite.toml` and Kossel declare the new observables.
- The reduced `kaolinite-aging-smoke.toml` supplies a finite damaged-site inventory and a reproducible qualitative H1 gate. The full write-up and plot are `docs/program/results/A5p0-aging-observables.md` and `docs/program/results/a5p0-aging-observables/aging-smoke.svg`.

Verification: `cargo test --workspace` = 96 passed / 0 failed; `cargo clippy --workspace --all-targets -- -D warnings` and `cargo fmt --all --check` pass. The paranoid smoke ran to its honest no-events stop at step 576; by step 40 its 33 fast sites were exhausted and total propensity had declined 286.96×. Existing seeded 20k kaolinite and Kossel parity gates remain byte-stable.
