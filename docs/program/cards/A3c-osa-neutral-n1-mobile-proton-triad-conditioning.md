# A3c-osa-neutral-n1-mobile-proton-triad-conditioning — condition the evidenced mobile triad

- status: done
- track: A (geochemistry)
- priority: P1
- machine: workstation
- depends: A3b-osa-neutral-n1-proton-microstate-stability verified inconclusive
- claimed-by: hermes-custom-build-001

## Objective

Determine whether the original Osa-neutral n=1 proton assignment can support a
B3LYP/def2-SVP/DF constrained stationary seed when all three protons now observed
to move are held at their original owners: `H50-O26`, `H52-O27`, and `H57-O29`.
This is a materially different conditioning mechanism prompted by A3b's new
`H52:O27->O32` evidence, not a replay of A3b's dual-owner budget. Do not release
constraints, launch n=2-4, search a saddle/barrier, write a store, or emit Petra
data in this card.

1. Start only from A3a's hash-pinned conditioned geometry
   `279d168cab2f604c5bb9c3eef83e04adbd0d66bff15a2a35237f750c1d087a47`.
2. Rehash A3a's closeout/conditioned receipts and A3b's independently verified
   terminal evidence before any calculator call.
3. Run one fresh B3LYP/def2-SVP/DF optimizer with the unchanged frozen shell and
   exactly the three owner-bond constraints above. Spend one predeclared maximum
   of 100 steps, with no continuation or retry; persist the raw endpoint before
   every structural or chemistry gate.
4. Classify exactly one terminal outcome: **triad-conditioned stationary seed**
   requires optimizer convergence, finite gated projected gradients, all original
   proton owners retained, constraint residuals <=`1e-4 A`, and frozen-shell /
   collision gates; otherwise report **triad-conditioning failure** with the exact
   owner/gradient/optimizer evidence. A later card, not this run, decides any
   release design.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one branch,
  branch-matched worktree, and PR.
- GPU-first under QI2; OMP/MKL/OpenBLAS threads <=16; nice and tee the long run.
- Preserve source/settings/seed/code fingerprints, one-budget reservation, raw
  endpoint, exact terminal stage, and independent verifier receipts externally.
- Do not weaken SCF, optimizer, ownership, frozen-shell, collision, finite
  energy/gradient, or constraint-residual gates.

## Acceptance

- CPU-only regressions prove the exact three constraints, fresh optimizer state,
  one-budget/no-retry behavior, raw-endpoint-before-gate persistence, and stale
  receipt rejection.
- The production endpoint is bound to all A3a/A3b source hashes, settings, seed,
  code revision, geometry hash, optimizer status, finite energy/gradient, owner
  map, constraint residuals, and terminal stage.
- An independent verifier rehashes the source and endpoint and recomputes proton
  ownership plus all terminal gates before any stationary-seed claim.
- Full QM tests, whole-tree Ruff check/format, CLI dry-run, and `git diff --check`
  pass. Update A3/PLAN/STATUS in the same PR. No A3b replay or downstream
  scientific output is emitted.

## Progress

- 2026-09-05 23:46 PDT (hermes-custom-build-001; profile=workstation) — The sole production budget and independent verifier completed cleanly at exact implementation head `065e845bd0d5174091561b03f0cdd188dbd2647d`. The optimizer retained `H50-O26`, `H52-O27`, and `H57-O29`, with no owner changes and a maximum constraint residual of `5.134241e-6 A`, but independent projected RMS/max gradients were `1.619628e-3/1.229259e-2 Eh/Bohr`, above the `3.0e-4/4.5e-4` stationary gates. The exact terminal classification is `triad-conditioning failure`; A3d now owns the evidence-only terminal-path adjudication.
- 2026-09-05 22:29 PDT (hermes-custom-build-001; profile=workstation) — Implemented the one-stage A3c route: exact hash/semantic binding to A3a's conditioned seed and A3b's independently verified inconclusive evidence, one exclusive 100-step B3LYP/def2-SVP/DF GPU optimizer budget with the H50-O26/H52-O27/H57-O29 owner triad, raw endpoint persistence before gates, code-revision binding, stale-terminal revocation, and a separate optimizer-free verifier. Eleven CPU-only regressions pass, including concurrent-launch exclusion, wrong-method forgery rejection, clean/pushed-code enforcement, stale-success revocation, and exact triad constraints. Production compute remains unstarted pending the full-suite and launch-envelope gates.
- 2026-09-05 20:47 PDT (hermes-custom-build-001; profile=workstation) — Filed from A3b's independently verified inconclusive endpoint: the dual-owner seed retained H50/H57 but transferred `H52:O27->O32` and remained nonstationary after its sole 100-step budget. This card adds that newly evidenced third constraint and authorizes one bounded conditioning run only.

## Result

- 2026-09-05 23:46 PDT (hermes-custom-build-001; profile=workstation) — DONE with the card's permitted independently verified terminal outcome: `triad-conditioning failure`. geomeTRIC reported a fresh-optimizer convergence and all original proton owners were retained, but the separately recomputed finite projected gradients failed stationarity by `5.40x` RMS and `27.32x` max. Raw endpoint SHA-256 is `831198e8343fcfdb74e20a639caee459c6c57416f18b914b6f94d715180e26b2`; gated endpoint SHA-256 is `9a171cd7657ef65a34baf382e03faf85254a90d5b7116255d385f09aefa2f9ec`; stage receipt SHA-256 is `2f4a95dd3511f7981f33475f443e75c2e9c761cfd321c6a232fcde4f2ff242c6`; candidate-terminal SHA-256 is `db1e8cc3730730ecf6886c9054f8c8f75c7ede57d2ea78740256b298adf4bc47`; independent verified-terminal SHA-256 is `fee49de7528aad175b185da24498234ecd34e36f3d650623d45a8909d40521cb`. Source/settings/seed/code identities rehashed, the one-budget/no-retry contract held, shared services were restored, and no release, n=2-4, saddle, barrier, store, or Petra output was emitted.