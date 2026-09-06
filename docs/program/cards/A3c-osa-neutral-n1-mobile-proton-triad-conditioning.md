# A3c-osa-neutral-n1-mobile-proton-triad-conditioning — condition the evidenced mobile triad

- status: active
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

- 2026-09-05 22:29 PDT (hermes-custom-build-001; profile=workstation) — Implemented the one-stage A3c route: exact hash/semantic binding to A3a's conditioned seed and A3b's independently verified inconclusive evidence, one exclusive 100-step B3LYP/def2-SVP/DF GPU optimizer budget with the H50-O26/H52-O27/H57-O29 owner triad, raw endpoint persistence before gates, code-revision binding, stale-terminal revocation, and a separate optimizer-free verifier. Eleven CPU-only regressions pass, including concurrent-launch exclusion, wrong-method forgery rejection, clean/pushed-code enforcement, stale-success revocation, and exact triad constraints. Production compute remains unstarted pending the full-suite and launch-envelope gates.
- 2026-09-05 20:47 PDT (hermes-custom-build-001; profile=workstation) — Filed from A3b's independently verified inconclusive endpoint: the dual-owner seed retained H50/H57 but transferred `H52:O27->O32` and remained nonstationary after its sole 100-step budget. This card adds that newly evidenced third constraint and authorizes one bounded conditioning run only.

## Result