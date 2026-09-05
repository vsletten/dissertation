# A3a-reactant-minimum-recovery — recover the first Osa-neutral reactant basin

- status: ready
- track: A (geochemistry)
- priority: P1
- machine: workstation
- depends: three hash-pinned A3 Osa-neutral n=1 failed receipts
- claimed-by:

## Objective

Recover exactly one scientifically valid production reactant minimum for
`osa-neutral-n1-s2-b3lyp-def2-svp`, then hand that checkpoint back to A3. The
current ladder driver lost its nonconverged B3LYP endpoint and its advisory
HF/STO-3G seed changed three termination-proton owners. Build one bounded,
microstate-preserving route; do not replay the current advisory seed, launch
n=2..4, or proceed to a saddle/barrier until the n=1 reactant basin closes.

1. Independently hash and preserve the three existing A3 evidence roots under
   `/mnt/data/vsletten/dissertation-data/task274-a3-barrier-ladder-20260905/`.
   Their terminal-receipt SHA-256 values are, in attempt order,
   `6e2bb923b461515295eae8f0bc357a6c13e9614904a3c0bc2833fe17272c4071`,
   `a25ba7611cd0d6d793e1a5fe291ce319c8a8cd825a14ad7e9bff23ccc48d9a66`,
   and `982ff963d9b2e421e0102881687771c085708e1ef2a938e7bff80bfa6b6092d9`.
   These owner-writable external files are preserved evidence, not immutable
   objects; every later use must rehash them against these pinned values.
2. Add a tested reactant-conditioning route that constrains each original O-H
   owner bond during only the cheap conditioning stage, preserves the exact
   frozen shell, atom order, charge, spin, and minimum-pair gate, then releases
   those conditioning constraints into a fresh B3LYP/def2-SVP/DF optimizer.
   Reject any unassigned, ambiguous, or changed oxygen-proton owner before the
   production PES is evaluated.
3. Make production optimization endpoint-preserving. On step exhaustion, write
   a hash/settings/input-bound nonconverged endpoint and terminal attempt
   receipt before raising. Permit exactly one predeclared continuation from
   that exact endpoint under a fresh optimizer/Hessian and unchanged production
   thresholds; never restart from the rejected advisory seed or silently add a
   third production budget.
4. Run only Osa-neutral n=1. If it converges, independently re-evaluate the final
   production gradient and apply exact proton-owner, frozen-shell,
   minimum-pair, finite-energy, and zero-significant-imaginary-mode gates before
   promoting `complex.xyz`. If it still fails, preserve a terminal failure
   receipt and file the next executable mechanism/convergence card. Either
   terminal outcome writes the atomic receipt at
   `/mnt/data/vsletten/dissertation-data/a3a-reactant-minimum-recovery/terminal-receipt.json`.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one branch,
  branch-matched worktree, and PR.
- GPU-first; OMP/MKL/OpenBLAS threads <=16; nice long CPU stages; tee every long
  run. Use the established QI2 GPU lease and bounded dead-man restoration.
- Existing A3 evidence is hash-pinned and must be revalidated before use. Do not
  weaken SCF, optimizer, basin,
  minimum, Hessian, or proton-ownership gates to rescue a number.
- Large run artifacts remain outside Git. Commit and push implementation and
  card checkpoints before the production run.

## Acceptance

- CPU-only regressions prove exact O-H owner constraints, owner-change /
  ambiguity / unassigned refusal, hash-pinned source evidence, settings/hash drift
  refusal, endpoint persistence before failure, one-continuation maximum, fresh
  optimizer state, and promotion only after convergence plus an independent
  final-gradient/minimum gate.
- The latest A3 advisory seed is mechanically rejected for exactly
  `H50:O26->O31`, `H52:O27->O32`, and `H57:O29->O20`; no production calculator is
  called after that rejection.
- The external root contains one atomic terminal receipt and a hash manifest.
  Success includes one provenance-bound Osa-neutral n=1 `complex.xyz` that
  passes every reactant-minimum gate. Failure names the exact exhausted stage
  and leaves no canonical complex, barrier, store, or Petra fragment.
- Full QM suite, whole-tree Ruff check/format, CLI dry-run, and `git diff --check`
  pass. A3/PLAN/STATUS are updated in the same PR; A3 remains blocked unless the
  reactant minimum is actually accepted.

## Progress

- 2026-09-05 10:19 PDT (hermes-custom-build-001; profile=workstation) — Pre-PR adversarial review found the generic Phase-2 publication path still trusted unsigned `complex.xyz`, failed to gate the production optimizer's returned geometry, omitted reactant-minimum and quick-IRC basin acceptance, and left stale `results.json`/`store.sqlite` canonical on a failed rerun. The code-only recovery slice now binds production checkpoints to exact input/settings/hash/geometry receipts, applies proton/frozen-shell/collision gates after production, rejects imaginary reactants and wrong IRC endpoints, and quarantines stale canonical outputs before live work. Thirty focused CPU-only tests plus Ruff/format/diff checks pass; no fourth calculation was launched.
- 2026-09-05 09:50 PDT (hermes-custom-build-001; profile=workstation) — Created from A3's third fail-closed Osa-neutral n=1 attempt. Independent review confirmed strict production geometry-optimization exhaustion, no valid scientific output, a discarded production endpoint, and three changed termination-proton owners in the advisory seed. This card is the only authorized next calculation; identical replay and n=2..4 remain prohibited.

## Result
