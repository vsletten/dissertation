# A3-barrier-ladder — connectivity × protonation barrier ladder (Phase 2)

- status: in-progress
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaigns; cluster-builder code is machine-any)
- depends: — (runs at the validated b3lyp/def2-svp/df tier; A2 re-tiers
  the finished ladder later; acid-state variants may reuse A1b's
  protonated-bridge machinery when it lands but are not blocked by it)
- claimed-by: fable (agents/A3-barrier-ladder, 2026-08-19)

## Objective
Execute qm/HANDOFF.md (Phase 2 edition — the full build plan lives
there): build the crystallographic cluster builder (terminated,
peripherally-frozen edge clusters cut from the kaolinite cell in
petra/examples/kaolinite.toml, per KMC site family and protonation
state), generalize the Phase-1 driver into a per-ladder-cell campaign
script, compute the connectivity- and protonation-resolved barriers for
the five site families, and emit them as petra `by_count`/`when` tables
with provenance. Closes CALCULATIONS.md rows CALC-002..005 from
`needed` to `computed`.

## Context
- qm/HANDOFF.md — the complete Phase-2 plan (read first)
- qm/SURVEY.md §6.2 (cluster norms), §7.3 (pKa table, kaolinite targets)
- qm/CALCULATIONS.md CALC-002..005 (the shopping list)
- Repo skills: quarry-campaign, petra-deck, pr-merge-gauntlet
- Reference result: si-neutral 113.05 kJ/mol (27.0 kcal/mol) in the
  claude/qm-phase1-scan-extend worktree runs/ dir
- .claude/learnings/ — the measured TS-hunt failure modes and gates

## Acceptance
- Cluster-builder unit gates green (CPU-only): taxonomy match against
  the deck cell, termination/charge bookkeeping, frozen shell
  peripheral, no collisions.
- Pilot cell (300s Si–O–Si, neutral, reference connectivity) reproduces
  the si-neutral barrier modulo a logged lattice-resistance shift.
- Per completed family: gated barriers in runs/ + store, an emitted
  petra fragment that compiles via petra-cli round-trip, and the
  CALCULATIONS.md row updated — one PR per family.
- Full fast test suite + ruff green on every PR.

## Progress
- 2026-08-21 — MECHANISM FINDING (pilot, oss-neutral-n4-s2): the embedded
  Si-O-Si neutral hydrolysis is SEQUENTIAL, not concerted like the free
  X&L dimer — Sella + PHVA verified a pentacoordinate-Si intermediate
  (r(Si-Ow)=1.80, r(Si-Obr)=1.78, proton mid-transfer at 1.26 A; zero
  significant imaginary modes) at ~+92 kJ/mol above the reactant
  complex. Lattice resistance stabilizes 5-coordinate Si. Per the
  handoff: both saddles recorded, rate-limiting TS2 (bridge rupture,
  crest near r(Si-Obr)=2.46) being converged with the directed
  Cartesian Sella (#40 machinery). Artifacts: intermediate.xyz +
  ts.rejected-slid-to-intermediate.xyz in the run dir. Ops lessons
  (GPU contention arbitration, memlock, cuTENSOR) in learnings +
  merged PR #55.
- 2026-08-19 — card created with the Phase-2 handoff rewrite (fable).
- 2026-08-19 — claimed; cluster builder (`qm/quarry/crystal.py`) +
  campaign driver (`qm/scripts/phase2_ladder.py`) built and gated
  (27 new CPU-only tests, full fast suite 117 green, ruff clean).
  DEVIATION: builder lives in a new module `quarry/crystal.py`, not as
  a clusters.py extension (clusters.py stays the hand-built benchmark
  set; the Cluster contract is shared). Pilot cell oss-neutral-n4-s2
  (Al3H25O29Si6, 63 atoms, 29 frozen heavies) launched on the GPU.
  Findings so far: (1) the deck cell faithfully imports the legacy
  data.cell crystallography, which is crude (bonds 1.38–2.88 A) — a
  logged systematic in every frozen-shell number (learnings/gotchas);
  (2) constructed protons are never frozen — only crystallographic
  heavy atoms; (3) frequencies() now does partial-Hessian analysis on
  frozen clusters (no trans/rot projection), so verify/IRC/thermo flow
  through the existing machinery unchanged.
