# A3-barrier-ladder — connectivity × protonation barrier ladder (Phase 2)

- status: active
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaigns; cluster-builder code is machine-any)
- depends: — (runs at the validated b3lyp/def2-svp/df tier; A2 re-tiers
  the finished ladder later; acid-state variants may reuse A1b's
  protonated-bridge machinery when it lands but are not blocked by it)
- claimed-by: hermes-custom-build-001

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

## Result (pilot cell — family campaigns remain)

**oss-neutral-n4-s2 @ b3lyp/def2-svp/df, PHVA thermochemistry:**
ΔG‡(298) = **205.7 kJ/mol** (49.2 kcal/mol), ΔH‡ = 187.0 kJ/mol,
ΔE_elec = 185.0 (vs complex) / 151.0 (vs fragments) kJ/mol, one
imaginary mode 128i cm⁻¹, Wigner κ = 1.02. Artifacts:
`qm/runs/phase2/oss-neutral-n4-s2-b3lyp-def2-svp/` (results.json,
store.sqlite, all stage checkpoints, teed log with all nine attempts).

**Mechanism**: sequential, not concerted — a PHVA-verified
pentacoordinate-Si intermediate (r(Si-Ow)=1.80, r(Si-Obr)=1.78, proton
mid-transfer; ~+92 kJ/mol above the complex) precedes rate-limiting
bridge rupture (TS2 at r(Si-Ow)=1.65, r(Si-Obr)=2.80). Quick-IRC
connects intermediate ↔ product. TS1 (attack) lies below TS2 and was
not separately converged (1-D scans cut its corner); not rate-limiting.

**Lattice resistance**: +92.7 kJ/mol vs the free-dimer 113.0 — above
Pelmenschikov's +20–70 window, but the absolute matches his ~205 kJ/mol
surface-embedded first-rupture ("self-healing") barrier almost exactly.
Suspected contributors to the high shift, all logged: the crude legacy
cell frozen at strained positions (learnings/gotchas), single explicit
water, pilot tier (A2 re-tiers). Trends (ΔΔEa by connectivity) are the
transferable payload, per the handoff.

**For the factory (campaign matrix per HANDOFF §3.4)**: the proven
recipe per cell is approach scan → proton-route fallback → product →
reverse crest scan → directed Cartesian Sella (reaction_path_vector
over the reactive core) → PHVA gates → thermo. Ops prerequisites on the
workstation: stop ollama + both ophir-email-pipeline services (dead-man
restarts), bounded memlock raise (prlimit 16 GiB), cuTENSOR preload —
all scripted in the run log and learnings. ~10 h GPU per cell realistic
including one saddle retry.

## Progress
- 2026-09-05 05:19 PDT (hermes-custom-build-001; profile=workstation) — The required-HF recovery failed honestly after 5,777.00 s with zero completed cells: constrained HF/STO-3G remained electronically stable but exhausted 100 geomeTRIC steps; over steps 91--100 its best RMS/max gradients were still 3.6x/6.7x above target and its worst were 17x/43x, so the endpoint is not a stationary point. No checkpoint or scientific output was promoted. DEVIATION: a cold scientific review selected the existing sibling-driver contract instead of an identical replay: persist the bounded HF endpoint only as an advisory seed, project the frozen shell back exactly, enforce atom/state/collision and settings/hash receipts, require one finite converged production gradient, and keep the production optimizer plus all downstream minimum/saddle/IRC gates unchanged. Seventy-three focused crystal/driver/pipeline tests, focused Ruff, Ruff format, and `git diff --check` pass.
- 2026-09-05 01:58 PDT (hermes-custom-build-001; profile=workstation) — The first Osa-neutral n=1 campaign terminated honestly after 707.50 s with zero completed cells: the raw 72-atom reactant complex improved from RMS/max gradient 0.07783/0.3742 to 0.04804/0.2275 Eh/Bohr, but its first geomeTRIC move reached 0.3404 A and the following B3LYP/def2-SVP DFRKS gradient exhausted the existing 150-cycle SCF bound. Constraints remained within 0.00110 A, no collision/OOM/lease failure occurred, and no trajectory geometry or barrier was promoted. The recovery adds checkpointed constrained HF/STO-3G relaxation of only the reactant complex before unchanged production optimization and scientific gates; the failed receipt remains preserved at `/mnt/data/vsletten/dissertation-data/task274-a3-barrier-ladder-20260905/terminal-receipt.json` (SHA-256 `6e2bb923b461515295eae8f0bc357a6c13e9614904a3c0bc2833fe17272c4071`).
- 2026-09-05 01:06 PDT (hermes-custom-build-001; profile=workstation) — Launched the bounded `osa-neutral` family campaign for connectivity rungs n=1..4 as transient user unit `task274-a3-osa-neutral.service` (PID 33659; `RuntimeMaxSec=44h`, cgroup kill, two-minute stop grace). The child driver acquired the canonical GPU lease as PID 33669, applied the 16 GiB CuPy ceiling and verified 16 GiB process memlock, then entered the n=1 reactant optimization at 89% GPU utilization / 7.0 GiB process VRAM. Durable external root: `/mnt/data/vsletten/dissertation-data/task274-a3-barrier-ladder-20260905/`; the atomic terminal signal is `terminal-receipt.json`. Until that receipt exists, every barrier and the original user outcome remain unverified.
- 2026-09-05 01:02 PDT (hermes-custom-build-001; profile=workstation) — Preflight completed: 54 focused crystal/driver/lease tests passed; direct builder inventory proved exact Osa rungs n=1..4 with 72/75/78/81 atoms, 20/23/26/29 frozen heavy atoms, neutral total charge, and 0.96 Å minimum pair distance; all three default family CLI dry-runs emitted geometry + metadata outside Git. Ruff 0.16.3 exposed one pre-existing 91-character line in `quarry/etiquette.py`; no source change or green Ruff claim was made.
- 2026-09-05 00:53 PDT (hermes-custom-build-001; profile=workstation) — Claimed `agents/A3-barrier-ladder` atomically from current `origin/main` (`59029f5e98b59ee74bfe1abbe1e34966827ac805`) after the upstream ByteQC QZ receipt released the GPU. Workstation preflight: 15-minute load 16.37 (<24), 30.0 GiB available, GPU at 716 MiB/24 GiB with no scientific compute process.
- 2026-08-22 — omnibus supervisor — restored `ready` after the pilot PR #67
  merged and released the claim; remaining family campaigns are unclaimed.
- 2026-08-21 — PILOT CELL DONE: ΔG‡ 205.7 kJ/mol, verified two-step
  mechanism, +92.7 lattice shift (see Result). Acceptance items 1–2 met
  (builder gates green; pilot reproduced modulo the logged shift, with
  the mechanism nuance recorded). Family campaigns + emission are
  factory work per Victor's 2026-08-20 handoff to mission-control.
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
