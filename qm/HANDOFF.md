# QM lab session handoff — Phase 2: the barrier ladder

*Written 2026-08-19 for a fresh session. Second edition — the Phase 0/1
handoff (build plan for the quarry platform itself) is in git history at
this path. Design authority remains `qm/SURVEY.md` (read §6–§8 at
minimum); the work ledger is `qm/CALCULATIONS.md` (CALC-002..005 are this
phase); the claimable card is `docs/program/cards/A3-barrier-ladder.md`
(claim per `docs/program/PROTOCOL.md` before starting).*

## 0. Mission

Compute the **connectivity- and protonation-resolved barrier ladder** for
the five KMC site families (100s Al, 200s Si, 300s Si–O–Si, 400s
Si–O–Al₂, 500s Al–OH–Al) on terminated, constrained cluster models cut
from kaolinite crystallography, and emit them as petra `by_count`/`when`
ΔEa tables with per-number provenance. This is the scientific payload the
whole platform was built for: the *non-flat environment tables the legacy
model never had* (its `data.rxn` gives every environment variant of a
reaction the same number).

## 1. State of the world (what already works — do not rebuild)

**Phase 0 + 1 are done and merged** (PRs #12–#19, more in flight):

- `quarry/` package: pipeline (PySCF energy/gradient/opt/freq, GPU via
  `DftSettings(use_gpu=True, density_fit=True)`), rates (quasi-RRHO →
  Eyring + Wigner/Eckart), emit (petra fragments + data.rxn), store
  (SQLite provenance), ts.py (Sella saddle search, verify, quick-IRC,
  auto-extending constrained scans with spike repair, CI-NEB guess).
- **First validated barrier**: si-neutral (H₂O + (HO)₃Si–O–Si(OH)₃,
  4-center concerted TS) — **ΔG‡(298) = 113.05 kJ/mol = 27.0 kcal/mol**
  at b3lyp/def2-svp/df vs Xiao & Lasaga ~29. TS: r(Si–Ow)=1.73,
  r(Si–Obr)=1.97, proton delivered, one imaginary mode (100i).
  Artifacts: `runs/phase1/si-neutral-b3lyp-def2-svp-flank/` in the
  `claude/qm-phase1-scan-extend` worktree (results.json, store.sqlite).
- **Measured hardware economics** (Phase 0): 4090 does analytic DF
  Hessians ~74× faster than the 32-thread CPU; gradients ~2–3 s at
  ~20 atoms/def2-svp. One barrier ≈ 45–90 min GPU wall. Frequency-heavy
  work ALWAYS runs on the GPU.
- **al-neutral is in flight, not yours**: fleet workers found it goes
  through a *sequential* mechanism (associative pentacoordinate
  intermediate, then bridge rupture — two saddles, rate-limiting step
  wins) — see TASK-168 in mission-control queue and dissertation PR #29.
  Do not duplicate; the mechanism lesson matters for the ladder though:
  **Si–O–Al chemistry may be stepwise where Si–O–Si is concerted.**
- Parallel cards on the board: `A1b-acid-mechanisms` closed the one-water
  al-acid path at 19.645 kcal/mol as a validated sequential mechanism, but
  the one-water si-acid protonated bridge is not an unconstrained minimum;
  it is BLOCKED on 3–6-water microsolvation (`qm/ACID_MECHANISMS.md`).
  `A2-production-energetics` remains
  (r²SCAN-3c / ωB97M-V/def2-TZVPD + SMD +
  DLPNO calibration — BLOCKED on Victor's ORCA install). Phase 2 runs at
  the validated b3lyp/def2-svp/df tier; A2 re-tiers the ladder later.
  Design every result file so re-tiering is a settings swap, not a redo.

**Read before running anything**: the three repo skills —
`quarry-campaign` (compute etiquette: tee everything, OMP cap 16,
GPU-first, checkpoint/resume, the failure signatures we already hit),
`pr-merge-gauntlet` (cloud Hermes owns the PR lifecycle now — open the
PR and stand down; `human-merge` label to hold), `petra-deck` (deck
authoring/validation for the emission step). Also
`.claude/learnings/` (debugging.md: the five TS-hunt traps, each
measured live; gotchas.md; performance.md).

## 2. What Phase 2 actually builds

Two genuinely new pieces, then a big campaign:

### 2a. Crystallographic cluster builder (`quarry/clusters.py` extension)

Phase 1 used hand-built dimers. Phase 2 cuts clusters from the real
mineral. **The crystallography is already machine-readable**:
`petra/examples/kaolinite.toml` `[cell]` carries the 26-site unit cell
(fractional coords, cell parameters, full bond list with `dcell` image
offsets, site kinds matching the five families) — parse it with stdlib
`tomllib`; cross-check against `legacy/cpp-model/data.cell`. Build:

- `from_deck_cell(deck_path, site_kind, ...)`: replicate the cell,
  pick a center site of the requested family at an edge, walk the bond
  graph to a cutoff (aim 8–20 tetrahedral/octahedral centers, 50–150
  atoms — SURVEY §6.2 norms), terminate dangling bonds with OH/H₂O
  (survey's termination rules), and mark the outer shell as
  `frozen_indices` (lattice resistance is physics, not a detail:
  Pelmenschikov's +5–16 kcal/mol vs free clusters). The frozen-atom
  contract already flows through optimize/scans/Sella/NEB untouched.
- Charge/protonation assignment per the pKa table (SURVEY §7.3):
  Si–OH ≈ 6.9, Al–OH₂OH ≈ 5.7, bridging Si–O ≈ 2, octahedral
  Al–OH₂ ≈ 10. For each family compute the states that matter across
  pH 2–10 (typically neutral + one protonated or deprotonated variant).
- Gates: stoichiometry/charge bookkeeping, no atom collisions, frozen
  shell actually peripheral, cluster matches the deck's site taxonomy.
  All testable CPU-only without DFT.

### 2b. Ladder campaign driver (`scripts/phase2_ladder.py`)

Generalize the Phase-1 driver from "one reaction" to "one (family,
connectivity n, protonation state) cell of the ladder":

- Reuse the whole proven stage machinery (pre-opt → opt → approach scan
  → proton scan → product construction → CI-NEB → Sella → in-channel +
  escaped-saddle gates → verify → quick-IRC → quasi-RRHO). It is all
  importable; the driver mostly supplies geometry and bookkeeping.
- Connectivity dimension: for the target bridge/site, vary the number of
  intact neighbors n by cutting clusters around sites with different
  coordination (edge vs kink vs terrace analogs). Each n is one ladder
  cell → one `dea[n]` entry.
- Every cell writes `runs/phase2/<family>-<state>-n<count>/` with
  checkpoints, results.json, store.sqlite — resumable, batchable.
  Expected failure modes and their signatures are in
  `.claude/learnings/debugging.md`; the gates abort loudly rather than
  report plausible-wrong numbers. Sequential (two-saddle) mechanisms
  are legitimate outcomes — record both saddles, report the
  rate-limiting ΔG‡ (the al-neutral lesson).

### 2c. Emission + ledger closure

- For each reaction family: baseline barrier at the reference
  connectivity → petra `rate = {eyring = ...}`; ΔEa(n) relative to
  baseline → `by_count = {dea = [...]}` tables; protonation variants →
  separate `[[reactions]]` entries gated by `when`/state vocabulary
  (state names documented in the kaolinite deck header). Emit with
  `quarry.emit` (provenance comments are mandatory — method, geometry
  hash, date); splice-test against the template and compile with
  petra-cli (the round-trip test pattern exists; `petra-deck` skill).
- Update `qm/CALCULATIONS.md` rows CALC-002..005 (`needed` → `computed`
  with the value/provenance pointer) in the same PR as each batch of
  numbers — invariant: merged means ledger updated.
- Regenerating legacy `data.rxn` stays Phase 3 (needs the full 28-slot
  table); don't block the ladder on it.

## 3. Suggested execution order

1. Claim `A3-barrier-ladder` per PROTOCOL (branch push = atomic claim),
   worktree `agents/A3-barrier-ladder`.
2. Build + gate the crystallographic cluster builder (CPU-only tests;
   this is the bulk of the new code). PR early if it stands alone.
3. Pilot ONE ladder cell end-to-end on the GPU: Si–O–Si bridge (300s),
   neutral, reference connectivity — should reproduce si-neutral's
   113 kJ/mol within the cluster-size effect (expect a shift UP from
   lattice resistance; that shift is itself a publishable observable —
   log it prominently).
4. Then batch: one family per overnight campaign (~5–10 cells each,
   ~1 h/cell GPU). Order by KMC impact: 400s Si–O–Al₂ (the kaolinite
   deck's biggest family) → 300s → 500s Al–OH–Al → attachment/
   detachment ladders (CALC-005, these are minima-only — no TS hunt,
   much cheaper) → 100s/200s adsorption sites.
5. Emit + ledger PR per completed family, not one giant PR at the end.

## 4. Hard rules (unchanged, learned the hard way)

- **Compute etiquette**: everything teed to files (never trust terminal
  scrollback), OMP_NUM_THREADS ≤ 16 for CPU stages, GPU for all DFT,
  no multi-hour full-CPU jobs without Victor's explicit go. The
  `quarry-campaign` skill is the checklist.
- **Worktree discipline**: one branch = one worktree = one PR; main
  worktree is read-only reference.
- **PR lifecycle**: open it, let cloud Hermes drive Sourcery/merge
  (`pr-merge-gauntlet` skill); after merge, sync local main.
- **Barriers cross the QM→KMC bridge, never pre-exponentiated rates**
  (petra's truncated R vs CODATA — gotchas.md).
- **DFT owns every number that enters the KMC**; xtb/MLIPs may only
  screen and guess (SURVEY §5).
- **A verified saddle is not necessarily the right saddle** — the gates
  exist because we produced three chemically-wrong-but-mechanically-
  valid saddles in one evening. Trust the gates; when one aborts, read
  the geometry before touching the code.
- Big outputs stay out of git (`runs/` is gitignored); derived numbers +
  provenance go in the store and the ledger.

## 5. Open questions a session may hit (decide, log, proceed)

- Cluster termination chemistry at Al edges (OH vs H₂O vs OH₂⁺) sets the
  cluster charge — follow the pKa table for the target pH regime and
  record the choice per cell in results.json.
- How many explicit waters at the reactive site: survey norm is 3–6 for
  water-mediated steps; Phase 1 used 1. Start with 1 for ladder
  *trends* (ΔΔEa by connectivity is more transferable than absolute
  heights), flag absolute numbers as needing the A2 tier + explicit
  microsolvation.
- If a ladder cell's mechanism differs qualitatively across n (concerted
  → stepwise), that's a finding, not a bug — record both and emit the
  rate-limiting barrier.
