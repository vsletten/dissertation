# A3-barrier-ladder — connectivity × protonation barrier ladder (Phase 2)

- status: blocked
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaigns; cluster-builder code is machine-any)
- depends: A3g-oaa-neutral-n2-proton-microstate-stability (READY)
- claimed-by:

## Constraints
Ready is **not** authority to continue the OSS-neutral (300s Si–O–Si)
series. A3f terminally rejected exact `oss-neutral-n1-s2`; OSS-neutral
n=2–4 has no accepted n=1 baseline and is **prohibited** (no A3e replay,
threshold relaxation, or surrogate barrier/store/Petra/CALCULATIONS
output). A3d already closed Osa-neutral n=1–4. A worker claiming this
parent must pick a **different** independently gated site-family /
protonation campaign already in this card's objective. HANDOFF §3 next
remaining family after the closed 400s and 300s neutral series is
**500s Al–OH–Al**, then attachment/detachment ladders and 100s/200s.
The crystallographic Oaa builder can realize only exact even coordination
rungs: n=2 at center 18, n=4 at center 23, and n=6 at center 18. Odd
requests alias lower even structures and are prohibited. The serial family
campaign starts at exact n=2 and must pass its reactant-minimum and normal
scientific gates before n=4, then n=6. The first n=2 attempt exhausted its
100-step HF/STO-3G advisory preoptimization and ended with ambiguous H52
ownership; no production optimization ran. Parent A3 is blocked on A3g's one
finite H52-O15 owner-constrained → released B3LYP test. The retained OSS-neutral
n=4 pilot below is historical/non-actionable.

## Objective
Execute qm/HANDOFF.md (Phase 2 edition — the full build plan lives
there): build the crystallographic cluster builder (terminated,
peripherally-frozen edge clusters cut from the kaolinite cell in
petra/examples/kaolinite.toml, per KMC site family and protonation
state), generalize the Phase-1 driver into a per-ladder-cell campaign
script, and compute connectivity- and protonation-resolved barriers for
remaining independently gated families. The original five-family list
is 100s Al, 200s Si, 300s Si–O–Si, 400s Si–O–Al₂, and 500s Al–OH–Al,
but this card **excludes** OSS-neutral n=1–4 (n=1 rejected; n=2–4
prohibited) and the closed Osa-neutral series. Claimable next work is
500s Al–OH–Al under its own reactant-minimum gate, not a prohibited
OSS n=2–4 route. Emit accepted remaining families as petra
`by_count`/`when` tables with provenance. Closes CALCULATIONS.md rows
CALC-002..005 from `needed` to `computed` only for families that
produce accepted reactant minima — never via OSS n=2–4.

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
- Historical OSS-neutral n=4 pilot thermochemistry is
  historical/non-actionable reconstruction evidence (see Result). It
  does **not** authorize OSS n=2–4, A3e replay, or treating 300s
  Si–O–Si neutral as the next claimable family.
- A worker claiming this ready parent must select a remaining
  independently gated family — not Osa-neutral n=1–4, not OSS-neutral
  n=1–4. Permitted next family: 500s Al–OH–Al (HANDOFF §3). That family
  uses only its exact crystallographic n=2/4/6 rungs (centers 18/23/18),
  rejects odd aliased requests, and requires its n=2 reactant minimum plus
  normal scientific gates before the serial campaign may advance.
- Per completed remaining family: gated barriers in runs/ + store, an
  emitted petra fragment that compiles via petra-cli round-trip, and the
  CALCULATIONS.md row updated — one PR per family.
- Full fast test suite + ruff green on every PR.

## Result (historical OSS-neutral n=4 pilot — non-actionable)

The OSS-neutral n=4 pilot material in this section is
**historical/non-actionable**. It is not live campaign authority, does
not keep the OSS-neutral family incomplete-but-claimable, and does
**not** permit OSS-neutral n=2–4.

- 2026-09-06 11:44 PDT (hermes-custom-build-001; profile=workstation) — A3f independently verified terminal rejection of the exact `oss-neutral-n1-s2` route. The owner-retaining constrained B3LYP endpoint passed RMS-gradient/residual/shell/collision/finite-energy gates but failed the immutable maximum projected-gradient gate (`4.836372724428273e-4 > 4.5e-4 Eh/Bohr`, `1.07474949432x`). Release, PHVA, and all downstream output remained absent; one 100-step budget per spent stage and zero retries are proven. This exact route is closed without a global no-basin claim: no A3e replay, threshold relaxation, OSS n=2-4, or surrogate output is authorized. Parent A3 is READY only for a different independently gated site-family/protonation campaign (permitted next: 500s Al–OH–Al). `docs/program/A3f-oss-neutral-n1-adjudication.md`.

**oss-neutral-n4-s2 @ b3lyp/def2-svp/df, PHVA thermochemistry (historical/non-actionable):**
ΔG‡(298) = **205.7 kJ/mol** (49.2 kcal/mol), ΔH‡ = 187.0 kJ/mol,
ΔE_elec = 185.0 (vs complex) / 151.0 (vs fragments) kJ/mol, one
imaginary mode 128i cm⁻¹, Wigner κ = 1.02. The historical ignored
artifact directory `qm/runs/phase2/oss-neutral-n4-s2-b3lyp-def2-svp/`
(results.json, store.sqlite, stage checkpoints, teed log) is **not
currently available** and is not banked source evidence. Reconstruct
the documented pilot under
`docs/program/cards/A2d-oss-neutral-n4-production-energetics.md`;
consume a later-discovered archive only after independent hash-bound
identity and provenance checks.

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

## Result — OSS-neutral family attempt (closed; n=2–4 not actionable)

Historical record of the closed OSS-neutral attempt. Do **not** read
"incomplete" dated bullets below as permission to finish the series via
OSS n=2–4. A3f closed the exact n=1 route; n=2–4 remains prohibited.

- 2026-09-06 09:02 PDT (hermes-custom-build-001; profile=workstation) — A3e completed with an independently verified `inconclusive terminal failure`. Its exact owner-constrained B3LYP endpoint retained `H35:O21`, passed RMS gradient/residual/shell/collision/finite-energy gates, but reproducibly failed the fixed maximum projected-gradient gate (`4.8363991055655475e-4 > 4.5e-4 Eh/Bohr`). Release, PHVA, and every downstream family stage correctly remained unspent. A3 stays BLOCKED without a barrier, store, Petra fragment, or CALCULATIONS closeout; READY A3f owns CPU-only terminal adjudication of this exhausted OSS-neutral n=1 route, and no A3e replay or OSS n=2-4 launch is authorized.
- 2026-09-06 07:06 PDT (hermes-custom-build-001; profile=workstation) — The bounded serial family supervisor worked as designed and stopped after `oss-neutral-n1-s2` failed its first chemistry gate. The exact `ec5bdb059911e32e4b790f849ff97595587c8519` execution completed zero cells: its HF/STO-3G advisory preoptimization converged at step 84, then strict ownership rejected `H35:O21->O14`. This was not a B3LYP production optimization and is not a production no-basin verdict. The rejected advisory endpoint was not persisted; `complex_preopt.xyz`, `complex.xyz`, `results.json`, and `store.sqlite` are absent, and n=2–4 never ran. Terminal receipt/log SHA-256 values are `c392dafe3c60f75687954295ae5b54da4c97a24ae1e54c9ec5b7cb8fdf9d9e53` / `b04a2991be5621925756222ae66bd41490c3152dc1ee442fcdfece590d346173`; shared services were restored and no A3 process or lease remains.
- 2026-09-06 07:06 PDT (hermes-custom-build-001; profile=workstation) — A3 remains BLOCKED without a barrier, store, Petra fragment, or CALCULATIONS closeout. READY `A3e-oss-neutral-n1-proton-microstate-stability` owns one materially different, finite owner-constrained → released B3LYP basin test with raw-endpoint persistence and independent verification. Blind replay, OSS n=2–4, and surrogate family emission remain prohibited until A3e produces an accepted n=1 reactant minimum or a verified terminal rejection.

## Result — Oaa-neutral family attempt (blocked on A3g)

- 2026-09-06 15:28 PDT (hermes-custom-build-001; profile=workstation) — **Independently verified computational failure at `oaa-neutral-n2-s2`, before production science.** The bounded family receipt completed with exit 1 after 2,153.50 seconds and zero completed cells. HF/STO-3G advisory preoptimization spent all 100 steps; the final logged RMS/max projected gradients were `8.942e-3/4.151e-2 Eh/Bohr`, and strict ownership rejected H52 as ambiguous between O15/O9 at `1.134/1.195 A` (margin `0.061 A`, below the fixed `0.15 A` gate). Because the raw endpoint was not persisted, this proves neither a completed H52 transfer nor a production-PES no-basin result. B3LYP production, PHVA, n=4/6, saddle, barrier, store, Petra fragment, and CALCULATIONS closeout did not run.
- 2026-09-06 15:28 PDT (hermes-custom-build-001; profile=workstation) — Receipt/progress/log/source/seed/metadata/launch/restoration hashes re-match the A3g evidence contract; exact source is `97cea5b585c95f18e499f53317ca8864a4c542d2`. The terminal and restoration units are inactive, the GPU lease is released, both email pipelines and both watchdog/deploy timers are active, and Honcho health passed. Independent evidence review found no inconsistency; a separate cold scientific review selected READY `A3g-oaa-neutral-n2-proton-microstate-stability` for one materially different, zero-retry H52-O15 owner-constrained → released B3LYP basin test with durable raw endpoints. Parent A3 remains BLOCKED; no identical replay or n=4/6 launch is authorized.
- 2026-09-06 16:00 PDT (hermes-custom-build-001; profile=workstation) — Final cold review found and the branch regression-closed five remaining false-green classes: stale approach/product reuse after TS incompatibility, quick-IRC frozen-coordinate drift, reduced lookalike Store schemas and inconsistent electronic barriers, bool/arbitrary-route acceptance, and missing physical all-rung Oaa builder coverage. Fresh final CPU-only verification passed `672` tests with no cache; whole-QM Ruff check passed, Ruff format reported all 77 files formatted, and `git diff --check` passed. The failed first `uv run pytest` attempt used an uninstalled fresh worktree environment and collected 31 import errors; the required main-venv invocation then passed the complete suite. No scientific calculator call was made during closeout.

## Progress
- 2026-09-06 15:28 PDT (hermes-custom-build-001; profile=workstation) — Reconciled the finite Oaa-neutral terminal receipt and two independent cold reviews. The supervisor correctly fail-closed at n=2 after advisory optimizer exhaustion plus ambiguous H52 ownership; no production or family output was emitted. READY A3g owns the only authorized Oaa continuation, with fixed gates and a terminal no-A3h default.
- 2026-09-06 16:00 PDT (hermes-custom-build-001; profile=workstation) — Final hardening and verification closed the last adversarial findings: incompatible TS checkpoints now quarantine transitive approach/product artifacts, quick-IRC checks finite/collision-safe endpoints and coordinate-level frozen-shell preservation, family validation requires the real Store schema and electronic-barrier consistency with canonical geometry hashes, booleans and unknown routes reject, and physical tests cover all Oaa n=2/4/6 centers/builds. Fresh full CPU-only QA is 672 passed plus whole-QM Ruff/format/diff gates. No production calculator call was made in closeout.
- 2026-09-06 14:02 PDT (hermes-custom-build-001; profile=workstation) — Claimed A3 for the independently gated neutral 500s Al-OH-Al family. A cold pre-launch audit blocked DFT after proving that requested Oaa n=1/3/5 silently alias even structures, default-center n=4 aliases n=2, terminal receipt ordering could contradict live progress, a junk SQLite file could pass as provenance, TS guesses were not bound to the optimized reactant/settings, and quick-IRC did not preserve charge/spin/frozen-shell or the structural bridge proton. The implementation now defines the exact Oaa domain and centers as n=2/4/6 at 18/23/18, fails closed on every alias, commits terminal receipts last, validates the read-only SQLite evidence graph, binds/quarantines TS-guess checkpoints, and gates endpoint identity plus structural bridge-H ownership. Physical dry-runs produce exact neutral 55/68/81-atom cells with 17/25/33 frozen atoms; odd n=1 rejects. Subsequent cold reviews regression-closed transitive stale-saddle and approach/product reuse, self-consistent lookalike stores, remote substrate-proton transfer, signal interruption at terminal commit, coordinate-level frozen-shell drift, bool/route false-greens, and missing all-rung builder coverage. Superseded test count: see the newer 16:00 line. No production DFT result is claimed.
- 2026-09-06 11:44 PDT (hermes-custom-build-001; profile=workstation) — A3f closed the exhausted OSS-neutral n=1 route by independent hash/arithmetic adjudication with no calculator calls. A3 returns to READY for a different family under its existing scientific gates; the rejected OSS route and its n=2-4 rungs remain prohibited.
- 2026-09-06 09:02 PDT (hermes-custom-build-001; profile=workstation) — A3e's independent verifier reproduced the constrained-production max-gradient failure with zero optimizer calls and confirmed release was correctly forbidden. Parent A3 is now blocked on READY A3f's one-pass terminal evidence adjudication; the fixed numerical gate will not be loosened post hoc and the exhausted A3e route may not replay.
- 2026-09-06 07:06 PDT (hermes-custom-build-001; profile=workstation) — Reconciled the exact failed family receipt and two independent cold reviews. Both confirmed the supervisor correctly fail-closed, and both rejected a terminal scientific inference from the discarded HF endpoint. The supervisor implementation stays in this slice; the parent is released as blocked on READY A3e rather than held active behind a dead campaign.
- 2026-09-06 03:23 PDT (hermes-custom-build-001; profile=workstation) — Atomically reclaimed `agents/A3-barrier-ladder` from exact `origin/main@5497d83da2d14dc7cd0dbb2679ce7acaeb2dd230` after A3d retired only the Osa-neutral series. The next permitted family is neutral 300s Si-O-Si (`oss`) across exact connectivity rungs n=1..4; CPU-only dry-runs produced distinct 37/54/60/63-atom complexes with neutral charge, peripheral 13/22/26/29-atom frozen shells, and no collisions. Added a bounded family supervisor that waits on the canonical GPU lease, pins clean local/remote source, runs explicit cells serially, and emits atomic progress/terminal receipts only after exact result identity, finite thermochemistry, route provenance, and store hashes pass. Full fast QA is 606 passed, 4 deselected; whole-QM Ruff check/format and `git diff --check` pass. The live A2b production campaign currently owns the GPU lease, so A3 will wait under a finite transient-unit ceiling rather than race it.
- 2026-09-06 01:21 PDT (hermes-custom-build-001; profile=workstation) — A3d independently verified terminal rejection of the original Osa-neutral n=1 assignment: geomeTRIC's constrained-run convergence did not satisfy the separately recomputed canonical-endpoint stationarity gate (`5.3988x/27.3169x` RMS/max over threshold). The exact Osa-neutral series is closed: no A3a/A3b/A3c replay, n=2–4 launch, or surrogate barrier/store/Petra value. A3 returns to READY for its other independently gated site-family/protonation campaigns; each still requires its own valid reactant minimum and normal scientific gates. Evidence: `docs/program/A3d-osa-neutral-n1-adjudication.md`.
- 2026-09-05 23:46 PDT (hermes-custom-build-001; profile=workstation) — A3c's one fresh three-owner B3LYP/def2-SVP/DF budget retained all original proton owners and met every `1e-4 A` constraint residual, but independently recomputed projected RMS/max gradients (`1.619628e-3/1.229259e-2 Eh/Bohr`) failed the stationary-seed gates despite geomeTRIC convergence. A3 remains blocked on READY A3d's evidence-only terminal-path adjudication; no A3a/A3b/A3c replay, n=2-4 launch, saddle, barrier, store, or Petra output is authorized.
- 2026-09-05 20:47 PDT (hermes-custom-build-001; profile=workstation) — A3b ended with an independently verified `inconclusive` outcome before either release route: its dual-owner production seed exhausted 100 steps, remained nonstationary, and exposed a third mobile termination proton `H52:O27->O32` while retaining the constrained H50/H57 bonds. READY A3c owns one materially different three-owner constrained-stationarity test. A3 remains blocked; no A3b replay, n=2-4 launch, saddle, barrier, store, or Petra output is authorized.
- 2026-09-05 17:37 PDT (hermes-custom-build-001; profile=workstation) — A3a terminated fail-closed after its B3LYP production optimizer converged but transferred `H50:O26->O31` and `H57:O29->O20`; no minimum, barrier, store, or Petra fragment was promoted. A3 now blocks on executable A3b's bounded two-order constraint-release test of whether the requested proton assignment is a production-level basin. No n=2–4 or fourth unconstrained replay is authorized.
- 2026-09-05 10:19 PDT (hermes-custom-build-001; profile=workstation) — Adversarial pre-PR review found four publication-path blockers beyond the failed advisory seed: unchecked production checkpoint reuse, no post-production microstate/shell gate, no reactant-minimum or quick-IRC basin gate, and stale canonical outputs surviving failed reruns. All four are now closed with strict/hash-bound checkpoint parsing, post-optimizer geometry acceptance, zero-significant-imaginary reactant plus hydrolysis-endpoint checks, stale-output quarantine, and focused adversarial regressions (30 passed). No fourth GPU attempt ran; A3 remains blocked on A3a's endpoint-preserving one-continuation recovery.
- 2026-09-05 09:50 PDT (hermes-custom-build-001; profile=workstation) — The one-shot advisory-seed recovery also failed honestly with zero completed cells after 13,781.27 s. Its HF/STO-3G endpoint passed finite production-gradient qualification, but the strict B3LYP/def2-SVP/DF reactant optimization exhausted 100 steps: final RMS/max gradients were 3.294e-4/1.274e-3 Eh/Bohr and RMS/max displacements were 3.651e-3/1.685e-2 A, all still above convergence. No `complex.xyz`, barrier, store, or Petra fragment was emitted. Independent read-only review also found the advisory seed changed termination-proton owners H50 O26→O31, H52 O27→O32, and H57 O29→O20. The driver now fails closed on exact oxygen-proton ownership (contract v2; regression-tested), and executable card `A3a-reactant-minimum-recovery` owns one microstate-preserving conditioning plus checkpointed production continuation. A3 remains blocked; no identical replay or n=2..4 launch is allowed.
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
