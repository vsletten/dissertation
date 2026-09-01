# Program status log

*Newest first. One line per merged unit of work (or notable event),
appended by the PR that did it: `date — actor — what + pointer`.
*Deviations get a `DEVIATION:` note; details live in the card.*

- 2026-09-01 — (hermes-8d9e30fee91548; profile=laptop) — C3 pit
  statistics DONE: 12,032 deterministic passive-film replicas, censor-aware
  Weibull/bootstrap and survival/hazard diagnostics, pit size/event-depth +
  fixed-occupancy clustering distributions, 157.0× nucleation-rate scaling,
  and a four-size half-initiation crossover point estimate land at
  `docs/program/results/C3-pit-statistics.md`.

- 2026-08-29 — fable — Fleet routing widened per Victor (POLICY v13 §10):
  workstation runs two concurrent lanes (GPU stays single-lease via QI2),
  macbots default-active for machine:any, cloud Hermes may hold one queue
  task alongside PR watch. Board: E3a re-scoped machine:any (classical
  potentials are CPU work), new C3-pit-statistics (READY, any) so
  laptop/cloud workers have C-track science. A8 enqueued into the
  workstation's second lane.

- 2026-08-28 — fable (hermes-custom-build-001; profile=workstation) —
  D2b done: explicit-surface rate campaign, 10/10 reactions clean on
  the 4090 (surface PFs + open-shell UCCSD(T)/CBS + (H₂O)₁,₂ site
  sampling). Gate verdict NO-GO with receipts — no D3b rate table
  (addition/H+CO 25–1600× under the Simons-corrected anchor at
  10–20 K; Eckart-vs-instanton is the isolated residual) — but the
  barrier tier validates (1w abstraction = S&K's surface 3030 K to
  2 K) and abstraction k(T) lands in-spread at every anchor. D2a's
  "4–5 orders" decomposed; much of it was the anchor itself
  (Simons+20 document the S&K fits as overestimates). New card
  D2c-instanton-tier (READY) flips the gate; D3b re-blocked on D2c.
  Results: docs/program/results/D2b-explicit-surface-rates.md.

- 2026-08-28 — (hermes-8d9e30fee91548; profile=laptop) — E2b runs seven
  32-replica muscovite volumes through 559,872 sites: the old-first-step
  staircase stabilizes at 49,152 sites and stays stable across two larger
  volumes; 559,872 is the completed host ceiling after the 1,327,104-site rung
  exits by `SIGKILL`. Full bootstrap bands, isotope-release distributions,
  exact-seed replay receipts, and scaling are in
  `docs/program/results/E2b-grain-size-sweep.md`.

- 2026-08-28 — (hermes-8d9e30fee91548; profile=laptop) — A5 Phase 1 runs
  4×64 finite-defect kaolinite-aging ensembles: all replicas are defect-free
  by step 90 at the latest, fresh/aged apparent-rate drops grow 36.36×→363.06×
  with initial damage density, spectra transiently broaden then down-shift
  three orders, and BET-like area rises while rate falls. The honest
  mechanism-tier verdict is partial plausibility (1.53–2.56 log units), not a
  calibrated field attribution. `docs/program/results/A5p1-aging-study.md`.

- 2026-08-27 — fable — Track E phase-3/paper-1 path carded: E3a
  classical-NEB barrier campaign (READY — the "uncomputed middle": Nteme
  replication gate, then dehydroxylate/defect-zone/trap/Xe barriers), E3b
  CP2K spot checks (blocked: E3a), E2b grain-size sweep (READY — closes
  E2's 768-site statistics deviation), E4 1998-comparison (blocked: E3a,
  E2b, A8; the paper-1 figure set). Also A5p1-aging-study (READY): M5 named it
  but it was never carded. Context: board starvation 08-25→08-27 was a
  feeder machine-parser bug, fixed in mission-control 387300f.

- 2026-08-27 — (hermes-8d9e30fee91548; profile=laptop) — C2 ships a
  strain-seeded passive-film CTMC deck plus a deterministic 7×32-replica study
  and three-size scan: initiation rises 0/32→32/32, mean pit count
  3.47→36.22, and mean largest cluster 1.25→7.28 across the potential sweep;
  descriptive induction plots, exact-deadline CSVs, and a metastable
  current-rate series are archived in `docs/program/results/C2-corrosion-deck.md`.

- 2026-08-27 — (hermes-custom-build-001; profile=workstation) — QI2 makes the RTX 4090 a single explicit lane: atomic process/PID lease, dead-PID stale-break receipts, fail-closed symlink-safe state, SIGTERM/fork hardening, shared `scripts/gpu_preflight.sh`, and a 16 GB CuPy pool cap preserving 6 GB for ollama. Phase 1/2 now refuse contention cleanly; 436 QM tests, Ruff/format, live 4090 preflight, separate-process busy exit, and signal release are green.

- 2026-08-26 — (hermes-custom-build-001; profile=workstation) — A2a's corrected full-system local dimer exhausted all 120 bounded translations inside the exact climb guard with strictly positive curvature (`+0.0274 eV/A²` minimum), finished wall-supported, and persisted no candidate or downstream number. Pushed head `95ea8df`; 425 QM tests and whole-QM Ruff/format pass; evidence manifest `9e163528bb43...` covers 478 files. A2a is scientifically blocked pending a new mechanism/path decision, not another local Sella/dimer retry.

- 2026-08-26 — hermes-workstation — A2a I↔P now has a tested neighbor-radius
  flat-bottom trust envelope, hard pre-PES escape guard, capped Sella step, and
  mandatory cuTENSOR GPU preflight. The 300-step cleavage localization stayed
  inside the local band neighborhood but cycled near `0.3--0.6 eV/A`, so no
  saddle or downstream number was emitted. Pushed head `9b9edb2`; evidence
  manifest `8542944ade6f...` (461 files). Next: bounded full-system local
  eigenvector-following/dimer localization, not another active-core retry.

- 2026-08-25 — hermes-workstation — A2a's internally conditioned production
  path now closes R↔I with a stable 714.83 cm^-1 reaction saddle and exact typed
  full IRC; I↔P pre-relax/climb/conditioning also converge, but its active-core
  Sella search escapes uphill into an SCF failure before a candidate is emitted.
  Commit `8ced110`; durable evidence manifest `82af53d69809...` (457 files).

- 2026-08-24 — hermes-workstation — A2 production infrastructure is pushed and
  fail-closed on the first real target: exact r2SCAN-3c D4+gCP optimization,
  same-surface finite-difference Hessians, repaired Sella/ASE full IRC, and
  verified Psi4/ByteQC calibration adapters. The si-neutral minima/TS have
  indices 0/0/1, but a real 260-step full IRC connects the banked saddle to the
  associative basin on both sides; directed Sella reconverges it. No production
  or CC barrier was fabricated. A2 is blocked on READY path-rebuild card A2a;
  evidence manifest `efdcfc48c533...` (47 files, 791,566 bytes).

- 2026-08-24 — hermes-workstation — A1g's hash-pinned bridge-side-H3O+ /
  distinct-neutral-attacker Si reactant is a zero-imaginary minimum, but its
  independent product converges outside the hydrolyzed family: Si--Ow 2.643 A,
  attacker H21 unassigned, basin `(False,False,True,1,7,0)`. This is a
  conclusive mechanism-v3/gate-v3 product rejection with no computational
  failure and no eligible Al/barrier/ordering. The only acid follow-through is
  A1i's one calibrated production-tier revisit after A2; no more SVP topology
  retries.

- 2026-08-24 — hermes-workstation — A1f's one four-water separated-donor /
  neutral-attacker Si path is a conclusive reactant-family rejection: the
  unconstrained B3LYP/def2-SVP/DF endpoint moved physical H16 from outer donor
  O15 to bridge-side relay O22, changing exact solvent occupancies from
  `(3,2,2,2)` to `(2,2,3,2)` while keeping the attacker neutral, bridge intact,
  and ownership valid. No computational stage failed; no Al/barrier/ordering
  was emitted. Mechanism-v2/gate-v2, 360-test proof, hash-bound archive, and
  READY A1g bridge-side-H3O+ follow-through land in the A1f card.

- 2026-08-24 — hermes-workstation — A1e post-merge activation closeout
  independently rehashed all 24 campaign evidence files, reconfirmed four
  typed Si endpoint rejections with zero failed/blocked/running paths and READY
  A1f follow-through, and passed 341 QM tests plus CLI/diff gates. A narrow
  follow-up normalizes the review-added regression test so whole-tree Ruff lint
  and formatting are both green on the final merged behavior.

- 2026-08-24 — hermes-workstation — A1e's exact concerted hydronium-as-
  nucleophile campaign is model-valid NO-GO: all four 3/4-water bridge-chain
  and compact-cycle Si paths are conclusive typed endpoint rejections, with
  zero optimizer/NEB/Hessian/IRC failures and therefore no eligible Al,
  barrier, or ordering. Mechanism-v1/gate-v1, full bidirectional Sella IRC,
  physical-H relay ownership, hash-bound receipts, and 338-test proof land with
  narrower READY `A1f-acid-neutral-water-attacker-relay`, which separates H3O+
  donor from the neutral-water nucleophile.

- 2026-08-24 — hermes-workstation — A1d closes A1c's sole inconclusive
  three-water seed and the finite pre-equilibrium: one hash-pinned fresh
  B3LYP/def2-SVP/DF continuation converged after 54 steps, then failed the exact
  protonated-bridge occupancy gate. Effective result is 16/16 conclusive Si
  rejections, so the pre-equilibrated acid bridge is model-valid NO-GO; no Al
  screen or unmatched ordering was emitted. A1c/A1d are done and executable
  `A1e-acid-concerted-hydronium-relay` is ready.

- 2026-08-24 — hermes-workstation — A1c's finite 16-seed Si conformer
  campaign is honestly BLOCKED, not falsely NO-GO: 15 B3LYP/def2-SVP/DF
  endpoints converged and failed exact protonated-bridge occupancy, while the
  three-water bridge-donor chain exhausted its 160-step bound. The hardened
  ensemble driver, hash-bound external archive, and executable A1d one-seed
  convergence closeout are in `A1c-acid-microsolvation-conformers`; no Al
  screen, barrier, or acid ordering was emitted.

- 2026-08-23 — hermes-workstation — A1b matched 3--6-water Si-acid
  deterministic-shell survey is honestly BLOCKED: all four converged
  B3LYP/def2-SVP/DF endpoints deprotonate Obr before any Hessian/TS work.
  Driver support, physical-H/shell identity gates, non-destructive reactant
  receipts, hashes, and executable follow-up card
  `A1c-acid-microsolvation-conformers` are in `qm/ACID_MECHANISMS.md` and
  the A1b card; no unmatched Al ordering was emitted.

- 2026-08-23 — fable — A2 unblocked: ORCA descoped from the platform
  (Victor's call, licensing). Calibration layer is now Psi4 1.11
  DLPNO-CCSD(T) + ByteQC canonical GPU CCSD(T), cross-validating; both
  installed and smoke-verified on the workstation (water dimer/cc-pVTZ:
  canonical CCSD(T) agrees across engines to 1e-8 Eh; DLPNO −0.53
  kJ/mol off canonical). Runbook `qm/CALIBRATION.md`; SURVEY §§2.3–2.4/§6.4
  amended; A2 card re-scoped → READY (feeder-eligible).

- 2026-08-22 — hermes-workstation — A1b partial closeout: validated one-water
  al-acid sequential hydrolysis at ΔG‡(298) = 82.193 kJ/mol (19.645 kcal/mol),
  with exact one-mode/quick-IRC gates for addition and cleavage. DEVIATION:
  one-water si-acid is not an unconstrained protonated-bridge minimum in gas
  phase or PCM, so the card is blocked on 3–6-water microsolvation rather than
  publishing a constrained fake barrier. `qm/ACID_MECHANISMS.md`.

- 2026-08-22 — omnibus supervisor — reconciled board metadata after merged
  work: C2 is ready after B3, A3 is ready for its remaining family campaigns,
  and the completed D2a card now has its missing PLAN.md row.

- 2026-08-22 — hermes-laptop — D3 DONE: a 24×24 open-system H₂ grain deck
  now exercises source deposition, thermal+tunneling diffusion, LH reaction,
  H/H₂ desorption, and seeded 20% quenched deep sites. The 16-replica sweep
  has the published finite window (efficiency 0.310 at 6 K, ~0.96 at 12–14 K,
  0.559 at 20 K, 0.055 at 26 K). DEVIATION: D2a surface rates are explicit
  NO-GO; D2b + blocked D3b cards preserve the CO follow-through honestly.
  `docs/program/results/D3-ice-mantle-deck.md`.

- 2026-08-21 — fable — A3 pilot cell DONE: embedded Si–O–Si neutral
  hydrolysis ΔG‡ = 205.7 kJ/mol (free dimer 113.0; +92.7 lattice shift,
  matching Pelmenschikov's ~205 embedded first-rupture) via a verified
  two-step mechanism (pentacoordinate-Si intermediate, rate-limiting
  bridge rupture, 128i, IRC-connected). CALC-002 neutral computed.
  Campaign recipe + workstation ops prerequisites in the card Result;
  family ladders are factory work per the 2026-08-20 handoff.

- 2026-08-22 — hermes-laptop — E2 DONE: scheduled full muscovite mechanism
  now distinguishes pristine/extended galleries, explicit locally driven
  delamination, ⁴⁰Ar/³⁹Ar/³⁶Ar, and octahedral recoil traps; deterministic
  release/age/contamination products and honest ±5 kcal proxy bands land at
  `docs/program/results/E2-muscovite-full-mechanism.md`. E1's 500/700 °C gates
  remain green.
- 2026-08-22 — hermes-laptop — B5 DONE: ordered piecewise-isothermal CTMC
  schedules now advance through exact wall-time boundaries, recompile full
  thermo/propensity tables + Fenwick trees, and run through native CLI and
  deterministic ensembles. Exact T1 replay, seeded T2 statistics, v1 20k-step
  parity, workspace tests, format, Clippy, and independent review are green.
- 2026-08-22 — hermes-laptop — A5 Phase 0 DONE: long-time rate spectra,
  projected-geometric versus BET-like surface accounting, and core-owned
  per-site exposure ages landed with deterministic alias-aware tracking. The
  finite-defect kaolinite smoke exhausts 33 fast sites by step 40 and cuts bulk
  propensity 286.96×; 96 tests, bitwise parity, Clippy, and paranoid replay are
  green. `docs/program/results/A5p0-aging-observables.md`.
- 2026-08-22 — hermes-laptop — B4 DONE: Rayon-parallel deterministic replicas,
  full distributions/means/bootstrap CIs, and declarative state-count, event-rate,
  rate-spectrum, cluster-size, and geometric/BET-like area outputs; Ising Binder
  dogfood + Kossel spectrum-broadening golden green, with all 88 Petra tests and
  Clippy clean. `petra/docs/OBSERVABLES.md` records the A5 exposure-age seam.
- 2026-08-21 — hermes-laptop — B3 DONE: Petra now runs synchronous CA,
  asynchronous Metropolis, and discrete-time PCA on schema-v2 grids; Conway,
  Ising (`Tc=2.273230` vs `2.269185`), and SIR (`R0*=1.017094`) analytic gates
  are green with CTMC bitwise parity preserved.
- 2026-08-21 — omnibus supervisor — reconciled merged-main board drift:
  B2 and QI3 are done; B3, B4, and B5 are now ready after B2 PR #33.
- 2026-08-20 — hermes x2 + fable — B2 DONE: petra engine behind the
  UpdateStrategy trait, schema v2 complete, bitwise parity green on
  three gates (incl. shipped-kossel full-lifetime). Two-worker relay +
  supervision closeout; 7 self-review findings all dispositioned and
  fixed. B3 unblocked.
- 2026-08-20 — hermes-laptop — QI1 campaign-driver etiquette done: all QM
  entry points now cap inherited numerical-library threads before heavy imports,
  apply best-effort niceness, and self-tee durable logs; 162 tests, Ruff, direct
  CLI smoke checks, and an unset-environment Phase-2 dry run are green.
- 2026-08-20 — hermes-cloud-default — QI3 Fenwick cancellation hardening:
  impossible zero totals now rebuild from authoritative positive leaves; a
  deterministic `1e12`/`1e-6 s⁻¹` regression, all 43 Petra tests, and an
  explicit 1000-step Kossel `--paranoid` run are green.

- 2026-08-20 — hermes-workstation — E1 exact-main lint closeout: restored the
  one-line `strain_gates.rs` identity-op fix that the prior squash merge omitted;
  all Rust/Python/lint gates plus both full trajectories and generated-product
  equivalence are green.
- 2026-08-20 — hermes-workstation — E1 Ruff-format closeout: reformatted the
  cloud review-fix commit's long SVG `sy()` expression so merged-main Ruff
  0.15.1 check/format gates are clean; no simulation semantics changed.
- 2026-08-20 — hermes-workstation — E1 merged-main closeout: removed four
  default-Ruff E731 violations from the muscovite SVG post-processor with
  typed local helpers, cleared current Clippy findings, and reverified both
  full trajectories plus byte-identical CSV/SVG/JSON products.
- 2026-08-20 — hermes-workstation — E1-muscovite-deck done: schema-v1
  500/700 °C paired-OH + divacancy decks, Hames–Bowring cylinder `D/a²`
  inversion, reproducible CSV/SVG/JSON products, and qualitative gates green
  (rise 2.13×/13.55×; fall 275×/32.3×). Result:
  `docs/program/results/E1-muscovite-phase1.md`. Added executable E2 and QI3
  cards for the full mechanism and Fenwick dynamic-range hardening.
- 2026-08-20 — fable — A3-barrier-ladder claimed + platform half landed:
  crystallographic cluster builder (`qm/quarry/crystal.py`, deck-cell
  parser + bond-valence termination/charge + n_intact pruning + frozen
  shell), PHVA frequencies for frozen clusters, and the per-cell campaign
  driver (`qm/scripts/phase2_ladder.py`). Pilot cell oss-neutral running
  on the GPU; barriers land in follow-up per-family PRs.
- 2026-08-20 — fable — Card QI1: campaign drivers must enforce their own
  thread caps + niceness (a live worker probe ran 22 threads past the
  prompt-level cap — prompts don't bind subprocesses).
- 2026-08-20 — fable — Board refresh: B1+D2a+A7 done, B2 READY (keystone
  merged), A3 active. PROTOCOL amendment: incremental bookkeeping
  (the D2a iteration-ceiling lesson). Hermes worker-profile max_turns
  raised 90 -> 1000 fleet-side (infra, not repo; the raise gives workers
  more headroom, but the PROTOCOL closeout rule still applies — budget
  the last ~15% of turns so a capped worker leaves legible state).
- 2026-08-20 — local-hermes + fable — D2a DONE: 6-reaction astro rate
  campaign on the 4090. Verdict: GO gas-phase / NO-GO surface-LH rates
  (4–5 orders low without explicit-surface treatment). First non-Claude
  runtime to work the board end-to-end; worker hit iteration ceiling at
  closeout, fable finished bookkeeping. results/D2a-rate-reproduction.md.
- 2026-08-18 — hermes-workstation — A7-kinetics-database done: schema v0,
  36 P&K silicate records / 74 mechanisms, stdlib validator, and a
  49-source 2004–present expansion inventory under `kinetics-db/`.

- 2026-08-19 — hermes cloud worker — B1-schema-rfc done:
  `petra/docs/RFC-001-DECK-V2.md` written (deck schema v2 + UpdateStrategy
  trait + v1 compat shim + determinism/RNG contract + B2/B3 test plan;
  corrosion + ice-mantle showcase fragments). B2 unblocked pending Victor's
  review ack of the RFC.
- 2026-08-18 — fable — Card A7-kinetics-database added (machine:any,
  macbot-shaped): schema + extraction groundwork for a machine-readable
  successor to Palandri & Kharaka 2004. Pointer task queued in
  mission-control for laptop workers.
- 2026-08-18 — fable — Program tracking system created (PLAN, PROTOCOL,
  STATUS, 10 seed cards, inbox). Board opens with B1-schema-rfc as the
  keystone READY card.
- 2026-08-18 — fable — Omnibus program plan + three scoping docs merged
  (PR #20); root README refreshed (PR #21).
- 2026-08-18 — fable — Phase 1 first barrier: si-neutral ΔG‡(298) =
  113.0 kJ/mol (27.0 kcal/mol) vs X&L ~29, full pipeline (scan → product
  construction → CI-NEB → Sella → gates → quasi-RRHO), results +
  provenance archived. PRs #16–#18. TASK-164 queued for al-neutral.
- 2026-08-18 — fleet worker — qm/tests/test_phase1_driver.py landed
  (driver test coverage, via mission-control drain).
- 2026-08-18 — fable — Phase 0 closed: quarry package (PR #12), smoke
  v2 (PR #14); measured 4090: 74× analytic Hessians vs 32-thread CPU,
  identical energies. GPU owns all frequency work.
- 2026-08-19 — fable — qm/HANDOFF.md rewritten as the Phase-2 edition
  (barrier-ladder build plan: crystallographic cluster builder from the
  kaolinite deck cell, per-cell campaign driver, by_count emission,
  CALC-002..005 closure); new READY card A3-barrier-ladder.
- 2026-08-19 — fable — Track E opened: Ar-in-muscovite scoping doc
  landed (docs/scoping/ar-muscovite.md — the 1998 Sletten & Onstott
  circle, Nteme 2022–24 barriers, the uncomputed middle, Villa
  arbitration). New cards: E1-muscovite-deck (READY, published numbers
  only) and B5-execution-schedule (blocked: B2 — piecewise-isothermal
  T(t), the one engine feature the flagship needs).
- 2026-08-20 — Fable — Autopilot engaged: board_feeder.py (mission-
  control, hourly) auto-enqueues READY/unblocked cards as queue pointer
  tasks; omnibus-supervisor Hermes cron (daily 07:50) does bookkeeping
  + stale claims + Telegram digest. PROTOCOL.md amended. Humans and
  Fable are now exception handlers, not schedulers.
- 2026-08-20 — fable — card QI2-gpu-lease (READY, P1): explicit
  single-lane GPU lease + 16 GB cupy cap + shared worker preflight,
  after reviewing overnight A3-pilot/ollama cohabitation (22 ollama
  loads, all-GPU, zero errors — but pilot peaked 21.7 GB; luck, not
  design) and TASK-168's per-tick hand-rolled preflights.
- 2026-08-20 — fable — PROTOCOL.md merge-doctrine paragraph aligned to
  POLICY v7 (open PR and STOP; watcher owns lifecycle; human-merge
  label) — the stale auto-merge/Sourcery wording was flagged by the
  cloud Hermes worker during its enrollment (it correctly refused to
  follow it). Companion fixes landed in `vsletten/mission-control`
  commit `a2edf91`: `factory/hermes/cron/queue-drain.md` now points to
  POLICY v7, and service activation moved to
  `factory/services/activate.sh` for cloud-safe installation.
- 2026-08-20 — fable — Victor recovered his complete 1999 dissertation
  archive from the NAS (thesis LaTeX, the original C-language model,
  88 timestamped MC campaign runs, 86 Gaussian B3LYP logs incl.
  oxalate/ligand systems). New card A8-thesis-archive-intake (READY,
  workstation): catalog, curate golden data into legacy/, build the
  1999-vs-2026 DFT barrier ledger.
