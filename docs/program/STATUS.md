# Program status log

*Newest first. One line per merged unit of work (or notable event),
appended by the PR that did it: `date — actor — what + pointer`.
*Deviations get a `DEVIATION:` note; details live in the card.*

- 2026-08-23 — fable — A2 unblocked: ORCA descoped from the platform
  (Victor's call, licensing). Calibration layer is now Psi4 1.11
  DLPNO-CCSD(T) + ByteQC canonical GPU CCSD(T), cross-validating; both
  installed and smoke-verified on the workstation (water dimer/cc-pVTZ:
  canonical CCSD(T) agrees across engines to 1e-8 Eh; DLPNO −0.53
  kJ/mol off canonical). Runbook `qm/CALIBRATION.md`; SURVEY §2.3/§6.4
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
