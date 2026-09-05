# A2-production-energetics — production-tier energies + open-source CC calibration

- status: blocked
- track: A (geochemistry)
- priority: P1
- machine: workstation
- depends: A2b-al-neutral-production-energetics, A2c-al-acid-production-energetics, A2d-oss-neutral-n4-production-energetics
- blocked-on: A2b, A2c, A2d
- claimed-by: hermes-custom-build-001


## Objective
Upgrade banked barriers to the survey protocol: r2SCAN-3c geometries
(needs pyscf-dispersion — already in the venv), wB97M-V/def2-TZVPD + SMD
single points on all stationary points, D4 dispersion; coupled-cluster
spot calibration on si-neutral (the "which functional do we trust"
half-day, SURVEY §6.4) via the **open-source calibration layer** —
Psi4 1.11 DLPNO-CCSD(T) cross-validated against ByteQC canonical GPU
CCSD(T), per `qm/CALIBRATION.md` (ORCA descoped 2026-08-23 on licensing
grounds). Full Gonzalez-Schlegel IRC on the soft-mode saddles replaces
quick-IRC. Report the barrier shifts SVP→production per reaction.

## Context
qm/SURVEY.md §6 (the tiered protocol) as amended 2026-08-23;
**qm/CALIBRATION.md — the calibration runbook: installed environments
(conda `calib-psi4`, `~/venvs/byteqc` + `~/opt/byteqc`), the focal-point
CBS recipe, and the DLPNO-vs-canonical acceptance gate**; quarry
pipeline supports solvent="smd" already. Banked stationary points exist for
si-neutral (A2a), al-neutral (TASK-168/PR #29), and al-acid one-water (A1b).
The `oss-neutral-n4` pilot result is documented by A3, but its ignored
stationary-point directory is no longer present and no hash-pinned external
archive was found on 2026-09-05; A2d must reconstruct it from the committed
builder and mechanism contract before re-tiering. Large single points that OOM
the 4090 may use rented H200/B200 capacity only through a separately authorized
spend gate (not B300 — FP64-gutted).

## Acceptance
- Production ΔG‡ table for every target named above, including the reconstructed
  embedded pilot, with method provenance in each run's `store.sqlite`.
- si-neutral calibration: CCSD(T)/CBS focal-point barrier from ByteQC
  (canonical) AND Psi4 (DLPNO, TightPNO), with the DLPNO−canonical
  delta reported (gate: ≤ 2 kJ/mol on the barrier); functional ranking
  stated (B3LYP-D4 vs r2SCAN-3c vs wB97M-V for Si–O chemistry).
- IRC-verified saddles (full Gonzalez-Schlegel) for the soft-mode TSs.
- Emitted petra fragments updated where re-tiered numbers change decks;
  CALCULATIONS.md rows updated to the new tier.

## Progress
- 2026-09-05 12:25 PDT (hermes-custom-build-001; profile=workstation) —
  Receipt-level inventory found that the completed A2a/si-neutral route is the
  only banked reaction satisfying A2's exact production, full-path, and
  provenance-store gates. The remaining targets are not a settings-only batch:
  TASK-168 al-neutral requires a barrierless-addition/one-saddle route contract;
  A1b al-acid requires two separately typed full IRCs and strict physical-proton
  ownership; the 63-atom embedded pilot's ignored stationary points are absent,
  so it requires deterministic reconstruction plus frozen-shell/PHVA semantics
  and may exceed the 4090 envelope. Launching them
  under the si-neutral-only CLI would either reject valid source roles or create
  scientifically false receipts. A2 is therefore decomposed into sequential,
  independently bounded cards A2b/A2c/A2d, while retaining this parent as the
  final table, functional-ranking, Petra-fragment, and CALCULATIONS closeout
  gate. No heavy calculation or shared-service mutation ran in this inventory
  slice.
- 2026-09-05 08:07 PDT (hermes-custom-build-001; profile=workstation) — **UNBLOCKED by A2a.** The exact neutral reactant now reaches the typed hydrolyzed product through a verified addition TS and a complete-grid barrierless cleavage shelf. Production wB97M-V ΔG‡ is `142.806783 kJ/mol`; the directive-adjusted canonical-TZ + TightPNO-CBS focal barrier is `132.960133 kJ/mol`, and the independent TZ DLPNO−canonical gate is `+0.245700 kJ/mol` (pass). The route, DFT values, provenance store, and focal receipt hashes are closed in the A2a card. A2 may resume the remaining banked-reaction production table and functional ranking without rebuilding si-neutral or launching canonical TS/QZ.
- 2026-08-24 22:53 PDT — bounded production verdict: **BLOCKED on an invalid
  banked si-neutral saddle, not on compute.** Exact r2SCAN-3c D4+gCP minima
  converged and central-difference composite Hessians give reactant/product/TS
  indices `0/0/1`; the TS imaginary mode is `79.185 cm^-1`. After repairing the
  Sella-2.5/ASE-3.29 zero-step compatibility defect, a real 260-step bilateral
  Gonzalez--Schlegel IRC terminates in associative basin `(True,True,True)` on
  both sides instead of reactant `(False,True,False)` and product
  `(True,False,True)`. An endpoint-mode directed Cartesian Sella retry changed
  the saddle only 0.0124 A RMS and reconverged the same basin. No wB97M-V,
  B3LYP-D4, or coupled-cluster barrier was computed on that wrong path.
  Executable READY card `A2a-si-neutral-production-path-rebuild` owns a persisted
  path/CI-NEB rebuild before A2 resumes. Evidence root has 47 files / 791,566
  bytes; manifest SHA-256 `efdcfc48c533dbadc1a30023f3a539f607711df0cd5cd354ff50fca3ba3556d1`
  and rehashes with zero mismatch. Both calibration adapters were independently
  smoke-run on water/cc-pVTZ: Psi4 1.11 TightPNO DLPNO-CCSD(T)
  `-76.3323355604 Eh`; ByteQC canonical CCSD(T) `-76.3325745245 Eh`.
- 2026-08-24 21:55 PDT — first production execution exposed and fixed two
  false-surface classes before any number could escape. The initial driver
  attached r2SCAN-3c D4+gCP to energies but PySCF geomeTRIC called an unpatched
  gradient factory; `078f060` now binds optimization to the full composite
  force. A force-converged r2SCAN-3c reactant/product then produced impossible
  11/7-mode indices through GPU4PySCF's analytic meta-GGA Hessian. Those receipts
  are preserved under the task207 evidence root, not relabeled. `f7c778a` adds
  a central-difference Hessian over the exact full composite gradient, tested
  against analytic HF/STO-3G frequencies; the resume-safe finite-difference
  si-neutral run is active from the converged v3 geometries. All implementation
  commits through `f7c778abf74995886aa69023cc606e21e83884fc` are pushed.
- 2026-08-24 20:58 PDT — claimed atomically by `hermes-workstation`; pushed
  `agents/A2-production-energetics` from exact `origin/main`
  `4d715c7597edd04ed1a0d180f02d5c600eb74128` and created the congruent
  dedicated worktree. Read-only bank inventory found four current production
  targets: free-dimer si-neutral, al-neutral, and al-acid plus the 66-atom
  embedded oss-neutral pilot. Their source stores/results remain byte-preserved
  in the historical worktrees/external archives. Environment probes verified
  PySCF 2.14.0 + GPU4PySCF 1.8.1, native def2-mTZVPP, `pyscf-dispersion`
  r2SCAN-3c D4+gCP primitives, Psi4 1.11, and ByteQC/PySCF 2.5.0. The first
  implementation slice is a tested, resume-safe production/calibration driver;
  no heavy calculation or shared-service mutation has started.
- 2026-08-18 — card created (fable).
- 2026-08-23 — fable — **descoped from ORCA** (Victor's call: licensing).
  Calibration layer replaced by Psi4 1.11 DLPNO-CCSD(T) + ByteQC
  canonical GPU CCSD(T), cross-validating. Both environments installed
  and smoke-verified on the workstation (water dimer/cc-pVTZ, frozen
  core): Psi4 canonical vs ByteQC canonical CCSD(T) agree to 1e-8 Eh
  (CPU vs GPU, independent codes); Psi4 DLPNO −0.53 kJ/mol off
  canonical. Runbook: qm/CALIBRATION.md. Card unblocked: status → ready.
