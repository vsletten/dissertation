# A2e-skala-functional-spike — Skala-1.1 single points vs the banked si-neutral CC focal barrier

- status: ready
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU single points; expected minutes, capped at the §12 4 h)
- depends: A2a-si-neutral-production-path-rebuild ✅ (banked structures + CC receipt)
- blocked-on: —
- claimed-by:

## Objective

Bounded survey-tier functional test, authorized by Victor on 2026-09-17
(motivation: `docs/scoping/msr-ai4science-scan.md` §2). Evaluate Skala-1.1
— Microsoft Research's neural exchange–correlation functional (MIT, native
PySCF/GPU4PySCF) — as single points on the four hash-bound A2a si-neutral
structures, and report the gas-phase electronic reactant → addition-TS
barrier against the banked canonical-TZ + TightPNO-CBS focal-point barrier
**132.960133 kJ/mol**. Deliver a receipt-bound comparison table and a
pre-declared verdict on whether Skala is a survey-tier single-point
functional candidate for the A-track. No geometry optimization, no
frequencies, no coupled cluster, no solvent-tier claims, no "production"
language.

## Context

- Design note: `docs/scoping/msr-ai4science-scan.md` (§2.2 verified
  environment fit, §2.3 known gaps, §2.4 anchor table, §2.5 this spike).
- Anchor: A2a card (`cards/A2a-si-neutral-production-path-rebuild.md`) and
  evidence root
  `/mnt/data/vsletten/dissertation-data/task208-a2a-path-rebuild-20260825/production-closeout/`.
  `store.sqlite` SHA-256
  `480cc244cf06bf7fd8edce86e4d57c3a8466e24e2793daf39be3fe27207dcd89`;
  `dft-summary.json` `b6f074c7351ad22b31b0340426e334b5777f04e081f5c8ceef61c0749b1959f8`;
  `cc-calibration/focal-point-summary.json`
  `2df64b9f1987f8e8924199bb51005156031c38e580bbcd958244a3956e340ff1`
  (all three rehashed equal on 2026-09-17).
- Structures (store table `structures`, all H8O8Si2, charge 0, spin 0):
  - id 1 `si-neutral-reactant` — geometry_hash
    `f976b1f07cfe0db48a593854fc5d9dd794b24a243d005c5c79c5368f9c635d40`
  - id 2 `si-neutral-intermediate` —
    `76ed3d0652cc93e6b024f7649102e7b53711d8054a0f0e6a3296b5540097df1e`
  - id 3 `si-neutral-addition-transition-state` —
    `b94fce89da1f03976e877447f154f527e7da73bd19659b7fc96a70733a849e92`
  - id 4 `si-neutral-released-product` —
    `dbbf58a07a595c0adfb97a3a36d29ca925b83abcf8b3a388dc5b427d3f9708e2`
- Gas-phase comparators on the same geometries (from `dft-summary.json` /
  `focal-point-summary.json`, kJ/mol, reactant → addition-TS): focal point
  132.960133; canonical CCSD(T)/cc-pVTZ 128.085164; TightPNO DLPNO/cc-pVTZ
  128.330864; r²SCAN-3c 125.742278; B3LYP/def2-SVP/DF 114.173360. The banked
  ωB97M-V (134.504248) and B3LYP-D4 (148.293616) rows are SMD(water) and are
  **not** gas-comparable — hence execution step 3d.
- Skala: https://github.com/microsoft/skala (README: install, `SkalaKS`,
  support matrix PySCF 2.14 / GPU4PySCF 1.8.1 / CUDA 12–13 / Python
  3.11–3.13). Environment fit verified 2026-09-17 by
  `uv pip install --dry-run skala-cuda12x` against quarry's exact pins
  (`pyscf==2.14.0`, `gpu4pyscf-cuda12x==1.8.1`, `cupy-cuda12x==14.1.1`):
  resolves with no conflicts, adding `skala==2026.9`,
  `skala-cuda12x==2026.9`, `torch==2.14.0`. torch is the only material new
  dependency (~3 GB; the qm venv has no torch today; `/mnt/data` has >1 TB
  free).
- API: `from skala.gpu4pyscf import SkalaKS; ks = SkalaKS(mol, xc="skala-1.1");
  ks.kernel()`; CPU path `from skala.pyscf import SkalaKS`.
- Run discipline: `.claude/skills/quarry-campaign/SKILL.md`, the QI2 GPU
  lease, `qm/quarry/etiquette.py`, `qm/quarry/store.py` (`Store` only for
  the new writable quarry store, same schema; never against A2a), fleet
  `POLICY.md` §12.

## Execution

1. **Environment.** Add an optional extra `skala = ["skala-cuda12x"]` to
   `qm/pyproject.toml`, run `uv lock`, and verify
   `env -u PYTHONPATH -u VIRTUAL_ENV uv run --frozen --extra gpu --extra skala
   python -c "import skala, torch, gpu4pyscf"` from `qm/`. Commit the lock.
   Do not change base dependencies; do not add torch to `dev`.
2. **Provenance in.** Rehash the A2a `store.sqlite` and refuse to proceed if
   it differs from the SHA-256 above. Open that file only as a SQLite URI
   with `mode=ro&immutable=1` (same pattern as
   `qm/quarry/calc005_store.py`:
   `sqlite3.connect(f"{path.as_uri()}?mode=ro&immutable=1", uri=True)`).
   Do **not** construct `Store(a2a_path)` from `qm/quarry/store.py` against
   A2a: `Store.__init__` opens a writable connection and runs the schema
   script, which can alter the evidence store before the hash-after check.
   Load structures 1–4 by `SELECT` on that read-only connection and verify
   each `geometry_hash` matches the list above. Write jobs/results only to
   the **new** quarry store from step 3.
3. **Single points** (GPU-first; each structure × each method is one job in
   a **new** quarry store under a new run dir on the data volume; record
   energy in Hartree, SCF converged flag, cycle count, wall time, device,
   package versions, and geometry hash per job):
   - a. `skala-1.1 / def2-tzvp`, gas phase, no dispersion.
   - b. `skala-1.1 / def2-tzvpd`, gas phase, no dispersion.
   - c. row (a) + D4 dispersion via `pyscf-dispersion` — a separately
     labelled row; Skala's docs are silent on dispersion and its training
     includes noncovalent sets, so do not fold D4 in silently and do not pick
     a winner.
   - d. `wb97m-v / def2-tzvpd`, **gas phase (no SMD)** — the
     solvent-consistent comparator the ledger lacks.
   Use `SkalaKS` default grids; no tuning. Run Skala without density fitting
   first (its docs do not mention DF); if memory forces DF, record it. If
   `SkalaKS` refuses any composition, record the failure receipt and skip —
   do not patch the package.
4. **Derive** ΔE‡ = E(TS) − E(reactant), ΔE(intermediate), ΔE(product) for
   every method row, and tabulate against the gas-phase comparators above.
   Report the SMD rows for context, flagged not comparable.
5. **Verdict** (pre-declared; gas-phase Skala rows vs the 132.960133 focal
   point): |Δ| ≤ 4.2 kJ/mol → *adopt as a survey-tier single-point
   functional candidate*; 4.2 < |Δ| ≤ 8.4 → *usable with stated
   uncertainty*; |Δ| > 8.4 → *reject for Si–O–Si hydrolysis*. State that one
   reaction is one data point, not a validation.
6. **Closeout.** Results doc `docs/program/results/A2e-skala-functional-spike.md`
   with hash-bound receipts (`results.json`, new `store.sqlite`, teed log);
   this card `done` with a `## Result`; STATUS line; PLAN row; add a Skala
   column/note to A2's Context (A2 stays blocked on Victor — this card
   neither closes nor unblocks it); one sentence in `qm/SURVEY.md` §6.4
   recording the outcome.

## Acceptance

- `skala` extra locked and importable in the qm venv; fast test suite +
  ruff green. Any new driver (`qm/scripts/skala_spike.py` or equivalent)
  carries a CPU-only tiny-molecule test guarded by
  `pytest.importorskip("skala")` so CI without torch stays green.
- 4 structures × 4 method rows = 16 converged single points with receipts,
  or an honest per-job failure receipt where a row could not run.
- Comparison table and pre-declared verdict in the results doc; no
  geometry changed; A2a store rehashes equal before and after; no
  "production" wording anywhere.
- Total workstation wall ≤ 4 h (POLICY §12); a run longer than 30 min goes
  under a bounded transient unit with `RuntimeMaxSec`, per A2b1's
  convention. Expected: minutes.

## Constraints

- POLICY §12: survey-tier only; no CC, no re-optimization, no solvent-tier
  claims, ≤ 4 h, never "production".
- Etiquette per the quarry-campaign skill: hold the QI2 GPU lease, tee the
  log, `OMP_NUM_THREADS ≤ 16`, prefix with `env -u PYTHONPATH -u VIRTUAL_ENV`.
- Skala's documented limits: molecular/gas-phase training; no documented
  solvation, dispersion, UKS or PBC support; ≤ ~180 atoms at TZVP (we run
  18). Respect them — record, don't work around.
- Do not hand-enqueue a mission-control pointer; the board feeder
  dispatches READY cards.
- A2a `store.sqlite` is immutable evidence: open it with
  `mode=ro&immutable=1`. `Store(path)` is writable and is only for the new
  spike store.

## Progress

- 2026-09-17 23:15 PDT — fable — card created from the MSR AI-for-Science
  scan (`docs/scoping/msr-ai4science-scan.md`); Victor authorized the spike
  2026-09-17. Environment resolution verified by dry-run only; A2a receipts
  rehashed equal; no compute run.
