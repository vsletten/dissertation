# QM lab session handoff — standing up the barrier pipeline

*Written 2026-08-18 for a fresh session. The design authority is
`qm/SURVEY.md` (read it in full — it is the survey AND the recommendation);
this handoff adds the build plan, house conventions, and what changed since
the survey was written (Petra exists now).*

## 0. Mission

Build the thin platform the survey's §8 says is "genuinely ours to build":
the orchestration layer that turns mineral crystallography + a reaction
taxonomy into **site-resolved activation energies and rate constants**,
computed on the recommended free stack, and emitted in the formats our KMC
consumes. Everything heavy (DFT engines, TS optimizers, thermochemistry) is
assembled, not written — see SURVEY.md §8's verdict and §2–§5 for exactly
which packages and why.

The stack, one line each (details + versions in SURVEY.md TL;DR table):
PySCF + GPU4PySCF (workhorse DFT, analytic Hessians, SMD), ORCA 6
(cross-check + DLPNO-CCSD(T) calibration; free but closed, registration,
**no redistribution** — never commit or containerize it), geomeTRIC + Sella
(TS/IRC), xtb + CREST (screening), quasi-RRHO + Eyring/tunneling (own the
~50-line module; GoodVibes/Arkane as audited references), ASE + plain
Python + SQLite as glue. No AiiDA/pyiron — one user, one box.

## 1. Hardware reality — read before running anything

SURVEY.md targets **Victor's workstation**: Ryzen 9 7950X (16C/32T), 64 GB
DDR5, RTX 4090 24 GB. Two consequences for a session:

1. **A cloud session container is not that machine.** Detect what you're on
   before making performance claims. All *code* (cluster builder, pipeline,
   parsers, rate module, tests) develops and tests fine anywhere on small
   molecules with CPU PySCF. The GPU smoke test and real benchmarks
   (SURVEY.md §9 Phase 0) only mean something on Victor's box — structure
   the work so he can run one script there and paste results back.
2. The 4090's FP64 is capped (~1.3 TFLOPS ≈ the CPU's peak): expect
   **2–5× (SCF/gradients) and 5–10× (Hessians)** over the 16 cores, not
   datacenter numbers. Frequency jobs dominate TST cost, so that is still
   the win that matters. Test the cuEST INT8-emulation path empirically
   (`nvidia-cuest-cu12`, already integrated in GPU4PySCF) — no published
   4090 numbers exist; ours would be genuinely informative. Details §3.

## 2. What changed since the survey: Petra

The survey predates Petra (`petra/` — read `petra/docs/DESIGN.md` §1–§4 for
orientation). Its §8.3 "bridge" targeted regenerating kmc-rs's `data.rxn`.
That is still worth doing (Phase 3 validation against the dissertation),
but **the primary output target is now a Petra deck fragment**: TOML
`[[reactions]]` entries with Arrhenius/Eyring parameters, plus `by_count`
ΔEa tables for connectivity dependence — i.e. the *non-flat environment
tables the legacy model never had*. The kaolinite deck's reaction taxonomy
(`petra/examples/kaolinite.toml`, states documented in its header) is the
exact shopping list of barriers: Si–O–Si, Si–O–Al (both sides), Al–OH–Al
hydrolysis ± protonation state, plus Si/Al attachment/detachment. Petra's
per-deck units (kcal/mol, kJ/mol, eV) mean no unit gymnastics at the
boundary — emit kJ/mol and say so.

## 3. Proposed shape (adjust freely, but start here)

New Python package in `qm/` — suggested name **`quarry`** (the place stone
comes from; petra is the rock). House rules from `.claude/CLAUDE.md`:
Python 3.13, **uv** for env/deps, ruff, tests for all new functions.

```
qm/
  SURVEY.md            (the design authority — do not fork its content)
  HANDOFF.md           (this file)
  pyproject.toml       (uv-managed; deps: pyscf, geometric, sella, ase,
                        goodvibes, numpy; gpu4pyscf-cuda12x + nvidia-cuest
                        as an optional [gpu] extra; xtb via conda-forge or
                        binary — document, don't vendor)
  quarry/
    clusters.py        cluster builder: terminated, constrained
                       aluminosilicate clusters per KMC site family
                       (100s–500s) and protonation state; peripheral-atom
                       constraints for lattice resistance (survey §6.2)
    pipeline.py        the tiered protocol as composable steps:
                       xtb screen → r²SCAN-3c (or B3LYP-D4/def2-TZVP)
                       opt+freq → geomeTRIC/Sella TS + IRC verify →
                       ωB97M-V/def2-TZVPD + SMD single points
    rates.py           quasi-RRHO ΔG‡ → Eyring k(T) with Wigner/Eckart
                       tunneling (the in-house ~50 lines, gated against
                       GoodVibes/Arkane numbers)
    emit.py            Petra deck fragments + legacy data.rxn, each number
                       with provenance (method, geometry hash, date)
    store.py           SQLite of structures/jobs/results (survey §4.3)
  scripts/
    phase0_smoke.py    env + engine sanity, CPU vs GPU timings table
                       (B3LYP/def2-TZVP E+grad+Hessian on Si(OH)₄·nH₂O) —
                       the script Victor runs on the 4090 box
  tests/               fast, CPU-only, tiny molecules (H₂O, Si(OH)₄ at
                       small basis): pipeline mechanics, rate math vs hand
                       Eyring values, emit round-trips into petra's parser
```

Testing philosophy mirrors petra's gates: the rate module is checked
against closed-form Eyring numbers to 1e-12; `emit.py` output must be
*compiled by petra-deck* in a test (add a dev workflow that runs
`cargo run -p petra-cli -- <emitted deck> --steps 0` or parse via a tiny
Rust test — cheapest is emitting into a template deck and invoking the
petra CLI, which exists in-repo).

## 4. Phase plan (from SURVEY.md §9, operationalized)

- **Phase 0 — bootstrap** (this session): package skeleton, uv env, CPU
  smoke tests green in-session; `scripts/phase0_smoke.py` ready for
  Victor's box (GPU optional-detect, prints the timing table). Deliver via
  PR (see §5).
- **Phase 1 — Xiao & Lasaga reproduction**: (HO)₃Si–O–Si(OH)₃ and
  (HO)₃Si–O–Al(OH)₃⁻ + H₂O/H₃O⁺ hydrolysis barriers; compare 1990s numbers
  (§7.1: ~24 kcal/mol acid, ~19 base) and modern values. In-session at
  modest basis; production numbers on the workstation. This validates every
  stack layer at once.
- **Phase 2 — the barrier ladder**: connectivity × protonation resolved
  barriers for the five KMC site families (survey §6.2/§7.3 for cluster
  norms and the pKa table that sets protonation states worth computing).
  ORCA DLPNO-CCSD(T) spot calibration (workstation only).
- **Phase 3 — close the loop**: emit the Petra kaolinite deck's rates
  (+ regenerate `data.rxn`), run both engines, compare against experimental
  Ea (kaolinite ~29 kJ/mol acidic, 33–51 alkaline; quartz 66–90 — §7.3).

## 5. House workflow (Victor's standing instructions)

- Branch `claude/<topic>-<suffix>`; every push ends in a **PR to `main`,
  merged** — enable auto-merge (branch protection wants 4 CodeQL checks
  that lag and read "expected"; don't hand-retry). After merge,
  **sync local `main`** (`git fetch origin main:main`). GitHub 500s
  intermittently; retry with backoff.
- ORCA: registration-gated, no redistribution — install instructions in a
  README, never the binary in git. Same for any gated ML weights (UMA).
  License flags table: SURVEY.md Appendix A.
- Big outputs (checkpoints, cube files) stay out of git; the SQLite store
  keeps derived numbers + provenance, structures as light XYZ.

## 6. Traps the survey already identified (don't rediscover)

- Psi4: no SMD + finite-difference DFT Hessians — third opinion only.
- MLIPs: fine for guesses/conformers; **DFT owns every number that enters
  KMC**; non-conservative variants have broken Hessians (§5).
- Pure continuum solvation is qualitatively wrong for water-mediated
  hydrolysis — cluster-continuum with 3–6 explicit waters (§6.3).
- Free clusters underestimate lattice resistance by up to tens of kJ/mol —
  constraint strategy is physics, not a detail (§6.2, Pelmenschikov).
- pysisyphus unmaintained (pin, don't build on); AIMNet2 has no Al;
  geomeTRIC's BSD has a no-AI-training-data clause (fine for our use).
