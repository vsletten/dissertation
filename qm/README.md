# quarry

The QM orchestration layer for the kaolin KMC models — the thin platform
[SURVEY.md](SURVEY.md) §8 says is "genuinely ours to build". It maps mineral
crystallography + a reaction taxonomy to **site-resolved activation energies
and rate constants** (transition state theory on terminated aluminosilicate
clusters), and emits them in the formats our KMC engines consume: petra deck
fragments (`[[reactions]]` TOML) and the legacy `data.rxn`.

Everything heavy is assembled, not written: PySCF (+ GPU4PySCF on the
workstation) for DFT, geomeTRIC/Sella for TS searches and IRC, xtb/CREST for
screening, Psi4 (DLPNO-CCSD(T)) + ByteQC (canonical GPU CCSD(T)) for
coupled-cluster calibration ([CALIBRATION.md](CALIBRATION.md); ORCA descoped
2026-08-23 on licensing grounds). Quarry owns
the cluster builder, the pipeline sequencing, the ~50-line quasi-RRHO/Eyring
rate module, the emitters, and the provenance store.

- Design authority: [SURVEY.md](SURVEY.md) (read in full)
- Build plan and house conventions: [HANDOFF.md](HANDOFF.md)

## Layout

| Module | Role |
|---|---|
| `quarry/clusters.py` | Terminated, constrained aluminosilicate clusters per KMC site family (100s–500s) and protonation state |
| `quarry/pipeline.py` | The tiered protocol as composable steps: xtb screen → DFT opt+freq → TS + IRC → high-level single points |
| `quarry/rates.py` | Quasi-RRHO ΔG‡ → Eyring k(T) with Wigner/Eckart tunneling (gated against closed-form values) |
| `quarry/emit.py` | Petra deck fragments + legacy `data.rxn`, every number with provenance |
| `quarry/store.py` | SQLite store of structures/jobs/results |
| `scripts/phase0_smoke.py` | Env + engine sanity and the CPU-vs-GPU timing table — run this on the 4090 box |

## Install

```bash
cd qm
uv sync --extra dev        # CPU stack + test tooling
uv run pytest              # fast, CPU-only gate (tiny molecules, small basis)
```

On the workstation (RTX 4090), add the GPU extra and run the smoke test:

```bash
uv sync --extra dev --extra gpu
uv run python scripts/phase0_smoke.py            # full timing table (DF on)
uv run python scripts/phase0_smoke.py --quick    # seconds-long sanity check
uv run python scripts/phase0_smoke.py --waters 24  # ~100-atom-regime probe
```

First measured table (2026-08-18, 15 atoms, **no** density fitting, cold
GPU): CPU/GPU energies bitwise-matched; Hessian 1.6× GPU, SCF/gradient
*slower* on GPU. That is the survey's predicted small-molecule worst case
plus untimed-warmup and missing-cuTENSOR artifacts — the defaults now warm
up the GPU, enable density fitting (the production path), and pull in
cuTENSOR, and bigger `--waters` probes the regime where the survey's
2–5×/5–10× claims live (details in `.claude/learnings/performance.md`).

### Not managed by uv (install separately, never commit)

- **Psi4** (calibration, CPU DLPNO-CCSD(T)) — conda env `calib-psi4`;
  **ByteQC** (calibration, GPU canonical CCSD(T)) — venv `~/venvs/byteqc`
  + source build at `~/opt/byteqc`. Both open source and rebuildable from
  public sources; full install runbook with the measured build gotchas in
  [CALIBRATION.md](CALIBRATION.md). (ORCA was descoped 2026-08-23.)
- **xtb / CREST** — `conda install -c conda-forge xtb crest`, or the static
  binaries from the Grimme-lab GitHub releases. Also optional-detected.
- Gated ML weights (e.g. UMA): same rule — document, don't vendor.

## Testing philosophy

Mirrors petra's gates:

- the rate module is checked against closed-form Eyring numbers to 1e-12;
- `emit.py` output must be *compiled by petra* in a test (`cargo run -p
  petra-cli -- <emitted deck> --steps 0`, marker `petra`, skipped when cargo
  is unavailable);
- everything in the default gate is CPU-only and tiny (H₂O, Si(OH)₄ at small
  basis) so it stays green in any session, not just on the workstation.

Big outputs (checkpoints, cube files) stay out of git; the SQLite store keeps
derived numbers + provenance, structures as light XYZ.
