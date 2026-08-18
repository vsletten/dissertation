# quarry

The QM orchestration layer for the kaolin KMC models — the thin platform
[SURVEY.md](SURVEY.md) §8 says is "genuinely ours to build". It maps mineral
crystallography + a reaction taxonomy to **site-resolved activation energies
and rate constants** (transition state theory on terminated aluminosilicate
clusters), and emits them in the formats our KMC engines consume: petra deck
fragments (`[[reactions]]` TOML) and the legacy `data.rxn`.

Everything heavy is assembled, not written: PySCF (+ GPU4PySCF on the
workstation) for DFT, geomeTRIC/Sella for TS searches and IRC, xtb/CREST for
screening, ORCA for cross-checks and DLPNO-CCSD(T) calibration. Quarry owns
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
uv run python scripts/phase0_smoke.py          # full timing table
uv run python scripts/phase0_smoke.py --quick  # seconds-long sanity check
```

### Not managed by uv (install separately, never commit)

- **ORCA 6** — registration-gated download from https://www.faccts.de/orca/.
  Free for academic/personal use, closed source, **no redistribution**: never
  commit the binary, never bake it into a container. Install anywhere on
  `PATH`; quarry treats it as optional and detects it at runtime.
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
