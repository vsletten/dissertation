---
name: quarry-campaign
description: Run a quarry QM campaign (barrier calculations, smoke tests) on the workstation — compute etiquette, logging, checkpoint/resume, failure signatures. Use before launching any qm/scripts/ campaign or long PySCF job.
---

# Running quarry campaigns

Canonical driver pattern: `qm/scripts/phase1_xiao_lasaga.py`. Design
authority `qm/SURVEY.md`; traps in `.claude/learnings/`.

## Non-negotiable run rules (history: a lost 88-minute run + a wedged box)

1. **Tee everything**: `... 2>&1 | tee <logfile>` — nothing may depend
   on terminal scrollback surviving.
2. **Thread cap**: `OMP_NUM_THREADS=16` (the box has 32; leave half).
3. **GPU-first**: pass `--gpu` on the workstation; frequency/Hessian
   jobs are 74× there. Only tiny test-gate molecules belong on CPU.
4. Long jobs run backgrounded/detached with a completion signal —
   never hold a terminal open as the only record.
5. No multi-hour full-CPU jobs without Victor's explicit go.
6. **Clear inherited Python state**: prefix project commands with
   `env -u PYTHONPATH -u VIRTUAL_ENV`; an outer Hermes venv otherwise imports
   the wrong NumPy/PySCF stack.
7. **Own the GPU before Hessians**: stop both email extraction services and
   Honcho's `honcho-deriver-1` / `honcho-api-1` containers, then unload
   `gemma4-deriver` and `bge-m3`. Honcho will immediately reload ~19 GB after
   `ollama stop` if its callers remain alive. Verify `ollama ps` is empty and
   establish a bounded dead-man restart before launching. Restart all four
   workloads at exfil.

```bash
cd qm && env -u PYTHONPATH -u VIRTUAL_ENV OMP_NUM_THREADS=16 \
  uv run --frozen --extra gpu --extra dev python scripts/phase1_xiao_lasaga.py \
  --reaction si-neutral --gpu 2>&1 | tee /path/to/run.log
```

## Checkpoint/resume semantics

Stage outputs are XYZ files in `runs/phase1/<reaction>-<xc>-<basis>-<approach>/`
(gitignored). Existing file = stage skipped on rerun. **If you change
scan/guess logic, delete the downstream checkpoints** (`ts_guess.xyz`,
`ts.xyz`, `product.xyz`) or you'll resume into stale geometry. Run dirs
are keyed by approach mode so geometries never mix.

## Failure signatures (all seen live; fixes in git history)

| Signature | Meaning | Response |
|---|---|---|
| SCF not converged at stage 1 | strained hand-built guess | pre-opt stage exists; check geometry sanity first |
| Scan monotonic to floor (`ScanNoMaximumError`) | driven coordinate doesn't cross the ridge | mechanism problem, not numerics — needs second coordinate / different attack geometry |
| Saddle "verified" but ΔG‡ implausibly small, soft imaginary mode | trivial rearrangement saddle | the structural gates abort these; if you see one pass, tighten gates, don't report it |
| `r(Si-Ow)` grew ≫ guess after Sella | escaped the channel | use CI-NEB guess (already in driver) |
| Free relax of "product" returns to reactant | intermediate isn't a minimum | construct the product deliberately (bond-breaking scan with new bonds held) |

**A verified saddle is not necessarily the right saddle.** Chemistry
gates (barrier magnitude vs literature, TS geometry) are part of
acceptance, always.

## Results discipline

`results.json` + `store.sqlite` per run dir; every number carries
method + geometry hash + date. Barriers cross tool boundaries — never
pre-exponentiated rates (petra's gas constant differs from CODATA by
0.17%; see learnings).
