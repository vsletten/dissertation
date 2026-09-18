# A2e — Skala-1.1 functional spike

## Decision

**Reject Skala-1.1 as a survey-tier single-point functional candidate for this
Si–O–Si hydrolysis route.** On the exact four A2a geometries, its gas-phase
reactant → addition-TS electronic barriers are 143.632106 kJ/mol with
`def2-TZVP` and 145.076794 kJ/mol with `def2-TZVPD`. Those differ from the
pre-declared 132.960133 kJ/mol focal-point comparator by 10.671973 and
12.116661 kJ/mol, respectively, beyond the 8.4 kJ/mol rejection boundary.

This is one reaction and one data point, not a general validation or rejection
of Skala. It does not authorize geometry optimization, frequencies, solvent
claims, or method selection for other A-track families.

## Receipt-bound comparison

All values are gas-phase electronic differences on unchanged A2a geometries.
`ΔE int` and `ΔE product` are relative to the reactant.

| Method row | ΔE‡ TS (kJ/mol) | ΔE int (kJ/mol) | ΔE product (kJ/mol) | |Δ from focal| (kJ/mol) | Pre-declared verdict |
|---|---:|---:|---:|---:|---|
| Skala-1.1 / def2-TZVP / no dispersion | 143.632106 | 121.108276 | 16.329596 | 10.671973 | reject for Si–O–Si hydrolysis |
| Skala-1.1 / def2-TZVPD / no dispersion | 145.076794 | 122.724602 | 17.533752 | 12.116661 | reject for Si–O–Si hydrolysis |
| Skala-1.1 / def2-TZVP + D4 sensitivity | 143.148973 | 120.530833 | 16.391493 | 10.188840 | reject for this route; sensitivity row only |
| ωB97M-V / def2-TZVPD / native VV10 | 125.214239 | 103.765396 | 17.242705 | 7.745894 | comparator only; Skala bands do not apply |

The D4 sensitivity row is deliberately separate. DFT-D4 provides no Skala
parameterization, so the receipt uses r²SCAN damping and is **not** a
parameterized Skala-D4 method. It does not select a winner; it only shows that
this sensitivity changes the barrier by about 0.48 kJ/mol and does not alter
the rejection decision.

Gas-phase banked comparators on the same geometries:

| Comparator | Barrier (kJ/mol) |
|---|---:|
| canonical-TZ + TightPNO-CBS focal point | 132.960133 |
| canonical CCSD(T)/cc-pVTZ | 128.085164 |
| TightPNO DLPNO-CCSD(T)/cc-pVTZ | 128.330864 |
| r²SCAN-3c | 125.742278 |
| B3LYP/def2-SVP/DF | 114.173360 |

The banked ωB97M-V and B3LYP-D4 SMD(water) barriers, 134.504248 and
148.293616 kJ/mol, are solvent-context rows and are not gas-comparable.

## Execution and integrity

- Run root:
  `/mnt/data/vsletten/dissertation-data/a2e-skala-functional-spike/20260918T001926-0700/`
- Device: `cuda:0` on `custom-build-001`; 16-thread caps; no density fitting;
  Skala default grids.
- Locked distributions: `skala==2026.9`, `skala-cuda12x==2026.9`,
  `pyscf==2.14.0`, `gpu4pyscf-cuda12x==1.8.1`, `torch==2.14.0`,
  `pyscf-dispersion==1.5.0`.
- Exact matrix: 4 structures × 4 rows = **16/16 converged single points**.
- Summed per-job wall: `1073.573485 s` (`17.892891 min`), below the 4 h
  survey envelope.
- A2a source store SHA-256 before and after:
  `480cc244cf06bf7fd8edce86e4d57c3a8466e24e2793daf39be3fe27207dcd89`.
  It was opened through an immutable, read-only SQLite URI; its bytes did not
  change.
- `results.json` SHA-256:
  `a058f1d702b7674c9e31b9c2bbe548447c7e05b731422d88590718df53e3736c`
- New `store.sqlite` SHA-256:
  `9f77e4a958fc72d62c1d06157767b7b5ae48aa77e58ffd23d0865ba3968dad30`
- Teed `run.log` SHA-256:
  `fc959c457ea63bf01fe424c87831514d86a0bc386813117f59de72b10f535479`
- `manifest.json` SHA-256:
  `58abcffff8b880f7de528151c3c35ef7b52adb448c331052d48a76736c07d053`
- The verifier reopened the new store read-only and reported 4 structures,
  16 jobs, 16 done jobs, 96 scalar results, byte-level receipt/result/store
  agreement, exact derived-row replay, and unchanged source-store bytes.

## Verification

From `qm/`:

```text
uv run --frozen --extra skala --extra dev ruff check scripts/skala_spike.py tests/test_skala_spike.py
All checks passed!

uv run --frozen --extra skala --extra dev pytest -q
935 passed, 468 warnings in 118.63s

uv run --frozen --extra skala python scripts/skala_spike.py ... verify ...
{"done_jobs": 16, "jobs": 16, "results": 96, "structures": 4, "verified": true, ...}
```

The CPU-only tiny-molecule regression executes `SkalaKS` without implicit
dispersion. Additional tests cover immutable-source loading, geometry-hash
refusal, exact 4 × 4 receipt derivation, honest failed-row typing, source hash
drift, and fail-closed disagreement among `results.json`, each immutable job
receipt, and the SQLite rows. One cold adversarial review reproduced every
number and verdict, then exposed that last verifier gap; the regression and
store-level cross-check close it.
