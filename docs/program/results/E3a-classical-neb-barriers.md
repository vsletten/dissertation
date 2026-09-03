# E3a classical NEB barriers — replication gate and input-contract stop

**State:** `incomplete-replication-hamiltonian`
**Run date:** 2026-09-02
**Operator:** `(hermes-custom-build-001; profile=workstation)`

## Verdict

Most of the Nteme et al. (2022) classical-potential setup is recoverable from
the open paper and its supplementary workbook. A real eight-image LAMMPS
climbing-image NEB **near-match** for published divacancy route 1 returned a
forward barrier of **68.956 kcal/mol**, versus **67.644 kcal/mol** for that
exact route in the workbook: **+1.312 kcal/mol (+1.94%)**. This passes the
numerical tolerance below, but it does **not** pass the scientific replication
gate: the source says it introduced remote Al-for-Si charge-compensation
substitutions for every defective route, while the workbook does not identify
those sites. This run therefore used a different, net -3 e Hamiltonian with
LAMMPS's Ewald neutralizing background. A close number is not permission to
call two Hamiltonians the same.

The requested new-barrier campaign is **not complete**. Neither the paper,
its supplement, nor the current E3a/scoping specification defines atomistic
endpoint structures for the dehydroxylate, extended-zone, or octahedral-trap
models. Producing numbers by inventing those structures during execution
would be false precision. The executable reconstruction is retained; the
campaign remains blocked pending an explicit atomistic input contract.

## Provisional numerical tolerance

The provisional tolerance is **±5 kcal/mol** for one route. This is justified by the
published route distributions rather than chosen after seeing the answer:
Nteme et al. report 69 ± 6 kcal/mol (1σ) across 432 Ar divacancy routes, and
the E2 model already uses ±5 kcal/mol sensitivity bands for uncomputed
proxies. The observed +1.312 kcal/mol difference is well inside both, but the
replication gate remains closed until the published charge-compensation
construction is recovered or a source-backed equivalent is specified.

The run reached a maximum per-atom force of
`5.08e-5 kcal mol^-1 Å^-1`; the LAMMPS whole-replica norm was
`0.0411 kcal mol^-1 Å^-1`, above the requested `0.01` norm after the bounded
13,000-step run. The barrier varied by only 0.0062 kcal/mol over the final 500
steps. It is therefore adequate for the ±5 kcal/mol replication gate but is
honestly labelled `incomplete-convergence` for any tighter numerical claim.

## Published-source reconstruction

Source: Nteme, J., Scaillet, S., Brault, P., & Tassan-Got, L. (2022),
“Atomistic simulations of 40Ar diffusion in muscovite,” *Geochimica et
Cosmochimica Acta* 331, 123–142,
[doi:10.1016/j.gca.2022.05.004](https://doi.org/10.1016/j.gca.2022.05.004).
The open manuscript was recovered from HAL record `insu-03668854`; the
Elsevier supplement is `1-s2.0-S0016703722002198-mmc1.zip`.

Reconstructed method:

- 2M1 muscovite, `6 × 3 × 1`, 1,512 atoms, periodic in all directions.
- Final 0.5 K / 0 kbar cell from the workbook:
  `a=30.900625 Å`, `b=26.747522 Å`, `c=19.549356 Å`,
  `α=89.862673°`, `β=95.576801°`, `γ=89.985803°`.
- ClayFF charges/Lennard-Jones parameters, White (1999) neutral-Ar
  parameters, Lorentz-Berthelot mixing, 10 Å real-space cutoff, Ewald
  `1e-4`, and harmonic O-H bonds.
- Eight images; `23 kcal mol^-1 Å^-2` spring; LAMMPS `quickmin`;
  climbing-image stage.
- Route: Ar at K-gallery site 3 hopping to vacancy site 4, with the second
  divacancy site at 88.

The paper states that the route defects were charge-compensated by substituting
Al for Si at tetrahedral sites several unit cells away from each migration path
(methods lines 291–294 in the open manuscript). The supplement gives neither
the selected sites nor per-route modified coordinates. The source workbook's
complete barrier catalogue was otherwise extracted without
rounding:

| Species / mechanism | n | mean ± sample σ (kcal/mol) | range |
|---|---:|---:|---:|
| Ar / monovacancy | 216 | 71.101 ± 5.433 | 58.749–84.401 |
| Ar / divacancy | 432 | 69.607 ± 6.488 | 51.612–88.448 |
| K / monovacancy | 216 | 78.607 ± 9.231 | 61.374–106.875 |
| K / divacancy | 216 | 66.367 ± 10.840 | 44.181–87.040 |

A pristine static-energy check returned `PE=-348835.31 kcal/mol`; adding the
classical 0.5 K kinetic expectation gives `-348833.06 kcal/mol`, only
1.68 kcal/mol (4.8 ppm) from the workbook's final reported total energy.
This independently checks the pristine force-field typing, bonds, cell, and
units. It does not validate the missing defective-route substitutions.

## Actual provisional NEB result

| Quantity | Value |
|---|---:|
| Published route-1 divacancy barrier | 67.644151 kcal/mol |
| Computed forward barrier | **68.955713 kcal/mol** |
| Computed reverse barrier | 69.402430 kcal/mol |
| Absolute / relative forward error | +1.311562 kcal/mol / +1.9389% |
| Images / spring | 8 / 23 kcal mol^-1 Å^-2 |
| Relax + climb bounds | 3,000 + 10,000 steps |
| Final max replica / atom force | 0.0411253 / 0.0000508 kcal mol^-1 Å^-1 |

Heavy evidence root:
`/mnt/data/vsletten/dissertation-data/e3a-classical-neb-barriers-20260902/`.
The route receipt is
`replication/divacancy-0001-converged/result-summary.json`; its screen-log
SHA-256 is
`155bb474fe40debaccc5097da3b735ecd2bc46e2875237973b7ee17e979480b7`.
The source workbook SHA-256 is
`292c28c2b182e45aa7fe6cc36a2e62779d7d8dbf0aa6e74525c0e898ee31237f`.

## Reproduction

```bash
ROOT=/mnt/data/vsletten/dissertation-data/e3a-classical-neb-barriers-20260902
uv run --with openpyxl python scripts/nteme_neb.py extract \
  --workbook "$ROOT/sources/nteme-2022-supplement/Results Of MD Simulations.xlsx" \
  --out "$ROOT/replication"

uv run --with openpyxl python scripts/nteme_neb.py prepare \
  --workbook "$ROOT/sources/nteme-2022-supplement/Results Of MD Simulations.xlsx" \
  --out "$ROOT/replication/divacancy-0001-converged" \
  --mechanism divacancy --route-number 1 \
  --relax-steps 3000 --climb-steps 10000 --ftol 0.01

cd "$ROOT/replication/divacancy-0001-converged"
OMP_NUM_THREADS=1 mpirun -np 8 lmp -partition 8x1 -in in.neb | tee screen.log
cd -
python3 scripts/nteme_neb.py summarize \
  --run-dir "$ROOT/replication/divacancy-0001-converged"
```

## Why the campaign stopped

The only atom-resolved supplement is the pristine 1,512-site structure and
Ar/K vacancy-hop catalogue. It contains no dehydroxylated structure, no
extended-defect coordinates, no stable octahedral-trap endpoint, and no Xe
model. The scoping document defines those as KMC state labels, not as atomic
structures. At minimum, the next executable contract must specify:

1. dehydroxylation stoichiometry, OH pair IDs, released species, residual-O
   charge/type, relaxation protocol, and collapsed-cell boundary condition;
2. an atomistic extended-zone geometry (strain/void dimensions, location,
   charge compensation, and endpoint sites);
3. a recoil-derived octahedral Ar trap minimum and a physically defined exit
   endpoint;
4. Xe Lennard-Jones parameters and validation source;
5. the per-route charge-compensation substitutions omitted from the published
   workbook.

Until those are explicit and source-backed, the E2 values remain labelled
proxies and no deck fragment is changed.
