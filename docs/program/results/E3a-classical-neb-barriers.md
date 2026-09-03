# E3a classical NEB barriers — bounded sensitivity campaign

**State:** `done-bounded-sensitivity`; all reported NEBs are
`incomplete-convergence`; octahedral trap is `incomplete-index-gate`
**Run date:** 2026-09-03
**Operator:** `(hermes-custom-build-001; profile=workstation)`

## Verdict

The source-constrained atomistic contract is executable. Seven model members were
relaxed under the reconstructed neutral ClayFF Hamiltonian. Six eight-image
LAMMPS climbing-image NEBs were run for 3,000 relaxation plus 10,000 climbing
steps; the octahedral-trap route was correctly withheld by its minimum/index
gate.

The neutral reconstructed replication returns **68.414811 kcal/mol** for Nteme
et al. (2022) divacancy route 1 versus **67.644151 kcal/mol** published:
**+0.770660 kcal/mol (+1.1393%)**, inside the predeclared ±5 kcal/mol tolerance.
It validates this source-constrained sensitivity Hamiltonian, not the unpublished
production structure: Nteme et al. did not disclose their route-specific remote
charge-compensation sites.

The campaign's principal bounded results are:

- post-reaction local-dehydroxylate Ar hop: **64.095991 kcal/mol**;
- extended-zone Ar hops across three geometries: **68.410690–68.411538
  kcal/mol**;
- Xe divacancy hop under the Bourg–Sposito/ClayFF screen: **90.744700
  kcal/mol**;
- octahedral-trap escape: **no barrier emitted**. The Ar seed displaced 2.316 Å
  out of its assigned cage during relaxation, and no full-system index-zero
  Hessian was established. The result remains `incomplete-index-gate`.
- delamination-adjacent gallery/surface release: **no barrier emitted**. The
  merged atomistic contract defines no delaminated slab or atom-resolved surface
  release endpoint, so this remains `incomplete-input-contract` rather than an
  invented path.

None of the six NEBs reached the requested whole-replica force norm of
0.01 kcal mol⁻¹ Å⁻¹. Final norms are 0.0272–0.0633, although maximum per-atom
forces are 0.00011–0.00022 and every barrier varies by at most 0.0101 kcal/mol
over the final 500 steps. The values below are therefore stable **bounded
sensitivity estimates**, not tightly converged barriers or confidence intervals.
No E2 proxy is silently overwritten.

## Replication gate

The numerical tolerance was declared before the new campaign values: **±5
kcal/mol** for one route. Nteme et al. report 69 ± 6 kcal/mol (1σ) across 432 Ar
divacancy routes, so ±5 is narrower than that route-family dispersion.

| Quantity | Value |
|---|---:|
| Published route-1 divacancy barrier | 67.644151 kcal/mol |
| Neutral reconstructed forward barrier | **68.414811 kcal/mol** |
| Neutral reconstructed reverse barrier | 68.790916 kcal/mol |
| Forward delta | +0.770660 kcal/mol (+1.1393%) |
| Images / spring | 8 / 23 kcal mol⁻¹ Å⁻² |
| Relax + climb bound | 3,000 + 10,000 steps |
| Final replica / atom force | 0.0500754 / 0.0001137 kcal mol⁻¹ Å⁻¹ |
| Final-500-step barrier span | 0.009419 kcal/mol |
| Numerical tolerance gate | PASS |
| Exact unpublished-input replication claim | NOT MADE |

The earlier net −3 e reconstruction produced 68.955713 kcal/mol but used an
Ewald neutralizing background. It is retained as historical diagnostic evidence,
not the campaign baseline. The present neutral model applies the contract's
charge-correct inverse operator: three remote Al2→Si substitutions, each paired
with four O3→O2 changes. Its provenance remains
`reconstructed-charge-compensation` because the paper does not identify the
production atom IDs.

## Barrier matrix and E2 proxy disposition

| Model / species | Geometry or contract | E2 proxy / source anchor | Forward barrier (kcal/mol) | Final replica force | Final-500 span | Disposition |
|---|---|---:|---:|---:|---:|---|
| reconstructed gallery / Ar | neutral route-1 divacancy | 67.644151 published | **68.414811** | 0.050075 | 0.009419 | numerical replication gate PASS; `incomplete-convergence` |
| local dehydroxylate / Ar | post-reaction topology; O73 residual, O76/H81/H84 removed | 57/60/69 published envelope | **64.095991** | 0.027248 | 0.010002 | deck-ready sensitivity replacement; −4.318820 vs reconstructed gallery; `incomplete-convergence` |
| extended zone small / Ar | a=6 Å, b=3 Å, opening=0.5 Å | 40 proxy | **68.411397** | 0.050048 | 0.008569 | proxy contradicted for this geometry; `incomplete-convergence` |
| extended zone middle / Ar | a=8 Å, b=5 Å, opening=1.0 Å | 40 proxy | **68.410690** | 0.049973 | 0.004555 | proxy contradicted for this geometry; `incomplete-convergence` |
| extended zone large / Ar | a=10 Å, b=7 Å, opening=1.5 Å | 40 proxy | **68.411538** | 0.050085 | 0.005221 | proxy contradicted for this geometry; `incomplete-convergence` |
| octahedral trap / Ar | vacant-M1 geometric seed → K675 gallery site | 64 proxy | — | — | — | retain proxy; `incomplete-index-gate` |
| delamination-adjacent / Ar | no atom-resolved endpoints in the merged contract | 58 ± 5 proxy | — | — | — | retain proxy; `incomplete-input-contract` |
| reconstructed gallery / Xe | Bourg–Sposito Xe LJ, Lorentz–Berthelot with ClayFF | no E2 Xe value | **90.744700** | 0.063274 | 0.009903 | screening only; `incomplete-convergence` |

The three extended-zone geometries span 0.000848 kcal/mol. Under this fixed-cell,
relaxed-endpoint family, opening the selected basal lens does **not** recover the
E2 deck's 40 kcal/mol proxy or its mapped 10⁶× enhancement: the middle member is
+28.410690 kcal/mol above that proxy and essentially equal to the reconstructed
gallery route. This rejects the proxy **for these three contract geometries**;
it does not reject all extended defects or measured fast pathways.

The local-dehydroxylate hop falls inside the 57/60/69 kcal/mol envelope but does
not model O–H bond breaking, H₂O formation, or the dehydroxylation reaction path.
It replaces only the post-reaction gallery-hop sensitivity value. Xe is
22.329889 kcal/mol above the reconstructed Ar route under the screening
parameterization; dry-muscovite Xe cross interactions remain uncalibrated.

## Deck-ready fragment

Use the following explicit phase-3 overrides only when the consuming deck can
retain provenance and typed convergence state:

```text
local_dehydroxylate_ar_hop:
  barrier_kcal_mol: 64.095991
  provenance: E3a/2026-09-03/reconstructed-charge-compensation
  state: incomplete-convergence

extended_zone_ar_hop_grid:
  - {radius_a_A: 6.0, radius_b_A: 3.0, opening_A: 0.5,
     barrier_kcal_mol: 68.411397, state: incomplete-convergence}
  - {radius_a_A: 8.0, radius_b_A: 5.0, opening_A: 1.0,
     barrier_kcal_mol: 68.410690, state: incomplete-convergence}
  - {radius_a_A: 10.0, radius_b_A: 7.0, opening_A: 1.5,
     barrier_kcal_mol: 68.411538, state: incomplete-convergence}

octahedral_trap_escape_ar:
  barrier_kcal_mol: 64.0
  provenance: E2-labelled-proxy-retained
  state: incomplete-index-gate

xenon_divacancy_screen:
  barrier_kcal_mol: 90.744700
  provenance: Bourg-Sposito-Xe/ClayFF-screen
  state: incomplete-convergence
```

The repository deck itself is not mutated by this PR: its existing schema cannot
carry the required `incomplete-*` state alongside a live barrier without making
the number look production-ready. E3b or a narrow deck-schema card must own that
promotion.

## Method and provenance

Source: Nteme, J., Scaillet, S., Brault, P., & Tassan-Got, L. (2022),
“Atomistic simulations of 40Ar diffusion in muscovite,” *Geochimica et
Cosmochimica Acta* 331, 123–142,
[doi:10.1016/j.gca.2022.05.004](https://doi.org/10.1016/j.gca.2022.05.004).
The source workbook SHA-256 is
`292c28c2b182e45aa7fe6cc36a2e62779d7d8dbf0aa6e74525c0e898ee31237f`.

Common campaign method:

- 2M1 muscovite, `6 × 3 × 1`, periodic fixed cell; source structure 1,512 atoms;
- ClayFF charges/Lennard-Jones, White (1999) neutral-Ar, explicit
  Lorentz–Berthelot unlike pairs, 10 Å real-space cutoff, Ewald `1e-4`, harmonic
  O–H bonds;
- endpoint FIRE minimization to requested max force 0.01 kcal mol⁻¹ Å⁻¹ with
  energy-tolerance early stopping disabled;
- eight images, 23 kcal mol⁻¹ Å⁻² spring, LAMMPS `quickmin`, 3,000 relaxation
  plus 10,000 climbing-image steps, one OpenMP thread per MPI rank;
- the trap candidate receives a 4.25 Å local-shell stage before full-cell
  release and is forbidden from NEB unless it stays in-cage and has an
  independently proven index-zero minimum.

No Hessian is inferred from force convergence. The trap result therefore remains
typed incomplete even though both endpoint force thresholds pass.

## Evidence and integrity

Canonical heavy evidence root:

`/mnt/data/vsletten/dissertation-data/e3a-classical-neb-campaign-grid-20260903/`

- 253 listed files (`manifest.json` is omitted from the listing by design);
- 296,414,715 bytes on disk including the unlisted `manifest.json`;
- every listed size and SHA-256 independently rechecked against disk;
- listed path set exactly matches every non-manifest file;
- `campaign-result.json` SHA-256:
  `2d4c13acad030b12955684af66e6abe91db6cc9a821d2ff7beec4da7eaa1f5a7`;
- `manifest.json` SHA-256:
  `81effdf35f1ce3600ee415c85b979da0b80669ab5d8156fb9aa56239fb1c1621`;
- campaign stderr is empty; the bounded command exited 0.

## Reproduction

```bash
ROOT=/mnt/data/vsletten/dissertation-data/e3a-classical-neb-campaign-grid-20260903
WORKBOOK=/mnt/data/vsletten/dissertation-data/e3a-classical-neb-barriers-20260902/sources/nteme-2022-supplement/Results\ Of\ MD\ Simulations.xlsx
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
nice -n 10 uv run --with openpyxl python scripts/e3a_campaign.py run \
  --workbook "$WORKBOOK" --out "$ROOT" \
  --min-ftol 0.01 --min-steps 5000 \
  --neb-ftol 0.01 --neb-relax-steps 3000 --neb-climb-steps 10000
```

The command refuses NEB step bounds that are not multiples of the 100-step
thermo interval, writes per-model inputs/logs/dumps/results, and seals all
non-manifest files into `manifest.json`.
