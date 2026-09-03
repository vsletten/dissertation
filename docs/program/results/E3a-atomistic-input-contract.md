# E3a atomistic input contract

**State:** executable input contract; barrier campaign not run

**Operator:** 2026-09-02 `(hermes-custom-build-001; profile=workstation)`
**Source workbook SHA-256:** `292c28c2b182e45aa7fe6cc36a2e62779d7d8dbf0aa6e74525c0e898ee31237f`

## Verdict

The four previously missing E3a environments now have deterministic, hash-pinned initial/endpoint builders, typed limitations, relaxation/acceptance gates, and LAMMPS run-0 smoke receipts for both ends. This closes the input-construction prerequisite; it does **not** claim that the emitted seeds are relaxed minima or supply migration barriers.

Three scientific limits remain deliberately visible:

1. ClayFF is non-reactive, so the dehydroxylation builder creates a post-reaction topology but cannot calculate the dehydroxylation reaction path.
2. Published microscopy does not define one atom-resolved extended-defect-zone structure, and the 6×3×1 cell is too small to contain an observed lens plus undisturbed host. The emitted structure is therefore one member of an explicitly parameterized local-opening sensitivity family, not an observed defect.
3. No stable octahedral Ar minimum is published for muscovite. The emitted cage-centroid seed is typed `incomplete-no-minimum`; E3a must prove an index-zero minimum before constructing an escape NEB.

## Shared source structure and Hamiltonian

Every model starts from the Nteme et al. 2M1 muscovite workbook: 6×3×1 periodic supercell, 1,512 atoms, `a=30.900625 Å`, `b=26.747522 Å`, `c=19.549356 Å`, `α=89.862673°`, `β=95.576801°`, and `γ=89.985803°`. The selected replication route is the first published Ar divacancy row: moving K site 3, vacancies 4 and 88, and published barrier 67.644151 kcal/mol.[1]

The builder retains the recovered ClayFF/LAMMPS contract: harmonic O–H bonds, 10 Å `lj/cut/coul/long`, Ewald `1e-4`, Lorentz–Berthelot mixing, and periodic `p p p`. Nteme et al. explicitly define σ as the LJ zero crossing, ε as the well depth, neutral Ar from White (1999), and Lorentz–Berthelot unlike-pair mixing.[1]

## Charge-compensation rule

Nteme et al. say the Ar/vacancy charge deficit was compensated by Al-for-Si substitution in tetrahedral sites “some unit cells away,” but neither the supplement nor workbook identifies those sites.[1] Under the printed ClayFF charges, literal Si→Al makes a K-vacancy cell *more* negative. The executable contract therefore uses the charge-correct inverse operator and records it as a **reconstructed sensitivity Hamiltonian**, not as Nteme's undisclosed production structure. Its atom-type changes follow the published ClayFF distinction between tetrahedral Al-bonded and Si-bonded oxygen environments.[1][11]

1. Compute the net charge from machine-readable per-atom charges after inserting the noble gas and vacancies.
2. For each missing K charge, choose the tetrahedral `Al2` site farthest from the route midpoint; ties use atom ID. Require selected sites to be at least 4.5 Å apart.
3. Change that `Al2` to `Si`, and retype its four nearest bonded `O3` atoms to `O2`.
4. Assert exactly four unique O neighbors inside 2.3 Å and assert final `|Q| ≤ 1e-8 e`.

With the recovered force-field values, each full operator contributes exactly +1 e: `(2.1000−1.5750) + 4×(−1.0500+1.16875) = +1.0000 e`. This preserves the paper's remote-compensator intent while refusing the sign error and omitted site IDs. Algebraic neutrality is verified; exact equivalence to Nteme's unpublished compensated structures is **not** claimed. Every resulting barrier must retain the `reconstructed-charge-compensation` provenance tag and cannot serve as a strict replication value.

For the selected divacancy route, the three deterministic compensators are `Al2` IDs **1273, 1189, 1281**, retyped with their recorded four-O neighborhoods. The trap model uses remote `Al2` ID **1443**.

## Model A — local dehydroxylate

### Exact transformation

Start from the neutral, compensated Ar divacancy route. Among all periodic O1–O1 pairs separated by at most 3.6 Å, choose the pair nearest the route midpoint; ties use O IDs. The selected hydroxyl oxygens are **73 and 76**, with hydrogens **81 and 84**.

Apply `2 OH(lattice) → H2O(removed) + O_residual(lattice)`:

- remove O1 ID **76** and H IDs **81, 84**;
- retain O ID **73** and retype it from `O1` (−0.95 e) to `O2` (−1.05 e);
- retain the periodic source cell for the seed;
- relax fixed-cell at 0 K, then allow `c`/tilts to relax at 0.5 K and 0 kbar while holding `a,b`.

The local operator is exactly charge-neutral. It produces 1,507 atoms after route vacancies and water removal. Before accepting a relaxed endpoint, require finite energy, minimum distance ≥0.70 Å, no retained water product, max force ≤0.01 kcal mol−1 Å−1, and two octahedral-Al neighbors around residual O73. High-temperature muscovite studies constrain the dehydroxylate as a distinct structure and do not justify calling this a “collapsed interlayer”; the reported crystallographic response includes expansion of `c` during dehydroxylate formation.[2][3]

**Scope limit:** this is a local post-reaction sensitivity topology. ClayFF cannot break the O–H bonds or establish a dehydroxylation barrier.

## Model B — extended-defect-zone sensitivity family

Microscopy shows (001)-parallel extended defects, including lenticular voids/basal partings, distributed through muscovite; their density correlates with high-diffusivity, Na-depleted zones.[4] Independent work supports stacking faults as an Ar fast-path mechanism.[5] Those observations define a defect *class*, not an atom-resolved 6×3×1 configuration.

The builder therefore seeds a periodic cosine-squared lenticular opening centered on the route midpoint. The committed member uses:

- in-plane semi-axes **8 Å × 5 Å**;
- peak opening **1.0 Å**;
- 229 moved atoms;
- periodic `p p p`, fixed source cell;
- route sites 3→4 with vacancy 88 and the same remote compensators.

The sensitivity grid is `a={6,8,10} Å`, `b={3,5,7} Å`, and opening `Δ={0.5,1.0,1.5} Å`. Those numbers are declared model parameters, not measurements. Relax each member with a restrained shell, release the interior, and reject a member if the void heals before NEB. Any result must be reported as cell-size-sensitive because the ~3.09×2.67×1.95 nm periodic cell cannot represent a ≥10 nm observed lens plus pristine matrix.[4]

## Model C — recoil-derived octahedral trap/escape

Recoil supplies a physically distinct intralayer 39Ar reservoir in muscovite experiments, while atomistic work on illite supports stable octahedral retention as an analog rather than a muscovite coordinate.[6][7] Nteme et al. find a small octahedral occupancy class but do not publish a stable octahedral muscovite structure.[1]

The exact seed is therefore typed **`incomplete-no-minimum`**:

- change the positive-(001)-normal gallery K selected by minimum angular departure from the 20° illite recoil analog (tie: distance, then atom ID), which resolves to atom **675**, into neutral Ar without changing atom count;
- move it to the geometric vacant-M1 cage centroid `(10.7175959581, 11.2801376226, 10.0589268001) Å`, computed from pinned oxygen IDs **743, 741, 708, 711, 628, 707**; validate that they are the six nearest oxygens, all within 2.30 Å, with the seventh beyond 2.50 Å;
- compensate the resulting −1 e with remote tetrahedral operator ID **1443**;
- define the explicit exit endpoint as the original K675 gallery coordinate `(11.6819095846, 12.6196715700, 15.0564655344) Å`;
- minimum-image seed-to-exit separation: **5.263045256 Å**.

Search the centroid plus Cartesian `{0, ±0.25 Å}` offsets. Freeze atoms farther than 4.25 Å from Ar, minimize first to 0.05 and then 0.01 kcal mol−1 Å−1, release the cell interior, and minimize again. Accept only a finite, index-zero endpoint with max force ≤0.01 kcal mol−1 Å−1, Ar remaining in the same octahedral cage, and minimum distance ≥0.70 Å. Otherwise preserve `incomplete-no-minimum` and run E3b's narrower endpoint-only periodic-DFT optimization. The ~55 kcal/mol illite analog must never be promoted as a computed muscovite barrier.[6]

## Model D — Xe divacancy screen

The Xe route uses the same sites and compensation as Model A, but replaces moving K3 with neutral Xe. The strongest directly compatible source found is the Bourg–Sposito neutral-Xe 12-6 LJ model used with ClayFF in clay-interlayer simulations:[8][9][10]

- `q_Xe = 0`;
- `σ_Xe = 4.00924 Å`;
- `ε_Xe = 0.0185900 eV = 0.4286955841712072 kcal/mol`;
- Lorentz–Berthelot mixing with ClayFF.

This is an **executable screening parameterization**, not a validated dry-muscovite Xe Hamiltonian. Any Xe barrier remains sensitivity-only until E3b or independent dry-muscovite adsorption/solubility data validate the cross interactions.

## Generated artifacts and verification

Generator: `scripts/e3a_input_contract.py`

Tests: `tests/test_e3a_input_contract.py`, `tests/test_nteme_neb.py`

Heavy artifact root: `/mnt/data/vsletten/dissertation-data/e3a-atomistic-input-contract-20260902/`
Receipt: 37 files, 989,628 bytes; `manifest.json` SHA-256 `b6b7d2503c84dfb9847365c000f44a6588c983eef76a779e162259519f8489d6`

Each model directory contains `structure.data`, `endpoint.data`, `in.static`, `in.endpoint.static`, per-end LAMMPS logs, and `manifest.json`. The root `manifest.json` aggregates provenance, atom/type counts, charge, periodic cell, minimum distance, static energy, force, and log hashes for both endpoint seeds.

| Model | Atoms | Q (e) | min distance (Å) | run-0 PE (kcal/mol) | run-0 max force | Status |
|---|---:|---:|---:|---:|---:|---|
| dehydroxylate-lattice | 1507 | 7.22e-15 | 1.02421 | −348211.59 | 104.99 | unrelaxed executable seed |
| extended-defect-zone | 1510 | 7.33e-15 | 1.01444 | −348033.78 | 212.17 | sensitivity-family seed |
| octahedral-trap-escape | 1512 | 7.77e-15 | 1.02421 | −348016.12 | 987.34 | `incomplete-no-minimum` |
| xenon-divacancy | 1510 | 7.33e-15 | 1.02421 | −348186.60 | 105.23 | screening seed |

The large run-0 forces are expected for imposed, unrelaxed defect seeds and are **not** convergence claims. They force the downstream campaign through each model's explicit relaxation/minimum gate.

Reproduction:

```bash
ROOT=/mnt/data/vsletten/dissertation-data/e3a-classical-neb-barriers-20260902
OUT=/mnt/data/vsletten/dissertation-data/e3a-atomistic-input-contract-20260902
uv run --with openpyxl python scripts/e3a_input_contract.py build \
  --workbook "$ROOT/sources/nteme-2022-supplement/Results Of MD Simulations.xlsx" \
  --out "$OUT"
python3 -m pytest tests/test_nteme_neb.py tests/test_e3a_input_contract.py -q
```

Observed result: **9 contract/replication tests passed**, including the real hash-pinned workbook integration; the full `qm` gate passed **468 tests with 1 skipped** via `uv run --extra dev python -m pytest -q`; all eight initial/endpoint LAMMPS run-0 jobs exited 0; every emitted endpoint was finite, periodic, charge-neutral to `1e-8 e`, and above the 0.70 Å collision floor.

## Sources

[1] https://insu.hal.science/insu-03668854/document
[2] http://minsocam.org/ammin/am72/am72_537.pdf
[3] https://doi.org/10.2465/ganko1941.69.381
[4] https://insu.hal.science/insu-03890387/document
[5] https://doi.org/10.1144/GSL.SP.2003.220.01.15
[6] https://doi.org/10.1016/j.gca.2015.03.005
[7] https://doi.org/10.1016/S0016-7037(97)00323-2
[8] https://doi.org/10.1016/j.gca.2008.02.012
[9] https://doi.org/10.1021/acs.jpcc.7b09768
[10] https://doi.org/10.1021/acs.jpcc.7b09768.s001
[11] https://doi.org/10.1016/j.jcis.2003.12.029 — Cygan et al. 2004 — Molecular models of hydroxide, oxyhydroxide, and clay phases
