# Al-neutral hydrolysis: sequential uphill-slide mechanism

Date: 2026-08-23

## Executive result

The free one-water Si–O–Al dimer at B3LYP/def2-SVP/DF resolves as a
**sequential associative mechanism with a quasi-barrierless uphill addition
slide**, followed by a verified bridge-cleavage saddle.

The reactant-referenced profile barrier is **134.811 kJ/mol = 32.221 kcal/mol**
at 298.15 K. This is **21.763 kJ/mol = 5.201 kcal/mol higher** than the
validated free-dimer `si-neutral` reference, 113.048 kJ/mol = 27.019 kcal/mol.
Xiao & Lasaga's predicted easier Si–O–Al ordering is therefore **not confirmed**
at this neutral model tier; it is not a test of their catalyzed ordering.

Absolute numbers remain flagged for the A2 re-tier. The mechanism and ordering
are the defensible Phase-1 result; this free-dimer value must not be substituted
for a lattice-embedded production barrier.

## Prior work reconciliation

### Victor's 1990s anionic dimer

The old and new campaigns use the **same terminated dimer core and electronic
state**: `(HO)3Si–O–Al(OH)3-`, formula `H6SiAlO7-`, with three terminal
hydroxyls on each center, one Si–O–Al bridge, total charge -1, and singlet spin.
TASK-168 constructed that core with `aluminosilicate_dimer()` and added one
neutral water, so its reacting complex is `H8SiAlO8-` and remains charge -1.
This identity follows both the current builder and the campaign's exact template,
not an inference from the result label. [TASK-168 code: `qm/quarry/clusters.py`
`aluminosilicate_dimer()` and `qm/scripts/task168_intermediate_sweep.py`
`run_sweep()`, which builds the template from that dimer plus one water.
The campaign tree is not in this repository; the durable `results.json`
SHA-256 is `94ecd8ea45dbea17743d7ce9ad354beb4ed150ff97d988714e6f0dd895ccb148`
in Durable receipts below.]

Victor's thesis used Gaussian 94/B3LYP, B3LYP/3-21G(3d,2p) geometries, and
6-31G(3d) single-point energies. Its dimer table reports **29.9267 kcal/mol**
for the direct net reaction from `H6SiAlO7- + H2O·H2O` to silicic acid,
`Al(OH)4-`, and a spectator water; the pre-adsorbed-water hydrolysis step is
**35.7769 kcal/mol**. [Victor Sletten, *Playing Dice with the Universe*
(~1999), ch. 4, “Computational Method” and “Reaction Model”; exact “Dimer
Hydrolysis” values from the thesis dimer table. Those LaTeX sources live in
the workstation thesis archive and are not yet ingested here; see
`docs/program/cards/A8-thesis-archive-intake.md`.]

The modern **32.221 kcal/mol** result is only **2.294 kcal/mol above** the old
29.9267 number and therefore corroborates the energy scale Victor found in the
1990s. The comparison is deliberately qualified: the thesis values are
zero-point-corrected reaction energies (`ΔE*`), whereas TASK-168 reports a
298.15 K quasi-RRHO activation free energy (`ΔG‡`); the old direct bookkeeping
also carries a second, spectator water while TASK-168 has one explicit water.
This is strong numerical continuity for the same anionic dimer core, not an
apples-to-apples reproduction or a claim of 2.294-kcal method accuracy.

### Catalyzed and embedded tiers

Xiao & Lasaga's ~24 kcal/mol acid and ~19 kcal/mol base anchors—and the
associated “Si–O–Al is easier than Si–O–Si” ordering—belong to catalyzed
channels. The neutral free-dimer result therefore **narrows rather than
refutes** that claim: neutral `al` is 5.201 kcal/mol higher than neutral `si`
here, while the matched microsolvated A1b acid pair is the direct test of the
catalyzed ordering. [Xiao & Lasaga, *Geochimica et Cosmochimica Acta* **58**
(1994) 5379–5400 and **60** (1996) 2283–2295; `SURVEY.md` §§7.1 and 9.]

The lattice-embedded ordering is a third, still-open tier. Periodic work reports
Si–O–Al **86 kJ/mol** versus Si–O–Si **169 kJ/mol** in metakaolinite, but a
free terminated dimer cannot confirm or overturn that lattice result. The
program must carry the question into the embedded barrier-ladder rungs.
[Izadifar et al., *Nanomaterials* **13** (2023) 1196 and *Nanoscale Advances*
**7** (2025); `SURVEY.md` §7.2.]

## Intermediate diagnostic

TASK-168 first ran a bounded eleven-seed proton/hydroxyl conformer sweep. Every
candidate was optimized at the production tier and rejected unless it retained
the exact associative signature `(Si–Ow, Si–Obr, H–Obr) = (true, true, true)`.

The original saved intermediate was not the lowest sampled conformer. Rotating
the terminal hydroxyl on O11 by 180 degrees produced the selected minimum:

| Quantity | Result |
|---|---:|
| Electronic energy | -1138.416589754644 Ha |
| Electronic energy vs reactant | +100.275 kJ/mol |
| Lowering vs reoptimized baseline intermediate | 9.551 kJ/mol = 2.283 kcal/mol |
| Quasi-RRHO `G(I) - G(R)` | +102.645 kJ/mol = 24.533 kcal/mol |
| Imaginary modes | 0 |

This triggered the approved endpoint-revision branch: both elementary bands
were rebuilt toward the lower intermediate instead of reusing the old
high-intermediate profile.

## Addition finding

The rebuilt seven-image R→I band used bounded MDMin pre-relaxation, an ODE
climbing guess, crash-persistent image checkpoints, and endpoint-directed
Cartesian Sella refinement. The band entered the associative basin, but the
refined stationary point had two soft imaginary modes, **10.18i and 30.32i
cm⁻¹**, and no accepted R↔I first-order reaction mode. It was rejected rather
than laundered into thermochemistry.

Per the pre-approved kill criterion, there was no attempt #5. The addition is
recorded as a **quasi-barrierless uphill slide into a high-lying associative
intermediate**. The previously observed ~150.471 kJ/mol electronic crest is
retained only as the bracketed lower-bound receipt requested in the campaign
decision; it is not called a transition state or a Gibbs barrier.

## Verified cleavage and profile

The new lower-I→P band was independently rebuilt and refined. Its saddle has
exactly one **84.42i cm⁻¹** mode, and 0.50 Å quick-IRC reaches the exact
associative-intermediate and hydrolyzed-product basins.

| Quantity | Result |
|---|---:|
| `G(I) - G(R)` | 102.645 kJ/mol |
| Cleavage local `ΔG‡ = G(TS_cleave) - G(I)` | 32.166 kJ/mol = 7.688 kcal/mol |
| Overall `ΔG‡ = G(TS_cleave) - G(R)` | **134.811 kJ/mol = 32.221 kcal/mol** |
| Overall `ΔH‡` | 129.646 kJ/mol |
| Wigner κ | 1.00692 |
| Overall `k(298 K)` | 1.5081×10⁻¹¹ s⁻¹ |

## Durable receipts

The ignored campaign directory is
`qm/runs/phase1/al-neutral-b3lyp-def2-svp-flank/`. Its durable copy is named in
the mission-control TASK-168 Result before worktree teardown.

| Artifact | SHA-256 |
|---|---|
| `results.json` | `94ecd8ea45dbea17743d7ce9ad354beb4ed150ff97d988714e6f0dd895ccb148` |
| `store.sqlite` | `d93aa48f4633b0f48d8e8a814246fca3adbc4aee179fbf3b835e91180f50dc88` |
| `task168-intermediate-sweep/manifest.json` | `8fe3c2143110d10b1e1477b6473829335d6d8913f8d4984e5c7b572330baf19f` |
| `task168-resolution-20260823.log` | `ab24ecc9144a7d7c3323471ccbae84316ffb692503207d23d47078e83baa5938` |
| `intermediate.task168-selected.xyz` | `663a0d43ce8e81a7cfb1096c167dfd8afaa7b52d7443e7fc5b842daed8b18f1e` |
| `task168-cleavage-lower-i-ts.xyz` | `83d6b2c2d9877fdf69105e0ac0cf0ecceb1fb88f7ae730807d3666bc879e8b83` |
| `task168-cleavage-lower-i-ts.frequency.npz` | `b718acfe7614d1b5ae3b921298af206f74a550f701ce675d841275bef1019c9b` |

SQLite `PRAGMA integrity_check` returned `ok`; the store contains three
structures, three completed frequency jobs, and three electronic results.
`run_status.json` is `completed` with the Victor-approved mechanistic-finding
closeout.
