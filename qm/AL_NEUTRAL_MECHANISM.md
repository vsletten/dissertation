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
at this model tier.

Absolute numbers remain flagged for the A2 re-tier. The mechanism and ordering
are the defensible Phase-1 result; this free-dimer value must not be substituted
for a lattice-embedded production barrier.

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
