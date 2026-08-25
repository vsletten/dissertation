# Phase-1 acid hydrolysis: verified result and model limit

Date: 2026-08-22

## Executive result

The one-water `al-acid` cluster closes as a **sequential associative acid
hydrolysis** mechanism at B3LYP/def2-SVP/DF:

1. residual water addition to the pre-protonated Si–O–Al complex;
2. cleavage of Si–Obr from the associative intermediate.

The validated profile barrier is **82.193 kJ/mol (19.645 kcal/mol)** at
298.15 K. It is lower than the validated `si-neutral` reference,
113.048 kJ/mol (27.019 kcal/mol), and close to Xiao & Lasaga's acid anchor
(~24 kcal/mol). Addition is rate limiting; cleavage is shallow once the
associative intermediate exists.

The matching one-water `si-acid` reactant is **not a minimum** at this model
tier. Unconstrained gas-phase B3LYP, a 4.5 Å residual-water placement, and
aqueous PCM all return the bridge proton to the residual water (H3O+).
SMD/GPU geometry optimization was numerically pathological and was stopped.
Therefore no honest Si–O–Si acid barrier or Si–O–Al/Si–O–Si acid ordering is
reported. The card is blocked on a microsolvated Si-acid model rather than
laundering a constrained non-minimum into thermochemistry.

## `al-acid` verified profile

Run directory (ignored campaign state):

`runs/phase1/al-acid-preprotonated-v2-b3lyp-def2-svp-flank/`

| Quantity | Result |
|---|---:|
| Mechanism | sequential-associative-acid |
| Addition ΔG‡(298 K) | 82.193 kJ/mol = 19.645 kcal/mol |
| Addition imaginary mode | 205.95i cm⁻¹ |
| Cleavage local ΔG‡(298 K) | 8.047 kJ/mol = 1.923 kcal/mol |
| Cleavage imaginary mode | 97.46i cm⁻¹ |
| Intermediate ΔG vs reactant | +62.814 kJ/mol |
| Overall profile ΔG‡ | 82.193 kJ/mol = 19.645 kcal/mol |
| Rate-limiting local step | addition |

Both saddles have exactly one imaginary mode. Quick-IRC endpoint signatures
are exact:

- addition: `(False, True, True, 1, 2, 0)` ↔ `(True, True, True, 1, 2, 0)`;
- cleavage: `(True, True, True, 1, 2, 0)` ↔ `(True, False, True, 1, 1, 1)`.

The tuple is `(Si–Ow bonded, Si–Obr bonded, opposite-center–Obr bonded,
bridge proton count, attacker proton count, terminal-framework proton count)`.
Every tracked acid proton is assigned uniquely to its nearest oxygen, so a
dissociated/missing proton cannot pass as the terminal rotamer. The addition
R/I endpoints also pass heavy reactive-core RMSD and same-method energy gates.
Reactant and associative-intermediate Hessians have zero imaginary modes.

Artifact receipts:

| Artifact | SHA-256 |
|---|---|
| `results.json` | `2c584be0340d17a2d9ee67049dced03d98515640985a12379ac755eb70a86ff2` |
| `store.sqlite` | `330d36b02f1c887d289715f9871afb6658cd8be746cc74c92cf17ecd345bcc3d` |
| `addition_scan/scan.json` | `28bac5ae20ae10f04cce986499f04682c1bb2ad324da65a9a94ec6286e914768` |
| `addition_scan_directed_ts.xyz` | `14256e743add5ede70f7c1f8ce1d842e8f27f4d5a04c482d9c79d401b4693175` |
| `cleavage_ts.xyz` | `e7ae4a246ef49d38c41d2b0236e898b88ecdf4bf9ffebb1e0fa30718f18dd68d` |

SQLite `PRAGMA integrity_check` returned `ok`; the store contains four
structures, four completed frequency jobs, and four electronic results.
`run_status.json` is `completed` and states that both saddle/IRC gates passed.

## Why the mechanism implementation changed

Live optimization exposed three false-green assumptions in the initial acid
route:

- Hydrolysis transfers one attacker proton to the leaving framework. The acid
  product is Si–OH, not Si–OH2; proton-count gates now encode that chemistry.
- The Al-acid path is sequential. The first validated saddle connects the
  associative intermediate to hydrolyzed product, not the initial reactant.
- Endpoint CI-NEB mixed unrelated terminal-OH conformers. The chemically
  relevant addition crest is instead located by a relaxed Si–Ow scan, then
  refined with an endpoint-directed Cartesian Sella search.

ODE climbing-image optimization replaced an MDMin climb that exhausted its
240-step bound. Long phase boundaries return CuPy pool blocks before Hessian
work. Rejected saddle checkpoints are quarantined instead of being reused on
the next run.

## `si-acid` negative result

The pre-equilibrium contract requires the optimized reactant signature
`(False, True, True, 1, 2, 0)`: protonated intact bridge plus residual H2O. The
one-water disilicate cluster instead optimized to
`(False, True, True, 0, 3, 0)`: intact neutral bridge plus H3O+.

Receipts:

| Probe | Outcome | Log SHA-256 |
|---|---|---|
| Gas-phase production geometry | bridge proton returned to water | `c0ab76b2af643f9fb961425bb9f58f049fba35c5e2f6dcfdc875e8543cba9f8a` |
| Residual water placed at 4.5 Å | same `(0, 3)` proton counts | `00d47740a2dda94869491308f3cb9d0e52aee1ed39c118b39631c0432d750753` |
| Aqueous PCM | same `(0, 3)` proton counts after 100 optimizer steps | `f7e192e4ae83611a7f446101f6ed258389c40d35997aa3a514a9ea57b454cc0c` |
| SMD/GPU probe | nonphysical energy descent / O(1) Ha/Bohr gradients; stopped | `233942dabcc318cc980ed38810397cde0ee1623b1b4006a5771d69c9c9a6d61d` |

The mechanism-v2 driver reverified this failure, wrote `run_status.json` as
`blocked`, and left no canonical `results.json`/`store.sqlite`; log SHA-256 is
`13d8c392cf991059d35274469b7db5148f0bab94db7b3c74b5aeabceb8ad9809`.

This is a model-validity failure, not a compute failure. Holding Obr–H by a
constraint would make the requested reactant exist but would invalidate the
unconstrained frequencies and ΔG‡ comparison.

## Matched microsolvation survey (3--6 waters)

Victor selected matched microsolvation for both dimers on 2026-08-23. The
continuation driver now supports 3--6 total explicit waters, isolates every
water count in its own checkpoint namespace, counts the attacker and each shell
water separately, checks every physical H atom (including framework OH), and
refuses all TS work until the unconstrained reactant passes connectivity and
zero-imaginary-mode gates.

The first bounded Si--O--Si survey tested one deterministic hydrogen-bond shell
seed at every requested water count. Every production B3LYP/def2-SVP/DF
geometry endpoint lost the bridge proton **before** a Hessian or TS search:

| total waters | production geometry | acid signature after optimization | disposition |
|---:|---|---|---|
| 3 | converged | `(False, True, True, 0, 7, 0)` | proton moved into the solvent network |
| 4 | converged after continuing the original 100-step endpoint | `(False, True, True, 0, 9, 0)` | residual water became H3O+ |
| 5 | converged | `(False, True, True, 0, 10, 1)` | proton moved to a terminal framework O |
| 6 | converged | `(False, True, True, 0, 12, 1)` | proton moved to a terminal framework O |

The signature is `(Si--Ow bonded, Si--Obr bonded, opposite-center--Obr bonded,
bridge mobile-H count, solvent mobile-H count, terminal-framework mobile-H
count)`. All four endpoints retain the intact Si--O--Si bridge but have zero
protons on Obr. The strict reactant gate requires `(False, True, True, 1, 2n,
0)` plus exact physical occupancies: attacking H2O, neutral two-H shell waters,
one H on Obr, and one baseline H on every terminal framework oxygen.

This is a bounded conformer survey, not a proof that no protonated-bridge local
minimum exists anywhere on the microsolvated potential-energy surface. It does
prove that merely adding 3--6 waters to the settled deterministic shell does
**not** rescue the pre-equilibrium assumed by this card. Running the matched Al
barrier would not produce a reportable pair while every Si count is red, so no
Al TS campaign was launched and no unsupported ordering was published.

Durable ignored-state receipts:

| waters | optimized `complex.xyz` SHA-256 | reactant-status SHA-256 | teed production log SHA-256 |
|---:|---|---|---|
| 3 | `785d6ea430f41e20ab1245c2539877779216f38a2a66fb72301274754812f971` | `c38114b0b2b18d66583e9c4156aa17be9a95713e6832cdf6d8399f7f26e1be51` | `cc9a702bd1ee449917a0441644f9a96c7b60a56dc305be181a581b2c0ebedeaa` |
| 4 | `e4dc46a999d3330a55574b3741e2e0b3b7483724a25c9a7300ed4da73b96a11d` | `6519337ce38565541af9efc4302da15e74e3d46ed0714b66eda9f4d148c0ac68` | refinement: `19f7ff2078d06d4302e0f60aa6c8aed64db1cb1462c72138a1c3d9542792a924` |
| 5 | `7f006d12573e64b1be7f27c223f56732b45cc507a74faaef04bb20ff678e326c` | `9f6312781701d6a2d0ac7d96f0cc529de31ddc75a59814daf261f4ce7dcf211c` | `c0cccbc3d96aa98aae84aab1bbdeb1f839e7f4efa3a01c794b03897b1fc7fe99` |
| 6 | `3fac28ca17c83c15a4161964de6d21d91de99ea97e5e2050350d05871563c95f` | `a238bc65ee782f28ae2e182231d9171057dd3a9b68c5a66fb53eb40871a8e5a9` | `dbc8d6bc1f3eedfbf4530ef847640009ce773d8f36e28522ec8917045c904ebf` |

## Finite conformer closeout (A1c/A1d)

A1c then ran four explicit proton-relay families at every water count 3--6.
Fifteen B3LYP/def2-SVP/DF endpoints converged and failed the same exact
physical-H/connectivity/solvent-occupancy basin; only
`si-acid:3w:bridge-donor-chain` exhausted its original 160-step bound. A1d
restarted from that exact archived production endpoint with a fresh geomeTRIC
state and one additional bound of 160 steps. It converged after 54 steps to
signature `(False, True, True, 0, 7, 0)`, again placing zero mobile H on Obr
instead of the required `(False, True, True, 1, 6, 0)`.

The hash-bound effective population is therefore exactly **16 rejected / 0
accepted / 0 failed / 0 blocked**. This is the predeclared finite-ensemble gate:
the pre-equilibrated protonated-bridge assumption is now a **model-valid NO-GO**
for Si--O--Si at this tier. No Al screen or unmatched acid ordering was emitted.
The immutable A1c manifest remains
`471a099c01c3fd8107409fd087e033bdd5da167cc6d41f09475a6a661cae5ac6`;
A1d's resolution manifest is
`6b02432f4920b2ff8406de5cd22d379ebc54ceed66800b1db8149efff7df62cf`
under `/mnt/data/vsletten/dissertation-data/task199-a1d-acid-3w-bridge-20260824/`.

## Concerted hydronium-as-nucleophile result (A1e)

A1e retired the pre-equilibrium and ran the predeclared finite concerted scope
at B3LYP/def2-SVP/DF: total waters 3 and 4 crossed with
`bridge-donor-chain` and `compact-cyclic-relay`. The reactant began with an
intact unprotonated Si--O--Si bridge and H3O+ in the network. The exact product
made that same hydronium oxygen the Si--Ow nucleophile while the physical-H
relay protonated Obr. Every endpoint was unconstrained in production.

All four Si paths are conclusive rejections, with no computational failures:

| path | terminal gate | exact result |
|---|---|---|
| 3w bridge-donor-chain | optimized product | `(False, False, True, 1, 5, 0)`; Si--Ow dissociated |
| 3w compact-cyclic-relay | optimized product | `(False, False, True, 1, 5, 0)`; Si--Ow dissociated |
| 4w bridge-donor-chain | optimized reactant | solvent occupancy `(2,2,2,3)` instead of H3O+ `(3,2,2,2)` |
| 4w compact-cyclic-relay | optimized reactant | `(False, True, True, 0, 8, 1)`; proton moved to framework O |

The required endpoint gates correctly stopped every path before CI-NEB. No
pre-climb band, saddle, IRC, Al screen, barrier, or ordering was fabricated from
a rejected basin. The implementation nevertheless lands the reusable,
regression-tested downstream contract: hash-bound finite driver, persisted
pre-climb/climb checkpoints, coupled relay/attack/cleavage one-mode gate, and
bounded full bidirectional Sella Gonzalez--Schlegel IRC.

Durable evidence root:
`/mnt/data/vsletten/dissertation-data/task200-a1e-acid-concerted-20260824/`.
Evidence manifest SHA-256 is
`3162e39f2f0f263444741b03e0a023b1209bf515a9f8291c330ca9995811e622`;
production log is
`627914bbb22ee655932508ea974f006a503389198f656b173fe225607a8ca784`;
run manifest is
`dfea131556fddf907ccabfc52945d390695f99d372cfdeac6859a4f2c0549982`.

## Separated donor / neutral attacker result (A1f)

A1f tested exactly one four-water
`separated-donor-neutral-attacker` Si--O--Si topology at
B3LYP/def2-SVP/DF. Atom order remained the A1e order, but the solvent roles were
outer H3O+ donor, measured-flank neutral H2O attacker, bridge-side relay water,
and spectator/solvator. The exact product preserved the attacker's two physical
H while moving donor H through the relay onto Obr; neither endpoint held a
production coordinate.

The Si reactant converged after 68 optimizer steps and passed finite-coordinate,
minimum-pair-distance (`0.964661 A`), connectivity, aggregate proton-count, and
unique-ownership gates. It conclusively failed only the exact role identity:
physical H16 moved from outer donor O15 to bridge-side relay O22. Solvent
occupancies changed from required `(3,2,2,2)` to `(2,2,3,2)` for
`(outer donor, neutral attacker, bridge-side relay, spectator)`. The neutral
attacker remained H2O, the Si--O--Si bridge remained intact, and no proton moved
to Obr or the framework.

This is a typed `rejected` terminal at `optimize-reactant`, not an optimizer,
NEB, Hessian, or IRC failure. The endpoint gate correctly stopped before path
finding. No Al path, barrier, or ordering was emitted. The migrated proton is a
new physical topology rather than permission to relabel A1f: executable A1g
tests bridge-side H3O+ with the distinct neutral attacker directly under a new
mechanism/gate identity.

Durable evidence root:
`/mnt/data/vsletten/dissertation-data/task201-a1f-acid-separated-20260824/`.
The evidence manifest covers six files / 36,212 bytes and has SHA-256
`940cd4b754f8b67244fd9e51812d4504a7872bb4cf1bb797bdc81776c7777d9c`;
production log is
`bcaa68f2c415653b235c61c9b537ea0fe8f8cae6f35911c3f500a9d40243073b`,
run manifest is
`3a72c6fdaef2ac8fe0730a91b44604e8b83dded4f9a5dfa7275818709dee9efb`,
and terminal is
`0b51c68336cf3dff739f51135664f6ac7918bd2374500f642c8f112bc1b5f095`.

## Bridge-side hydronium / neutral attacker result (A1g)

A1g promoted the exact topology selected by A1f rather than relabelling its
rejected endpoint. The mechanism-v3/gate-v3 Si run began from A1f's archived
optimized reactant SHA-256
`3a5107648bd44cbba15dd7482b4a9e8b3c6434f049442453744737e2c32aba75`.
Atom order and physical-H identity stayed fixed: migrated H16 made bridge-side
O22 hydronium, measured-flank O19/H20/H21 remained the distinct neutral
attacker, and bridge-facing H23 was the proton assigned to Obr in the product.

The source reactant reconverged in six optimizer steps. It retained intact
Si--O--Si, exact solvent occupancies `(2,2,3,2)`, unique H ownership, a neutral
attacker, and zero imaginary modes. The independent hydrolyzed product then
converged in 89 steps but left the required family:

- Si--Ow remained nonbonded at `2.642649 A`;
- Si--Obr cleaved to `3.535115 A`;
- attacker H21 became unassigned to any oxygen;
- exact basin was `(False,False,True,1,7,0)`, not required
  `(True,False,True,1,8,0)`;
- minimum pair distance remained collision-safe at `0.964884 A`.

This is a conclusive `rejected` terminal at `optimize-product`, not a compute
failure. The endpoint gate stopped before CI-NEB, saddle, IRC, Al, barrier, or
ordering. At B3LYP/def2-SVP/DF, accepting the A1f proton migration still does
not stabilize an atom-matched hydrolyzed Si product with a distinct neutral
attacker.

Durable evidence root:
`/mnt/data/vsletten/dissertation-data/task202-a1g-acid-bridge-side-20260824/`.
Evidence manifest SHA-256 is
`76d1d0275a773d7941a855250696cedc78a169ad68d58305195aa9df05d4564f`;
production log is
`7c6ab526c88b86a9f7de4031b2beedbd3a2f4b8c3e8276d8a05927b0ae133a46`;
run manifest is
`27f7a7040cb5f087d9aee7126676e25bed309eed243fd024a1579aa1f946fba2`;
terminal is
`cf74851972575f2a734d03b38c95fea85bc49acf67598739d3f01aea191c7834`.

## Required follow-up

Do not retry A1e's hydronium-as-nucleophile topologies, A1f's outer-donor role,
A1g's B3LYP/def2-SVP topology, or weaken any physical-H/family gate. The only
authorized acid-topology follow-through is board card
`A1i-acid-production-tier-bridge-side-revisit`, blocked until A2 settles the
calibrated production geometry/solvation tier. It retests this exact hash-bound
topology once; endpoint rejection there closes the acid topology program.
Only a fully accepted Si path may trigger exact matched Al. Never hold Obr--H,
Si--Ow, or Si--Obr to manufacture thermochemistry.

Reproduction commands must clear inherited environments:

```bash
env -u PYTHONPATH -u VIRTUAL_ENV \
  uv run --frozen --extra gpu --extra dev \
  python scripts/bridge_side_acid_relay.py --gpu \
  --run-dir <durable-run-dir> --threads 16 --nice 10 \
  --log <durable-log> \
  --source-reactant <a1f-reactant-optimized.xyz> \
  --source-reactant-sha256 \
  3a5107648bd44cbba15dd7482b4a9e8b3c6434f049442453744737e2c32aba75
```

The RTX 4090 must be isolated from both email extraction and Honcho memory
inference during GPU4PySCF work; see `.claude/learnings/gotchas.md`.
