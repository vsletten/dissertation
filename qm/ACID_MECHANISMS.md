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

## Required follow-up

Build a **3–6 explicit-water microsolvated Si-acid cluster**, as already
recommended by `HANDOFF.md`, and prove the protonated-bridge reactant is an
unconstrained minimum before any TS search. Then run the same validated
addition/cleavage gates for both Si–O–Si and Si–O–Al under like-for-like
microsolvation. Only that paired result can close the acid ordering claim.

Reproduction commands must clear inherited environments:

```bash
env -u PYTHONPATH -u VIRTUAL_ENV \
  uv run --frozen --extra gpu --extra dev \
  python scripts/phase1_xiao_lasaga.py --reaction al-acid --gpu \
  --threads 16 --nice 10 --log <durable-log>
```

The RTX 4090 must be isolated from both email extraction and Honcho memory
inference during GPU4PySCF Hessians; see `.claude/learnings/gotchas.md`.
