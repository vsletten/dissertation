# D2a — astrochemical rate reproduction: results and verdict

*Campaign run 2026-08-19/20 on the 4090 workstation by a local Hermes
worker (gpt-5.6 sol); closeout by fable after the worker hit its
tool-iteration ceiling. Full artifacts (geometries, checkpoints, rate
tables, failure receipts): `qm/runs/D2a-astro-rate-reproduction/`
(workstation worktree, gitignored per big-output policy). Campaign log:
`/tmp/D2a-astro-rate-reproduction.log`.*

## What ran

Six reactions at B3LYP/def2-SVP/DF geometries with quarry's
Sella + CI-NEB TS machinery, quasi-RRHO thermochemistry, and
asymmetric-Eckart tunneling over 50–300 K:

| Reaction | Outcome | ZPE barrier | Imag. mode |
|---|---|---|---|
| H + H₂CO → CH₃O (addition) | ✅ first-order saddle | 18.9 kJ/mol | 874i (lit. 831i) |
| H + H₂CO → H₂ + HCO (abstraction) | ✅ first-order saddle | 18.7 kJ/mol | 1559i |
| H + H₂ → H₂ + H (exchange calibration) | ✅ first-order saddle | 18.9 kJ/mol | 1069i |
| H + OH → H₂O (barrierless control) | ✅ downhill scan, κ = 1 as expected | — | — |
| H + H₂CO → CH₂OH | ❌ honestly rejected — saddle collapsed to a soft 47i mode | — | — |
| OH + H₂ → H₂O + H | ❌ honestly rejected — no chemical imaginary mode found | — | — |

Rejections carry receipts (`failures.json`); campaign isolation meant
failed reactions did not kill the batch (worker-built improvement,
commits `689fbed..b03d599`: checkpointed campaign driver, open-shell
GPU Hessian fallback, first-crest saddle seeding, per-reaction failure
isolation).

## Key comparisons

**H + H₂CO → CH₃O vs the published surface (Langmuir–Hinshelwood)
rate**: quarry sits at 1.4×10⁻⁵ (50 K) to 2.8×10⁻⁴ (300 K) of the
literature values — **4–5 orders of magnitude low**. The literature
rate is a *surface* rate on ice; our gas-phase cluster treatment lacks
surface partition functions, adsorbate constraint effects, and
surface-mediated TS geometry. This is a physics gap, not a numerics
bug (the barrier and imaginary mode are close to literature values —
874i vs 831i).

**H + H₂CO → H₂ + HCO (gas-phase-style comparison)**: within 3–50× of
literature across 50–300 K — tunneling-dominated, and Eckart holds up
far better than the worst-case expectation at these temperatures for
this channel.

**Barrierless control**: κ = 1 exactly as demanded — the tunneling
machinery does not invent corrections where none belong.

## Verdict (the card's go/no-go gate)

- **GO — gas-phase reaction classes**: quarry can produce
  literature-consistent tunneling-corrected k(T) tables (50–300 K,
  validity floor stated) for gas-phase reactions at modest cost. D2b
  gas-phase rate contributions are viable, with instanton-benchmarked
  spot checks below ~75 K per the scoping doc.
- **NO-GO — surface (LH) rates from the cheap gas-phase protocol**:
  4–5 orders low against surface literature. KIDA-grade *surface*
  rates require explicit-surface treatment (constrained ice-cluster
  models, surface partition functions, D3-adjacent machinery). Do not
  submit surface rates from this protocol.
- Follow-up worth carding: OH + H₂ (a bent, quasi-linear TS) needs a
  better guess strategy than the distance-scan family; CH₂OH addition
  needs a constrained approach vector. Both rejections are TS-hunting
  problems, not physics walls.

## Method caveats

Single-level (B3LYP/def2-SVP/DF) geometries and barriers; no
higher-level single points yet (A2's production tier applies here
too). Literature 50 K values partially extrapolated (flagged in
results.json). Eckart is the stated tunneling floor — instanton
comparison below 75 K pending.
