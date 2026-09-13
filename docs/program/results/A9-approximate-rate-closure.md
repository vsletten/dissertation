# A9 approximate rate closure — preliminary campaign evidence

**Status:** INCOMPLETE / NOT ACCEPTED
**Progress:** 2026-09-13 16:29 PDT — (hermes-custom-build-001; profile=workstation)

## What completed

The A9 implementation now supplies:

- a 298 K survey-tier kaolinite deck with 22 reaction-rate mappings and explicit `computed`, `literature`, or `heuristic` provenance;
- fail-closed validation of temperature, reaction-family coverage, observables, eight unique seeds, and the declared dilute reservoir;
- a bounded resumable runner for one nominal scenario plus four perturbations for each of seven families;
- observed gross/net event accounting, exact event-to-moles and lattice-area conversions, population and `results.dat`-equivalent series, zero-event Poisson bounds, and a separately labeled integrated-CTMC-propensity estimator;
- raw/derived hash verification and sabotage regressions.

A real campaign completed all 29 scenarios x 8 fixed-seed replicas x 20,000 events: **232/232 Petra runs**, with no incomplete or absorbed trajectories. The whole campaign took `56.36 s`; every child had a 600-second timeout. Because the ensemble stayed far below the card's 30-minute threshold, no `systemd-run` wrapper was required.

## Why this is not the A9 result

Independent scientific review found three material acceptance failures. No downstream card may consume the provisional ranking or rate comparison.

### 1. Steady state is not established for the response being ranked

The current event gate types zero-release trajectories as `steady-zero` even when the integrated desorption propensity is still evolving. One reviewed case (`al-o-al-analogue__prefactor-x0.1`, seed `90401`) changes Si-desorption propensity by roughly `5.03e4` between steps 8,000 and 20,000. The population diagnostics also label 22/232 combined-population traces and 32/232 Si-population traces `evolving`.

Completion of 20,000 events is therefore not proof of the card's required steady state. The next pass must add an explicit stationarity gate for the declared response, use adaptive/longer finite trajectories, and refuse to rank a scenario until that gate passes.

### 2. The sensitivity table is censored, not a complete ranking

The current table omits undefined perturbation responses when computing a family's maximum and then assigns all seven ordinary ranks. Adsorption `Ea + 3` and Al-O-Al `Ea - 3` have nonpositive combined propensity and undefined log responses. Those omitted values can change the apparent top order.

[`sensitivity-ranking.csv`](a9-approximate-rate-closure/derived/sensitivity-ranking.csv) is retained only as preliminary diagnostic evidence. The analyzer must represent censored families as unrankable or with defensible bounds; `sensitivity_complete` must not mean merely that seven names and ordinal rows exist.

### 3. The provisional lab-rate comparison used the wrong physical observable

All eight nominal runs observed zero desorption and net Si uptake. Independent event tracing found that transitions into desorption-eligible `Si.oh4` arose from `adsorb-si`, not hydrolysis of original lattice Si. The current expected desorption propensity can therefore measure possible return of newly adsorbed solution Si rather than mineral dissolution.

No laboratory dissolution-rate gap is accepted from this campaign. The next deck/estimator must distinguish original lattice cations from reservoir-derived cations—by a tracer/origin state contract or an adsorption-free numerical-sink design that Petra can validate—before comparing to the pH 3-5 Palandri-Kharaka/Brantley ledger.

## Reservoir and provenance

The current deck fixes `T=298 K`, `activity(Si)=activity(Al)=1e-12`, and `mu(Si)=mu(Al)=-1 kcal/mol`, giving an effective consumed-cation factor of `1.847676567496443e-13`. The generated provenance table records those exact choices, derivation, and limitation. They are an open-flow numerical product sink, **not** a calibrated pH-specific chemical potential; Petra has no H+/OH- coupling in this deck.

Direct survey anchors are `27.019 kcal/mol` from `qm/SURVEY.md` CALC-002 and `32.221 kcal/mol` from CALC-003, both B3LYP/def2-SVP/DF activation free energies—not r2SCAN-3c. That correction is encoded in the deck/registry and is not an acceptance blocker.

## Preserved evidence

External raw root:

```text
/mnt/data/vsletten/run-outputs/dissertation/A9-approximate-rate-closure-20260913/raw
```

- raw files: `1,191`
- raw bytes: `348,737,069`
- source deck SHA-256: `169b0b4cea9dafb00e6becfcf2ea36ca5cfaef3cc772b46e2b6c8402b13e86d6`
- Petra binary SHA-256: `bbedccd9b53beef113ec73cbb169fc3917c132eae93c3e51ba4bd5b45ae44043`
- raw manifest SHA-256: `30db1753f5da41149ae886318512fd1fbef096e1f1ce514e548c294535e46536`
- raw checkpoint SHA-256: `48c46e3d5aa594b903e4c11c60fd85e72dee433cffe8663eb31d419723cae2a7`
- preliminary verification SHA-256: `4fa7f908a82b895dc903cce20efe853fcece71221db3a6f4e63308530e7e8c94`

The committed evidence is split so the strict verifier sees only its exact derived inventory:

- [`derived/`](a9-approximate-rate-closure/derived/) — nine byte-verified generated products;
- [`run-receipts/`](a9-approximate-rate-closure/run-receipts/) — manifest, checkpoint, and execution logs.

The preliminary bundle currently reproduces byte-for-byte:

```bash
python3 petra/scripts/approximate_rate_closure.py verify \
  /mnt/data/vsletten/run-outputs/dissertation/A9-approximate-rate-closure-20260913/raw \
  docs/program/results/a9-approximate-rate-closure/derived
```

Expected terminal line:

```text
verified: scenarios=29 replicas=8 outcome=no-dissolution acceptance_passed=False
```

That verifies artifact integrity and the typed negative event outcome. It does **not** satisfy the missing scientific gates above.

## Continuation contract

1. Introduce a physically explicit lattice-origin dissolution observable; prevent adsorbed reservoir cations from masquerading as dissolved mineral.
2. Gate both population and declared-response stationarity, then run an adaptive but finite same-seed ensemble until every scenario is steady or explicitly censored.
3. Replace forced ordinal ranks with censor-aware ordering/bounds and require scientific completeness, not row-count completeness.
4. Only then compute the laboratory order-of-magnitude gap, select A2/A3 targets, update program statuses, and open the single A9 PR.

No QM, GPU lease, service mutation, or credential work occurred in this preliminary pass.
