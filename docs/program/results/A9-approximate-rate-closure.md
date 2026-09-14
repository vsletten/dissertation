# A9 — Approximate Kaolinite Rate Closure

**Result:** 2026-09-13 20:52 PDT (hermes-macbot-zero; profile=laptop) — **DONE AS AN HONEST PARTIAL / NEGATIVE SURVEY RESULT BY SCOPE RULING; NOT CALIBRATED KINETICS**. The historical 29-scenario × 8-seed campaign observed zero original-lattice Si and Al releases, but every sensitivity family is censored, ranked-response stationarity is not established under the corrected gates, and no ordinal ranking exists.

## Recommendation

Do **not** spend the next QM increment on an A9 barrier family from this result. The data can rank no family, show no family irrelevant at this tier, and support no final A2/A3 top-k; the temporary QM target set is empty. The three review findings are executable follow-up cards: `A9b-mechanism-reachability`, `A9b-reservoir-origin-contract`, and blocked `A9b-sensitivity-ranking` after both merge.

This is useful negative evidence and closes A9's current branch/PR, not a calibrated absolute-rate closure. The next ranking must first establish a finite reachability/sampling method and a defensible pH 3–5 reservoir/origin contract.

## Rationale

The corrected observable follows cation lineage by lattice site. An Al or Si release contributes to the physical dissolution numerator only if that site still contains its original lattice cation. Adsorbed reservoir cations are never relabeled as lattice material. Expected CTMC desorption propensity is used only for trajectories with zero adsorption events; otherwise it is censored as origin-contaminated.

The campaign imposed a far-from-equilibrium numerical product sink (`activity(Al) = activity(Si) = 1e-30`) and then verified zero adsorption events. That makes every realized occupied cation site origin-safe. Even under that favorable open-flow condition, no sampled observation or recorded event transition populated a state with nonzero Si or Al desorption propensity. That is finite non-observation, not a topology proof.

## Evidence

### Campaign coverage and verification

- 29 scenarios: nominal plus four perturbations for each of seven reaction families.
- 8 fixed seeds per scenario: 232 trajectories total.
- 200,000 events per trajectory: 46,400,000 KMC events total.
- 21 equal-cadence samples per trajectory (`report_every = 10,000`).
- 298.0 K; Arrhenius prefactor `1e13 s^-1`; barriers and provenance are recorded in `derived/provenance-conversions.csv`.
- Verifier result: `outcome=no-dissolution`, `sensitivity_complete=True`, `ordinal_ranking_complete=False`, `acceptance_passed=False`.
- Gate census: 232 `steady-zero`; 0 `steady-positive`; 0 `nonsteady`; 0 `absorbed`; 0 `incomplete`.
- Origin audit across all trajectories: 0 Si adsorption, 0 Al adsorption, 0 original-lattice Si release, 0 original-lattice Al release.
- Raw campaign: `/Volumes/DATA/hermes/run-outputs/TASK-309-A9-20260913-1740` (2.720523 GiB), with 232 hash-addressed receipts.

### Absolute rate result

| Quantity | Result |
|---|---:|
| Nominal original-lattice Si release rate | `0 mol m^-2 s^-1` |
| Nominal original-lattice Al release rate | `0 mol m^-2 s^-1` |
| Pooled nominal 95% zero-event upper bound, each species | `7.3623390865168e-12 mol m^-2 s^-1` |
| Nominal expected origin-safe Si desorption propensity flux | `0 mol m^-2 s^-1` |
| Nominal expected origin-safe Al desorption propensity flux | `0 mol m^-2 s^-1` |
| Si:Al release stoichiometry | undefined (0:0) |

The zero-event upper bound is a detection limit, not an estimated positive dissolution rate. The two zero expected-propensity fluxes are historical diagnostics from the origin-safe bundle, not validated absolute rates: ranked-response stationarity was not established under the corrected complete-event/state-distribution gates, and the corrected real trajectory is `mechanism-unsampled`.

### Comparison to the literature/lab discussion

Using the kaolinite reference law in `kinetics-db/minerals/kaolinite.toml`, the acid-plus-neutral 298 K laboratory rate is approximately:

- pH 3: `8.8925e-14 mol m^-2 s^-1`
- pH 4: `6.9889e-14 mol m^-2 s^-1`
- pH 5: `6.6708e-14 mol m^-2 s^-1`

The A9 pooled zero-event upper bound is 82.8–110.4 times those values (1.92–2.04 log10 units higher). Therefore the finite campaign does **not** falsify laboratory rates: its event-level detection limit is too loose. Its sampled CTMC desorption hazard also remained zero, but that is a mechanism-sampling diagnostic rather than an absolute-rate comparison. The laboratory gap is not calibrated because the `1e-30` dissolved-cation sink is not a defensible pH 3–5 reservoir and the historical bundle cannot prove the corrected origin/reachability contract.

The result also fails to reproduce the 1999 golden runs' nonzero cation-population evolution. A9's original-lattice Si and Al inventories remained unchanged in every scenario. The current absolute-rate deck is thus not yet phenomenologically equivalent to the legacy model.

### Sensitivity result

All seven families are **censored / unrankable**:

1. `siloxane-neutral`
2. `sioal-si-neutral`
3. `sioal-al-neutral`
4. `connectivity-ladder`
5. `al-o-al-analogue`
6. `adsorption`
7. `cation-desorption`

For every ±3 kcal/mol and ×/÷10 perturbation, the historical combined origin-safe Si+Al propensity response was zero in all eight same-seed pairs. Consequently `delta log10(rate)` is undefined, and assigning an ordinal rank would manufacture information. `sensitivity_complete=True` means only that every family/perturbation received a typed estimate-or-censor classification; `ordinal_ranking_complete=False` remains the scientific verdict.

**Ranking verdict:** this data can rank **none** of the seven families and cannot demonstrate **any** family irrelevant at this tier. There is no defensible top-k. The temporary QM target set is empty because the evidence is unresolved, not because seven families were shown insensitive.

## Trust and method corrections

This continuation closes two implementation defects in the preliminary branch but does not close the scientific card:

1. **Origin-safe observable.** Event replay tracks each cation site's lattice/reservoir lineage. Gross `desorb-*` counts remain diagnostic only; physical rates and stoichiometry use original-lattice release.
2. **Explicit convergence gates.** Every trajectory must have 21 samples, four tail blocks, stable Si and Al populations, and stationary observed-release plus sampled-propensity responses. A stationary complete zero is a typed `steady-zero` outcome with a Poisson upper bound; it is not silently treated as a positive estimate.

Generated output hashing now follows the declared file contract rather than filesystem enumeration, so macOS AppleDouble sidecars on external storage cannot poison declared hashes. Verification still rejects undeclared files in the committed derived bundle.

### Diagnostic-gate correction (implementation verified; campaign not rerun here)

The runner now requires one event row for every contiguous step from 1 through the deck's configured final step, replays the complete state-count vector from the step-0 populations row, validates each event's old state and every cadence population vector, and binds population, observable, and event times exactly at every sampled step. Recomputed artifact hashes therefore cannot launder a dropped event into an analysis. It records separate lattice/reservoir lineage diagnostics for event-level `Si.oh4` and `Al.l6` eligibility. A trajectory with zero nonfrozen original-lattice target eligibility is now typed `mechanism-unsampled` / `unresolved`; it cannot pass as steady zero. Population stationarity now covers every state fraction within every kind, including oxygen and empty states, rather than only conserved Si/Al totals.

When `snapshot.pgif.json` asserts `pgif: 1`, its production node/edge/meta structure, coordinates, categorical state/type/kind mappings, frozen flags, and per-row consistency are validated before its per-site state vector is cross-checked against the step-0 populations row. Only nonfrozen Si/Al centers count toward initial desorption-target eligibility. If the artifact is absent or does not assert PGIF, the step-0 populations row is authoritative: per-site state and nonfrozen eligibility are inferred only when that site first appears in a fired event, while aggregate old-state availability and every emitted cadence vector are still checked exactly. That fallback cannot independently prove the unobserved step-0 site-to-state assignment; the trajectory diagnostic records this limitation instead of inventing one. A malformed document that asserts `"pgif": 1` is rejected rather than downgraded to fallback.

The existing 232-run derived bundle predates this gate schema and was not regenerated in this correction slice because its raw campaign is retained at the external macOS path named above, not in this worktree. Its historical `steady-zero` census must therefore not be treated as acceptance evidence under the corrected runner. A focused real 200,000-event Petra contract run exercises the new path and is classified `mechanism-unsampled`, with both target maxima zero and the PGIF snapshot cross-check verified.

**DEVIATION:** The deck uses B3LYP/def2-SVP/DF survey anchors where the card text labels the 27.0/32.2 kcal/mol values r2SCAN-3c. The provenance table records the actual method and explicitly rejects the card's stale method label.

**DEVIATION:** The `1e-30` dissolved-cation activity is a numerical open-flow sink, not a realistic pH 3–5 solution model. It makes origin-safe attribution possible without engine-level particle identity but precludes accepting the adsorption-family sensitivity requested by the card.

## Tradeoffs

- The `1e-30` product activity is an explicit numerical open-flow sink, not a measured pH-dependent activity. It isolates dissolution from reservoir re-entry but does not constitute calibrated solution chemistry.
- The finite 200,000-event cap is enough to diagnose a mechanism-unsampled recrossing trap under the corrected gate, not to classify physical dissolution as stationary zero; the historical observed-event upper bound is therefore non-acceptance telemetry rather than a laboratory-rate estimate.
- A longer run could tighten the Poisson bound, but it would not resolve the more important zero-hazard topology problem unless the mechanism begins reaching desorption-eligible states.
- Survey barriers mix computed, literature-analogue, and heuristic values. The result must not be presented as production kinetics.

## Follow-up ownership

A9 closes with this preserved partial result. Its three unresolved review findings are executable cards rather than more gates on this branch:

1. READY `A9b-mechanism-reachability` proves lattice-cation paths to release and selects a finite rare-event/stiffness method without manufacturing rates.
2. READY `A9b-reservoir-origin-contract` replaces the `1e-30` numerical sink with a cited pH 3–5 activity/origin model and adversarially proves reservoir re-entry cannot count as lattice dissolution.
3. BLOCKED `A9b-sensitivity-ranking` reruns the 29 × 8 design after both contracts merge and owns the final censor-aware ranking or independently verified NO-GO.

Until the ranking card closes, A2/A3 should spend zero higher-level QM slots on A9 families. The historical all-censored result is evidence, not a ranking.

## Artifacts

- `docs/program/results/a9-approximate-rate-closure/derived/ensemble-rates.csv`
- `docs/program/results/a9-approximate-rate-closure/derived/per-replica-rates.csv`
- `docs/program/results/a9-approximate-rate-closure/derived/sensitivity-ranking.csv`
- `docs/program/results/a9-approximate-rate-closure/derived/steady-state-gates.json`
- `docs/program/results/a9-approximate-rate-closure/derived/state-populations.csv`
- `docs/program/results/a9-approximate-rate-closure/derived/stoichiometry.csv`
- `docs/program/results/a9-approximate-rate-closure/derived/results-dat-equivalent.csv`
- `docs/program/results/a9-approximate-rate-closure/derived/provenance-conversions.csv`
- `docs/program/results/a9-approximate-rate-closure/derived/verification.json`
- `docs/program/results/a9-approximate-rate-closure/run-receipts/manifest.json`
- `docs/program/results/a9-approximate-rate-closure/run-receipts/checkpoint.json`
