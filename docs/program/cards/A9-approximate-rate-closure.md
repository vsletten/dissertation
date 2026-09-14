# A9-approximate-rate-closure — run the kaolinite deck with approximate absolute rates

- status: done
- track: A (geochemistry)
- priority: P0
- machine: any
- depends: —
- claimed-by: hermes-macbot-zero

## Objective

Objective. The deck (`petra/examples/kaolinite.toml`) still carries the 1999
relative (k, dE) pairs at T = 8000 K; A5p1 ran it with a fixed 1000x
defect/terrace contrast. Close the loop once with real-shaped rates:

1. Build `petra/examples/kaolinite-approx.toml`: set its execution temperature
   to 298 K (do not copy the legacy `temperature = 8000.0` knob) and add a
   validation/acceptance check that rejects any other temperature; replace
   every relative pair with an `eyring` (or `arrhenius`, A = 1e13 s^-1) rate
   at 298 K from the best AVAILABLE number, each tagged by provenance class
   in a comment — `computed` (si-neutral 27.0 kcal/mol, al-neutral 32.2; both
   r2SCAN-3c survey tier, from qm/), `literature` (Xiao & Lasaga 1994/1996
   acid ~24 / base ~19; Pelmenschikov 2000/2001; Criscenti 2006; Nangia &
   Garrison 2008; Liu & Ruiz Pestana 2024 Q1/Q2/Q3 = 54/71/81 kJ/mol for the
   connectivity ladder shape), `heuristic` (legacy 6/12 kcal-per-bucket
   attach/detach and legacy ratios rescaled to an absolute anchor). Realistic
   delta-mu for far-from-equilibrium pH 3-5 dissolution. Every choice goes in
   a provenance table in the results doc.
2. Run it: >= 8 replicas to steady state; define and implement the
   lattice-to-m^2 conversion (`surface_area.geometric` is Å² from the cell
   matrix; 1 Å² = 1e-20 m²) and the event-to-moles conversion (one dissolved
   Si or Al atom = 1/N_A mol); require both quantities in the provenance
   table before accepting the result; emit dissolution rates (mol Si and Al
   per m^2 per s), Si:Al stoichiometry, site-population evolution, and a
   results.dat-equivalent series.
3. Compare: (a) qualitatively to the 1999 golden runs (A8 archive) — same
   phenomenology or not; (b) to measured far-from-equilibrium kaolinite lab
   rates from the ledger in `docs/scoping/field-lab-discrepancy.md` /
   Brantley's review — state the order-of-magnitude gap honestly.
4. Sensitivity: one-at-a-time +/-3 kcal/mol on each barrier family (and x/÷10
   on prefactors), same ensemble; rank families by |delta log10 rate|.

## Acceptance

**2026-09-13 scope ruling for this closeout:** the original calibrated-rate and
ordinal-ranking gate below is superseded as A9's terminal done criterion. Ship
one review-ready PR containing the 298 K deck/runner, the preserved zero-release
campaign and conversions, an explicit censored/unrankable verdict with
stationarity and reservoir/origin caveats, and executable follow-up cards for
each unresolved finding. No new campaign or review gate belongs on this branch.

Historical acceptance (not scientifically satisfied; transferred to the A9b
follow-ups): `docs/program/results/A9-approximate-rate-closure.md` with the
provenance table (including lattice-to-m^2 and event-to-moles conversions),
rate + stoichiometry with ensemble bands, the lab-rate comparison, and the
ranked sensitivity table; deck + runner script committed; ONE PR. Deck
execution temperature must be 298 K and the runner/acceptance check must
reject any other temperature. No QM runs on this card. Verdict names the
top-k barrier families A2/A3 may spend GPU on and the ones that are
irrelevant at this tier. Bounded CPU compute (`systemd-run --user -p
RuntimeMaxSec`) with a `file:` receipt if any ensemble exceeds 30 min.

## Progress

- 2026-09-13 20:52 PDT (hermes-macbot-zero; profile=laptop) — DONE BY SCOPE RULING / REVIEW-READY CLOSEOUT: no further runner hardening or campaign replay was performed. The branch ships the observed zero original-lattice Si/Al release, zero historical origin-safe expected propensity flux, pooled detection bounds, censored/unrankable seven-family sensitivity table, and explicit caveats: ranked-response stationarity was not established under the corrected gates, the old `1e-30` sink is not a realistic pH 3–5 reservoir, and no final top-k or irrelevant-family ranking is supported. Follow-up cards `A9b-mechanism-reachability`, `A9b-reservoir-origin-contract`, and blocked `A9b-sensitivity-ranking` own those gaps. A9 itself is closed as an honest partial/negative result so the work reaches review instead of accumulating more gates.

- 2026-09-13 20:26 PDT (hermes-custom-build-001; profile=workstation) — CONTINUATION: the pushed diagnostic-gate slice now requires a contiguous complete event stream, validates production-shaped PGIF state/kind/type/frozen/edge data, anchors lattice/reservoir lineage to the initial snapshot, excludes frozen targets, detects transient event-level `Si.oh4`/`Al.l6` eligibility, cross-checks event/population/observable timestamps, gates the complete per-kind state distribution, types unsampled mechanisms as unresolved across per-replica and ensemble products, exits nonzero on failed campaign acceptance, and enforces the aggregate 16-thread cap. Independent cold review's five concrete false-green findings were regression-closed. Fresh QA: 72 Python tests, focused Ruff/format, full Petra workspace tests/check, Rust format, and diff check pass; a real 200,000-event run is correctly rejected as `mechanism-unsampled`. Card acceptance still requires a validated rare-event/stiffness method plus realistic pH 3–5 reservoir/origin contract and a regenerated 29 x 8 ranking; no PR, QM, GPU, live service, database, or credential action occurred.

- 2026-09-13 18:54 PDT (hermes-macbot-one; profile=laptop) — CONTINUATION: corrected origin lineage, response gates, censor-aware ranking, macOS sidecar-safe declared hashing, and the full 29-scenario × 8-seed × 200,000-event campaign are preserved. Verification reproduces `no-dissolution`, `sensitivity_complete=True`, `ordinal_ranking_complete=False`, `acceptance_passed=False`; all seven families are censored and there is no defensible top-k. Cold review also found that 10,000-event sampling cannot prove transient eligibility absence, total-cation population gates do not establish state-distribution stationarity, and the `1e-30` sink does not satisfy realistic pH 3–5 adsorption sensitivity. DEVIATION: the actual 27.0/32.2 kcal/mol survey anchors are B3LYP/def2-SVP/DF, not the card's stale r2SCAN-3c label. Raw data and 232 receipts remain at `/Volumes/DATA/hermes/run-outputs/TASK-309-A9-20260913-1740`. Resume from pushed branch `agents/A9-approximate-rate-closure` with event-level eligibility/reachability, state-distribution convergence, a defensible reservoir/origin contract, then rerun before completion. No QM or GPU work ran.

- 2026-09-13 16:29 PDT (hermes-custom-build-001; profile=workstation) — CONTINUATION: the implementation, 232-run preliminary campaign, strict raw/derived verifier, and condition provenance are preserved, but independent stage-1 review rejected completion. The current 20,000-event windows do not establish stationarity of the ranked propensity response; two perturbations are censored yet receive forced ordinal ranks; and nominal desorption eligibility came from adsorbed solution Si rather than proven lattice-Si hydrolysis, invalidating the provisional laboratory-rate comparison. Resume on this branch/worktree with an origin-aware dissolution observable, response/population stationarity gates, adaptive finite runs, and censor-aware sensitivity ordering. Evidence and exact hashes: `docs/program/results/A9-approximate-rate-closure.md`. No QM or GPU work ran.

- 2026-09-13 12:04 PDT (hermes-custom-build-001; profile=workstation) — Card created verbatim from Victor's TASK-274 ruling. A9 is READY at P0 and machine-any; the board feeder owns its mission-control pointer after this PR merges.

## Result

2026-09-13 20:52 PDT (hermes-macbot-zero; profile=laptop) — **DONE as an honest partial/negative survey result under the 20:50 ship-it ruling; not accepted as calibrated kinetics.** The committed deck is fixed at 298 K and carries source-classed approximate barriers. The historical 29-scenario × 8-seed campaign completed 46.4 million events and observed zero original-lattice Si release and zero original-lattice Al release; its pooled zero-event upper bound is `7.3623390865168e-12 mol m^-2 s^-1` per species. Historical origin-safe expected Si/Al propensity fluxes were also zero, but the corrected complete-event/PGIF/lineage gate classifies a real 200,000-event trajectory `mechanism-unsampled`, so ranked-response stationarity is not established and the old bundle is not acceptance evidence under the final schema.

All seven barrier families are censored/unrankable. The data can rank **no** family, demonstrate **no** family irrelevant at this tier, and support **no** final A2/A3 QM top-k; the temporary QM target set is empty. The laboratory comparison is likewise caveated because the `1e-30` dissolved-cation activity was an origin-isolating numerical sink, not a realistic pH 3–5 reservoir. The preserved result and generated bundle are under `docs/program/results/A9-approximate-rate-closure.md` and `docs/program/results/a9-approximate-rate-closure/`.

Three executable follow-ups own every unresolved review finding: READY `A9b-mechanism-reachability`, READY `A9b-reservoir-origin-contract`, and BLOCKED `A9b-sensitivity-ranking` after both merge. A3, A2, and A2b now block on the ranking follow-up rather than on closed A9.
