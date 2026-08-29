# D2b-explicit-surface-rates — surface-valid astrochemical LH rates

- status: done
- track: D (astrochemistry)
- priority: P1
- machine: workstation
- depends: D2a
- claimed-by: fable (hermes-custom-build-001; profile=workstation)

## Objective

Resolve D2a's explicit NO-GO for importing cheap gas-phase rates into a grain
surface deck. Build and validate the smallest defensible explicit-surface rate
protocol (constrained ice cluster or equivalent adsorbate treatment, surface
partition functions, and site sampling) for H + CO and the first formaldehyde
addition/abstraction branches. Compare against published instanton/experimental
surface rates and state a validity floor.

## Acceptance

- At least H + CO and H + H2CO surface channels have provenance-carrying
  barriers/rates and quantitative literature comparisons at 3+ temperatures.
- The 4–5-order D2a discrepancy is either closed to within the documented
  literature spread or the protocol remains an explicit NO-GO with receipts.
- A machine-readable rate table safe for D3b exists only if the gate passes.
- All long runs are bounded, GPU-first, thread-capped, and teed to durable logs.

## Progress

- 2026-08-22 — created by D3 after refusing to treat D2a gas-phase rates as
  Langmuir–Hinshelwood surface rates.
- 2026-08-27 — fable (hermes-custom-build-001; profile=workstation) — claimed;
  protocol designed and driver landed. Key findings during design, from D2a
  artifacts re-rated in place: (1) the field's implicit-surface partition
  functions (vibration-only, Meisner/Lamberts/Kästner convention) move the
  4–5-order addition gap by only ~2× — the dominant term is *barrier error*
  at pwb6k/def2-SVP (addition 2273 K vs lit 1900 K → 10⁴ slow; abstraction
  2249 K vs 3030 K → 10² fast); (2) the literature itself spans ~4 orders at
  60 K — Simons, Lamberts & Cuppen 2020 (A&A 634 A52) document that the Song
  & Kästner instanton fits D2a compared against are overestimates (DFT
  barriers low vs CCSD(T)-F12, worst for abstraction) and adopt plateau rates
  of 2·10⁵ s⁻¹ (H+CO) / 1.5·10⁵ s⁻¹ (H+H2CO channels) in their
  experiment-benchmarked KMC. D2a's gas-phase 60 K addition rate is within
  ~5× of that adopted value. Protocol therefore: rung A surface PFs +
  rung B UHF-UCCSD(T)/CBS barrier correction (A2 engines are closed-shell
  only, so open-shell CC added in-driver) + rung C explicit (H2O)₁,₂ clusters
  with orientation site sampling. Driver:
  `qm/scripts/surface_rate_protocol.py` (10 reactions: 3 gas + 7 wet sites);
  anchors embedded with provenance: Song & Kästner 2017 fits, Simons 2020
  adopted plateaus, Fuchs 2009 effective barriers (12–16.5 K, own
  convention), Andersson 2011 T_c = 79 K.
- 2026-08-27 (evening) — campaign launched on the 4090 (driver commit
  `394fa00`): GPU-first, threads=16, teed to
  `qm/runs/D2b-campaign{,-stdout}.log` (worktree
  `agents/D2b-explicit-surface-rates`, gitignored). GPU etiquette
  applied per quarry-campaign skill: both email pipeline services +
  honcho deriver/api stopped, `ollama ps` verified empty, 12 h
  dead-man restart armed (`d2b-deadman-restart.timer`, fires 08:30).
  qm suite green at launch (455 passed, 1 skipped); CPU smoke of the
  driver at STO-3G exercised scan→Sella→IRC→gates→failure-isolation
  end-to-end.
- 2026-08-27 (night) — pass 1: 5/10 clean; four campaign-hardening
  fixes landed from live failures (`d1d7cb6` shallow-crest fallback,
  `f60c810` name-insensitive CC receipt hashes, `c311c3e`/`335e9c9`
  tiered spectator-mode gates for wet sites, `87ea5eb` wet-seed
  pre-optimization). Catch-up pass running detached (pid 3010250 on
  the 4090): 9/10 reactions done or resumed clean as of 22:52; final
  reaction h-h2co-h2-hco-2w scanning from its pre-optimized seed.
  REMAINING FOR CLOSEOUT once `complete:` appears in
  `qm/runs/D2b-campaign-stdout.log`: analyze
  `qm/runs/D2b-explicit-surface-rates/results.json` against the gate
  (rungs A/B/C decomposition vs Song&Kästner fits / Simons plateaus /
  Fuchs effective barriers), write
  `docs/program/results/D2b-explicit-surface-rates.md` with the
  go/no-go verdict, emit the D3b rate table ONLY if the gate passes,
  card→done, STATUS.md line, open PR and STOP. Then restore GPU
  workloads (email pipelines + honcho deriver/api — dead-man
  `d2b-deadman-restart.timer` fires 08:30 regardless) and stop the
  `d2b` monitor.
- 2026-08-28 — fable (hermes-custom-build-001; profile=workstation) —
  campaign complete: 10/10 reactions, 0 failures (23:05 the night
  before). GPU workloads restored (dead-man revived the email
  pipelines at 08:30; honcho containers restarted by hand at
  closeout). Verdict written: docs/program/results/
  D2b-explicit-surface-rates.md. Card closed.

## Result

**Split verdict, receipts in the results doc.** NO-GO for the D3b
rate table: at 10–20 K the addition and H+CO channels sit 25–1600×
below the Simons-corrected literature anchor with a ~60× site spread,
and the implied branching contradicts the experiment-benchmarked
network — per acceptance, no machine-readable table ships. GO for the
protocol tier: the 1-water CC abstraction barrier reproduces Song &
Kästner's published surface value to 2 K (3028 vs 3030 K), the
abstraction channel's k(T) lands inside the documented spread at every
anchor, and D2a's "4–5 orders" is fully decomposed (≈2× partition
functions + ≤10² barrier error + up to 10³·⁵ *anchor overestimate*
documented by Simons+20 + 1.5–3 orders Eckart-vs-instanton in deep
tunneling, independently confirmed by the T_c mismatch on H+CO: 149 K
at PWB6K width vs Andersson's 79 K). Validity floor stated: ≥100 K
≲1 order; 59–100 K ≲2 orders; <59 K not publishable as absolute
surface rates without instanton-class tunneling. Follow-up:
D2c-instanton-tier (READY) is the single upgrade that flips the gate.
Surprises for the next worker: the anchor itself carried much of
D2a's discrepancy; water moves barriers ±450 K by site; four driver
fixes landed from live failures (shallow-crest fallback, name-blind
geometry hashes, tiered spectator-mode gates, wet-seed pre-opt).
