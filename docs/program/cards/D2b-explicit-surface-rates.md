# D2b-explicit-surface-rates — surface-valid astrochemical LH rates

- status: in-progress
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
