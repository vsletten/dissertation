# D2b-explicit-surface-rates — surface-valid astrochemical LH rates

- status: ready
- track: D (astrochemistry)
- priority: P1
- machine: workstation
- depends: D2a
- claimed-by:

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
