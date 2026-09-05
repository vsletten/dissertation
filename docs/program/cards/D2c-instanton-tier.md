# D2c-instanton-tier — deep-tunneling rates that flip the D2b gate

- status: ready
- track: D (astrochemistry)
- priority: P2
- machine: workstation
- depends: D2b
- claimed-by:

## Objective

D2b's NO-GO residual is isolated to tunneling method: asymmetric
Eckart is 1.5–3 orders low against instanton-class anchors below
~59 K, because it inherits the barrier *width* from a single imaginary
frequency (and PWB6K's width is demonstrably too stiff — T_c 149 K vs
Andersson's 79 K for H+CO). Add one instanton-grade tunneling tier to
quarry: small-curvature tunneling (SCT) along the existing
`full_irc`/`quick_irc` path machinery, or a ring-polymer instanton on
the driver's stationary points. Benchmark reaction-by-reaction against
Song & Kästner 2017 (H2CO+H channels) and Andersson 2011 / Simons 2020
(H+CO), reusing the D2b campaign's checkpointed geometries and CC
barriers (`qm/runs/D2b-explicit-surface-rates/`, campaign worktree).

## Acceptance

- A tunneling correction beyond Eckart implemented against the
  existing IRC machinery, unit-tested, with the same receipt
  discipline as D2b.
- H+CO and both H2CO channels within the documented literature spread
  at 12–20 K (the Simons plateau ±1 order), with the branching ratio
  defensible against the experiment-benchmarked network.
- The D2b validity floor restated with the new tier; if the gate now
  passes, emit the machine-readable D3b rate table D2b withheld —
  that table unblocks D3b.
- Bounded, GPU-first where applicable, thread-capped, teed logs.

## Progress

- 2026-08-28 — created by D2b closeout: the campaign validated the
  barrier tier (1w abstraction barrier = S&K's surface value to 2 K)
  and named Eckart-vs-instanton as the sole remaining gap.
