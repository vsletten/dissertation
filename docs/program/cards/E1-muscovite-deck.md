# E1-muscovite-deck — Ar release from muscovite, phase 1 (published numbers)

- status: ready
- track: E (muscovite / the perfect circle)
- priority: P1
- machine: any (KMC is CPU-light; no QM in this phase)
- depends: — (phase 1 is isothermal; staircase spectra later need B5)
- claimed-by:

## Objective
Phase 1 of docs/scoping/ar-muscovite.md §6: the minimal muscovite deck
using published barriers only. Single-gallery 2M₁ reduction — K/vacancy/
Ar40 interlayer gallery with divacancy hops (Em ≈ 66 kcal/mol, Nteme
et al. 2022), paired OH sites with bulk dehydroxylation kinetics
(by_count steric table raising Ea with neighboring O_residual —
Guggenheim mechanism), open c-axis boundary with fast Ar desorption at
the surface. Isothermal runs at 500 °C and 700 °C (the Sletten &
Onstott 1998 in-vacuo conditions). Post-process the event log into
cumulative-release curves and D/a²(t) via the same Fickian inversion
the 1998 paper used.

Read the scoping doc first — it is the design authority for site
kinds, rules, and barrier provenance (§4 table).

## Acceptance
- Deck under petra/decks/ validates and runs on current main (v1
  schema is fine if v2 is not yet merged; note which).
- Release curves at both temperatures show the two-stage non-Fickian
  shape; D/a²(t) inversion shows the rise-then-fall signature
  (scoping §5 claim 1 — qualitative gate, not a fit).
- Results doc in docs/program/results/E1-muscovite-phase1.md with
  plots/data pointers and an honest gap list for phase 2.

## Progress
- 2026-08-19 — card created (fable); scoping doc landed same PR.
