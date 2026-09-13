# A9-approximate-rate-closure — run the kaolinite deck with approximate absolute rates

- status: ready
- track: A (geochemistry)
- priority: P0
- machine: any
- depends: —
- claimed-by:

## Objective

Objective. The deck (`petra/examples/kaolinite.toml`) still carries the 1999
relative (k, dE) pairs at T = 8000 K; A5p1 ran it with a fixed 1000x
defect/terrace contrast. Close the loop once with real-shaped rates:

1. Build `petra/examples/kaolinite-approx.toml`: replace every relative pair
   with an `eyring` (or `arrhenius`, A = 1e13 s^-1) rate at 298 K from the best
   AVAILABLE number, each tagged by provenance class in a comment —
   `computed` (si-neutral 27.0 kcal/mol, al-neutral 32.2; both r2SCAN-3c
   survey tier, from qm/), `literature` (Xiao & Lasaga 1994/1996 acid ~24 /
   base ~19; Pelmenschikov 2000/2001; Criscenti 2006; Nangia & Garrison 2008;
   Liu & Ruiz Pestana 2024 Q1/Q2/Q3 = 54/71/81 kJ/mol for the connectivity
   ladder shape), `heuristic` (legacy 6/12 kcal-per-bucket attach/detach and
   legacy ratios rescaled to an absolute anchor). Realistic delta-mu for
   far-from-equilibrium pH 3-5 dissolution. Every choice goes in a provenance
   table in the results doc.
2. Run it: >= 8 replicas to steady state; emit dissolution rate (mol Si and
   Al per m^2 per s via A5p0's surface accounting), Si:Al stoichiometry,
   site-population evolution, and a results.dat-equivalent series.
3. Compare: (a) qualitatively to the 1999 golden runs (A8 archive) — same
   phenomenology or not; (b) to measured far-from-equilibrium kaolinite lab
   rates from the ledger in `docs/scoping/field-lab-discrepancy.md` /
   Brantley's review — state the order-of-magnitude gap honestly.
4. Sensitivity: one-at-a-time +/-3 kcal/mol on each barrier family (and x/÷10
   on prefactors), same ensemble; rank families by |delta log10 rate|.

## Acceptance

Acceptance. `docs/program/results/A9-approximate-rate-closure.md` with the
provenance table, rate + stoichiometry with ensemble bands, the lab-rate
comparison, and the ranked sensitivity table; deck + runner script
committed; ONE PR. No QM runs on this card. Verdict names the top-k barrier
families A2/A3 may spend GPU on and the ones that are irrelevant at this
tier. Bounded CPU compute (`systemd-run --user -p RuntimeMaxSec`) with a
`file:` receipt if any ensemble exceeds 30 min.

## Progress

- 2026-09-13 12:04 PDT (hermes-custom-build-001; profile=workstation) — Card created verbatim from Victor's TASK-274 ruling. A9 is READY at P0 and machine-any; the board feeder owns its mission-control pointer after this PR merges.
