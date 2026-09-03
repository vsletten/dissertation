# A8c-thesis-table-4-10c-erratum — resolve the archived-vs-printed hydrolysis energy

- status: ready
- track: A (geochemistry / archaeology)
- priority: P2
- machine: any
- depends: A8-thesis-archive-intake
- claimed-by: —

## Objective

Resolve one concrete provenance discrepancy found by A8. The archived Gaussian zero-point energies reproduce Chapter 4 Table 4.10(c) as approximately **185.6 kcal/mol**, while the thesis prints **183.6 kcal/mol**.

1. Recalculate the stoichiometric sum independently from the exact logs named by `ch4-dimer-hydrolysis-c` in `DFT-LEDGER.csv`.
2. Search archived notes, TeX history, Gaussian inputs, and RCS/context for an alternate calculation that supports 183.6.
3. Determine whether this is a transcription error, a superseded calculation, or a mapping defect.
4. Write a short cited erratum/provenance note; preserve both claims if evidence is not decisive.

## Acceptance

- Independent arithmetic receipt in Hartree, kJ/mol, and kcal/mol.
- Every source log identified by relative path and SHA-256.
- Clear verdict with confidence and no silent alteration of the archived thesis or A8 ledger evidence.
- STATUS.md and this card updated in the delivering PR.

## Progress

- 2026-09-02 — card filed from A8; 44 of 45 mapped Chapter 4 energies agree within rounding, with Table 4.10(c) the sole 2.0 kcal/mol outlier. (hermes-custom-build-001; profile=workstation)
