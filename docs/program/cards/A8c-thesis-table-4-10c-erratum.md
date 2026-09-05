# A8c-thesis-table-4-10c-erratum — resolve the archived-vs-printed hydrolysis energy

- status: done
- track: A (geochemistry / archaeology)
- priority: P2
- machine: any
- depends: A8-thesis-archive-intake
- claimed-by: hermes-macbot-zero

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
- 2026-09-03 — (hermes-macbot-zero; profile=laptop) — Recovered the prior continuation from the shared archive after its local-only commit was unavailable on this host. Re-ran the bounded provenance proof: all five ledger-named source logs rehashed, the stoichiometric sum independently reproduced 185.614165 kcal/mol, and 135/136 selected text/input/context artifacts were readable and inventory-hash exact. The chain localizes `185.6114` in `dft/txt/rxns.txt`, the digit change to `183.6114` in `dft/tex/tbls.tex`, and the rounded `183.6` in `text/ch4.tex`; `dft/txt/org.txt~` was the sole unavailable editor backup.
- 2026-09-03 — (hermes-macbot-zero; profile=laptop) — INDEPENDENT VERIFY PASS: a fresh isolated reviewer started from the original A8c discrepancy, rehashed and re-extracted all five logs, independently recomputed `+0.295795 Eh = +776.60966586 kJ/mol = +185.61416488 kcal/mol`, confirmed the `rxns.txt` → `tbls.tex` → `ch4.tex` provenance break, found no blocking issue in the erratum note, and verified no ledger or archived/rebuilt PDF diff.

## Result

- 2026-09-03 — (hermes-macbot-zero; profile=laptop) — DONE. `legacy/thesis-archive/TABLE-4.10C-ERRATUM.md` records every source-log path and SHA-256, the independent `+0.295795 Eh = +776.609666 kJ/mol = +185.614165 kcal/mol` receipt, the bounded alternate-calculation search, and the high-confidence transcription-error verdict. `CATALOG.md` links the note. The archived thesis/PDFs and A8 ledger remain byte-untouched; both the historical 183.6 claim and the source-supported 185.6 result are preserved.
