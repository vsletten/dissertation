# Table 4.10(c) hydrolysis-energy erratum

## Verdict

**High-confidence transcription error.** The archived Gaussian frequency logs and the archived reaction worksheet independently give **185.6114 kcal/mol** for Chapter 4 Table 4.10(c). The value first changes to **183.6114 kcal/mol** in the hand-authored DFT table source and is then rounded to **183.6 kcal/mol** in the dissertation source and printed thesis. No alternate Gaussian calculation supporting 183.6 was found in the available archive evidence.

This note records provenance only. It does **not** alter the archived dissertation, its rebuilt or archived PDFs, the source archive, or `DFT-LEDGER.csv`.

## Independent arithmetic receipt

The ledger reaction is:

```text
h7sialo7-h2o + h2o -> h4sio4 + h4alo4 + h3o
```

The final normally terminated `Sum of electronic and zero-point Energies` value from each ledger-named source log is:

| sign | source log (relative to the recovered thesis archive) | SHA-256 | electronic + ZPE (Eh) |
|---:|---|---|---:|
| − | `dft/results/h7sialo7+h2o/sialhw-b3lyp.freq.log` | `6d8c437403df279ccd6f1bb0540c55828a9a6c29ab85dba378525fe75fe708ff` | −1139.390135 |
| − | `dft/results/h2o/h2o.freq.log` | `bd166f0407a42766664caab1d8a9a8dd8eff2186d3e76de3adb6d6a98a4c3ffd` | −76.387600 |
| + | `dft/results/h4sio4/h4sio4.freq.log` | `da17b1385c4706cf92cb8f2767d77d90738c9dbbaa0fc0f0c3d14fc9d3ddd5ff` | −592.911551 |
| + | `dft/results/h4alo4/h4alo4.freq.log` | `ac21baefe1d92578f63f5db1cb76623c2f684830127c2bf45c0d8820969cb634` | −545.917743 |
| + | `dft/results/h3o/h3o.freq.log` | `6a576cc077e057613cba25c06b6e430962ba36a08244df36802307d1d22b7d31` | −76.652646 |

```text
products  = -592.911551 - 545.917743 - 76.652646
          = -1215.481940 Eh
reactants = -1139.390135 - 76.387600
          = -1215.777735 Eh
ΔE        = products - reactants
          = +0.295795 Eh
          = +776.609666 kJ/mol   (2625.4996394799 kJ mol⁻¹ Eh⁻¹)
          = +185.614165 kcal/mol (627.5094740631 kcal mol⁻¹ Eh⁻¹)
```

Rounded as the dissertation tables are rounded, the result is **185.6 kcal/mol**, matching the curated `ch4-dimer-hydrolysis-c` row in `DFT-LEDGER.csv`.

## Provenance chain

The recovered archive preserves both the correct worksheet value and the later printed value:

1. `dft/txt/zpe.txt` contains the same five zero-point-corrected molecular energies as the source logs.
2. `dft/txt/rxns.txt:10` records the reaction as `189.7729` kcal/mol electronic and **`185.6114` kcal/mol electronic + ZPE**.
3. `dft/tex/tbls.tex:28` changes only this table entry to **`183.6114`**. Its generated `tbls.dvi` and `tbls.ps` retain the same number.
4. `text/ch4.tex:955` prints **`183.6`** in the table, and `text/ch4.tex:1009` repeats **`183.6`** in the discussion. The archived thesis and clean rebuild therefore faithfully reproduce the source's later value; they are not the origin of the discrepancy.

The MATLAB utilities establish the intended arithmetic but not a preserved execution history: `dft/txt/CalcZpe.m` reads `zpe.txt` and calculates products minus reactants, while `CalcRxn.m` reads `nrg.txt` for the uncorrected electronic calculation. The archive contains no transcript proving which utility invocation wrote `rxns.txt`, so this note does not claim one.

## Alternate-calculation search and limits

A bounded provenance scan selected 136 archive entries covering:

- all Gaussian `.com` input decks under `dft/results/`;
- `dft/scripts/`, `dft/tex/`, `dft/txt/`, and `dft/work/` context;
- text-like dissertation and miscellaneous-source files.

Of these, 135 were readable and matched the SHA-256 recorded in `ARCHIVE-INVENTORY.csv`. The only unavailable entry was the editor-backup path `dft/txt/org.txt~`; its live `org.txt` counterpart was present and hash-matched. The archive inventory contains no DFT RCS file. The current Git history begins with the modern curation, so it adds no 1999 revision trail.

Within that evidence set, the exact reaction's correct `185.6114` value appears in `rxns.txt`; the incorrect `183.6114` appears in `tbls.tex` and its generated outputs; and the rounded `183.6` appears in `ch4.tex`. No alternate input deck, worksheet row, energy list, or note supports 183.6 as a separately calculated result. Because one irrelevant editor backup is unavailable and no DFT RCS history survives, the verdict is **high confidence**, not absolute proof of the typist's action.
