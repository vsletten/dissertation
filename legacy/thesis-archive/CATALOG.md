# 1999 thesis archive catalog

## Provenance and immutability

The source is Victor Sletten's recovered dissertation archive for *Playing Dice with the Universe: Monte Carlo and Quantum Mechanics in Geochemistry* (circa 1999). Per board card A8, it was recovered from NAS on 2026-08-20 and is treated as **read-only** at:

```text
/home/vsletten/Documents/personal/science/thesis
```

The intake never wrote beneath that path. Before and after curation, the following deterministic metadata manifest was generated:

```bash
LC_ALL=C find "$SRC" -printf '%P\t%y\t%s\t%T@\n' | LC_ALL=C sort
```

| check | value |
|---|---|
| pre-intake manifest SHA-256 | `edf7dea5a690823c9a131d0e28cac5ca68634a0538b3df424dd8c91235198fdb` |
| post-intake manifest SHA-256 | `edf7dea5a690823c9a131d0e28cac5ca68634a0538b3df424dd8c91235198fdb` |
| byte comparison | identical (`cmp` exit 0) |
| source extent | 1,426 files; 155 directories; 559,833,403 bytes |
| source mtime range | 2007-01-13 09:10:42 through 2007-01-13 11:46:56 (local archive timestamps) |

`ARCHIVE-INVENTORY.csv` is the exhaustive file-level inventory: relative path, category, byte count, UTC mtime, SHA-256, and annotation for all 1,426 files. It makes later bit-rot or provenance checks mechanical without committing the 559 MB source tree.

## Archive map

| subtree | files | directories | bytes | annotation |
|---|---:|---:|---:|---|
| `dft/` | 528 | 50 | 513,862,586 | Gaussian inputs, outputs, checkpoints, geometries, wavefunctions, rendered structures, scripts, and chapter support material |
| `mc/` | 820 | 98 | 32,653,764 | Original RCS-tracked C KMC model, input templates, scripts/utilities, and 88 simulation run directories |
| `misc/` | 14 | 0 | 241,977 | Abstracts, status notes, and historical TeX/DVI/PostScript summaries |
| `ref/` | 8 | 0 | 2,496,608 | Literature-search notes and archived reference documents |
| `text/` | 56 | 2 | 10,578,468 | LaTeX dissertation source, bibliography, generated artifacts, figures by absolute reference, MATLAB helper, and a plain C source snapshot |

### DFT subtree

- `dft/results/`: 44 named chemical-system directories, 494 files, and 83 Gaussian `.log` files. Most systems carry some combination of `.com`, `.log`, `.chk`, `.xyz`, `.wfn`, `.eps`, `.jpg`, frequency jobs, and a historical `done` sentinel.
- `dft/old-results/`: four prior-result files retained only in the full inventory.
- `dft/scripts/`, `dft/work/`: four files each, supporting the 1999 Gaussian workflow.
- `dft/tex/`, `dft/txt/`: ten and twelve files respectively, supporting the dissertation narrative and tabulation.
- Binary checkpoints, wavefunctions, rendered PostScript, and the full DFT tree are deliberately **not** copied into Git. Their paths and hashes remain in `ARCHIVE-INVENTORY.csv`.
- `DFT-SOURCE-LOGS.csv` is the exhaustive parser coverage edge: one row for every one of the 83 Gaussian result logs, including seven without a normally terminated energy-bearing job.
- `DFT-LEDGER.csv` contains 44 source-calculation rows plus 45 Chapter 4 reaction-energy mappings. `DFT-LEDGER.md` documents selection and status semantics.

**Scientific finding:** this archive is not a transition-state/barrier archive. Chapter 4 explicitly says the study focused on total reaction energies rather than transition states. No source log can be defensibly assigned as a TS. Barrier columns therefore remain blank and the relevant rows are `unmatched`; no negative frequency was promoted into a TS by filename or frequency count alone.

### Monte Carlo subtree

- `mc/model/`: 72 files. Its root contains the active 1999 Makefile, 28 C headers/sources, four `data.*` inputs, and five generated/runtime artifacts. `RCS/` contains 30 revision-control archives; `input/` contains the four prototype inputs (`prot.cell`, `prot.lattice`, `prot.rxn`, `prot.sim`).
- `text/src/`: a distinct 29-file source snapshot printed through Appendix A. Many files differ from the active `mc/model/` working tree, so both snapshots are curated and never conflated.
- `mc/results/`: 88 run directories across `hotrox/` and `jasper/`. Each is keyed by a Unix epoch-like run ID. `MC-RUN-INVENTORY.csv` records every run's inferred UTC start, thermodynamic inputs, simulation length, seed, completion indicators, output line counts, and bytes.
- Two archived run directories lack one or both `surfAl.out`/`surfSi.out`; this is represented by `required_files_present=false` rather than silently repaired.
- `mc/scripts/`, `mc/util/`, `mc/tex/`, and `mc/txt/` retain historical batch, analysis, model-documentation, and notes material in the source inventory.

### Dissertation text

| source | content |
|---|---|
| `ch1.tex` | Introduction |
| `ch2.tex` | Description of the Monte Carlo Model |
| `ch3.tex` | Monte Carlo Simulation Results |
| `ch4.tex` | Organic Acids and Aluminosilicates I: Reactions |
| `ch5.tex` | Organic Acids and Aluminosilicates II: Bonding |
| `appA.tex` | Monte Carlo Model Source Code |
| `appB.tex` | Minimum Energy Geometries |

The source also contains `thesis.tex`, `thesis.bib`, old generated TeX artifacts, and a MATLAB plotting helper. Figure references use the original absolute `/home/accts/vsletten/thesis` prefix; the clean builder rewrites that prefix only in a disposable copy.

## Curated payload

### C model

`c-model/` contains:

- `source/`: byte-for-byte copies of the active `mc/model/` Makefile, 28 C headers/sources, and four `data.*` inputs (33 files total).
- `appendix-source/`: byte-for-byte copies of the 29-file `text/src/` snapshot used for Appendix A. Differences from `source/` are preserved as historical evidence.
- `RCS/`: byte-for-byte copies of the 30 original RCS archives from `mc/model/RCS/`, preserving revision history and commit-era metadata.
- `input/`: byte-for-byte copies of the four prototype input files from `mc/model/input/`.

No encoding conversion was necessary. The inventory and parser read historical text as ISO-8859-1/Latin-1, which is byte-preserving for this material. Original hashes remain in `ARCHIVE-INVENTORY.csv`.

### Representative KMC golden runs

Only small textual inputs/results were copied; each run's `plot.ps` was excluded. All selected runs have the complete required seven-file input/output set and 5,000 `results.dat` records.

| curated run | UTC run ID | T (K) | Δμ Si | Δμ Al | seed | historical `done` | why selected |
|---|---|---:|---:|---:|---:|---|---|
| `golden-runs/hotrox/935077498` | 1999-08-19 15:44:58 | 5000 | 0 | 0 | -30 | yes | canonical neutral baseline and the only selected run with an explicit writer completion sentinel |
| `golden-runs/hotrox/936930575` | 1999-09-10 02:29:35 | 5000 | -1 | -1 | -2 | no | low-chemical-potential endpoint; pairs with the identical Jasper condition for cross-host determinism checks |
| `golden-runs/hotrox/937172019` | 1999-09-12 21:33:39 | 6000 | -1 | -1 | -2 | no | temperature-sensitivity case at otherwise identical low chemical potentials |
| `golden-runs/jasper/933892971` | 1999-08-05 22:42:51 | 5000 | -1 | -1 | -2 | no | cross-host duplicate of the Hotrox low-potential condition |
| `golden-runs/jasper/935835145` | 1999-08-28 10:12:25 | 5000 | +1 | +1 | -2 | no | high-chemical-potential endpoint, bounding the selected sweep |

The `done` sentinel was inconsistently used in 1999. Absence is not treated as proof of truncation: the four sentinel-free selected runs each contain all required files, 5,000 result rows, and both surface outputs. `MC-RUN-INVENTORY.csv` preserves the distinction.

### Gaussian ledger

Regenerate and validate:

```bash
python3 legacy/thesis-archive/tools/parse_gaussian_archive.py \
  --source /home/vsletten/Documents/personal/science/thesis \
  --out-dir legacy/thesis-archive
python3 legacy/thesis-archive/tools/parse_gaussian_archive.py \
  --source /home/vsletten/Documents/personal/science/thesis \
  --out-dir legacy/thesis-archive --check
```

The parser:

1. enumerates every `dft/results/**/*.log` file;
2. treats explicit Gaussian Link1 sections independently;
3. selects the last normally terminated energy-bearing job in each log;
4. prefers a normally terminated frequency/ZPE log as each system's primary result;
5. records all source logs, including failures, in `DFT-SOURCE-LOGS.csv`;
6. reproduces the 45 reaction energies printed in Chapter 4 Tables 4.4, 4.5, 4.7, 4.8, and 4.10–4.14; and
7. refuses to calculate a barrier without a defensible reactant/TS/product assignment.

The generated reaction energies reproduce the printed rounded values within 2.01 kcal/mol. Forty-four of forty-five agree within rounding tolerance; the archived files reproduce Table 4.10(c) as 185.6 kcal/mol while the thesis prints 183.6 kcal/mol. The ledger preserves and flags the discrepancy instead of coercing the archive to the publication.

Adversarial checks are encoded in the parser validator: exact fixed schema; enumerated status/treatment values; exhaustive source-log row count; no `ok` row without a TS and forward barrier; deterministic regeneration; exact 45-reaction map; and bounded comparison against the independent thesis tables.

### Thesis build

The repository retains both provenance and reproducibility artifacts:

| artifact | pages | SHA-256 | note |
|---|---:|---|---|
| `text/Playing-Dice-with-the-Universe.archived-2007.pdf` | 134 | `9d591ca1210e9c9b3752c33156ac6049aaa6b5daa451d35fa64e91eec6a36fa4` | byte-for-byte archived PDF; Acrobat Distiller 7.0.5 metadata |
| `text/Playing-Dice-with-the-Universe.rebuilt-2026.pdf` | 133 | `624b2a40a70d898cedbe6bfd5fa0fc0e9185f47dcc64a28425847951c54880b0` | clean build from the archived source and referenced figures |

The clean rebuild intentionally starts without archived `.aux`, `.bbl`, `.toc`, `.dvi`, `.ps`, or `.pdf` outputs. It uses the committed `doublespace.sty` compatibility shim backed by modern `setspace`, and rewrites old absolute figure prefixes only inside a temporary directory.

```bash
python3 legacy/thesis-archive/tools/build_thesis.py \
  --source /home/vsletten/Documents/personal/science/thesis \
  --output legacy/thesis-archive/text/Playing-Dice-with-the-Universe.rebuilt-2026.pdf
```

Verified toolchain:

- Python 3.11.15
- GNU findutils 4.9.0
- Latexmk 4.83
- pdfTeX 3.141592653-2.6-1.40.25 (TeX Live 2023/Debian)
- BibTeX 0.99d
- dvips 2023.1
- GPL Ghostscript 10.02.1

The build completed through four LaTeX passes, BibTeX, DVI, and PDF generation. Final output is US-letter PDF 1.7. Remaining diagnostics are historical overfull-box warnings and Helvetica mapping warnings from two EPS figures; there are no fatal errors or unresolved final bibliography pass.

The one-page difference between the archived Distiller PDF and the clean rebuild is retained as provenance evidence, not hidden. The source tree does not contain enough build-environment state to assert byte- or pagination-equivalence across the 2007 and 2026 TeX stacks.

## Reproduction checks

```bash
python3 legacy/thesis-archive/tools/inventory_archive.py \
  --source /home/vsletten/Documents/personal/science/thesis \
  --out-dir legacy/thesis-archive --check

python3 legacy/thesis-archive/tools/parse_gaussian_archive.py \
  --source /home/vsletten/Documents/personal/science/thesis \
  --out-dir legacy/thesis-archive --check
```

`inventory_archive.py` validates all 1,426 source files and all 88 KMC run directories. The two `--check` modes rebuild into temporary directories and byte-compare generated CSV/Markdown outputs against the committed versions.

## Exclusions

The following remain in the read-only archive and are represented by path/hash in `ARCHIVE-INVENTORY.csv`, but are not copied into Git:

- Gaussian `.chk`, `.wfn`, `.grd`, and full result directories;
- bulk EPS/JPEG/PostScript renderings;
- all 88 complete KMC run trees and their plots;
- historical generated TeX intermediates; and
- reference binaries/tarballs.

This keeps the curated tree reviewable while preserving a complete, independently hashable index of the source.
