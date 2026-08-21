# A8-thesis-archive-intake — catalog and curate the 1999 dissertation archive

- status: ready
- track: A (geochemistry)
- priority: P1
- machine: workstation (the archive lives only at
  /home/vsletten/Documents/personal/science/thesis — 538 MB, 1,426
  files, recovered from Victor's NAS 2026-08-20; READ-ONLY source)
- depends: —
- claimed-by:

## What the archive is
Victor's complete PhD archive: *"Playing Dice with the Universe: Monte
Carlo and Quantum Mechanics in Geochemistry"* (Sletten, ~1999).
Layout:

- `text/` — full thesis LaTeX (ch1–ch5, appA/appB, figures, matlab)
- `mc/model/` — the ORIGINAL C-language KMC model (pre-C++; the
  missing generation 0 of legacy/cpp-model) + its input decks and an
  `8000.out`
- `mc/results/hotrox/` (65 runs) + `mc/results/jasper/` (23 runs) —
  complete July–August 1999 MC campaigns, each timestamped dir
  self-contained: data.{cell,lattice,rxn,sim} + results.dat +
  surfSi.out/surfAl.out + plot.ps
- `dft/results/` — the original Gaussian B3LYP calculations (86 .log
  files, with .com inputs and freq jobs): Si/Al cluster + water
  hydrolysis systems (sialsuw…) AND oxalate/organic-ligand systems
  (c2o4, h2c2o4, h2c3o4, strained variants…)
- `ref/`, `misc/` — reference library and sundries

## Objective
Make the archive durably useful to the program without bloating the
repo:

1. **Catalog**: `legacy/thesis-archive/CATALOG.md` in-repo — full
   annotated inventory (what each dir/run is, dates, sizes, formats),
   plus a provenance note (NAS recovery date, source path). The
   archive itself stays on the workstation; only curated small
   artifacts enter git.
2. **Golden data extraction**: copy IN (small, texty) the original C
   model source (`mc/model/*.c,*.h,Makefile,data.*`) and 3–5
   representative 1999 run dirs (inputs + results.dat + surf*.out —
   no .ps) under `legacy/thesis-archive/`, chosen to span temperature/
   chemistry settings. These become future parity/behavior targets for
   kmc-rs and petra (follow-up card, not this one).
3. **Barrier ledger**: parse every dft/results/*/ Gaussian log:
   system, method/basis, final energy, frequencies if present, and the
   derived barriers where reactant/TS/product triplets exist. Emit
   `legacy/thesis-archive/DFT-LEDGER.md` (+ CSV). This is the
   1999-vs-2026 comparison base for A3/A7 — same author, same
   reactions, 27 years of method development apart. Flag the oxalate/
   ligand systems distinctly: they are a candidate future track
   (ligand-promoted dissolution, ERW-relevant).
4. **Thesis text**: verify the LaTeX builds (or convert cleanly);
   place a built PDF (or the .tex tree if small) under
   `legacy/thesis-archive/text/`. The thesis is the definitive spec of
   the 1999 model's intent — cite it in CATALOG.md by chapter.
5. **Follow-up cards** for what this exposes (at minimum: the 1999-run
   parity gate; the 1999-vs-2026 barrier comparison writeup; a
   possible ligand-promotion track scoping) — file them, don't do them.

## Constraints
- The NAS archive directory is READ-ONLY: never modify, move, or
  delete anything under it. Copy out only.
- Repo big-output policy applies: nothing bulky/binary into git
  (no .chk, no .ps, no full 88-run trees). Curate hard.
- Text encodings may be 1999-era (latin-1); convert deliberately.

## Acceptance
- CATALOG.md + DFT-LEDGER.md (+ CSV) merged; curated artifacts in
  legacy/thesis-archive/ per above; archive source untouched
  (verifiable: no mtime changes under the source dir).
- Follow-up cards filed.

## Progress
- 2026-08-20 — card created (fable), night the archive came off the NAS.
