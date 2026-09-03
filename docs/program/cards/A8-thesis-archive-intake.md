# A8-thesis-archive-intake — catalog and curate the 1999 dissertation archive

- status: done
- track: A (geochemistry)
- priority: P1
- machine: workstation (the archive lives only at
  /home/vsletten/Documents/personal/science/thesis — 538 MB, 1,426
  files, recovered from Victor's NAS 2026-08-20; READ-ONLY source)
- depends: —
- claimed-by: hermes-custom-build-001

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
3. **Barrier ledger**: parse every `dft/results/*/` Gaussian log:
   system, method/basis, final energy, frequencies if present, and the
   derived barriers where reactant/TS/product triplets exist. Emit
   `legacy/thesis-archive/DFT-LEDGER.md` plus a companion CSV. Use one
   CSV row per reaction candidate (including unmatched calculations,
   with unavailable fields left empty) and this fixed schema:
   `system_id`, `reaction_id`, `ligand_class`, `ref_log`, `ts_log`,
   `product_log`, `method`, `basis`, `energy_treatment`,
   `e_ref_hartree`, `e_ts_hartree`, `e_product_hartree`,
   `barrier_forward_kj_mol`, `barrier_reverse_kj_mol`,
   `reaction_energy_kj_mol`, `imaginary_frequencies_cm_1`,
   `normal_termination`, `parse_status`, `notes`. Source energies are
   hartree and derived energies are kJ mol⁻¹; never emit a unitless
   energy. `energy_treatment` must state `electronic` or
   `electronic+ZPE`, and a row must not mix treatments. The Markdown
   ledger explains triplet assignments and every non-`ok` parse status.
   Encode imaginary frequencies as semicolon-delimited numeric values
   in ascending order, without unit suffixes (the column unit is cm⁻¹).
   An empty field means no imaginary mode was found; `parse_status` is
   one of `ok`, `unmatched`, `truncated`, `ambiguous`, or `parse_error`
   and distinguishes a clean empty result from incomplete parsing.
   Encode `normal_termination` as `true` or `false`.
   This is the 1999-vs-2026 comparison base for A3/A7 — same author,
   same reactions, 27 years of method development apart. Flag the
   oxalate/ligand systems distinctly: they are a candidate future track
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
- Before copying, write a sorted source metadata manifest outside the
  archive; repeat it after all work and require a clean diff:
  ```sh
  SRC=/home/vsletten/Documents/personal/science/thesis
  LC_ALL=C find "$SRC" -printf '%P\t%y\t%s\t%T@\n' |
    LC_ALL=C sort > /tmp/thesis-archive.before
  # Perform read-only discovery and copy-out work here.
  LC_ALL=C find "$SRC" -printf '%P\t%y\t%s\t%T@\n' |
    LC_ALL=C sort > /tmp/thesis-archive.after
  diff -u /tmp/thesis-archive.before /tmp/thesis-archive.after
  ```
  The `SRC` value is the normative source on this workstation; the
  `/tmp` output names are examples but must remain outside the archive.
  Run both manifests on that workstation with GNU `find`, preserve the
  raw `%T@` precision (no rounding), and compare them byte-for-byte;
  cross-platform manifest equality is not an acceptance requirement.
  Store the manifest SHA-256 values and clean-diff result in
  `CATALOG.md`. Do not store snapshots under `$SRC`.
- Repo big-output policy applies: nothing bulky/binary into git
  (no .chk, no .ps, no full 88-run trees). Curate hard.
- Make parsing reproducible: commit the parser, record its invocation
  and the output of `python3 --version` and `find --version` in
  `CATALOG.md`, process each Gaussian Link1 section independently,
  select only normally terminated final jobs, and report
  truncated/ambiguous logs instead of guessing assignments.
- Record the successful thesis build command and versions of the TeX
  engine and `latexmk` in `CATALOG.md`. Prefer `latexmk -pdf` with the
  archive's existing engine assumptions; document any compatibility
  flags or unavailable legacy packages.
- Text encodings may be 1999-era. Detect and record the source encoding,
  preserve the original file hash, and convert only copied artifacts to
  UTF-8 (for known Latin-1, `iconv -f ISO-8859-1 -t UTF-8`). Never
  rewrite files in place or silently replace undecodable bytes.

## Acceptance
- CATALOG.md + DFT-LEDGER.md (+ schema-compliant CSV) merged; curated
  artifacts in legacy/thesis-archive/ per above; archive source
  untouched (verifiable by identical pre/post metadata manifests and
  their recorded SHA-256 values).
- Follow-up cards filed.

## Progress
- 2026-09-02 — cataloged and independently rehashed all 1,426 files; curated both C source snapshots, 30 RCS archives, prototype inputs, and five complete 1999 KMC runs spanning hosts, temperature, and chemical potential. (hermes-custom-build-001; profile=workstation)
- 2026-09-02 — parsed all 83 Gaussian result logs into exhaustive source coverage plus a fixed-schema 89-row ledger. Forty-five Chapter 4 thermochemical reactions reproduce the thesis tables; no defensible transition-state calculation exists, so barriers remain explicitly unmatched. (hermes-custom-build-001; profile=workstation)
- 2026-09-02 — rebuilt the thesis cleanly through LaTeX/BibTeX/DVI/PDF, retained the original 134-page PDF beside the 133-page clean rebuild, and filed A8a/A8b/A8c/A8d follow-ups. Independent acceptance verifier passed seven evidence groups. (hermes-custom-build-001; profile=workstation)
- 2026-08-20 — card created (fable), night the archive came off the NAS.

## Result

- `legacy/thesis-archive/CATALOG.md` documents provenance, the complete archive map, curation choices, five golden runs, DFT semantics, thesis build, toolchain, and exclusions.
- `ARCHIVE-INVENTORY.csv` hashes all 1,426 files; `MC-RUN-INVENTORY.csv` catalogs all 88 KMC runs. Reproducible generators and independent verifier are committed under `tools/`.
- `DFT-SOURCE-LOGS.csv` covers all 83 Gaussian result logs exactly once. `DFT-LEDGER.csv` has the required 19-column schema, 44 unmatched source calculations, and 45 mapped Chapter 4 reaction energies. No reactant/TS/product triplet exists in the archive; no barrier was fabricated.
- Pre/post GNU-find metadata manifests are byte-identical at SHA-256 `edf7dea5a690823c9a131d0e28cac5ca68634a0538b3df424dd8c91235198fdb`.
- Independent closeout receipt: 1,426 files rehashed; DFT coverage/schema/arithmetic/TS-premise checked; 96 model/input and 36 golden-run copies byte-compared; both PDFs parsed and text-probed; seven checks passed.
- Follow-ups filed: A8a KMC conformance, A8b modern TS/barrier bridge, A8c Table 4.10(c) erratum, and A8d ligand-promotion scoping.
