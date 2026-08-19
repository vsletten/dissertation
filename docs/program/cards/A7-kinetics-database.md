# A7-kinetics-database — machine-readable dissolution-kinetics compilation, v0

- status: ready
- track: A (geochemistry, community-infrastructure lane)
- priority: P2
- machine: any (laptop-shaped by design: reading, extraction, structuring — no GPU, no long compute)
- depends: —
- claimed-by: (see agents/A7-kinetics-database branch)

## Objective

Groundwork for a modern successor to Palandri & Kharaka (2004), *A
compilation of rate parameters of water-mineral interaction kinetics for
application to geochemical modeling* (USGS Open-File Report 2004-1068) —
still the field's standard rate compilation, 20+ years stale, and not
machine-readable. This card is v0: schema + the core extraction.

1. **Schema v0**: define the record format (TOML, one file per mineral)
   in `kinetics-db/SCHEMA.md`. Required fields: mineral name + formula;
   per-mechanism (acid/neutral/base) log k25, Ea, reaction orders (n_H+,
   n_OH-, others); rate-law form; surface-area normalization basis
   (BET vs geometric — record which, this is the classic comparability
   trap); provenance (source, page/table, year); uncertainty fields
   (even if initially empty — the schema must carry them so later UQ
   work has a home); free-text caveats.
2. **Extract P&K's silicate core**: the phyllosilicates (kaolinite
   first), feldspars, quartz/silica polymorphs, and olivine/pyroxene
   groups — target ≥25 minerals — into `kinetics-db/minerals/*.toml`
   per the schema, each value carrying its P&K table/page provenance.
   (Full P&K coverage + carbonates → follow-up card.)
3. **Source inventory 2004→present**: `kinetics-db/SOURCES.md` — a
   candidate list (~30–50 entries, DOIs) of primary rate studies and
   compilations published since P&K, one line each on what they'd add.
   Include the obvious anchors: Ganor/Cama kaolinite work, Brantley
   group compilations, Hermanska et al. 2022-23 (the recent
   thermodynamic-consistency compilation — check scope overlap
   honestly), ERW-driven olivine/wollastonite studies.
4. **A tiny validator**: `kinetics-db/validate.py` (stdlib tomllib
   only) that checks every mineral file parses and carries required
   fields; wire as a runnable check documented in SCHEMA.md.

## Context

- docs/OMNIBUS.md Track A; the ideation lineage: this is the
  "community infrastructure" play — unglamorous, citable, and the
  natural calibration substrate for the field-lab discrepancy flagship
  (docs/scoping/field-lab-discrepancy.md needs exactly these rates for
  its validation gates).
- P&K 2004-1068 is a public USGS document (pubs.usgs.gov). Extraction
  is transcription with provenance, not judgment — flag ambiguities in
  caveats rather than resolving them silently.

## Acceptance

- `kinetics-db/SCHEMA.md`, `SOURCES.md`, `validate.py` exist on main;
  ≥25 mineral TOML files parse cleanly under the validator; every numeric
  value has provenance; kaolinite record cross-checked against the
  experimental Ea range quoted in qm/SURVEY.md §7.3 (29 kJ/mol acidic /
  33–51 alkaline) with any discrepancy noted in caveats.
- STATUS.md line + this card updated in the delivering PR (protocol).

## Progress

- 2026-08-18 — card created (fable), macbot-targeted by design.
