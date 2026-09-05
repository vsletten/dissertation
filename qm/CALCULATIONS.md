# Calculation inventory — numbers the KMC decks are waiting on

A living ledger connecting Petra decks to the quarry pipeline: every
physical number a deck carries as a placeholder, wants at higher quality,
or already has from the literature. **Add an entry whenever deck work
invents a parameter; update its status when quarry (or a paper) delivers.**
Keep entries small — provenance and detailed job specs live in quarry's
SQLite store once a calculation runs; this file is the shopping list and
its state.

Statuses: `needed` (no number yet — deck uses a placeholder or the deck is
blocked), `literature` (a published value is in use or ready to use),
`queued` (quarry job defined), `computed` (quarry produced it; provenance
in the store), `landed` (the number is in a merged deck).

| id | quantity | consumer | status | current value / source |
|---|---|---|---|---|
| CALC-001 | Kaolinite interlayer H-bond energy (per OH···O bond, inner-surface OH → basal siloxane O) | `petra/examples/kaolinite-multilayer.toml` — `e_hb` in the two desorption `per_match` modifiers | **needed** | placeholder 5.0 kcal/mol (generic OH···O; deck header flags it) |
| CALC-002 | Si–O–Si hydrolysis barrier, neutral/acid/base (Xiao & Lasaga reproduction) | kaolinite decks R0/R1; quarry Phase 1 validation | **computed** (neutral free dimer production tier + embedded pilot; acid/base pending) | A2a accepted sequential free-dimer path: addition TS then barrierless cleavage shelf; wB97M-V/def2-TZVPD+SMD **ΔG‡ 142.807 kJ/mol** (electronic 134.504), r2SCAN-3c electronic 125.742, B3LYP-D4 electronic 148.294, directive-adjusted canonical-TZ + TightPNO-CBS focal point **132.960 kJ/mol** with DLPNO−canonical TZ gate **+0.246 kJ/mol** (pass). Provenance: external `production-closeout/store.sqlite` and focal summary, exact hashes in A2a card. Embedded n4 cell remains **205.7 kJ/mol** (sequential, TS2 rate-limiting, +92.7 lattice shift) @ b3lyp/def2-svp/df in `qm/runs/phase2/oss-neutral-n4-s2-b3lyp-def2-svp/store.sqlite` |
| CALC-003 | Si–O–Al hydrolysis barrier, both protonation sides | kaolinite decks R2–R13 family | **computed** (free-dimer neutral + one-water Al-acid; matched microsolvated acid pair pending) | al-neutral sequential uphill-slide profile **134.811 kJ/mol (32.221 kcal/mol)**, cleavage TS rate-limiting. Same `H6SiAlO7-` / `(HO)3Si–O–Al(OH)3-` core and charge as Victor's 1990s calculation: 29.9267 kcal/mol direct `ΔE*` (35.7769 pre-adsorbed step), numerically corroborative but not the same observable. Neutral non-confirmation narrows rather than refutes Xiao–Lasaga's catalyzed ordering; matched A1b acids test it directly; Izadifar's embedded 86-vs-169 kJ/mol ordering remains open for embedded rungs. B3LYP/def2-SVP/DF; `AL_NEUTRAL_MECHANISM.md` + `ACID_MECHANISMS.md`; A2 re-tiers |
| CALC-004 | Al–OH–Al hydrolysis barrier ± protonation state | kaolinite decks R14/R15 | needed | legacy relative rates |
| CALC-005 | Si and Al attachment/detachment energetics vs coordination (the `by_count` ladders) | kaolinite decks desorb/adsorb modifiers | needed | legacy 6/12 kcal-per-bucket heuristic tables |
| CALC-006 | Quartz Q1/Q2/Q3 per-connectivity dissolution barriers | planned quartz framework deck (handoff deck #3) | **literature** | Liu & Ruiz Pestana 2024: 54/71/81 kJ/mol — drops straight into a `by_count` table |
| CALC-007 | Gibbsite step-edge vs terrace Al detachment barriers | planned gibbsite deck (handoff deck #1) | literature | Guo et al. gibbsite step-edge barriers — extract numbers when the deck is authored |
| CALC-008 | Kaolinite elastic inputs for strain fields: μ (and anisotropy check), ν | `[[defects]]` sections of future kaolinite defect decks | needed | μ ≈ 20 GPa isotropic assumed (STRAIN.md §6 open question) |

## Entry notes

### CALC-001 — interlayer H-bond energy

The multilayer deck models the anhydrous interlayer with one H-bond per
inner-surface hydroxyl (6 per cell, donor→acceptor assignments in the deck
header, O···O 2.67–3.40 Å) and charges `e_hb` per intact bond against the
removal (desorption) rates, Blum–Lasaga style (all of it in the forward
barrier). Wanted from quarry: the per-bond interlayer binding energy —
e.g. periodic DFT (r²SCAN or B3LYP-D4) bilayer cohesion divided per
H-bond, or a constrained cluster pair. Total interlayer cohesion
literature is ~20–30 kJ/mol per formula unit (~2 H-bonds), so the 5.0
kcal/mol placeholder is the right order but not a number to publish on.
At the legacy T = 8000 K knob it shifts rates by exp(5/15.9) ≈ 1.4× per
bond; at real temperatures it is decisive (exp(5/0.6) ≈ 4000×) — which is
why the real number matters before any real-T multilayer run.

### CALC-002..005 — the barrier ladder

The kaolinite deck's environment machinery is validated but runs on the
legacy flat tables (every variant of a reaction shares one number —
DESIGN.md §10.2). These entries are the qm/HANDOFF.md Phase 1–2 targets:
connectivity- and protonation-resolved barriers for the five site
families, emitted as Petra `by_count`/`when` tables with provenance.

### CALC-003 — prior-work reconciliation

- **Charge and termination match:** TASK-168 uses exactly
  `(HO)3Si–O–Al(OH)3-` (`H6SiAlO7-`, charge -1, singlet), then adds one neutral
  attacking water. That is the same terminated anionic dimer core documented in
  Victor's thesis; the old direct reaction additionally carried a spectator
  water through `H2O·H2O` bookkeeping.
- **Victor's number:** the thesis reports 29.9267 kcal/mol for the direct net
  neutral-water `ΔE*` and 35.7769 kcal/mol for the pre-adsorbed-water hydrolysis
  step (B3LYP/3-21G(3d,2p) geometries, 6-31G(3d) energies). TASK-168's 32.221
  kcal/mol is a 298.15 K `ΔG‡`, so the 2.294-kcal numerical agreement
  corroborates scale, not equality of observables. Sources: Victor Sletten
  dissertation ch. 4 (“Computational Method” / “Reaction Model”) and the
  thesis dimer table (workstation archive; not yet ingested — see
  `docs/program/cards/A8-thesis-archive-intake.md`); current provenance is
  in `AL_NEUTRAL_MECHANISM.md` and the TASK-168 durable receipt hashes there.
- **Three-tier interpretation:** neutral free-dimer `al > si` does not directly
  test Xiao & Lasaga's catalyzed ~24/~19 kcal/mol channels; matched A1b acids do.
  Izadifar's embedded Si–O–Al 86 versus Si–O–Si 169 kJ/mol ordering remains an
  embedded-rung question. Citations: Xiao & Lasaga, *GCA* 58 (1994) 5379–5400
  and 60 (1996) 2283–2295; Izadifar et al., *Nanomaterials* 13 (2023) 1196 and
  *Nanoscale Advances* 7 (2025); `SURVEY.md` §§7.1–7.2.
