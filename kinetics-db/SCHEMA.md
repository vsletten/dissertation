# Kinetics database schema v0

This directory is a machine-readable, deliberately conservative transcription layer for mineral dissolution kinetics. Version 0 starts with the silicate core of Palandri & Kharaka (2004), USGS Open-File Report 2004-1068. It preserves omissions and ambiguities instead of filling them with guessed values.

Primary compilation: <https://pubs.usgs.gov/of/2004/1068/pdf/OFR_2004_1068.pdf>

## Layout

- `minerals/<slug>.toml` — one mineral or compositional mineral group per file.
- `SOURCES.md` — post-2004 expansion inventory; it is not imported by the validator.
- `validate.py` — Python 3.11+ stdlib validator (`tomllib` only).

Run from the repository root:

```bash
python3 kinetics-db/validate.py
```

The command fails on malformed TOML, missing required fields, unsupported mechanism kinds, non-finite numeric parameters, missing uncertainty placeholders, broken provenance pointers, or fewer than 25 mineral records.

## Record contract

```toml
schema_version = 0
name = "Kaolinite"
formula = "Al2Si2O5(OH)4"
mineral_group = "clay"
rate_law = "Palandri-Kharaka Eq. 6 (...)"
surface_area_basis = "unspecified" # BET | geometric | unspecified
caveats = ["Free-text limitations or transcription choices."]

[source]
id = "pk2004-table-29-p39"
citation = "..."
url = "https://..."
year = 2004
page = 39
# Printed report page, not PDF object index.
table = 29
notes = "..."

[[mechanism]]
kind = "acid" # acid | neutral | base | other
log_k25_mol_m2_s = -11.31
activation_energy_kj_mol = 65.9
reaction_orders = { H_plus = 0.777 }
source_ref = "pk2004-table-29-p39"

# Uncertainty fields are present even when the source reports none. Empty arrays
# mean “not reported”; they are not zero uncertainty.
uncertainty_status = "not_reported"
uncertainty_log_k25 = []
uncertainty_activation_energy_kj_mol = []
uncertainty_reaction_orders = []
uncertainty_notes = []
```

### Required top-level fields

| Field | Type | Meaning |
|---|---|---|
| `schema_version` | integer | Must be `0`. |
| `name` | non-empty string | Display mineral name. |
| `formula` | non-empty string | Formula or honest compositional-group formula. |
| `mineral_group` | non-empty string | Broad class used for inventory/filtering. |
| `rate_law` | non-empty string | Parameterization contract. |
| `surface_area_basis` | enum | `BET`, `geometric`, or `unspecified`. `unspecified` is required when the summary table does not identify the underlying normalization; never guess. |
| `caveats` | array of strings | Ambiguities, cutoffs, pooled datasets, formula scope, or known comparability traps. |
| `source` | table | Provenance anchor for every numeric parameter in the record. |
| `mechanism` | non-empty array of tables | Mechanism-specific kinetic parameters. |

### Mechanism semantics

Palandri & Kharaka Eq. 6 expresses a surface-normalized far-from-equilibrium rate as a sum of Arrhenius mechanisms:

`r = Σ k25_i exp[-Ea_i/R (1/T - 1/298.15)] Π a_j^(n_ij) (1 - Ω^(p_i))^(q_i)`

- `log_k25_mol_m2_s` is base-10 `log(k)` at 25 °C and pH 0, in `mol m^-2 s^-1`, as printed in the source table.
- `activation_energy_kj_mol` is `Ea` in `kJ mol^-1`.
- `reaction_orders` maps activity names to exponents. P&K prints the base mechanism using a **negative order in `H_plus`**, not a positive order in OH-. The record preserves that convention.
- `source_ref` must equal `[source].id`. This is the mechanical guarantee that each numeric kinetic value has an explicit page/table provenance.
- Empty `reaction_orders` means the source table did not print an order or the mechanism is pH-independent. It does not mean a measured order of zero.
- `uncertainty_* = []` means the compilation did not report a machine-readable uncertainty for that quantity. Future records may carry one numeric standard uncertainty in each list; the schema avoids fake zeroes.

### Surface-area normalization

BET versus geometric area is a major comparability trap. P&K explicitly presents both for quartz, so `quartz.toml` selects and labels the BET row. Most summary tables do not state a single basis after combining primary studies; those records say `unspecified` and carry a caveat where needed. Version 0 does **not** silently convert between bases.

## Transcription policy

1. Preserve source precision; do not add significant figures.
2. Use the printed report page and table in `[source]`.
3. Never infer a missing mechanism, reaction order, uncertainty, or area basis.
4. Put report-imposed neutral cutoffs and pooled/interpolated parameters in `caveats`.
5. Formula strings may describe a solid-solution family. They are identifiers, not stoichiometric inputs to a solver.
6. Numeric kinetic values must be backed by `source_ref`; the validator rejects dangling provenance.
7. A correction changes both the TOML and its source/caveat receipt in one commit.

## Kaolinite cross-check

`kaolinite.toml` records P&K Table 29 exactly: acidic `Ea = 65.9 kJ/mol`, neutral `22.2`, and base `17.9`. `qm/SURVEY.md` section 7.3 separately quotes experimental targets of about `29 kJ/mol` acidic and `33–51 kJ/mol` alkaline. Those ranges do not agree with the P&K fitted acid/base values. The TOML caveat states the discrepancy rather than selecting a preferred number: the P&K mechanism-specific pooled fit and individual apparent experimental activation energies are not interchangeable.
