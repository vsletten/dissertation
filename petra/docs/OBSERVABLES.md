# Ensembles and observables

Petra's B4 framework implements RFC-001 §6.3 in the `petra-observables`
crate. The input deck remains the scientific document: execution selects the
replica count and seed policy, while `[[observables.series]]` selects what is
sampled. The CLI writes long-form `observables.csv`, final-replica
`ensemble.csv`, and `ensemble-summary.csv` with the full distributions,
means, and deterministic percentile-bootstrap 95% confidence intervals.

## Execution and reproducibility

`run_ensemble` uses Rayon only across replicas. A replica is still one serial
engine and one serial RNG stream. Indexed parallel collection restores
replica-index order before aggregation, so thread scheduling cannot alter
output. Replica seeds are derived from the deck's pinned `increment` or
SplitMix64 `hash` policy. A CLI `--seed` override becomes the ensemble's base
seed without changing that policy.

Every run records the initial sample, cadence samples, and the final state even
when the stopping step falls between cadence boundaries. Bootstrap resampling
uses a separate caller-supplied seed and cannot consume simulation randomness.

## Implemented declarations

```toml
[observables]
report_every = 200

[[observables.series]]
kind = "state_counts"

[[observables.series]]
kind = "event_rates"

[[observables.series]]
kind = "rate_spectra"

[[observables.series]]
kind = "cluster_sizes"
state = "M_site.empty"       # qualify as Kind.state

[[observables.series]]
kind = "surface_area"
state = "M_site.occupied"
axis = 2
```

- `state_counts`: one value per dense state ID.
- `event_rates`: the instantaneous total CTMC propensity per rule.
- `rate_spectra`: the distribution of positive instantaneous total
  propensities over sites (the Fischer/Luttge observable).
- `cluster_sizes`: connected-component sizes for a required state filter.
- `surface_area`: two deliberately separate measures. Projected geometric
  area counts occupied surface-normal columns once; the BET-like proxy counts
  every exposed solid site, including internal pit/roughness surfaces. A site
  is exposed when it lies on the declared open-axis face or touches a site
  outside the selected solid-state set.
- `exposure_age`: one value per currently exposed solid site, in deterministic
  site order. The value is physical/strategy time since that site most recently
  became exposed.

State filters compile to dense IDs before execution and may name either one
qualified `Kind.state` or a multi-state `@alias`. Observation is read-only,
uses deterministic site/adjacency order, and never draws from an RNG.

## Validation gates

The 2-D Ising Binder-cumulant crossing now runs through the parallel ensemble
and declarative `state_counts` framework. The Kossel etch-pit deck declares the
first full observable set; its rate-spectrum log-width broadens from the
pinned initial golden value as the pit nucleates. Both are CI tests.

## A5 exposure-age implementation

A5 Phase 0 needs the age distribution of *currently exposed* sites, not merely
a timestamp inferred by a reporter. Petra implements that in two layers:

1. **Core-owned transition metadata.** `Lattice::exposed_since` stores a
   per-site optional timestamp. During the same deterministic commit that
   applies state changes and refreshes dirty sites, the engine recomputes
   exposure for the changed site and its immediate neighbors: it timestamps a
   newly exposed site at the committing event's time, preserves the timestamp
   while continuously exposed, and clears it on burial/removal. Exposure
   history is therefore replayable trajectory state, not observer-local state.
2. **A read-only observable.** `kind = "exposure_age"` accepts the same solid
   `state`/alias and `axis` parameters as `surface_area`. At cadence it exports
   `current_time - exposed_since` for every currently exposed site. Existing
   ensemble aggregation can consume that distribution without touching the
   simulation RNG.

The update preserves changed-site and adjacency insertion order and uses
marking arrays for deterministic deduplication. CTMC reports physical time;
discrete strategies report their RFC-defined strategy time (`1`, `1/N`, etc.).
The Kossel and production kaolinite examples both declare the observable; the
finite-inventory smoke and plot are documented in
`docs/program/results/A5p0-aging-observables.md`.
