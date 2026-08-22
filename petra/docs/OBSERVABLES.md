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
- `surface_area`: two deliberately separate measures: projected geometric
  area from exposed sites and the BET-like exposed-site count proxy. A site is
  exposed when it lies on the declared open-axis face or touches a site outside
  the selected solid state.

State filters compile to dense IDs before execution. Observation is read-only,
uses deterministic site/adjacency order, and never draws from an RNG.

## Validation gates

The 2-D Ising Binder-cumulant crossing now runs through the parallel ensemble
and declarative `state_counts` framework. The Kossel etch-pit deck declares the
first full observable set; its rate-spectrum log-width broadens from the
pinned initial golden value as the pit nucleates. Both are CI tests.

## A5 exposure-age attachment

A5 Phase 0 needs the age distribution of *currently exposed* sites, not merely
a timestamp in a reporter. It should attach in two layers:

1. **Core-owned transition metadata.** Add a per-site `exposed_since: Option<f64>`
   (or strategy-time equivalent) beside lattice state. During the same
   deterministic commit that applies state changes and refreshes dirty sites,
   recompute exposure for the changed neighborhood: set the timestamp on a
   newly exposed site, preserve it while continuously exposed, and clear it on
   burial/removal. This makes exposure history part of replayable trajectory
   state rather than observer-local inference.
2. **A read-only observable.** Add `kind = "exposure_age"` with the same solid
   `state` and `axis` parameters as `surface_area`. At cadence it exports
   `current_time - exposed_since` for exposed sites. Ensemble aggregation then
   uses the existing distribution/mean/bootstrap machinery unchanged.

No extra RNG stream is needed. The update must preserve the engine's current
changed-site insertion order and expand the deterministic dirty neighborhood;
using an unordered set would violate the trajectory contract. CTMC reports
physical time, while discrete strategies report their RFC-defined strategy
time (`1`, `1/N`, etc.).
