# Scientific domain of the kaolinite KMC model

Victor's kinetic Monte Carlo model is the Monte Carlo half of *Playing Dice
with the Universe: Monte Carlo and Quantum Mechanics in Geochemistry*. It models
atomic-scale growth and dissolution of kaolinite at a mineral-water interface:
the microscopic weathering reactions underlying larger geochemical cycles.

## State and events

The model uses a rigid lattice. Each site has a finite discrete chemical state;
the numeric site codes are species/coordination states rather than arbitrary
implementation labels. The principal bridging-site families are Si-O-Si,
Al-OH-Al, and interlayer Si-O-Al2. Elementary events include hydrolysis and its
reverse condensation, plus adsorption and desorption. Hydrolysis removes or
changes bonded material; reverse channels represent growth.

At each step the model enumerates valid local events, selects one in proportion
to its rate, updates the affected site neighborhood, and advances time with the
N-fold-way expression

```text
dt = -ln(xi) / sum(k_i)
```

The interlayer Si-O-Al2 site has five hydrolysis pathways. Together with the
other bridge-site pathways, the dissertation's table defines fourteen elementary
reactions before directionality and ancillary adsorption/desorption expand the
runtime reaction list.

## Why the hydrolysis and adsorption tables are flat

The shipped hydrolysis and adsorption tables in `data.rxn` are flat because
barriers from the dissertation's other half: ab-initio electronic-structure
calculations plus transition-state theory, were not available when the model
was written. Desorption retains distinct environment-dependent rate buckets
across its environment variants; the flatness of the hydrolysis and adsorption
tables does not mean that local-environment machinery has no physical purpose.

Any future environment-dependent rate deck must preserve the dissertation's
thermodynamic constraints, including microscopic reversibility and consistency
between alternate pathway sums. Derived rates should carry provenance to the
calculation, temperature, standard state, and unit convention used to obtain
them.

## Interpretation of historical repairs

Repairs that restore a reverse condensation channel change scientifically
observable growth behavior, not merely code coverage. Likewise, restoring
variable random seeds enables ensemble sampling needed for rate and activation
energy studies; a single fixed trajectory cannot establish those distributions.

For behavior-as-implemented, state encodings, exact reaction mappings, and parity
claims, use the living code and archaeology documents. This note supplies the
physical interpretation and must not override executable evidence.
