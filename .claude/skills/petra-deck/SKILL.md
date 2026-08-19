---
name: petra-deck
description: Author or modify a petra deck (TOML lattice-KMC model) and validate it — required sections, units, gotchas, the compile/run loop. Use when creating decks, emitting deck fragments from quarry, or debugging deck compile errors.
---

# Petra deck authoring

Full reference: `petra/docs/DESIGN.md`. Working examples:
`petra/examples/kossel.toml` (tutorial), `kaolinite.toml` (the real
thing), `kossel-etchpit.toml` (defects).

## Required sections (learned the hard way)

A deck does NOT compile without ALL of: `[deck]`, `[cell]` (+ sites,
bonds), `[[species]]`, `[[kinds]]` (+ states), `[lattice]`, `[thermo]`,
and **`[simulation]`** (steps/seed/report_every — required even if you
only want compile-validation; its absence is the classic
`missing field 'simulation'` error).

## Units

`[deck] units = "kcal/mol" | "kJ/mol" | "eV"` — applies to every
energy-valued field (ea, dh, ds-per-K, dea tables, mu). Default
kcal/mol. quarry emits kJ/mol and its splice guard refuses templates
that don't declare it. Temperatures always Kelvin, prefactors 1/s.

## Rates and modifiers

```toml
rate = { arrhenius = { prefactor = 1.0e13, ea = 41.84 } }
rate = { eyring = { dh = 60.0, ds = -0.02 } }        # ΔH‡, ΔS‡/K

[[reactions.modifiers]]
select = { distance = 1, state = ["occupied"] }
by_count = { dea = [0.0, 8.4, 16.8] }   # ΔEa by match count; last extends
```

Exactly one of constant/arrhenius/eyring; exactly one of
per_match/by_count/when per modifier. Selectors use `Kind.state`
qualified names or `@group` aliases.

## Validate/run loop

```bash
cd petra
cargo run -p petra-cli -- path/to/deck.toml --steps 0        # compile check
cargo run -p petra-cli -- deck.toml --steps 10000 --seed 42 --out /tmp/run
#   add --paranoid for differential invariant checks on new decks
#   add --ensemble N for seed sweeps
cargo test                                                    # includes analytic + parity gates
```

Determinism contract: same deck + seed ⇒ identical trajectory. If your
change breaks a golden gate, the change is wrong until proven otherwise.

## Gotchas

- State declaration order is the `shift` ladder — reordering states
  silently changes `shift = 1` semantics.
- Engine runtime is canonical kcal/mol internally; deck units convert
  once at compile. Don't "help" by pre-converting.
- Frozen boundary sites are excluded from selectors (not fatal aborts —
  deliberate divergence from legacy).
