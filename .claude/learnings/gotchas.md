# Gotchas

### Petra decks require a [simulation] section to compile (2026-08-17)
`petra-cli` rejects a deck missing `[simulation]` (steps/seed/report_every)
with "missing field `simulation`" even when you only want `--steps 0`
validation. Any deck template used for emit round-trips (qm/tests/data/
template_deck.toml) must carry one.

### pyscf thermo silently applies the rotational symmetry number (2026-08-17)
`pyscf.hessian.thermo.thermo` detects the molecular point group and divides
the rotational partition function by sigma (2 for water). Comparing an
independent thermochemistry implementation with sigma=1 shows an entropy
discrepancy of exactly R·ln(sigma) — 5.76 J/(mol K) for water. Pass the
matching `symmetry_number` when cross-checking (qm/tests/test_pipeline.py).

### petra's gas constant differs from CODATA by 0.17% (2026-08-17)
petra-core rate.rs uses the dissertation's truncated R = 1.987e-3
kcal/(mol K); CODATA is 1.9872e-3. At a 10 kcal/mol barrier and 298 K the
rate differs by ~0.17%. Negligible against DFT barrier error, but it means
the QM→KMC bridge should always transmit barriers/ΔH‡/ΔS‡ and let the
engine exponentiate — never pre-computed rate constants
(documented in qm/tests/test_rates.py::test_petra_gate_value).

### The kaolinite cell crystallography is legacy-crude (2026-08-19)
The deck cell (petra/examples/kaolinite.toml [cell]) faithfully imports the
dissertation-era data.cell, whose bond lengths run 1.38-2.88 A (real
kaolinite: Si-O ~1.6, Al-O ~1.9). Fine for lattice KMC (topology only), but
QM clusters cut from it start strained and the *frozen shell keeps those
positions forever* — a systematic in every lattice-resistance number until
someone re-derives the cell from a modern refinement (Bish 1993). The free
region re-optimizes, so trends (ΔΔEa by connectivity) are safer than
absolute barriers. Never freeze constructed protons: H positions are
guesses, not crystallography (quarry/crystal.py freezes heavy atoms only).

### Extreme lumped rates can erase slow Fenwick totals (2026-08-20)
A phase-1 muscovite deck used `1e12 s^-1` for effectively instantaneous
surface release beside `~1e-6 s^-1` migration/dehydroxylation events. Removing
the fast event required subtracting nearly equal f64 Fenwick sums; after many
updates the tree reported total rate zero while slow site events still
existed. Use a rate that is merely fast relative to the load-bearing process
(e.g. `1 s^-1` here), keep the dynamic range bounded, and run `--paranoid`.
An engine hardening follow-up should rebuild/compensate the total when this
`events exist but all rates are zero` condition occurs.

### The 4090 is contested: extraction LLM vs QM campaigns (2026-08-20)
The email-poc extraction daemons keep gemma4-deriver (~14 GB) resident on
the workstation GPU with a self-refreshing keep-alive. When it is loaded,
gpu4pyscf spills the DF integral tensor to *pinned* host memory, which
collides with the 7.7 GB memlock ulimit → cudaErrorInvalidValue at
~60 atoms/def2-svp (raising the ulimit is the wrong fix on a box that
swap-spiraled in July: 7 GB unswappable RAM). Arbitration recipe for
long GPU campaigns: `systemctl --user stop ophir-email-pipeline.service
ophir-email-pipeline-hotmail.service`, run with a trap that restarts
them on exit + a `systemd-run --user --on-active=10h` dead-man restart,
and a VRAM-below-threshold gate before launch (ingest daemons keep
running; extraction catches up afterward).

Honcho's `honcho-deriver-1` and `honcho-api-1` containers also call local
Ollama (`gemma4-deriver` + `bge-m3`) and can reload ~19 GB immediately after
`ollama stop`. Stop those two containers for Hessian campaigns, unload both
models, verify `ollama ps` is empty, and restart the containers at exfil.

### One-water Si-acid is not a protonated-bridge minimum (2026-08-22)
At B3LYP/def2-SVP/DF, the neutral disilicate + one H3O+ model returns the
bridge proton to residual water in gas phase, with water moved to 4.5 A, and
under aqueous PCM. Do not hold Obr-H and publish constrained thermochemistry.
Use 3–6 explicit waters and first prove an unconstrained protonated-bridge
minimum. The anionic Si-O-Al cluster *does* stabilize this state and closes as
a sequential addition/cleavage mechanism; receipts are in
`qm/ACID_MECHANISMS.md`.
