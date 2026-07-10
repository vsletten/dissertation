# Dynamics parity oracle (M6)

`traj_seed0_full.txt` — a per-step trace of the **C++ reference model**
(`dissertation/main/model`, read-only) running the golden inputs under the
legacy fixed seed (spec B2), 20,000 steps. It is the oracle the M6 parity gate
(`crates/mckaol-cli/tests/parity_m6.rs`) checks the Rust port against,
step-by-step and bitwise.

## Format

One line per MC step, whitespace-separated:

```
step  site  rxn  0xDT_BITS  occupied  0xSTATE_HASH
```

- `site`, `rxn` — the event the C++ selected and applied that step.
- `DT_BITS` — the raw 32-bit pattern of the C++ `float dt` (exact, not rounded).
- `occupied` — count of occupied non-EDGE sites after the step
  (`state % 100 > 0 && state != 9`).
- `STATE_HASH` — FNV-1a-64 over every site's 4-byte (little-endian) state
  code, in site order. **Non-canonical basis**: the capture harness used
  `1469598103934665603` (a digit-truncated FNV basis); the Rust checker
  matches that exact constant. A hash is only a fingerprint, so the basis need
  only agree between capture and check, which it does.

## How it was captured

The read-only C++ model was rebuilt verbatim in scratch (`g++ 13.3.0
-std=c++11 -O3 -ffast-math`, the golden toolchain) with a single added trace
`fprintf` in the step loop — no change to the model logic. Proof the
instrumentation is inert: the same build reproduces the golden
`start.msi`/`end.msi` SHA-256 hashes byte-for-byte (mission-control
`projects/kmc/golden/manifest.md`).

The RNG reference values in `crates/kmc-engine/src/rng.rs`'s test were captured
the same way (a `main` calling `ran2()` 20 times from seed 0).
