# Learnings Index

- [gotchas.md](gotchas.md) — Petra decks need `[simulation]` even for validation; pyscf thermo auto-applies rotational symmetry numbers (R·ln σ entropy trap); petra's truncated gas constant vs CODATA (send barriers over the bridge, never rates); legacy-crude kaolinite cell geometry is a frozen-shell systematic (freeze heavies only, never constructed H)
- [performance.md](performance.md) — first 4090-vs-CPU DFT timings: small jobs favor CPU; GPU wins need density fitting + warmup + cuTENSOR + ~100 atoms
- [debugging.md](debugging.md) — trivial-saddle trap: verified ≠ right saddle; interior scan maximum required before Sella; CuPy pool fragmentation OOMs multi-hour GPU runs (trim at stage boundaries; pipefail with tee)
