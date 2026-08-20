# Learnings Index

- [gotchas.md](gotchas.md) — Petra decks need `[simulation]` even for validation; pyscf thermo auto-applies rotational symmetry numbers (R·ln σ entropy trap); petra's truncated gas constant vs CODATA (send barriers over the bridge, never rates)
- [performance.md](performance.md) — first 4090-vs-CPU DFT timings: small jobs favor CPU; GPU wins need density fitting + warmup + cuTENSOR + ~100 atoms
- [debugging.md](debugging.md) — trivial-saddle trap: verified ≠ right saddle; interior scan maximum required before Sella; GPU4PySCF 1.8.1 DF-UKS Hessian assertion needs a bounded tiny-system CPU retry
