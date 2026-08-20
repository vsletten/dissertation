# QI1-driver-enforced-etiquette — campaign drivers enforce their own compute caps

- status: ready
- track: A (quarry infrastructure)
- priority: P2
- machine: any (code + tests; no campaigns needed)
- depends: —
- claimed-by: (see agents/QI1-driver-enforced-etiquette branch)

## Objective

Compute etiquette (OMP_NUM_THREADS ≤ 16, background niceness, teed
logs) currently lives in prompts, skills, and docs — and 2026-08-20
proved prompts don't bind subprocesses: a fleet worker's CPU probe ran
~22 threads at nice 0 on the shared workstation (harmless, but only
because the box was lightly loaded). Move enforcement into the code:

1. A `quarry.etiquette` module: `enforce(threads=16, niceness=10)` that
   hard-sets OMP/MKL/OPENBLAS thread env vars BEFORE numpy/pyscf import
   and applies os.nice() — called at the top of every campaign
   entry-point (phase1_xiao_lasaga.py, phase2_ladder.py, phase0_smoke,
   and any future scripts/), overridable only by explicit CLI flag
   (--threads, --nice), never by environment inheritance.
2. Logging: entry-points refuse to run without a --log/tee path OR
   self-tee to a timestamped file under runs/ and print the path first.
3. Tests: enforce() sets the env vars when unset and clamps when set
   higher; entry-points call it before heavy imports (import-order
   test); CLI overrides work.

## Acceptance

- All qm/scripts entry-points call enforce() before heavy imports;
  tests green; ruff clean; a deliberately-unset environment run shows
  ≤16 threads (document the check in the PR).

## Progress

- 2026-08-20 — card created (fable) after the live 22-thread probe.
