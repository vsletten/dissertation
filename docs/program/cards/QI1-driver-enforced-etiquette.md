# QI1-driver-enforced-etiquette — campaign drivers enforce their own compute caps

- status: done
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
   hard-sets thread env vars BEFORE numpy/pyscf import and applies
   os.nice() — called at the top of every campaign entry-point
   (phase1_xiao_lasaga.py, phase2_ladder.py, phase0_smoke, and any
   future scripts/), overridable only by explicit CLI flag (--threads,
   --nice), never by environment inheritance.
2. Logging: entry-points refuse to run without a --log/tee path OR
   self-tee to a timestamped file under runs/ and print the path first.
3. Tests: enforce() sets the env vars when unset, clamps when set
   higher, and leaves them untouched when at/below the cap; entry-points
   call it before heavy imports (import-order test); CLI overrides work
   (including `--threads` raising the cap and negative `--nice` being
   rejected for non-root); `os.nice()` failure warns and continues.

## Spec clarifications

- **Thread env vars managed** — `enforce()` sets exactly
  `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, and `OPENBLAS_NUM_THREADS`.
  Clamping is one-directional: if a var is unset it is set to `threads`;
  if already set **above** the cap it is clamped down to `threads`; if
  already set **at or below** the cap it is left untouched (never
  raised). This makes the "≤ cap, never inflate" contract unambiguous.
- **`os.nice()` failure** — on platforms where `os.nice()` is unavailable
  (non-POSIX) or raising niceness fails (e.g. containerized without
  CAP_SYS_NICE), `enforce()` must **warn and continue** (log a
  one-line warning to stderr), not fail hard. Thread caps still apply;
  only the niceness best-effort degrades.
- **CLI override semantics** — `--threads` may **raise** the cap above 16
  (explicit operator intent wins over the default), but may not lower it
  below 1. `--nice` accepts any value `os.nice()` will accept; negative
  (higher-priority) values are **rejected** unless the process is
  running as root, and root-only values are otherwise refused with a
  clear error. When multiple sources conflict (CLI flag vs. env var vs.
  default), precedence is: explicit CLI flag > default; environment
  inheritance is **never** a source (env vars are overwritten, not
  read).

## Acceptance

- All qm/scripts entry-points call enforce() before heavy imports;
  tests green; ruff clean; a deliberately-unset environment run shows
  ≤16 threads (document the check in the PR).

## Progress

- 2026-08-20 22:55 UTC — Done. Full QM test/lint gates and direct CLI smoke checks are green; implementation, card, board, and program status are ready in the single card PR.
- 2026-08-20 22:51 UTC — Implementation and verification checkpoint: added `quarry.etiquette`, wired all four campaign entry points before heavy imports, and added CLI/logging gates. `uv run --frozen pytest -q` passed (162 tests), `ruff check .` passed, and an unset-environment Phase-2 dry run reported all three managed thread variables at 16.
- 2026-08-20 — card created (fable) after the live 22-thread probe.

## Result

Implemented `qm/quarry/etiquette.py` with a default 16-thread cap across
`OMP_NUM_THREADS`, `MKL_NUM_THREADS`, and `OPENBLAS_NUM_THREADS`, best-effort
niceness, non-root negative-nice rejection, and automatic stdout/stderr teeing
to an explicit `--log` path or a timestamped file under `qm/runs/`. All four
campaign entry points bootstrap it before NumPy/PySCF/quarry heavy imports and
accept explicit `--threads`, `--nice`, and `--log` overrides. The entry-point
inventory test makes this convention fail closed for newly added scripts.

Verification:

- `uv run --frozen pytest -q` — 162 passed.
- `uv run --frozen ruff check .` — clean.
- Direct `--help --nice 0 --log ...` smoke checks passed for all four entry points.
- With all three managed environment variables deliberately unset,
  `phase2_ladder.py --family oss --dry-run --threads 16 --nice 0` printed
  `OMP_NUM_THREADS=16 MKL_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16` before heavy
  imports and completed the no-DFT geometry/metadata run.
