# 1999 C model conformance harness

This harness turns the curated C source and five archived runs into an executable,
machine-readable oracle for the Rust port. It preserves the model as implemented,
including the inactive diffusion slots 24–27.

## Canonical pinned run

From the repository root, with Docker available:

```bash
mkdir -p /tmp/a8a-conformance
rm -rf /tmp/a8a-conformance/runs

docker build --platform linux/amd64 \
  -f legacy/thesis-archive/c-model/Dockerfile.conformance \
  -t dissertation-a8a-conformance \
  legacy/thesis-archive

docker run --rm --platform linux/amd64 \
  -v /tmp/a8a-conformance:/out \
  dissertation-a8a-conformance
```

The Dockerfile pins the complete multi-architecture image index by digest and fixes
the selected platform to `linux/amd64`. Its restricted build context excludes the
worktree's host-only `.git` pointer and unrelated files. The harness additionally
refuses to run unless compiler identity, architecture, all curated source bytes,
and all 20 fixture input files match `conformance-toolchain.json`. Use a new reviewed
lock if any changes; do not pass `--allow-compiler-drift` for the canonical oracle.

Outputs:

- `/tmp/a8a-conformance/conformance-report.json` — machine-readable oracle;
- `/tmp/a8a-conformance/conformance-report.md` — human summary;
- `/tmp/a8a-conformance/runs/` — disposable replay artifacts.

The five 5,000,000-step replays are CPU-bound. `--jobs 4` is the maximum allowed by
fleet policy and avoids an unbounded worker pool.

## Explicit host-comparison run

A non-canonical compiler may be useful to measure architecture/compiler drift. It
must be labeled, never mistaken for the pinned oracle:

```bash
python3 legacy/thesis-archive/tools/run_c_model_conformance.py \
  --source legacy/thesis-archive/c-model/source \
  --fixtures legacy/thesis-archive/golden-runs \
  --report legacy/thesis-archive/c-model/conformance-report.json \
  --markdown legacy/thesis-archive/c-model/conformance-report.md \
  --compiler clang \
  --allow-compiler-drift \
  --jobs 4 \
  --timeout-seconds 2400
```

## Comparison contract

For each `results.dat`, `surfAl.out`, and `surfSi.out`, the report records both
SHA-256 hashes, row/column shape, parse status, exact matching prefix, first
mismatch row, and maximum numeric delta.

- `byte_parity`: exact SHA-256 equality.
- `compiler_prng_drift_candidate`: only for an explicitly allowed non-canonical
  compiler/architecture comparison whose locked source and inputs match, every
  output has the historical shape and row count, and at least ten initial
  `results.dat` rows are exact before the stochastic trajectory diverges. The
  mismatch stays visible and the label remains a candidate, not a fact.
  Drift candidates never pass the canonical conformance gate; only five exact
  `byte_parity` outcomes do.
- `behavioral_mismatch`: nonzero/timeout, missing or malformed output, changed row
  count, or divergence beginning at the first trajectory row.

The two independent 1999 hosts' identical-input fixtures are also required to have
identical archived hashes. A sabotage gate perturbs a plausible value after eleven
matching rows and must remain `numeric_divergence`, never promoted to drift. Both
controls and the diffusion guard are hard success gates.

## Tests

```bash
python3 -m unittest legacy/thesis-archive/tools/tests/test_c_model_conformance.py -v
```

The compatibility history and unmodified failure are frozen in
`UNMODIFIED-BUILD-FAILURE.md`.
