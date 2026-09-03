# 1999 C model conformance harness

This harness turns the curated C source and five archived runs into an executable,
machine-readable oracle for the Rust port. It preserves the model as implemented,
including the inactive diffusion slots 24–27.

## Canonical pinned run

From the repository root, with Docker available:

```bash
mkdir -p /tmp/a8a-conformance
rm -rf /tmp/a8a-conformance/runs

docker build \
  -f legacy/thesis-archive/c-model/Dockerfile.conformance \
  -t dissertation-a8a-conformance .

docker run --rm \
  -v /tmp/a8a-conformance:/out \
  dissertation-a8a-conformance
```

The Dockerfile pins the complete multi-architecture image index by digest. The
harness additionally refuses to run unless `gcc --version` exactly matches the
identity in `conformance-toolchain.json`. Use a new reviewed lock if either changes;
do not silently pass `--allow-compiler-drift` for the canonical oracle.

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
- `compiler_prng_drift`: every output has the historical shape and row count,
  `results.dat` starts with at least one exact row, then the stochastic trajectory
  diverges. The mismatch is retained; no tolerance normalizes it away.
- `behavioral_mismatch`: nonzero/timeout, missing or malformed output, changed row
  count, or divergence beginning at the first trajectory row.

The two independent 1999 hosts' identical-input fixtures are also required to have
identical archived hashes. A sabotage unit test perturbs row 1 and must be classified
as `behavioral_mismatch`.

## Tests

```bash
python3 -m unittest legacy/thesis-archive/tools/tests/test_c_model_conformance.py -v
```

The compatibility history and unmodified failure are frozen in
`UNMODIFIED-BUILD-FAILURE.md`.
