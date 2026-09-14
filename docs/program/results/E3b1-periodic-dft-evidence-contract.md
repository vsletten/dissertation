# E3b1 periodic-DFT evidence contract

**Status:** complete machinery/evidence contract; no E3b DFT barrier is claimed
**Operator:** `(hermes-custom-build-001; profile=workstation)`
**Evidence root:** `/mnt/data/vsletten/dissertation-data/e3b1-periodic-dft-evidence-contract-20260914`
**Implementation head:** `1baf3dc766fd1cad7c1bc0947c4d404f630d3184`

## Verdict

E3b1 is complete. The E3b calibration route now fails closed unless CP2K and
LAMMPS results are independently parsed from immutable, hash-bound native raw
outputs. Prepared inputs execute only from a copied runtime tree, all three E3a
classical anchors remain explicitly `incomplete-convergence`, and the
five-coordinate-Al requirement is an atom-identity/coordination/path-interaction
criterion rather than a midpoint-distance proxy.

This card does **not** satisfy E3b's scientific acceptance. The parent remains
at **0/3 DFT barriers checked**. It is now unblocked and ready to run those
bounded calculations with the corrected evidence machinery.

## Prepared models and source gates

Fresh preparation from the final implementation produced the same valid wrapped
2×2×1 route crop at window `(5, 2)`:

| model | atoms | E3a reference |
|---|---:|---|
| reconstructed-replication | 334 | `incomplete-convergence` |
| dehydroxylate-screen | 331 | `incomplete-convergence` |
| stoichiometric-screen | 334 | `incomplete-convergence` |

The preparation and manifest are retained under `prepared-reviewed/`. The final
receipt proves that the exact input and every CP2K data/geometry dependency used
by the completed timing run are byte-identical to the final reviewed
preparation.

## Evidence contract closed

- CP2K endpoint and BAND evidence is parsed from copied native outputs with
  immutable size/SHA-256 records. BAND collection requires force-evaluation
  accounting, optimization completion, climbing-image evidence, nontrivial
  image evolution, and a profile-derived barrier.
- Matched LAMMPS endpoint/NEB evidence is likewise copied, hash-bound, and
  independently reparsed. Manually asserted barriers or convergence booleans
  cannot emit a calibration number.
- Cross-model calibration is allowed only when native evidence is complete,
  atom/model identity matches, and both solver-specific gates pass. The three
  current E3a references remain non-promotable.
- The dehydroxylate criterion proves Al atoms 10 and 11 change from six- to
  five-coordinate when hydroxyl O 76 is removed and requires the explicit
  eight-image noble-gas path to interact with that chemical environment.
- Prepared inputs remain immutable: runtime files live outside the manifest,
  and smoke execution/removal does not invalidate a later smoke or analysis.

## Completed reconstructed-Ar timing probe

A real 334-atom reconstructed-Ar PBE-D3 `ENERGY_FORCE` probe completed in a
bounded transient systemd unit using pinned CP2K 2024.3 Docker image
`cp2k/cp2k@sha256:979e011c2ea15a2a56eae16b5f3a2d214142bb74fe644bfb302d366dddc2a443`.

| field | value |
|---|---|
| status | `converged` |
| wall time | 780.6223025830113 s |
| energy | -3796.7518690224792 Eh |
| resources | 4 MPI × 4 OpenMP = 16 CPU threads; 48 GiB; nice 10 |
| unit bound | 3700 s |
| cleanup | disposable runtime removed; CID file absent; no image-matching container remains |

The systemd journal records the transient unit start. Because the unit used
collect semantics, its runtime properties were garbage-collected after success;
the atomic timing receipt retains the declared bound/resources, and the native
CP2K receipt proves normal termination, SCF convergence, output hash, elapsed
time, method identity, and cleanup. This is planning evidence only—not a BAND
lower bound or a calibration value.

Key durable receipts:

- `final-evidence.json` — SHA-256
  `cc199ab545fa29e9d6ceb5c24db84f9229f0d0208f12008631dda561f1d721e3`
- `evidence-manifest.json` — SHA-256
  `fe8e6cac49343bb661e989983689c937b751ca563e993c01342e73761eaf7a2b`
- `timing-probe-receipt.json`
- `runtime-reconstructed-ar/smoke-result.json`
- `runtime-reconstructed-ar/smoke.stdout.log`
- `prepared-reviewed/preparation.json`
- `prepared-reviewed/manifest.json`

The evidence manifest covers 146 retained files totaling 6,613,627 bytes.

## Verification

Final post-review checks at implementation head
`1baf3dc766fd1cad7c1bc0947c4d404f630d3184`:

```text
uv run --with openpyxl --with pytest pytest tests/test_e3b_periodic_dft.py -q
30 passed

uv run --with ruff ruff check scripts/e3b_periodic_dft.py tests/test_e3b_periodic_dft.py
All checks passed

uv run --with ruff ruff format --check scripts/e3b_periodic_dft.py tests/test_e3b_periodic_dft.py
2 files already formatted

python3 -m py_compile scripts/e3b_periodic_dft.py tests/test_e3b_periodic_dft.py
git diff --check
passed
```

The single allowed cold review found six correctness/safety defects in the first
implementation. All were fixed before final verification: native CP2K barrier
parsing, climbing-image proof, LAMMPS structure mapping, swap-resistant Al/O
identity, stricter E3a convergence evidence, and process-group timeout cleanup.
No second cold review or expensive BAND run was performed.
