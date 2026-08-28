# QI2-gpu-lease — single-lane GPU lease + VRAM headroom for cohabiting services

- status: done
- track: A (quality/infrastructure)
- priority: P1
- machine: any (implementation + mocked tests; live verification needs workstation)
- depends: —
- claimed-by: hermes-custom-build-001

## Objective
The RTX 4090 is shared: quarry DFT campaigns (10–22 GB, spiky), ollama
(email-poc pipeline + other local-model users, small models loaded on
demand), and occasional container jobs. Overnight 2026-08-19/20 the A3
pilot and ollama coexisted only by luck (pilot peaks hit 21.7 GB; ollama
models happened to be small), and queue workers serialized via ad-hoc
hand-rolled nvidia-smi preflights re-implemented per task (see
mission-control TASK-168 progress log). Make the lane explicit:

1. `qm/quarry/etiquette.py` (or extend QI1's module if it lands first):
   `acquire_gpu(owner, expected_gb, ttl_hours)` — file lease under
   `~/.local/state/gpu-lease/lease.json` (owner, pid, started, ttl,
   expected_gb). Stale detection: lease pid dead → break with a dated
   note. Refuse to acquire when a live lease exists; drivers exit with
   a clear "GPU lane busy (owner=X)" message instead of contending.
2. VRAM headroom: leased DFT processes set a CuPy memory-pool limit
   (default 16 GB, `--gpu-mem-gb` override) so ollama always has ≥6 GB
   to land models in. Document the floor in qm/README.md.
3. Shared preflight for queue workers: `scripts/gpu_preflight.sh`
   (repo) that checks the lease + nvidia-smi and exits 0/1 — replaces
   each worker's hand-rolled check; reference it from the TASK-168
   family and future campaign cards.
4. Wire `phase2_ladder.py` and `phase1_xiao_lasaga.py` to acquire the
   lease at start and release on exit (context manager; release on
   SIGTERM too).

## Acceptance
- Unit tests (mocked filesystem/pid) for acquire/refuse/stale-break;
  ruff + suite green.
- Both campaign drivers take/release the lease; running a second
  driver while one holds it exits non-zero with the busy message.
- CuPy pool limit verifiably applied (test may mock CuPy).
- qm/README.md documents the lane rules + ollama floor.

## Progress
- 2026-08-27T12:33:25-07:00 — (hermes-custom-build-001; profile=workstation) — DONE. Added the atomic process-bound GPU lease, dead-PID-only stale break receipts, fail-closed state-file handling, shared live preflight, 16 GB CuPy ceiling with a hard 6 GB ollama floor, and Phase 1/2 driver wiring. Full QM tests, Ruff, formatting, shell syntax, live `nvidia-smi`, separate-process contention, and SIGTERM release proofs are green.
- 2026-08-20 — card created (fable) after reviewing the overnight
  A3-pilot/ollama cohabitation and the TASK-168 preflight churn.

## Result

`quarry.etiquette` now owns a serialized, mode-0600 GPU lease record with
atomic replacement, exact owner/PID/start identity, declared TTL and expected
VRAM, and dead-PID-only stale breaking with UTC receipts. Live PIDs are never
evicted on TTL alone. Lease/lock/stale-receipt paths reject symlinks and special
files. Forked children cannot release a parent's lease; SIGTERM/SIGINT release
blocks reentrant delivery around the file lock and then preserves the original
signal behavior.

`phase1_xiao_lasaga.py` and `phase2_ladder.py` acquire the lane before GPU
imports, set CuPy's default pool to 16 GB, allow an explicit maximum-18-GB
`--gpu-mem-gb` override, and leave at least 6 GB nominal headroom for ollama.
A conflicting real Phase-2 process exits 1 with only `GPU lane busy
(owner=...)`; termination of the holder removes the lease. The executable
`scripts/gpu_preflight.sh` checks both the shared lease and live `nvidia-smi`.
The README and TASK-168 mechanism receipt point future campaigns at that shared
preflight.

Verification: **436 passed / 1 skipped** in the complete QM suite; whole-QM
Ruff and formatting passed; shell syntax and `git diff --check` passed. The
workstation preflight reported the live RTX 4090 and healthy NVIDIA driver path.
Sixteen focused tests cover acquire/refuse/stale-break, exact release, mocked
CuPy limits, separate-process CLI contention, fork ownership, real SIGTERM
reentrancy, and symlink fail-closed behavior. Independent adversarial review
found the fork-release, signal-lock reentrancy, state-symlink, and noisy CLI
failure classes; all four were fixed and regression-tested before closeout.
