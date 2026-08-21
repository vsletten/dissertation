# QI2-gpu-lease — single-lane GPU lease + VRAM headroom for cohabiting services

- status: ready
- track: A (quality/infrastructure)
- priority: P1
- machine: any (implementation + mocked tests; live verification needs workstation)
- depends: —
- claimed-by:

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
- 2026-08-20 — card created (fable) after reviewing the overnight
  A3-pilot/ollama cohabitation and the TASK-168 preflight churn.
