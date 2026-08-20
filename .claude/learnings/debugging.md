# Debugging

### Saddle searches slide to trivial saddles from uncrossed-ridge guesses (2026-08-18)
First live hydrolysis TS hunt: the r(Si-Ow) scan ended (1.9 A) with energy
still rising, its endpoint became the Sella guess, and the search relaxed
back to a 208i cm^-1 water-wag saddle — mechanically valid (1 imaginary
mode, IRC converges) but chemically wrong (dG‡ 2.5 vs ~120 kJ/mol lit.).
Diagnosis signature: the saddle's driven coordinate far past the guess and
both quick-IRC ends identical to the reactant complex. Fixes in quarry:
scan_to_maximum (extend until interior max) + escaped-channel guard in
the phase-1 driver. A verified saddle is not necessarily the right saddle.

### CuPy pool fragmentation OOMs multi-hour GPU campaigns (2026-08-19)
First Phase-2 pilot: stage 1 (2 h of geomeTRIC steps on a 66-atom cluster)
ran fine at ~7 GB, then the first stage-2 DF-K gradient died with
cudaErrorMemoryAllocation on the 24 GB card — the CuPy memory pool had
fragmented across thousands of allocations. Fix: free_all_blocks() on the
default + pinned pools at every stage boundary and scan point
(phase2_ladder.trim_gpu_pool). Also: `cmd | tee log` reports tee's exit
code — a crashed driver looks like exit 0 unless `set -o pipefail`.
