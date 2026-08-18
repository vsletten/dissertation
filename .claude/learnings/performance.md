# Performance

### RTX 4090 vs 16C CPU, B3LYP/def2-TZVP, first real numbers (2026-08-18)
Si(OH)4·2H2O (15 atoms), no density fitting, cold GPU, quarry phase0_smoke:
CPU E/grad/Hess = 22.4/7.1/316.5 s; GPU = 54.3/14.0/197.6 s (0.41×/0.51×/1.60×);
energies matched to 8 decimals. Three artifacts depressed the GPU numbers:
(1) one-time CUDA/JIT warmup landed inside the first timed SCF, (2) no RI —
GPU4PySCF's published wins assume density fitting, (3) cuTENSOR missing, so
gpu4pyscf warned and fell back to CuPy contractions. smoke v2 fixes all
three (warmup pass, DF default, cutensor-cu12 in the [gpu] extra) and adds
--waters up to 24 to probe the 100+-atom regime; expect the survey's 2–5×
(SCF/grad) and 5–10× (Hessian) only there. Small jobs should just run CPU.
