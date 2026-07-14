# Movie-frame golden reference (`step{5000,10000,15000}.msi`)

The primary golden capture (TASK-004, `mission-control/projects/kmc/golden/`)
ran with `msteps = 1000000 > nsteps = 20000`, so the C++'s movie-frame path
(`mckaol.cpp`: `if (msteps && i && i % msteps == 0) writeMSI("step<i>")`)
**never fired** — the M7 movie writer would otherwise be gated by nothing.

These three frames were captured for TASK-019 from a scratch rebuild of the
read-only C++ model (`dissertation/main/model`, sources copied verbatim, no
edits — same recipe as TASK-004/TASK-018):

- toolchain: g++ (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0, stock `Makefile`
- inputs: the golden `data.cell` / `data.lattice` / `data.rxn`, plus
  `../../inputs/data.sim.msteps5000` — identical to the golden `data.sim`
  except `msteps = 5000`
- provenance check: the run's `start.msi` and `end.msi` hashed identical to
  the primary golden capture (`68bc498e…` / `c044cc3a…`), and its
  `step0.dat`/`step19000.dat`/`end.dat`/`start.xyz`/`surf*.out` byte-matched —
  `msteps` only adds snapshots; the trajectory is untouched (output never
  feeds back into the RNG or the state). So these frames are renderings of
  the *same* golden trajectory at steps 5000/10000/15000.

SHA-256:

```
5a569833f66b977d8dc1fe2495ec2e0e12f947d4f75631a6c9493b9e87600021  step5000.msi
83899d48d1b7d5e0d1e9c30fb2225b455d5b3a04a83e052db3ff37d11de9935a  step10000.msi
1986ad66cd52f4270b00ddb49e25283076b7a6b84382a740f3c61acc59f5b5a4  step15000.msi
```

Gated by `crates/mckaol-cli/tests/golden_m7.rs::movie_frames_match_the_cpp`.
