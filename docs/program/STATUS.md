# Program status log

*Newest first. One line per merged unit of work (or notable event),
appended by the PR that did it: `date — actor — what + pointer`.
*Deviations get a `DEVIATION:` note; details live in the card.*

- 2026-09-05 03:06 PDT (hermes-custom-build-001; profile=workstation) — A2a's remaining TightPNO matrix (TS/QZ, reactant/TZ, TS/TZ) is now supervised by `a2a-task220-psi4-remaining.service` (`1e0b30a25020425f8fe574deea951f08`) at pinned head `63d604d`; it waits at most eight hours for the active A3 recovery job, then applies load/memory/swap gates before the sequential Psi4 run. The continuation is capped at 30 h, 42 GiB RAM, 8 GiB swap, 16 cores, and nice 10; launcher SHA-256 `4a04f02e45bae12ab4fee79703da51e3ca20491a6c30783fa3dd539f8c111910`. Completion status file is `production-closeout/cc-calibration/psi4-remaining-20260905.status.json`.

- 2026-09-05 02:53 PDT (hermes-custom-build-001; profile=workstation) — A2a's reactant/cc-pVQZ TightPNO DLPNO-CCSD(T) completed in `6732.640926 s` at `-1184.6871167652455 Eh`; against the accepted same-geometry ByteQC canonical total `-1184.7350551217924 Eh`, DLPNO−canonical is `+125.86213783105343 kJ/mol`. This large absolute offset is recorded as method provenance, not misreported as the barrier gate; exact receipt/comparison hashes, independent arithmetic, identity checks, and complete service restoration are green in the A2a card. The approved next computation is TS/QZ DLPNO, never another canonical QZ slow-mode run.

- 2026-09-05 00:15 PDT (hermes-custom-build-001; profile=workstation) — A2a reactant/cc-pVQZ TightPNO DLPNO-CCSD(T) is durably running in finite system unit `a2a-task220-psi4-reactant-qz.service` (invocation `802711d984f643cdad51568286026c23`) from exact pushed head `b494bc5`; Psi4 1.11 converged DF-RHF and entered pair-domain work. A 25-hour dead-man and exit trap restore the intentionally isolated email/Honcho/watchdog/deploy-watch services. No DLPNO energy or canonical comparison is claimed before the atomic completion receipt appears; details and exact hashes are in the A2a card.

- 2026-09-04 21:36 PDT (hermes-custom-build-001; profile=workstation) — TASK-249 closed the A2a canonical reactant/cc-pVQZ bottleneck: the atomic CCSD(T) receipt completed after `105.570297 h` with triples `-0.0901434944124883 Eh` and total `-1184.7350551217924 Eh` (receipt SHA `f4aec84e...`). Exact source/input/restart identity, receipt arithmetic, 75,902-row resource telemetry, `32` focused CC tests, Ruff, diff checks, and full runtime restoration pass. A2a returns to four TightPNO receipts, the reactant/QZ DLPNO-vs-canonical delta, directive-adjusted barrier/focal-point/docs/store, and its single final PR; Victor's 2026-09-04 routing directive forbids any additional canonical CCSD(T)/QZ slow-mode campaign and sends TS/QZ through DLPNO.

- 2026-08-31 09:49 PDT (hermes-custom-build-001; profile=workstation) — A2a canonical reactant/QZ advanced to cumulative converged CCSD cycle 28 (`E_corr=-2.391227800954846 Eh`; checkpoint SHA `26f73c2a...`; marker SHA `834ea1bb...`), then an eighth exact 120-minute service again expired in slow-mode canonical `(T)` without a triples value, QZ receipt, or accepted partial total. Five-second cgroup sampling recorded 1,440 samples, `32.301662/34.147346 GiB` maximum current/peak memory, and `0.796658 GiB` maximum swap under unchanged 40 GiB/8 GiB/120m bounds. Runtime/services/QI2/GPU are restored. Resume from cycle 28 or a separately verified triples-capable route; six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-31 07:05 PDT (hermes-custom-build-001; profile=workstation) — A2a canonical reactant/QZ advanced to cumulative converged CCSD cycle 27 (`E_corr=-2.391227800875518 Eh`; checkpoint SHA `08f002f1...`; marker SHA `6d6bf46b...`), then a seventh exact 120-minute service again expired in slow-mode canonical `(T)` without a triples value, QZ receipt, or accepted partial total. Five-second cgroup sampling recorded 1,440 samples, `34.109318/34.113285 GiB` maximum current/peak memory, and zero swap under unchanged 40 GiB/8 GiB/120m bounds. Runtime/services/QI2/GPU are restored. Resume from cycle 27 or a separately verified triples-capable route; six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-31 04:11 PDT (hermes-custom-build-001; profile=workstation) — A2a canonical reactant/QZ advanced to cumulative converged CCSD cycle 26 (`E_corr=-2.391227800811323 Eh`; checkpoint SHA `1e4269b3...`; marker SHA `036c6ba1...`), then a sixth exact 120-minute service again expired in slow-mode canonical `(T)` without a triples value, QZ receipt, or accepted partial total. Five-second cgroup sampling recorded 1,440 samples, `33.930679/33.942219 GiB` maximum current/peak memory, and zero swap under unchanged 40 GiB/8 GiB/120m bounds. Runtime/services/QI2/GPU are restored. Resume from cycle 26 or a separately verified triples-capable route; six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-31 01:28 PDT (hermes-custom-build-001; profile=workstation) — A2a canonical reactant/QZ advanced to cumulative converged CCSD cycle 25 (`E_corr=-2.391227800548428 Eh`; checkpoint SHA `ac98ed4a...`; marker SHA `12ec08b5...`), then a fifth exact 120-minute service again expired in slow-mode canonical `(T)` without a triples value, QZ receipt, or accepted partial total. Five-second cgroup sampling recorded 1,440 samples, `29.753899/31.599800 GiB` maximum current/peak memory, and zero swap under unchanged 40 GiB/8 GiB/120m bounds. Runtime/services/QI2/GPU are restored. Resume from cycle 25 or a separately verified triples-capable route; six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-30 22:50 PDT (hermes-custom-build-001; profile=workstation) — A2a canonical reactant/QZ advanced to cumulative converged CCSD cycle 24 (`E_corr=-2.391227800306946 Eh`; checkpoint SHA `f8c56fea...`; marker SHA `4a8b2f84...`), then a fourth exact 120-minute service again expired in slow-mode canonical `(T)` without a triples value, QZ receipt, or accepted partial total. Five-second cgroup sampling recorded 1,440 samples, `32.121014/32.223732 GiB` maximum current/peak memory, and zero swap under unchanged 40 GiB/8 GiB/120m bounds. Runtime/services/QI2/GPU are restored. Resume from cycle 24 or a separately verified triples-capable route; six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-30 20:12 PDT (hermes-custom-build-001; profile=workstation) — A2a canonical reactant/QZ advanced to cumulative converged CCSD cycle 23 (`E_corr=-2.39122779940603 Eh`; checkpoint SHA `4d16b00d...`; marker SHA `c329d6a2...`) and then spent the full exact 120-minute service in slow-mode canonical `(T)` without a triples value, QZ receipt, or accepted partial total. Five-second cgroup sampling recorded 1,440 samples, `34.196129/34.200188 GiB` maximum current/peak memory, and zero swap under unchanged 40 GiB/8 GiB/120m bounds. Runtime/services/QI2/GPU are restored. Resume from the cycle-23 marker; six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-30 17:35 PDT (hermes-custom-build-001; profile=workstation) — A2a canonical reactant/QZ revalidated its bound convergence state and advanced to durable cumulative cycle 22 (`E_corr=-2.3912277985117507 Eh`; checkpoint SHA `d06324dd...`; convergence marker SHA `f7911bdc...`). Canonical `(T)` then consumed the full exact 120-minute service and timed out without a triples value, QZ JSON receipt, or accepted partial total. Five-second cgroup sampling recorded 1,440 samples, `32.974773/33.996983 GiB` maximum current/peak memory, and zero swap under unchanged 40 GiB/8 GiB/120m bounds. Runtime/services/QI2/GPU are restored. Resume from the cycle-22 marker; six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-30 14:58 PDT (hermes-custom-build-001; profile=workstation) — A2a canonical reactant/QZ CCSD converged durably at cumulative cycle 21 (`E_corr=-2.391227795318941 Eh`; checkpoint SHA `2a4ae4a2...`; bound convergence marker SHA `17123a41...`). Canonical triples then consumed the remaining exact 120-minute service window and timed out without a triples value, QZ JSON receipt, or accepted partial total. Five-second cgroup sampling recorded 1,440 samples, `29.684486/29.711086 GiB` maximum current/peak memory, and zero swap under unchanged 40 GiB/8 GiB/120m bounds. Runtime/services/QI2/GPU are restored. Resume from the convergence marker so the full next finite window goes to triples; six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-30 12:14 PDT (hermes-custom-build-001; profile=workstation) — A2a canonical reactant/QZ advanced from cumulative cycle 19 to crash-durable cycle 20 (`E_corr=-2.3912277920195737 Eh`; checkpoint SHA `6fa2f2d2...`) and deliberately refused before cycle 21 because `|dE|=1.17314944e-8 Eh` remains narrowly above the unchanged `1e-8 Eh` convergence threshold. Five-second cgroup sampling recorded 385 samples, `32.733540/33.998081 GiB` maximum current/peak memory, and zero swap under unchanged 40 GiB/8 GiB/120m bounds. No convergence marker, triples, QZ receipt, or partial accepted energy exists; runtime/services/QI2/GPU are restored. Resume exact cycle 20; six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-30 05:18 PDT (hermes-custom-build-001; profile=workstation) — A2a canonical reactant/QZ advanced from cumulative cycle 13 to crash-durable cycle 14 (`E_corr=-2.3912273260470793 Eh`; checkpoint SHA `447c6fd0...`) and deliberately refused before cycle 15. Independent five-second cgroup sampling recorded 386 samples, `34.076252/34.080013 GiB` maximum current/peak memory, and zero swap under unchanged 40 GiB/8 GiB/120m bounds. No convergence marker, triples, QZ receipt, or partial accepted energy exists; runtime/services/QI2/GPU are restored. Resume exact cycle 14; six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-30 04:08 PDT (hermes-custom-build-001; profile=workstation) — A2a canonical reactant/QZ advanced from cumulative cycle 12 to crash-durable cycle 13 (`E_corr=-2.3912265627460307 Eh`; checkpoint SHA `f029b3bc...`) and deliberately refused before cycle 14. Independent five-second cgroup sampling recorded 386 samples, `29.926502/30.762199 GiB` maximum current/peak memory, and zero swap under unchanged 40 GiB/8 GiB/120m bounds. No convergence marker, triples, QZ receipt, or partial accepted energy exists; runtime/services/QI2/GPU are restored. Resume exact cycle 13; six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-30 02:57 PDT (hermes-custom-build-001; profile=workstation) — A2a canonical reactant/QZ advanced from cumulative cycle 11 to crash-durable cycle 12 (`E_corr=-2.391225913040409 Eh`; checkpoint SHA `a4107fbc...`) and deliberately refused before cycle 13. Independent five-second cgroup sampling recorded 386 samples, `34.171780/34.176086 GiB` maximum current/peak memory, and zero swap under unchanged 40 GiB/8 GiB/120m bounds. No convergence marker, triples, QZ receipt, or partial accepted energy exists; runtime/services/QI2/GPU are restored. Resume exact cycle 12; six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-30 01:48 PDT (hermes-custom-build-001; profile=workstation) — A2a canonical reactant/QZ advanced from cumulative cycle 10 to crash-durable cycle 11 (`E_corr=-2.3912224736503305 Eh`; checkpoint SHA `9697cdf1...`) and deliberately refused before cycle 12. Independent five-second cgroup sampling recorded 414 samples, `29.707680/29.707832 GiB` maximum current/peak memory, and zero swap under unchanged 40 GiB/8 GiB/120m bounds. No convergence marker, triples, QZ receipt, or partial accepted energy exists; runtime/services/QI2/GPU are restored. Resume exact cycle 11; six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-30 00:35 PDT (hermes-custom-build-001; profile=workstation) — A2a canonical reactant/QZ advanced from cumulative cycle 9 to crash-durable cycle 10 (`E_corr=-2.3912200137100634 Eh`; checkpoint SHA `11474fe9...`) and deliberately refused before cycle 11. Independent five-second sampling of the exact system-service cgroup resolved the prior bogus terminal accounting receipt: 387 samples measured `32.450630 GiB` maximum current, `34.296589 GiB` cgroup peak, and zero swap under the unchanged 40 GiB/8 GiB/120m bounds. No convergence marker, triples, QZ receipt, or partial accepted energy exists; runtime/services/QI2/GPU are restored. Resume exact cycle 10; six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-29 23:24 PDT (hermes-custom-build-001; profile=workstation) — A2a canonical reactant/QZ advanced from exact cumulative cycle 8 to crash-durable cycle 9 (`E_corr=-2.3912034084094502 Eh`; checkpoint SHA `f35747a5...`) and deliberately refused before cycle 10. No convergence marker, triples, QZ receipt, or partial accepted energy exists; runtime/services/QI2/GPU are restored. The system-slice launch retained explicit 40 GiB RAM / 8 GiB swap / 120m bounds, but `systemd-run` reported an anomalous `256.0K` peak despite journal evidence placing Python in the intended cgroup, so resolve or independently measure memory accounting before resuming exact cycle 9. Six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-29 20:54 PDT (hermes-custom-build-001; profile=workstation) — A2a canonical reactant/QZ advanced from exact cumulative cycle 6 to crash-durable cycle 7 (`E_corr=-2.3911073991528538 Eh`; checkpoint SHA `ef553056...`) under a preflight-verified system-slice envelope that escaped the user-slice `systemd-oomd` failure mode while retaining finite 40 GiB RAM / 8 GiB swap / 120m bounds. The deliberate one-checkpoint stop finished in 32m19s at 34.2 GiB peak with zero swap and no oomd event; no convergence marker, triples, QZ receipt, or partial accepted energy exists. Both email pipelines/watchdog, Honcho health, QI2/dead-man, idle GPU, and clean worktree are restored. Resume exact cycle 7 under the same envelope; six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-29 19:43 PDT (hermes-custom-build-001; profile=workstation) — A2a canonical reactant/QZ cycle 7 was OOM-killed during its first resumed CCSD update at the exact 34.0 GiB cgroup ceiling, before any cycle-7 checkpoint or accepted energy. The byte-identical cycle-6 checkpoint remains authoritative (`E_corr=-2.3910759996876836 Eh`; SHA `c7d6decc...`); explicit cleanup restored both email pipelines/watchdog, Honcho health, QI2/dead-man, and the idle GPU. Resume exact cycle 6 only under a separately verified host-memory-safe envelope; six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-29 18:36 PDT (hermes-custom-build-001; profile=workstation) — A2a canonical reactant/QZ advanced from exact schema-v2 cycle 5 to crash-durable cycle 6 (`E_corr=-2.3910759996876836 Eh`; checkpoint SHA `c7d6decc...`) and deliberately refused before cycle 7. The finite 37m10s service peaked at 33.3 GiB plus 1.3 GiB swap; no convergence marker, triples, QZ receipt, or partial accepted energy exists. Both email pipelines/watchdog, Honcho health, QI2/dead-man, idle GPU, and clean worktree are restored. Resume exact cycle 6; six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-29 17:21 PDT (hermes-custom-build-001; profile=workstation) — A2a's physically verified one-checkpoint guard advanced canonical reactant/QZ from exact cycle 4 to cycle 5 (`E_corr=-2.3905665287527667 Eh`; checkpoint SHA `2c8a45f1...`) and refused before cycle 6. The deliberate continuation boundary used 27.8 GiB peak versus the prior 35.8 GiB OOM, emitted no convergence marker/triples/QZ receipt/partial energy, and restored both email pipelines/watchdog, Honcho health, QI2, and the idle GPU. Pushed implementation head `0d6fa69`; `488 passed, 1 skipped` plus whole-tree Ruff/format are green. Resume exact cycle 5; six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-29 16:02 PDT (hermes-custom-build-001; profile=workstation) — A2a's canonical reactant/QZ restart advanced durably from cycle 3 to cycle 4 (`E_corr=-2.390489726825349 Eh`; 4.59 GB checkpoint file SHA `4d6b0bb2...`). `systemd-oomd` then killed the finite service during cycle 5 at a measured 35.8 GiB peak; no cycle-5 state, convergence, triples, receipt, or partial energy was accepted. Explicit cleanup restored both email pipelines/watchdog, Honcho health, the QI2 lane, and the idle GPU. Resume exact cycle 4 under a one-checkpoint launch bound or separately verified memory-safe route; six receipts plus focal-point/docs/PR closeout remain.

- 2026-08-29 12:20 PDT (hermes-custom-build-001; profile=workstation) — A2a's canonical reactant/QZ continuation is now physically restartable: schema v2 atomically binds exact RHF state plus t1/t2, rejects legacy v1, validates energy/Fock physics and all checksums, and always reruns convergence. A bounded fresh launch persisted cycle 1 (`E_corr=-2.337220110032075 Eh`, manifest `efa84b8f...`); a second launch restored it without SCF replay and entered cycle 2. Pushed implementation head `134def0`; 29 actual-ByteQC tests and the complete `485 passed, 1 skipped` QM suite are green. No QZ receipt or partial energy is accepted; six receipts plus focal-point/docs/PR closeout remain, and all services/GPU are restored.

- 2026-08-29 07:33 PDT (hermes-custom-build-001; profile=workstation) — A2a now has fingerprint-bound, atomic ByteQC CCSD amplitude restart at pushed head `78ab878`: every completed post-DIIS cycle persists checksummed t1/t2 plus exact geometry/engine/driver/orbital and convergence provenance, while a bound convergence marker resumes directly into triples after post-CCSD interruption. Adversarial timing/collision/terminal-cycle defects are closed; 16 focused tests, 473 complete QM tests, and whole-tree Ruff/format/diff gates pass. Next is the bounded reactant/cc-pVQZ continuation, then canonical TS/QZ, four TightPNO receipts, focal-point gate, final docs, and one PR.

- 2026-08-29 06:17 PDT (hermes-custom-build-001; profile=workstation) — A2a's deliberately extended 180-minute canonical reactant/cc-pVQZ run completed four CCSD cycles (`E_corr=-2.39119870365418 Eh`) but did not converge or emit an atomic receipt before the finite cgroup stopped. Full restarts are now exhausted: the next continuation must add fingerprint-bound atomic amplitude restart/checkpointing before another QZ launch. All shared services and the GPU lane are restored; six CC receipts plus focal-point/docs/PR closeout remain.

- 2026-08-29 02:38 PDT (hermes-custom-build-001; profile=workstation) — A2a's corrected canonical reactant/cc-pVQZ restart again crossed the frozen-core transform ceiling, converged RHF/MP2, and reached CCSD cycle 2 before its finite 115-minute cgroup terminated. No receipt or partial energy was accepted; the identical 115-minute replay is exhausted, and the next continuation needs a longer finite GPU envelope or verified restart-capable route. All shared services, the QI2 lane, and the clean worktree are restored.

- 2026-08-29 00:03 PDT (hermes-custom-build-001; profile=workstation) — A2a's first canonical reactant/cc-pVQZ attempt exposed the cc-pVTZ full-auxiliary-block workaround as a 7.73 GB QZ OOM. Pushed head `28d801d` restores bounded DF blocks with a numerically verified frozen-core-safe transform, scoped monkeypatch, and source-hashed receipt identity; 466 QM tests and whole-tree Ruff gates pass. A 40-minute physical proof crossed the former OOM through converged RHF and into CCSD before deliberate stop, emitted no partial receipt, and restored services/QI2. Six CC receipts plus the focal-point/docs/PR closeout remain.

- 2026-08-28 22:37 PDT (hermes-custom-build-001; profile=workstation) — A2a's exact addition-TS/cc-pVTZ ByteQC canonical CCSD(T) receipt completed in `4477.6269 s`; with the accepted reactant/TZ receipt it gives a finite-basis canonical electronic barrier of `128.085164 kJ/mol`. Two of eight CC receipts are now complete; six receipts, the focal-point gate, final docs/store, and one PR remain. Runtime/services and the QI2 lane are restored.

- 2026-08-28 15:29 PDT (hermes-custom-build-001; profile=workstation) — A2a's barrierless I↔P follow-through is closed: fresh unconstrained release from exact product cell `(7,7)` moves `0.2510 A`, drops `0.9169 kJ/mol`, preserves exact typed P, and has zero finite-difference imaginary modes (lowest real `32.64 cm^-1`). Adversarial review rejected and quarantined the earlier `(2,7)` release with one `30.88 cm^-1` imaginary mode. Pushed head `be96273`; 457 QM tests and whole-tree lint/format pass; evidence v7 is 982 files / 10,300,380 bytes (`a22dbd1eb93ab21c...`). Remaining: production/method-shift/CC energetics and final one-PR documentation closeout using the accepted R↔I saddle as route TS.

- 2026-08-28 14:05 PDT (hermes-custom-build-001; profile=workstation) — A2a's 9x9 production coupled-coordinate scan is complete at 81/81 finite, residual-gated cells and decisively classifies I↔P as a barrierless shelf: the physical minimax path's boundary bottleneck is already `-44.250652 kJ/mol` relative to I, with no proton-first minimum or interior crest above 2 kJ/mol. Evidence manifest v6 covers 969 files / 10,019,351 bytes (`a30e22780fb64044...`); runtime and the QI2 lane are restored. Next is the required fresh downhill typed-product proof, then production/CC energetics using the accepted R↔I saddle; no cleavage saddle or final barrier exists yet.

- 2026-08-28 12:06 PDT (hermes-custom-build-001; profile=workstation) — A2a's production coupled scan is crash-durable at 65/81 exact cells. New cells `(6,6)` through `(6,8)` and `(7,8)` through `(7,7)` are finite, residual-gated hydrolyzed product with cleaved Si--Obr and water-owned H16; the bounded cgroup stopped during unaccepted `(7,6)`, and all shared services plus the QI2 lane are restored. The incomplete surface has no classifier, saddle, IRC, production/CC energy, barrier, or PR; resume at `(7,6)`.

- 2026-08-28 11:09 PDT (hermes-custom-build-001; profile=workstation) — A2a's production coupled scan is crash-durable at 60/81 exact cells. New cells `(6,1)` through `(6,5)` are finite, residual-gated hydrolyzed product with cleaved Si--Obr and water-owned H16; the bounded cgroup stopped during unaccepted `(6,6)`, and all shared services plus the QI2 lane are restored. The incomplete surface has no classifier, saddle, IRC, production/CC energy, barrier, or PR; resume at `(6,6)`.

- 2026-08-28 08:25 PDT (hermes-custom-build-001; profile=workstation) — A2a's production coupled scan is crash-durable at 46/81 exact cells. New cells `(4,4)` through `(4,8)` and `(5,8)` are finite, residual-gated hydrolyzed product with cleaved Si--Obr and water-owned H16; the bounded cgroup stopped during unaccepted `(5,7)`, and all shared services plus the QI2 lane are restored. The incomplete surface has no classifier, saddle, IRC, production/CC energy, barrier, or PR; resume at `(5,7)`.

- 2026-08-28 07:29 PDT (hermes-custom-build-001; profile=workstation) — A2a's production coupled scan is crash-durable at 40/81 exact cells. New cell `(3,0)` and cells `(4,0)` through `(4,3)` are finite, residual-gated hydrolyzed product with cleaved Si--Obr and water-owned H16; the bounded cgroup stopped during unaccepted `(4,4)`, and all shared services plus the QI2 lane are restored. The incomplete surface has no classifier, saddle, IRC, production/CC energy, barrier, or PR; resume at `(4,4)`.

- 2026-08-28 06:35 PDT (hermes-custom-build-001; profile=workstation) — A2a's production coupled scan is crash-durable at 35/81 exact cells. Five new cells through `(3,1)` remain finite, residual-gated hydrolyzed product with cleaved Si--Obr and water-owned H16; the bounded cgroup stopped during unaccepted `(3,0)`, and all shared services plus the QI2 lane are restored. The incomplete surface has no classifier, saddle, IRC, production/CC energy, barrier, or PR; resume at `(3,0)`.

- 2026-08-28 05:40 PDT (hermes-custom-build-001; profile=workstation) — A2a's production coupled scan is crash-durable at 30/81 exact cells. Six new cells through `(3,6)` remain finite, residual-gated hydrolyzed product with cleaved Si--Obr and water-owned H16; the bounded cgroup stopped during unaccepted `(3,5)`, and all shared services plus the QI2 lane are restored. The incomplete surface has no classifier, saddle, IRC, production/CC energy, barrier, or PR; resume at `(3,5)`.

- 2026-08-28 04:45 PDT (hermes-custom-build-001; profile=workstation) — A2a's production coupled scan is crash-durable at 24/81 exact cells. New cells `(2,0)` through `(2,5)` are finite, residual-gated, and typed as hydrolyzed product with cleaved Si--Obr while H16 remains water-owned; energies descend to `-67.750969 kJ/mol` relative to `(0,0)`. The bounded cgroup stopped during unaccepted `(2,6)`, and all shared services plus the QI2 lane are restored. The incomplete surface has no classifier, saddle, IRC, production/CC energy, barrier, or PR; resume at `(2,6)`.

- 2026-08-28 03:52 PDT (hermes-custom-build-001; profile=workstation) — A2a's production coupled scan is crash-durable at 18/81 exact cells and has completed both initial Si--Obr rows. All accepted cells remain finite and pass typed-topology, H16-water-ownership, intact-Si--Obr, and `1e-4 A` residual gates. The bounded systemd cgroup stopped during unaccepted `(2,0)`, cleanup marked the status interrupted, and all shared services plus the QI2 GPU lane are restored. No classifier, saddle, IRC, production/CC energy, barrier, or PR exists; resume at `(2,0)`.

- 2026-08-28 00:58 PDT (hermes-custom-build-001; profile=workstation) — A2a's production coupled scan is crash-durable at 3/81 exact cells; all three have finite r2SCAN-3c energies, typed associative topology, and distance residuals below `4.35e-6 A`. The bounded worker stopped before classification after catching the email watchdog defeating GPU isolation; the external launcher now stops/restores that timer and resumes at `(0,3)`. No saddle, IRC, production/CC energy, barrier, or PR exists yet.

- 2026-08-27 23:57 PDT (hermes-custom-build-001; profile=workstation) — A2a's commissioned coupled-coordinate cleavage redesign now has a tested, resume-safe 9x9 r2SCAN-3c scan framework with fixed H16 mapping, atomic per-cell geometry/energy/topology receipts, and a topology-backed three-outcome classifier. Real I/P mapping and the complete 451-test QM suite pass; the production grid and classifier-selected saddle/shelf follow-through remain unrun and no barrier was emitted.

- 2026-08-28 — (hermes-8d9e30fee91548; profile=laptop) — E2b runs seven
  32-replica muscovite volumes through 559,872 sites: the old-first-step
  staircase stabilizes at 49,152 sites and stays stable across two larger
  volumes; 559,872 is the completed host ceiling after the 1,327,104-site rung
  exits by `SIGKILL`. Full bootstrap bands, isotope-release distributions,
  exact-seed replay receipts, and scaling are in
  `docs/program/results/E2b-grain-size-sweep.md`.

- 2026-08-28 — (hermes-8d9e30fee91548; profile=laptop) — A5 Phase 1 runs
  4×64 finite-defect kaolinite-aging ensembles: all replicas are defect-free
  by step 90 at the latest, fresh/aged apparent-rate drops grow 36.36×→363.06×
  with initial damage density, spectra transiently broaden then down-shift
  three orders, and BET-like area rises while rate falls. The honest
  mechanism-tier verdict is partial plausibility (1.53–2.56 log units), not a
  calibrated field attribution. `docs/program/results/A5p1-aging-study.md`.

- 2026-08-27 — fable — Track E phase-3/paper-1 path carded: E3a
  classical-NEB barrier campaign (READY — the "uncomputed middle": Nteme
  replication gate, then dehydroxylate/defect-zone/trap/Xe barriers), E3b
  CP2K spot checks (blocked: E3a), E2b grain-size sweep (READY — closes
  E2's 768-site statistics deviation), E4 1998-comparison (blocked: E3a,
  E2b, A8; the paper-1 figure set). Also A5p1-aging-study (READY): M5 named it
  but it was never carded. Context: board starvation 08-25→08-27 was a
  feeder machine-parser bug, fixed in mission-control 387300f.

- 2026-08-27 — (hermes-8d9e30fee91548; profile=laptop) — C2 ships a
  strain-seeded passive-film CTMC deck plus a deterministic 7×32-replica study
  and three-size scan: initiation rises 0/32→32/32, mean pit count
  3.47→36.22, and mean largest cluster 1.25→7.28 across the potential sweep;
  descriptive induction plots, exact-deadline CSVs, and a metastable
  current-rate series are archived in `docs/program/results/C2-corrosion-deck.md`.

- 2026-08-27 — (hermes-custom-build-001; profile=workstation) — QI2 makes the RTX 4090 a single explicit lane: atomic process/PID lease, dead-PID stale-break receipts, fail-closed symlink-safe state, SIGTERM/fork hardening, shared `scripts/gpu_preflight.sh`, and a 16 GB CuPy pool cap preserving 6 GB for ollama. Phase 1/2 now refuse contention cleanly; 436 QM tests, Ruff/format, live 4090 preflight, separate-process busy exit, and signal release are green.

- 2026-08-26 — (hermes-custom-build-001; profile=workstation) — A2a's corrected full-system local dimer exhausted all 120 bounded translations inside the exact climb guard with strictly positive curvature (`+0.0274 eV/A²` minimum), finished wall-supported, and persisted no candidate or downstream number. Pushed head `95ea8df`; 425 QM tests and whole-QM Ruff/format pass; evidence manifest `9e163528bb43...` covers 478 files. A2a is scientifically blocked pending a new mechanism/path decision, not another local Sella/dimer retry.

- 2026-08-26 — hermes-workstation — A2a I↔P now has a tested neighbor-radius
  flat-bottom trust envelope, hard pre-PES escape guard, capped Sella step, and
  mandatory cuTENSOR GPU preflight. The 300-step cleavage localization stayed
  inside the local band neighborhood but cycled near `0.3--0.6 eV/A`, so no
  saddle or downstream number was emitted. Pushed head `9b9edb2`; evidence
  manifest `8542944ade6f...` (461 files). Next: bounded full-system local
  eigenvector-following/dimer localization, not another active-core retry.

- 2026-08-25 — hermes-workstation — A2a's internally conditioned production
  path now closes R↔I with a stable 714.83 cm^-1 reaction saddle and exact typed
  full IRC; I↔P pre-relax/climb/conditioning also converge, but its active-core
  Sella search escapes uphill into an SCF failure before a candidate is emitted.
  Commit `8ced110`; durable evidence manifest `82af53d69809...` (457 files).

- 2026-08-24 — hermes-workstation — A2 production infrastructure is pushed and
  fail-closed on the first real target: exact r2SCAN-3c D4+gCP optimization,
  same-surface finite-difference Hessians, repaired Sella/ASE full IRC, and
  verified Psi4/ByteQC calibration adapters. The si-neutral minima/TS have
  indices 0/0/1, but a real 260-step full IRC connects the banked saddle to the
  associative basin on both sides; directed Sella reconverges it. No production
  or CC barrier was fabricated. A2 is blocked on READY path-rebuild card A2a;
  evidence manifest `efdcfc48c533...` (47 files, 791,566 bytes).

- 2026-08-24 — hermes-workstation — A1g's hash-pinned bridge-side-H3O+ /
  distinct-neutral-attacker Si reactant is a zero-imaginary minimum, but its
  independent product converges outside the hydrolyzed family: Si--Ow 2.643 A,
  attacker H21 unassigned, basin `(False,False,True,1,7,0)`. This is a
  conclusive mechanism-v3/gate-v3 product rejection with no computational
  failure and no eligible Al/barrier/ordering. The only acid follow-through is
  A1i's one calibrated production-tier revisit after A2; no more SVP topology
  retries.

- 2026-08-24 — hermes-workstation — A1f's one four-water separated-donor /
  neutral-attacker Si path is a conclusive reactant-family rejection: the
  unconstrained B3LYP/def2-SVP/DF endpoint moved physical H16 from outer donor
  O15 to bridge-side relay O22, changing exact solvent occupancies from
  `(3,2,2,2)` to `(2,2,3,2)` while keeping the attacker neutral, bridge intact,
  and ownership valid. No computational stage failed; no Al/barrier/ordering
  was emitted. Mechanism-v2/gate-v2, 360-test proof, hash-bound archive, and
  READY A1g bridge-side-H3O+ follow-through land in the A1f card.

- 2026-08-24 — hermes-workstation — A1e post-merge activation closeout
  independently rehashed all 24 campaign evidence files, reconfirmed four
  typed Si endpoint rejections with zero failed/blocked/running paths and READY
  A1f follow-through, and passed 341 QM tests plus CLI/diff gates. A narrow
  follow-up normalizes the review-added regression test so whole-tree Ruff lint
  and formatting are both green on the final merged behavior.

- 2026-08-24 — hermes-workstation — A1e's exact concerted hydronium-as-
  nucleophile campaign is model-valid NO-GO: all four 3/4-water bridge-chain
  and compact-cycle Si paths are conclusive typed endpoint rejections, with
  zero optimizer/NEB/Hessian/IRC failures and therefore no eligible Al,
  barrier, or ordering. Mechanism-v1/gate-v1, full bidirectional Sella IRC,
  physical-H relay ownership, hash-bound receipts, and 338-test proof land with
  narrower READY `A1f-acid-neutral-water-attacker-relay`, which separates H3O+
  donor from the neutral-water nucleophile.

- 2026-08-24 — hermes-workstation — A1d closes A1c's sole inconclusive
  three-water seed and the finite pre-equilibrium: one hash-pinned fresh
  B3LYP/def2-SVP/DF continuation converged after 54 steps, then failed the exact
  protonated-bridge occupancy gate. Effective result is 16/16 conclusive Si
  rejections, so the pre-equilibrated acid bridge is model-valid NO-GO; no Al
  screen or unmatched ordering was emitted. A1c/A1d are done and executable
  `A1e-acid-concerted-hydronium-relay` is ready.

- 2026-08-24 — hermes-workstation — A1c's finite 16-seed Si conformer
  campaign is honestly BLOCKED, not falsely NO-GO: 15 B3LYP/def2-SVP/DF
  endpoints converged and failed exact protonated-bridge occupancy, while the
  three-water bridge-donor chain exhausted its 160-step bound. The hardened
  ensemble driver, hash-bound external archive, and executable A1d one-seed
  convergence closeout are in `A1c-acid-microsolvation-conformers`; no Al
  screen, barrier, or acid ordering was emitted.

- 2026-08-23 — hermes-workstation — A1b matched 3--6-water Si-acid
  deterministic-shell survey is honestly BLOCKED: all four converged
  B3LYP/def2-SVP/DF endpoints deprotonate Obr before any Hessian/TS work.
  Driver support, physical-H/shell identity gates, non-destructive reactant
  receipts, hashes, and executable follow-up card
  `A1c-acid-microsolvation-conformers` are in `qm/ACID_MECHANISMS.md` and
  the A1b card; no unmatched Al ordering was emitted.

- 2026-08-23 — fable — A2 unblocked: ORCA descoped from the platform
  (Victor's call, licensing). Calibration layer is now Psi4 1.11
  DLPNO-CCSD(T) + ByteQC canonical GPU CCSD(T), cross-validating; both
  installed and smoke-verified on the workstation (water dimer/cc-pVTZ:
  canonical CCSD(T) agrees across engines to 1e-8 Eh; DLPNO −0.53
  kJ/mol off canonical). Runbook `qm/CALIBRATION.md`; SURVEY §§2.3–2.4/§6.4
  amended; A2 card re-scoped → READY (feeder-eligible).

- 2026-08-22 — hermes-workstation — A1b partial closeout: validated one-water
  al-acid sequential hydrolysis at ΔG‡(298) = 82.193 kJ/mol (19.645 kcal/mol),
  with exact one-mode/quick-IRC gates for addition and cleavage. DEVIATION:
  one-water si-acid is not an unconstrained protonated-bridge minimum in gas
  phase or PCM, so the card is blocked on 3–6-water microsolvation rather than
  publishing a constrained fake barrier. `qm/ACID_MECHANISMS.md`.

- 2026-08-22 — omnibus supervisor — reconciled board metadata after merged
  work: C2 is ready after B3, A3 is ready for its remaining family campaigns,
  and the completed D2a card now has its missing PLAN.md row.

- 2026-08-22 — hermes-laptop — D3 DONE: a 24×24 open-system H₂ grain deck
  now exercises source deposition, thermal+tunneling diffusion, LH reaction,
  H/H₂ desorption, and seeded 20% quenched deep sites. The 16-replica sweep
  has the published finite window (efficiency 0.310 at 6 K, ~0.96 at 12–14 K,
  0.559 at 20 K, 0.055 at 26 K). DEVIATION: D2a surface rates are explicit
  NO-GO; D2b + blocked D3b cards preserve the CO follow-through honestly.
  `docs/program/results/D3-ice-mantle-deck.md`.

- 2026-08-21 — fable — A3 pilot cell DONE: embedded Si–O–Si neutral
  hydrolysis ΔG‡ = 205.7 kJ/mol (free dimer 113.0; +92.7 lattice shift,
  matching Pelmenschikov's ~205 embedded first-rupture) via a verified
  two-step mechanism (pentacoordinate-Si intermediate, rate-limiting
  bridge rupture, 128i, IRC-connected). CALC-002 neutral computed.
  Campaign recipe + workstation ops prerequisites in the card Result;
  family ladders are factory work per the 2026-08-20 handoff.

- 2026-08-22 — hermes-laptop — E2 DONE: scheduled full muscovite mechanism
  now distinguishes pristine/extended galleries, explicit locally driven
  delamination, ⁴⁰Ar/³⁹Ar/³⁶Ar, and octahedral recoil traps; deterministic
  release/age/contamination products and honest ±5 kcal proxy bands land at
  `docs/program/results/E2-muscovite-full-mechanism.md`. E1's 500/700 °C gates
  remain green.
- 2026-08-22 — hermes-laptop — B5 DONE: ordered piecewise-isothermal CTMC
  schedules now advance through exact wall-time boundaries, recompile full
  thermo/propensity tables + Fenwick trees, and run through native CLI and
  deterministic ensembles. Exact T1 replay, seeded T2 statistics, v1 20k-step
  parity, workspace tests, format, Clippy, and independent review are green.
- 2026-08-22 — hermes-laptop — A5 Phase 0 DONE: long-time rate spectra,
  projected-geometric versus BET-like surface accounting, and core-owned
  per-site exposure ages landed with deterministic alias-aware tracking. The
  finite-defect kaolinite smoke exhausts 33 fast sites by step 40 and cuts bulk
  propensity 286.96×; 96 tests, bitwise parity, Clippy, and paranoid replay are
  green. `docs/program/results/A5p0-aging-observables.md`.
- 2026-08-22 — hermes-laptop — B4 DONE: Rayon-parallel deterministic replicas,
  full distributions/means/bootstrap CIs, and declarative state-count, event-rate,
  rate-spectrum, cluster-size, and geometric/BET-like area outputs; Ising Binder
  dogfood + Kossel spectrum-broadening golden green, with all 88 Petra tests and
  Clippy clean. `petra/docs/OBSERVABLES.md` records the A5 exposure-age seam.
- 2026-08-21 — hermes-laptop — B3 DONE: Petra now runs synchronous CA,
  asynchronous Metropolis, and discrete-time PCA on schema-v2 grids; Conway,
  Ising (`Tc=2.273230` vs `2.269185`), and SIR (`R0*=1.017094`) analytic gates
  are green with CTMC bitwise parity preserved.
- 2026-08-21 — omnibus supervisor — reconciled merged-main board drift:
  B2 and QI3 are done; B3, B4, and B5 are now ready after B2 PR #33.
- 2026-08-20 — hermes x2 + fable — B2 DONE: petra engine behind the
  UpdateStrategy trait, schema v2 complete, bitwise parity green on
  three gates (incl. shipped-kossel full-lifetime). Two-worker relay +
  supervision closeout; 7 self-review findings all dispositioned and
  fixed. B3 unblocked.
- 2026-08-20 — hermes-laptop — QI1 campaign-driver etiquette done: all QM
  entry points now cap inherited numerical-library threads before heavy imports,
  apply best-effort niceness, and self-tee durable logs; 162 tests, Ruff, direct
  CLI smoke checks, and an unset-environment Phase-2 dry run are green.
- 2026-08-20 — hermes-cloud-default — QI3 Fenwick cancellation hardening:
  impossible zero totals now rebuild from authoritative positive leaves; a
  deterministic `1e12`/`1e-6 s⁻¹` regression, all 43 Petra tests, and an
  explicit 1000-step Kossel `--paranoid` run are green.

- 2026-08-20 — hermes-workstation — E1 exact-main lint closeout: restored the
  one-line `strain_gates.rs` identity-op fix that the prior squash merge omitted;
  all Rust/Python/lint gates plus both full trajectories and generated-product
  equivalence are green.
- 2026-08-20 — hermes-workstation — E1 Ruff-format closeout: reformatted the
  cloud review-fix commit's long SVG `sy()` expression so merged-main Ruff
  0.15.1 check/format gates are clean; no simulation semantics changed.
- 2026-08-20 — hermes-workstation — E1 merged-main closeout: removed four
  default-Ruff E731 violations from the muscovite SVG post-processor with
  typed local helpers, cleared current Clippy findings, and reverified both
  full trajectories plus byte-identical CSV/SVG/JSON products.
- 2026-08-20 — hermes-workstation — E1-muscovite-deck done: schema-v1
  500/700 °C paired-OH + divacancy decks, Hames–Bowring cylinder `D/a²`
  inversion, reproducible CSV/SVG/JSON products, and qualitative gates green
  (rise 2.13×/13.55×; fall 275×/32.3×). Result:
  `docs/program/results/E1-muscovite-phase1.md`. Added executable E2 and QI3
  cards for the full mechanism and Fenwick dynamic-range hardening.
- 2026-08-20 — fable — A3-barrier-ladder claimed + platform half landed:
  crystallographic cluster builder (`qm/quarry/crystal.py`, deck-cell
  parser + bond-valence termination/charge + n_intact pruning + frozen
  shell), PHVA frequencies for frozen clusters, and the per-cell campaign
  driver (`qm/scripts/phase2_ladder.py`). Pilot cell oss-neutral running
  on the GPU; barriers land in follow-up per-family PRs.
- 2026-08-20 — fable — Card QI1: campaign drivers must enforce their own
  thread caps + niceness (a live worker probe ran 22 threads past the
  prompt-level cap — prompts don't bind subprocesses).
- 2026-08-20 — fable — Board refresh: B1+D2a+A7 done, B2 READY (keystone
  merged), A3 active. PROTOCOL amendment: incremental bookkeeping
  (the D2a iteration-ceiling lesson). Hermes worker-profile max_turns
  raised 90 -> 1000 fleet-side (infra, not repo; the raise gives workers
  more headroom, but the PROTOCOL closeout rule still applies — budget
  the last ~15% of turns so a capped worker leaves legible state).
- 2026-08-20 — local-hermes + fable — D2a DONE: 6-reaction astro rate
  campaign on the 4090. Verdict: GO gas-phase / NO-GO surface-LH rates
  (4–5 orders low without explicit-surface treatment). First non-Claude
  runtime to work the board end-to-end; worker hit iteration ceiling at
  closeout, fable finished bookkeeping. results/D2a-rate-reproduction.md.
- 2026-08-18 — hermes-workstation — A7-kinetics-database done: schema v0,
  36 P&K silicate records / 74 mechanisms, stdlib validator, and a
  49-source 2004–present expansion inventory under `kinetics-db/`.

- 2026-08-19 — hermes cloud worker — B1-schema-rfc done:
  `petra/docs/RFC-001-DECK-V2.md` written (deck schema v2 + UpdateStrategy
  trait + v1 compat shim + determinism/RNG contract + B2/B3 test plan;
  corrosion + ice-mantle showcase fragments). B2 unblocked pending Victor's
  review ack of the RFC.
- 2026-08-18 — fable — Card A7-kinetics-database added (machine:any,
  macbot-shaped): schema + extraction groundwork for a machine-readable
  successor to Palandri & Kharaka 2004. Pointer task queued in
  mission-control for laptop workers.
- 2026-08-18 — fable — Program tracking system created (PLAN, PROTOCOL,
  STATUS, 10 seed cards, inbox). Board opens with B1-schema-rfc as the
  keystone READY card.
- 2026-08-18 — fable — Omnibus program plan + three scoping docs merged
  (PR #20); root README refreshed (PR #21).
- 2026-08-18 — fable — Phase 1 first barrier: si-neutral ΔG‡(298) =
  113.0 kJ/mol (27.0 kcal/mol) vs X&L ~29, full pipeline (scan → product
  construction → CI-NEB → Sella → gates → quasi-RRHO), results +
  provenance archived. PRs #16–#18. TASK-164 queued for al-neutral.
- 2026-08-18 — fleet worker — qm/tests/test_phase1_driver.py landed
  (driver test coverage, via mission-control drain).
- 2026-08-18 — fable — Phase 0 closed: quarry package (PR #12), smoke
  v2 (PR #14); measured 4090: 74× analytic Hessians vs 32-thread CPU,
  identical energies. GPU owns all frequency work.
- 2026-08-19 — fable — qm/HANDOFF.md rewritten as the Phase-2 edition
  (barrier-ladder build plan: crystallographic cluster builder from the
  kaolinite deck cell, per-cell campaign driver, by_count emission,
  CALC-002..005 closure); new READY card A3-barrier-ladder.
- 2026-08-19 — fable — Track E opened: Ar-in-muscovite scoping doc
  landed (docs/scoping/ar-muscovite.md — the 1998 Sletten & Onstott
  circle, Nteme 2022–24 barriers, the uncomputed middle, Villa
  arbitration). New cards: E1-muscovite-deck (READY, published numbers
  only) and B5-execution-schedule (blocked: B2 — piecewise-isothermal
  T(t), the one engine feature the flagship needs).
- 2026-08-20 — Fable — Autopilot engaged: board_feeder.py (mission-
  control, hourly) auto-enqueues READY/unblocked cards as queue pointer
  tasks; omnibus-supervisor Hermes cron (daily 07:50) does bookkeeping
  + stale claims + Telegram digest. PROTOCOL.md amended. Humans and
  Fable are now exception handlers, not schedulers.
- 2026-08-20 — fable — card QI2-gpu-lease (READY, P1): explicit
  single-lane GPU lease + 16 GB cupy cap + shared worker preflight,
  after reviewing overnight A3-pilot/ollama cohabitation (22 ollama
  loads, all-GPU, zero errors — but pilot peaked 21.7 GB; luck, not
  design) and TASK-168's per-tick hand-rolled preflights.
- 2026-08-20 — fable — PROTOCOL.md merge-doctrine paragraph aligned to
  POLICY v7 (open PR and STOP; watcher owns lifecycle; human-merge
  label) — the stale auto-merge/Sourcery wording was flagged by the
  cloud Hermes worker during its enrollment (it correctly refused to
  follow it). Companion fixes landed in `vsletten/mission-control`
  commit `a2edf91`: `factory/hermes/cron/queue-drain.md` now points to
  POLICY v7, and service activation moved to
  `factory/services/activate.sh` for cloud-safe installation.
- 2026-08-20 — fable — Victor recovered his complete 1999 dissertation
  archive from the NAS (thesis LaTeX, the original C-language model,
  88 timestamped MC campaign runs, 86 Gaussian B3LYP logs incl.
  oxalate/ligand systems). New card A8-thesis-archive-intake (READY,
  workstation): catalog, curate golden data into legacy/, build the
  1999-vs-2026 DFT barrier ledger.
