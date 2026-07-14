# TASK-004 Golden Capture Manifest

## Toolchain
- compiler: g++ (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0
- build command: `make clean && make -j1`
- binary: `.scratch/task-004/dissertation-model/mckaol`

## Run configuration
- working directory: `.scratch/task-004/dissertation-model`
- command: `./mckaol > run.log 2>&1`
- source-input mode: copied source tree from `dissertation/main/model`, no source edits
- data files used: `data.sim`, `data.rxn`, `data.cell`, `data.lattice`
- `data.sim` values used in this run:
  - `nsteps=20000`
  - `wsteps=1000`
  - `msteps=1000000`
  - `seed=-2` (upstream bug swallows seed)
  - `drawbonds=1`
- rationale: truncated run for runtime control; deterministic behavior confirmed by duplicate execution with identical hashes.

## Inputs
- `projects/kmc/golden/inputs/data.sim`
- `projects/kmc/golden/inputs/data.cell`
- `projects/kmc/golden/inputs/data.lattice`
- `projects/kmc/golden/inputs/data.rxn`
- preserved presets: `data.lattice.100x10`, `data.lattice.50x10`, `data.lattice.500x10`

## Inputs SHA-256
a7db7cb729ecfc74c9e7006b398a2d071d26f12984f46c524d6f118e4c139dfe  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/inputs/data.sim
4a93dde3656a8ca9f05c5be9fbd01db8ca5b73d2463beae628e5f0aaf335c362  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/inputs/data.cell
27cc032af9d4df38d6da37b7966e43e0c08257f9e181eccc13aeab87aff472b8  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/inputs/data.lattice
964805da201d936c9028a4c102a03a64cc27d4ce1cda1c45d08bd2d559cd2899  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/inputs/data.rxn
15e6e39a45f8e34ca9a7bd23d056771dd048e298c194a9bdac76ec9c693a1ad3  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/inputs/data.lattice.100x10
1dead3b6ab7757b4a0c47006febd56b01adde684c7798ea69f8b593067eb1169  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/inputs/data.lattice.50x10
a7cc26f8847af4de073f73be2b8ca3ca414d47fe02f447c85be7479cf6ceaccb  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/inputs/data.lattice.500x10

## Outputs
- outputs/start.msi
- outputs/start.xyz
- outputs/end.msi
- outputs/end.dat
- outputs/surfSi.out
- outputs/surfAl.out
- outputs/step*.dat (all 20 files present)
- Note: no output/end.xyz file is produced by this code path in this run.

## Output SHA-256
68bc498e8e3d107e175f29b5e39d509d93c15bd34731ff76eca4fc300a12295a  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/start.msi
d2b4af42866cfb993a2be302c241371de981e723a8c119956870e59dc2026821  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/start.xyz
c044cc3af83089574f672544b5fd5a5b03da0b016fa6115634ba8d552fe648b0  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/end.msi
ef7e435eb91ea9814677731490ce1d7e5e1bb034ec542cc6025d59c6a8ebf00d  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/end.dat
d11b6aeeefc4e3c8f1d564a2392b919ce38fbab0e7f103bd072de6f03740145f  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/surfSi.out
ecb69c9da80b51cb69a5a66b665e53a7dcf9ba5ba5a1074fbd752e14d11a2996  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/surfAl.out
84824a70de544e51c7614244a137ad78c7f0dcbafad777d85c479f7303a045e0  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step0.dat
5518c4d30040f81f5760859a9adaedbbf93fe31d0dd08adec226b599110bcd2d  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step1000.dat
24e59b955bf2f5a091e2d7caa9b5882349749107d86610f3b3545411a9b37b31  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step10000.dat
b94c343d6afb588b81aeeccf28f1ac9d04e62828cf97a21dde1fda2acd8982aa  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step11000.dat
8de9a5b1d330c77c36512f7b1ea7824573bfa8108b6812d4907e5c88325fc40b  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step12000.dat
64f69c5d0214f4ff7725c109945609212ebbc03eb8b264ea25e1f835ef927002  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step13000.dat
bedb16c39561c2e23c082824e7782dd18fb489cda958d54ff60fc846dabbe9fe  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step14000.dat
957b2532cc0f7d064d64be93ff57ef38d3bf98b8689d728d1baaf91627b05cf9  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step15000.dat
504f99a34947b74083a69ca378ab9c7a969910455b8e3fad42aa6a578b243b7e  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step16000.dat
cd84c7cbee2e216c854b5967f95f3b3a3b5d5412a9d8736a8b130461c5b53ad0  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step17000.dat
bead949eac87374f535e691dd0ed836526de8e18a8dc66269490be870423a053  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step18000.dat
35ce7e87137dc185881177ea0233b553b94a9fb3b2f2205fa6fec72acdb96a3f  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step19000.dat
8fadbb7d73fdfdf1a027aced4c47e496e84e232236157024d6f1a81cd6608630  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step2000.dat
7a38a71ef73ab388000f0b2b17315de3e5763536e63644c688e39bfe5c0b6f6e  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step3000.dat
cdc80f3f5593c2860f25ecfb8565eb6efedc5b8312adb18352c4e99fc87edda8  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step4000.dat
70e166d8995ff1ea75f8652cb2ea34abdeaebbdc6434c85200d95e0a4e449147  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step5000.dat
5c0505b3cfafae9d9b7e6c48730feae1bbe39073742f1570da9eabb5217e85b6  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step6000.dat
37978341dd629259462f2e5d40bb4b3228cf8ae522cd1c5b40e363539e7d9481  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step7000.dat
35f8ec41a286ce6bcf50b7bf02ec3e232fd04a728ec5113de236b427ac74ff0c  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step8000.dat
c8a655429974503117119624d6946c3b1b3cf56f95d1031cf38a97018c6639aa  /mnt/data/vsletten/src/vsletten/mission-control/main/projects/kmc/golden/outputs/step9000.dat

## Step file selection policy
- Total step data files: 20
- Policy (first 10 + last 10 + every 100th) is equivalent to including all files at this scale.

## Reproducibility
- Ran a second full run in the same scratch tree with identical config.
- Byte comparison (`sha256sum` over tracked outputs) was identical (`reproducible`).
