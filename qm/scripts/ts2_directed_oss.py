#!/usr/bin/env python
"""One-off TS2 recovery for the oss-neutral pilot (receipt in run log).

Sella (internal, undirected) slid from the rupture crest into the
pentacoordinate intermediate twice. Rebuild the crest by constrained
optimization at r(Si-Obr)=2.46 with the new bonds held, then run the
directed Cartesian Sella (#40 machinery: reaction_path_vector over the
reactive core) and leave ts.xyz for the campaign driver to verify.
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    bootstrap_cli(
        "ts2_directed_oss",
        default_run_root=Path(__file__).resolve().parent.parent / "runs",
    )


import numpy as np  # noqa: E402

from quarry.clusters import (
    water,  # noqa: E402
)
from quarry.crystal import attack_complex, from_deck_cell  # noqa: E402
from quarry.pipeline import DftSettings  # noqa: E402
from quarry.ts import constrained_scan, find_ts, reaction_path_vector  # noqa: E402
from scripts.phase2_ladder import (  # noqa: E402
    load_xyz,
    log,
    preload_cutensor,
    save_xyz,
)

RUN = (
    Path(__file__).resolve().parent.parent
    / "runs/phase2/oss-neutral-n4-s2-b3lyp-def2-svp"
)
BR, SI, OW = 0, 1, 63
ACTIVE = [0, 1, 63, 64, 65]  # bridge O, attacked Si, water O + H's

def _cli():
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--nice", type=int, default=10, help="process niceness increment")
    ap.add_argument("--log", help="tee output to this path (default: qm/runs)")
    return ap.parse_args()


ARGS = _cli()
preload_cutensor()
settings = DftSettings(xc="b3lyp", basis="def2-svp", density_fit=True, use_gpu=True)

deck = Path(__file__).resolve().parent.parent.parent / "petra/examples/kaolinite.toml"
cc = from_deck_cell(deck, "Oss", metal_shells=2)
template, _ = attack_complex(cc, water())

intermediate = load_xyz(RUN / "intermediate.xyz", template)
product = load_xyz(RUN / "product.xyz", template)

crest_path = RUN / "ts2_crest.xyz"
if crest_path.exists():
    log("resume: ts2_crest.xyz exists")
    crest = load_xyz(crest_path, template)
else:
    log("rebuilding TS2 crest: constrained opt at r(Si-Obr)=2.46, new bonds held")
    pts = constrained_scan(
        product, settings, atom_i=SI, atom_j=BR, distances_a=[2.46],
        fixed_distances=[(SI, OW, 1.70), (BR, 64, 0.98)],
    )
    crest = pts[0][2]
    save_xyz(crest, crest_path)

log(
    "directed Cartesian Sella from the crest "
    "(mode: intermediate -> product, reactive core)"
)
mode = reaction_path_vector(intermediate, product, active_indices=ACTIVE)
ts = find_ts(
    crest, settings,
    initial_mode=mode, internal=False,
    trajectory=str(RUN / "sella.ts2.traj"),
)
save_xyz(ts, RUN / "ts.xyz")
r_ow = float(np.linalg.norm(ts.coords[SI] - ts.coords[OW]))
r_br = float(np.linalg.norm(ts.coords[SI] - ts.coords[BR]))
log(f"TS2 saddle: r(Si-Ow)={r_ow:.2f} A  r(Si-Obr)={r_br:.2f} A  -> ts.xyz written")
log("resume phase2_ladder for verification + thermochemistry")
