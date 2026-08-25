#!/usr/bin/env python
"""Run A1g's bridge-side-hydronium / neutral-attacker acid campaign.

The finite production scope is exactly one four-water
``bridge-side-hydronium-neutral-attacker`` Si path selected by A1f's exact
unconstrained reactant. The migrated H makes the bridge-side water H3O+; its
bridge-facing physical H then transfers directly to Obr while a distinct
neutral H2O attacks Si. A matched Al path is
legal only after Si passes every endpoint, CI-NEB, saddle, coupled-mode, and
full-IRC gate. Every stage is bounded and hash-journaled.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

_ETIQUETTE_SESSION = None
if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    _ETIQUETTE_SESSION = bootstrap_cli(
        "bridge_side_acid_relay",
        default_run_root=Path(__file__).resolve().parent.parent / "runs",
    )

import numpy as np  # noqa: E402

from quarry.clusters import (  # noqa: E402
    BRIDGE_SIDE_ACID_FAMILY,
    BRIDGE_SIDE_ACID_WATER_COUNT,
    BridgeSideHydroniumEndpoints,
    Cluster,
    bridge_side_hydronium_endpoints,
)
from quarry.pipeline import DftSettings, FrequencyResult  # noqa: E402
from scripts import concerted_acid_relay as shared  # noqa: E402
from scripts import phase1_xiao_lasaga as phase1  # noqa: E402

SCHEMA_VERSION = shared.SCHEMA_VERSION
MECHANISM_VERSION = 3
GATE_VERSION = 3
MECHANISM_NAME = "bridge-side-hydronium-neutral-water-attacker-direct-transfer"
MIN_PAIR_DISTANCE_A = shared.MIN_PAIR_DISTANCE_A
MODE_COMPONENT_MIN = shared.MODE_COMPONENT_MIN
TERMINAL_STATUSES = shared.TERMINAL_STATUSES
Bounds = shared.Bounds

A1F_SOURCE_REACTANT = Path(
    "/mnt/data/vsletten/dissertation-data/task201-a1f-acid-separated-20260824/"
    "acid-separated-relay-v2-g2-b3lyp-def2-svp/si/4w/"
    "separated-donor-neutral-attacker/reactant-optimized.xyz"
)
A1F_SOURCE_REACTANT_SHA256 = (
    "3a5107648bd44cbba15dd7482b4a9e8b3c6434f049442453744737e2c32aba75"
)


def proton_indices(ends: BridgeSideHydroniumEndpoints, n_water: int) -> tuple[int, ...]:
    if n_water != BRIDGE_SIDE_ACID_WATER_COUNT:
        raise ValueError("bridge-side acid relay requires exactly 4 waters")
    expected_count = 2 * n_water + 1
    if len(ends.solvent_h_indices) != expected_count:
        raise ValueError(
            "bridge-side endpoint solvent-H map drifted from A1f atom order"
        )
    return ends.solvent_h_indices


def expected_basin(
    n_water: int, endpoint: str
) -> tuple[bool, bool, bool, int, int, int]:
    if n_water != BRIDGE_SIDE_ACID_WATER_COUNT:
        raise ValueError("bridge-side acid relay requires exactly 4 waters")
    if endpoint == "reactant":
        return (False, True, True, 0, 2 * n_water + 1, 0)
    if endpoint == "product":
        return (True, False, True, 1, 2 * n_water, 0)
    raise ValueError(f"unknown endpoint {endpoint}")


def role_occupancy_reason(
    cluster: Cluster,
    ends: BridgeSideHydroniumEndpoints,
    *,
    endpoint: str,
) -> str | None:
    all_hydrogens = tuple(
        index for index, symbol in enumerate(cluster.symbols) if symbol == "H"
    )
    if reason := phase1.acid_hydrogen_ownership_reason(cluster, all_hydrogens):
        return reason
    assignments = phase1._all_hydrogen_assignments(cluster)
    actual_solvent = tuple(
        assignments.count(oxygen) for oxygen in ends.solvent_oxygen_indices
    )
    expected_solvent = (2, 2, 3, 2) if endpoint == "reactant" else (2, 2, 2, 2)
    if actual_solvent != expected_solvent:
        return (
            f"bridge-side solvent occupancies {actual_solvent} != {expected_solvent} "
            "(outer solvator, neutral attacker, bridge-side hydronium, spectator)"
        )
    expected_bridge = 0 if endpoint == "reactant" else 1
    actual_bridge = assignments.count(phase1.BR_INDEX)
    if actual_bridge != expected_bridge:
        return f"bridge-side bridge occupancy {actual_bridge} != {expected_bridge}"
    solvent_set = set(ends.solvent_oxygen_indices)
    framework_oxygens = [
        index
        for index, symbol in enumerate(cluster.symbols)
        if symbol == "O" and index not in {phase1.BR_INDEX, *solvent_set}
    ]
    actual_framework = tuple(
        sorted(assignments.count(oxygen) for oxygen in framework_oxygens)
    )
    expected_framework = tuple(1 for _ in framework_oxygens)
    if actual_framework != expected_framework:
        return (
            f"bridge-side framework occupancies {actual_framework} "
            f"!= {expected_framework}"
        )
    return None


def endpoint_gate_reason(
    cluster: Cluster,
    ends: BridgeSideHydroniumEndpoints,
    *,
    n_water: int,
    endpoint: str,
) -> str | None:
    if ends.hydronium_oxygen_index == ends.ow_index:
        return "bridge-side hydronium and neutral attacker oxygen must be distinct"
    if cluster.symbols != ends.reactant.symbols or not np.all(
        np.isfinite(cluster.coords)
    ):
        return f"{endpoint} atom mapping or coordinates are invalid"
    minimum_distance = shared.minimum_pair_distance(cluster)
    if minimum_distance < MIN_PAIR_DISTANCE_A:
        return (
            f"{endpoint} minimum pair distance {minimum_distance:.3f} A is below gate"
        )
    protons = proton_indices(ends, n_water)
    signature = phase1.acid_basin_signature(
        cluster,
        ends.ow_index,
        protons,
        solvent_oxygen_indices=ends.solvent_oxygen_indices,
    )
    expected = expected_basin(n_water, endpoint)
    if signature != expected:
        return f"{endpoint} basin {signature} != {expected}"
    return role_occupancy_reason(cluster, ends, endpoint=endpoint)


def coupled_mode_components(
    ts: Cluster,
    frequency: FrequencyResult,
    ends: BridgeSideHydroniumEndpoints,
) -> dict[str, Any]:
    if frequency.n_imaginary != 1 or frequency.imaginary_mode is None:
        return {"accepted": False, "reason": "expected exactly one imaginary mode"}
    mode = np.asarray(frequency.imaginary_mode, dtype=float)
    if mode.shape != ts.coords.shape or not np.all(np.isfinite(mode)):
        return {"accepted": False, "reason": "imaginary mode is invalid"}
    norm = float(np.linalg.norm(mode))
    if norm < 1.0e-12:
        return {"accepted": False, "reason": "imaginary mode is zero"}
    mode = mode / norm
    cleavage = shared._distance_derivative(ts, mode, 1, 0)
    if cleavage < 0.0:
        mode = -mode
        cleavage = -cleavage
    components = {
        "si_obr_cleavage": cleavage,
        "si_ow_attack": shared._distance_derivative(ts, mode, 1, ends.ow_index),
        "hydronium_h_release": shared._distance_derivative(
            ts, mode, ends.hydronium_oxygen_index, ends.transferred_h_index
        ),
        "hydronium_to_obr": shared._distance_derivative(
            ts, mode, phase1.BR_INDEX, ends.transferred_h_index
        ),
    }
    accepted = (
        components["si_obr_cleavage"] >= MODE_COMPONENT_MIN
        and components["si_ow_attack"] <= -MODE_COMPONENT_MIN
        and components["hydronium_h_release"] >= MODE_COMPONENT_MIN
        and components["hydronium_to_obr"] <= -MODE_COMPONENT_MIN
    )
    reason = (
        None
        if accepted
        else (
            "imaginary mode does not couple direct hydronium transfer, "
            "attack, and cleavage"
        )
    )
    return {"accepted": accepted, "reason": reason, "components": components}


def irc_channel_reason(
    backward: Cluster,
    forward: Cluster,
    ends: BridgeSideHydroniumEndpoints,
    *,
    n_water: int,
) -> str | None:
    matched: list[str] = []
    for endpoint in (backward, forward):
        ownership = phase1.acid_hydrogen_ownership_reason(
            endpoint,
            tuple(
                index for index, symbol in enumerate(endpoint.symbols) if symbol == "H"
            ),
        )
        if ownership:
            return f"IRC endpoint {ownership}"
        label = None
        fallback = None
        for candidate in ("reactant", "product"):
            reason = endpoint_gate_reason(
                endpoint,
                ends,
                n_water=n_water,
                endpoint=candidate,
            )
            if reason is None:
                label = candidate
                break
            if fallback is None or "basin" not in reason:
                fallback = reason
        if label is None:
            return fallback
        matched.append(label)
    if set(matched) != {"reactant", "product"}:
        return f"full IRC endpoint basins {sorted(matched)} != ['product', 'reactant']"
    return None


def load_source_reactant(
    path: Path,
    expected_sha256: str,
    ends: BridgeSideHydroniumEndpoints,
) -> BridgeSideHydroniumEndpoints:
    if not path.is_file():
        raise ValueError(f"A1f source reactant does not exist: {path}")
    actual = shared.sha256_path(path)
    if actual != expected_sha256:
        raise ValueError(
            f"A1f source reactant SHA-256 {actual} != expected {expected_sha256}"
        )
    source = shared.load_xyz_strict(path, ends.reactant)
    loaded = replace(ends, reactant=replace(ends.reactant, coords=source.coords))
    if reason := endpoint_gate_reason(
        loaded.reactant,
        loaded,
        n_water=BRIDGE_SIDE_ACID_WATER_COUNT,
        endpoint="reactant",
    ):
        raise ValueError(f"A1f source reactant fails A1g identity: {reason}")
    return loaded


def run_path(
    run_dir: Path,
    settings: DftSettings,
    bounds: Bounds,
    *,
    model: str,
    n_water: int = BRIDGE_SIDE_ACID_WATER_COUNT,
    family: str = BRIDGE_SIDE_ACID_FAMILY,
    source_reactant: Path | None = None,
    source_reactant_sha256: str | None = None,
) -> dict[str, Any]:
    if (source_reactant is None) != (source_reactant_sha256 is None):
        raise ValueError("source reactant path and SHA-256 must be supplied together")
    source_ends = None
    if source_reactant is not None and source_reactant_sha256 is not None:
        template = bridge_side_hydronium_endpoints(
            shared.SI_MODELS[model](), n_water=n_water, family=family
        )
        source_ends = load_source_reactant(
            source_reactant, source_reactant_sha256, template
        )

    def endpoint_builder(
        dimer: Cluster,
        *,
        n_water: int = n_water,
        family: str = family,
    ) -> BridgeSideHydroniumEndpoints:
        ends = bridge_side_hydronium_endpoints(dimer, n_water=n_water, family=family)
        if source_ends is None:
            return ends
        if ends.reactant.symbols != source_ends.reactant.symbols:
            raise ValueError(
                "A1f source reactant is not atom-matched to requested model"
            )
        return replace(
            ends,
            reactant=replace(ends.reactant, coords=source_ends.reactant.coords.copy()),
        )

    return shared.run_path(
        run_dir,
        settings,
        bounds,
        model=model,
        n_water=n_water,
        family=family,
        mechanism_version=MECHANISM_VERSION,
        gate_version=GATE_VERSION,
        mechanism_name=MECHANISM_NAME,
        endpoint_builder=endpoint_builder,
        endpoint_gate=endpoint_gate_reason,
        mode_gate=coupled_mode_components,
        irc_gate=irc_channel_reason,
        identity_extra={"source_reactant_sha256": source_reactant_sha256},
    )


def _expected_manifest(
    settings: DftSettings,
    bounds: Bounds,
    log_path: str,
    *,
    source_reactant: Path | None,
    source_reactant_sha256: str | None,
) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "mechanism_version": MECHANISM_VERSION,
        "gate_version": GATE_VERSION,
        "mechanism": MECHANISM_NAME,
        "settings": asdict(settings),
        "bounds": asdict(bounds),
        "water_count": BRIDGE_SIDE_ACID_WATER_COUNT,
        "family": BRIDGE_SIDE_ACID_FAMILY,
        "source_reactant": str(source_reactant) if source_reactant else None,
        "source_reactant_sha256": source_reactant_sha256,
        "paths": {},
        "log": log_path,
        "summary": None,
    }


def _record_path(
    manifest: dict[str, Any],
    run_dir: Path,
    record: dict[str, Any],
) -> None:
    path_dir = shared.path_directory(
        run_dir,
        record["model"],
        record["water_count"],
        record["family"],
    )
    terminal = path_dir / "terminal.json"
    manifest["paths"][record["key"]] = {
        "terminal": str(terminal.relative_to(run_dir)),
        "sha256": shared.sha256_path(terminal),
        "status": record["status"],
    }


def run_campaign(
    run_dir: Path,
    settings: DftSettings,
    bounds: Bounds,
    *,
    log_path: str,
    source_reactant: Path | None,
    source_reactant_sha256: str | None,
) -> dict[str, Any]:
    bounds.validate()
    if (source_reactant is None) != (source_reactant_sha256 is None):
        raise ValueError("source reactant path and SHA-256 must be supplied together")
    run_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = run_dir / "manifest.json"
    expected = _expected_manifest(
        settings,
        bounds,
        log_path,
        source_reactant=source_reactant,
        source_reactant_sha256=source_reactant_sha256,
    )
    identity_fields = (
        "schema_version",
        "mechanism_version",
        "gate_version",
        "mechanism",
        "settings",
        "bounds",
        "water_count",
        "family",
        "source_reactant_sha256",
    )
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text())
        for field in identity_fields:
            if manifest.get(field) != expected[field]:
                raise RuntimeError(f"existing A1g manifest has different {field}")
        manifest["source_reactant"] = expected["source_reactant"]
        manifest["log"] = log_path
    else:
        manifest = expected
        shared.atomic_json(manifest_path, manifest)

    si_record = run_path(
        run_dir,
        settings,
        bounds,
        model="si",
        source_reactant=source_reactant,
        source_reactant_sha256=source_reactant_sha256,
    )
    _record_path(manifest, run_dir, si_record)
    shared.atomic_json(manifest_path, manifest)

    al_record = None
    matched = None
    if si_record["status"] == "accepted":
        al_record = run_path(
            run_dir,
            settings,
            bounds,
            model="al",
            source_reactant=None,
            source_reactant_sha256=None,
        )
        _record_path(manifest, run_dir, al_record)
        shared.atomic_json(manifest_path, manifest)
        if al_record["status"] == "accepted":
            si_barrier = float(si_record["barrier_kj_mol"])
            al_barrier = float(al_record["barrier_kj_mol"])
            matched = {
                "water_count": BRIDGE_SIDE_ACID_WATER_COUNT,
                "family": BRIDGE_SIDE_ACID_FAMILY,
                "si_terminal": manifest["paths"][si_record["key"]],
                "al_terminal": manifest["paths"][al_record["key"]],
                "si_barrier_kj_mol": si_barrier,
                "al_barrier_kj_mol": al_barrier,
                "ordering": (
                    "Si-O-Al < Si-O-Si"
                    if al_barrier < si_barrier
                    else "Si-O-Si <= Si-O-Al"
                ),
            }

    if si_record["status"] in {"failed", "blocked"}:
        verdict = "incomplete-si-campaign"
    elif si_record["status"] == "rejected":
        verdict = "si-path-conclusive-rejection"
    elif al_record is None or al_record["status"] in {"failed", "blocked"}:
        verdict = "incomplete-matched-al-campaign"
    elif al_record["status"] == "rejected":
        verdict = "matched-al-conclusive-rejection"
    else:
        verdict = "matched-si-al-barrier-ready"
    manifest["summary"] = {
        "verdict": verdict,
        "si_status": si_record["status"],
        "al_status": al_record["status"] if al_record else None,
        "matched": matched,
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    shared.atomic_json(manifest_path, manifest)
    print("A1G_SUMMARY " + json.dumps(manifest["summary"], sort_keys=True), flush=True)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--nice", type=int, default=10)
    parser.add_argument("--log")
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--xc", default="b3lyp")
    parser.add_argument("--basis", default="def2-svp")
    parser.add_argument("--run-dir", type=Path)
    parser.add_argument("--source-reactant", type=Path, default=A1F_SOURCE_REACTANT)
    parser.add_argument("--source-reactant-sha256", default=A1F_SOURCE_REACTANT_SHA256)
    parser.add_argument("--endpoint-steps", type=int, default=160)
    parser.add_argument("--neb-images", type=int, default=7)
    parser.add_argument("--neb-pre-steps", type=int, default=100)
    parser.add_argument("--neb-climb-steps", type=int, default=200)
    parser.add_argument("--saddle-steps", type=int, default=300)
    parser.add_argument("--irc-steps", type=int, default=400)
    args = parser.parse_args()
    if _ETIQUETTE_SESSION is None:
        raise RuntimeError("campaign etiquette was not initialized")
    settings = DftSettings(
        xc=args.xc,
        basis=args.basis,
        density_fit=True,
        use_gpu=args.gpu,
    )
    bounds = Bounds(
        endpoint_steps=args.endpoint_steps,
        neb_images=args.neb_images,
        neb_pre_steps=args.neb_pre_steps,
        neb_climb_steps=args.neb_climb_steps,
        saddle_steps=args.saddle_steps,
        irc_steps=args.irc_steps,
    )
    method = f"{args.xc}-{args.basis}".lower().replace("/", "-")
    run_dir = args.run_dir or (
        Path(__file__).resolve().parent.parent
        / "runs"
        / "phase1"
        / f"acid-bridge-side-relay-v{MECHANISM_VERSION}-g{GATE_VERSION}-{method}"
    )
    manifest = run_campaign(
        run_dir,
        settings,
        bounds,
        log_path=str(_ETIQUETTE_SESSION.log_path),
        source_reactant=args.source_reactant,
        source_reactant_sha256=args.source_reactant_sha256,
    )
    return 1 if manifest["summary"]["verdict"].startswith("incomplete-") else 0


if __name__ == "__main__":
    raise SystemExit(main())
