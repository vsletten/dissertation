#!/usr/bin/env python
"""Run A1f's one donor/attacker-separated acid-relay campaign.

The finite production scope is exactly one four-water
``separated-donor-neutral-attacker`` Si path.  H3O+ donates through one relay
water while a distinct neutral water attacks Si; a matched Al path is legal
only after the Si path passes every endpoint, CI-NEB, saddle, coupled-mode, and
full-IRC gate.  Every stage is bounded and hash-journaled.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

_ETIQUETTE_SESSION = None
if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    _ETIQUETTE_SESSION = bootstrap_cli(
        "separated_acid_relay",
        default_run_root=Path(__file__).resolve().parent.parent / "runs",
    )

import numpy as np  # noqa: E402

from quarry.clusters import (  # noqa: E402
    SEPARATED_ACID_FAMILY,
    SEPARATED_ACID_WATER_COUNT,
    Cluster,
    SeparatedRelayEndpoints,
    separated_acid_relay_endpoints,
)
from quarry.pipeline import DftSettings, FrequencyResult  # noqa: E402
from scripts import concerted_acid_relay as shared  # noqa: E402
from scripts import phase1_xiao_lasaga as phase1  # noqa: E402

SCHEMA_VERSION = shared.SCHEMA_VERSION
MECHANISM_VERSION = 2
GATE_VERSION = 2
MECHANISM_NAME = "separated-hydronium-donor-neutral-water-attacker-relay"
MIN_PAIR_DISTANCE_A = shared.MIN_PAIR_DISTANCE_A
MODE_COMPONENT_MIN = shared.MODE_COMPONENT_MIN
TERMINAL_STATUSES = shared.TERMINAL_STATUSES
Bounds = shared.Bounds


def proton_indices(ends: SeparatedRelayEndpoints, n_water: int) -> tuple[int, ...]:
    if n_water != SEPARATED_ACID_WATER_COUNT:
        raise ValueError("separated acid relay requires exactly 4 waters")
    expected_count = 2 * n_water + 1
    if len(ends.solvent_h_indices) != expected_count:
        raise ValueError("separated endpoint solvent-H map drifted from A1e atom order")
    return ends.solvent_h_indices


def expected_basin(
    n_water: int, endpoint: str
) -> tuple[bool, bool, bool, int, int, int]:
    if n_water != SEPARATED_ACID_WATER_COUNT:
        raise ValueError("separated acid relay requires exactly 4 waters")
    if endpoint == "reactant":
        return (False, True, True, 0, 2 * n_water + 1, 0)
    if endpoint == "product":
        return (True, False, True, 1, 2 * n_water, 0)
    raise ValueError(f"unknown endpoint {endpoint}")


def role_occupancy_reason(
    cluster: Cluster,
    ends: SeparatedRelayEndpoints,
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
    expected_solvent = (3, 2, 2, 2) if endpoint == "reactant" else (2, 2, 2, 2)
    if actual_solvent != expected_solvent:
        return (
            f"separated solvent occupancies {actual_solvent} != {expected_solvent} "
            "(donor, neutral attacker, relay, spectator)"
        )
    expected_bridge = 0 if endpoint == "reactant" else 1
    actual_bridge = assignments.count(phase1.BR_INDEX)
    if actual_bridge != expected_bridge:
        return f"separated bridge occupancy {actual_bridge} != {expected_bridge}"
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
            f"separated framework occupancies {actual_framework} "
            f"!= {expected_framework}"
        )
    return None


def endpoint_gate_reason(
    cluster: Cluster,
    ends: SeparatedRelayEndpoints,
    *,
    n_water: int,
    endpoint: str,
) -> str | None:
    if ends.donor_oxygen_index == ends.ow_index:
        return "acid donor and neutral attacker oxygen must be distinct"
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
    ends: SeparatedRelayEndpoints,
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
        "donor_h_release": shared._distance_derivative(
            ts, mode, ends.donor_oxygen_index, ends.donor_h_index
        ),
        "donor_h_to_relay": shared._distance_derivative(
            ts, mode, ends.relay_oxygen_index, ends.donor_h_index
        ),
        "relay_h_release": shared._distance_derivative(
            ts, mode, ends.relay_oxygen_index, ends.transferred_h_index
        ),
        "relay_to_obr": shared._distance_derivative(
            ts, mode, 0, ends.transferred_h_index
        ),
    }
    accepted = (
        components["si_obr_cleavage"] >= MODE_COMPONENT_MIN
        and components["si_ow_attack"] <= -MODE_COMPONENT_MIN
        and components["donor_h_release"] >= MODE_COMPONENT_MIN
        and components["donor_h_to_relay"] <= -MODE_COMPONENT_MIN
        and components["relay_h_release"] >= MODE_COMPONENT_MIN
        and components["relay_to_obr"] <= -MODE_COMPONENT_MIN
    )
    reason = (
        None
        if accepted
        else "imaginary mode does not couple donor/relay transfer, attack, and cleavage"
    )
    return {"accepted": accepted, "reason": reason, "components": components}


def irc_channel_reason(
    backward: Cluster,
    forward: Cluster,
    ends: SeparatedRelayEndpoints,
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


def run_path(
    run_dir: Path,
    settings: DftSettings,
    bounds: Bounds,
    *,
    model: str,
    n_water: int = SEPARATED_ACID_WATER_COUNT,
    family: str = SEPARATED_ACID_FAMILY,
) -> dict[str, Any]:
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
        endpoint_builder=separated_acid_relay_endpoints,
        endpoint_gate=endpoint_gate_reason,
        mode_gate=coupled_mode_components,
        irc_gate=irc_channel_reason,
    )


def _expected_manifest(
    settings: DftSettings, bounds: Bounds, log_path: str
) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "mechanism_version": MECHANISM_VERSION,
        "gate_version": GATE_VERSION,
        "mechanism": MECHANISM_NAME,
        "settings": asdict(settings),
        "bounds": asdict(bounds),
        "water_count": SEPARATED_ACID_WATER_COUNT,
        "family": SEPARATED_ACID_FAMILY,
        "log": log_path,
        "paths": {},
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
) -> dict[str, Any]:
    bounds.validate()
    run_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = run_dir / "manifest.json"
    expected = _expected_manifest(settings, bounds, log_path)
    identity_fields = (
        "schema_version",
        "mechanism_version",
        "gate_version",
        "mechanism",
        "settings",
        "bounds",
        "water_count",
        "family",
    )
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text())
        for field in identity_fields:
            if manifest.get(field) != expected[field]:
                raise RuntimeError(f"existing A1f manifest has different {field}")
        manifest["log"] = log_path
    else:
        manifest = expected
        shared.atomic_json(manifest_path, manifest)

    si_record = run_path(run_dir, settings, bounds, model="si")
    _record_path(manifest, run_dir, si_record)
    shared.atomic_json(manifest_path, manifest)

    al_record = None
    matched = None
    if si_record["status"] == "accepted":
        al_record = run_path(run_dir, settings, bounds, model="al")
        _record_path(manifest, run_dir, al_record)
        shared.atomic_json(manifest_path, manifest)
        if al_record["status"] == "accepted":
            si_barrier = float(si_record["barrier_kj_mol"])
            al_barrier = float(al_record["barrier_kj_mol"])
            matched = {
                "water_count": SEPARATED_ACID_WATER_COUNT,
                "family": SEPARATED_ACID_FAMILY,
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
    print("A1F_SUMMARY " + json.dumps(manifest["summary"], sort_keys=True), flush=True)
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
        / f"acid-separated-relay-v{MECHANISM_VERSION}-g{GATE_VERSION}-{method}"
    )
    manifest = run_campaign(
        run_dir,
        settings,
        bounds,
        log_path=str(_ETIQUETTE_SESSION.log_path),
    )
    return 1 if manifest["summary"]["verdict"].startswith("incomplete-") else 0


if __name__ == "__main__":
    raise SystemExit(main())
