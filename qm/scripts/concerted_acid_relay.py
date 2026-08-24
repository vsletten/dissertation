#!/usr/bin/env python
"""Run A1e's finite concerted hydronium/proton-relay campaign.

The production scope is exactly four Si paths: 3/4 waters crossed with
bridge-donor-chain/compact-cyclic-relay. A matched Al path is legal only after
its exact Si path passes endpoint, CI-NEB, saddle, coupled-mode, and full-IRC
gates. Every stage is bounded and hash-journaled.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from collections.abc import Callable
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

_ETIQUETTE_SESSION = None
if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    _ETIQUETTE_SESSION = bootstrap_cli(
        "concerted_acid_relay",
        default_run_root=Path(__file__).resolve().parent.parent / "runs",
    )

import numpy as np  # noqa: E402

from quarry.clusters import (  # noqa: E402
    CONCERTED_ACID_FAMILIES,
    CONCERTED_ACID_WATER_COUNTS,
    Cluster,
    ConcertedRelayEndpoints,
    aluminosilicate_dimer,
    concerted_acid_relay_endpoints,
    disilicate,
)
from quarry.pipeline import (  # noqa: E402
    DftSettings,
    FrequencyResult,
    frequencies,
    frequency_geometry_fingerprint,
    frequency_settings_fingerprint,
    optimize_bounded,
)
from quarry.ts import (  # noqa: E402
    find_ts,
    full_irc,
    neb_ts_guess,
    reaction_path_vector,
)
from scripts import phase1_xiao_lasaga as phase1  # noqa: E402

SCHEMA_VERSION = 1
MECHANISM_VERSION = 1
GATE_VERSION = 1
MIN_PAIR_DISTANCE_A = 0.75
MODE_COMPONENT_MIN = 0.02
TERMINAL_STATUSES = frozenset({"accepted", "rejected", "failed", "blocked"})
SI_MODELS = {"si": disilicate, "al": aluminosilicate_dimer}


@dataclass(frozen=True)
class Bounds:
    endpoint_steps: int = 160
    neb_images: int = 7
    neb_pre_steps: int = 100
    neb_climb_steps: int = 200
    saddle_steps: int = 300
    irc_steps: int = 400

    def validate(self) -> None:
        if not 1 <= self.endpoint_steps <= 240:
            raise ValueError("endpoint_steps must be in 1..240")
        if not 5 <= self.neb_images <= 11 or self.neb_images % 2 == 0:
            raise ValueError("neb_images must be odd and in 5..11")
        for name in ("neb_pre_steps", "neb_climb_steps", "saddle_steps", "irc_steps"):
            if not 1 <= getattr(self, name) <= 600:
                raise ValueError(f"{name} must be in 1..600")


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save_xyz(cluster: Cluster, path: Path) -> None:
    path.write_text(cluster.to_xyz())


def load_xyz_strict(path: Path, template: Cluster) -> Cluster:
    lines = path.read_text().splitlines()
    if len(lines) != len(template.symbols) + 2:
        raise ValueError(f"{path.name} atom count/line count drift")
    if int(lines[0]) != len(template.symbols):
        raise ValueError(f"{path.name} atom count drift")
    symbols: list[str] = []
    coords: list[list[float]] = []
    for line in lines[2:]:
        fields = line.split()
        symbols.append(fields[0])
        coords.append([float(value) for value in fields[1:4]])
    if symbols != template.symbols:
        raise ValueError(f"{path.name} atom order drift")
    array = np.asarray(coords, dtype=float)
    if array.shape != template.coords.shape or not np.all(np.isfinite(array)):
        raise ValueError(f"{path.name} coordinates are invalid")
    return replace(template, coords=array)


def save_frequency(path: Path, result: FrequencyResult) -> None:
    temporary = path.with_name(f".{path.name}.{time.time_ns()}.tmp.npz")
    np.savez(
        temporary,
        frequencies_cm=result.frequencies_cm,
        imaginary_cm=result.imaginary_cm,
        electronic_hartree=np.array(result.electronic_hartree),
        molar_mass_kg=np.array(result.molar_mass_kg),
        rotational_temperatures_k=np.array(result.rotational_temperatures_k or ()),
        linear=np.array(result.linear),
        imaginary_mode=(
            result.imaginary_mode
            if result.imaginary_mode is not None
            else np.empty((0, 3))
        ),
        imaginary_modes=(
            result.imaginary_modes
            if result.imaginary_modes is not None
            else np.empty((0, 0, 3))
        ),
        geometry_fingerprint=np.array(result.geometry_fingerprint),
        settings_fingerprint=np.array(result.settings_fingerprint),
    )
    temporary.replace(path)


def cached_frequency(
    path: Path,
    cluster: Cluster,
    settings: DftSettings,
    *,
    frequency_fn: Callable[[Cluster, DftSettings], FrequencyResult] = frequencies,
) -> FrequencyResult:
    if path.is_file():
        with np.load(path, allow_pickle=False) as data:
            geometry = str(data["geometry_fingerprint"].item())
            method = str(data["settings_fingerprint"].item())
            if geometry == frequency_geometry_fingerprint(cluster) and method == (
                frequency_settings_fingerprint(settings)
            ):
                rotational = tuple(
                    float(value) for value in data["rotational_temperatures_k"]
                )
                imaginary_mode = np.asarray(data["imaginary_mode"])
                imaginary_modes = np.asarray(data["imaginary_modes"])
                return FrequencyResult(
                    frequencies_cm=np.asarray(data["frequencies_cm"]),
                    imaginary_cm=np.asarray(data["imaginary_cm"]),
                    electronic_hartree=float(data["electronic_hartree"]),
                    molar_mass_kg=float(data["molar_mass_kg"]),
                    rotational_temperatures_k=rotational or None,
                    linear=bool(data["linear"]),
                    imaginary_mode=imaginary_mode if imaginary_mode.size else None,
                    imaginary_modes=imaginary_modes if imaginary_modes.size else None,
                    geometry_fingerprint=geometry,
                    settings_fingerprint=method,
                )
        path.replace(path.with_name(f"{path.stem}.stale-{time.time_ns()}{path.suffix}"))
    result = frequency_fn(cluster, settings)
    save_frequency(path, result)
    return result


def minimum_pair_distance(cluster: Cluster) -> float:
    delta = cluster.coords[:, None, :] - cluster.coords[None, :, :]
    distances = np.linalg.norm(delta, axis=2)
    distances[np.eye(len(cluster.symbols), dtype=bool)] = np.inf
    return float(np.min(distances))


def proton_indices(ends: ConcertedRelayEndpoints, n_water: int) -> tuple[int, ...]:
    protons, solvents = phase1.acid_mobile_indices(ends.ow_index, n_water)
    if solvents != ends.solvent_oxygen_indices:
        raise ValueError("concerted endpoint solvent map drifted from A1c atom order")
    return protons


def expected_basin(
    n_water: int, endpoint: str
) -> tuple[bool, bool, bool, int, int, int]:
    if endpoint == "reactant":
        return (False, True, True, 0, 2 * n_water + 1, 0)
    if endpoint == "product":
        return (True, False, True, 1, 2 * n_water, 0)
    raise ValueError(f"unknown endpoint {endpoint}")


def endpoint_gate_reason(
    cluster: Cluster,
    ends: ConcertedRelayEndpoints,
    *,
    n_water: int,
    endpoint: str,
) -> str | None:
    if cluster.symbols != ends.reactant.symbols or not np.all(
        np.isfinite(cluster.coords)
    ):
        return f"{endpoint} atom mapping or coordinates are invalid"
    minimum_distance = minimum_pair_distance(cluster)
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
    return phase1.acid_occupancy_reason(
        cluster,
        ends.ow_index,
        solvent_oxygen_indices=ends.solvent_oxygen_indices,
        attacker_h_count=3 if endpoint == "reactant" else 2,
        bridge_h_count=0 if endpoint == "reactant" else 1,
        extra_framework_h_count=0,
    )


def _distance_derivative(
    cluster: Cluster, mode: np.ndarray, atom_i: int, atom_j: int
) -> float:
    vector = cluster.coords[atom_i] - cluster.coords[atom_j]
    distance = float(np.linalg.norm(vector))
    if distance < 1.0e-12:
        raise ValueError("coupled-mode distance has coincident atoms")
    return float(np.dot(vector / distance, mode[atom_i] - mode[atom_j]))


def coupled_mode_components(
    ts: Cluster,
    frequency: FrequencyResult,
    ends: ConcertedRelayEndpoints,
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
    cleavage = _distance_derivative(ts, mode, 1, 0)
    if cleavage < 0.0:
        mode = -mode
        cleavage = -cleavage
    attack = _distance_derivative(ts, mode, 1, ends.ow_index)
    proton = _distance_derivative(ts, mode, 0, ends.transferred_h_index)
    components = {
        "si_obr_cleavage": cleavage,
        "si_ow_attack": attack,
        "relay_to_obr": proton,
    }
    accepted = (
        cleavage >= MODE_COMPONENT_MIN
        and attack <= -MODE_COMPONENT_MIN
        and proton <= -MODE_COMPONENT_MIN
    )
    reason = (
        None
        if accepted
        else "imaginary mode does not couple relay, attack, and cleavage"
    )
    return {"accepted": accepted, "reason": reason, "components": components}


def irc_channel_reason(
    backward: Cluster,
    forward: Cluster,
    ends: ConcertedRelayEndpoints,
    *,
    n_water: int,
) -> str | None:
    actual: list[tuple[bool, bool, bool, int, int, int]] = []
    protons = proton_indices(ends, n_water)
    for endpoint in (backward, forward):
        ownership = phase1.acid_hydrogen_ownership_reason(
            endpoint,
            tuple(
                index for index, symbol in enumerate(endpoint.symbols) if symbol == "H"
            ),
        )
        if ownership:
            return f"IRC endpoint {ownership}"
        actual.append(
            phase1.acid_basin_signature(
                endpoint,
                ends.ow_index,
                protons,
                solvent_oxygen_indices=ends.solvent_oxygen_indices,
            )
        )
    expected = {expected_basin(n_water, "reactant"), expected_basin(n_water, "product")}
    if set(actual) != expected:
        return f"full IRC endpoint basins {sorted(actual)} != {sorted(expected)}"
    return None


def path_key(model: str, n_water: int, family: str) -> str:
    return f"{model}:{n_water}w:{family}"


def path_directory(run_dir: Path, model: str, n_water: int, family: str) -> Path:
    return run_dir / model / f"{n_water}w" / family


def _artifact_hashes(path_dir: Path) -> dict[str, str]:
    artifacts: dict[str, str] = {}
    for path in sorted(path_dir.rglob("*")):
        if path.is_file() and path.name != "terminal.json":
            artifacts[str(path.relative_to(path_dir))] = sha256_path(path)
    return artifacts


def terminal_record_is_reusable(
    path_dir: Path,
    *,
    key: str,
    settings: DftSettings,
    bounds: Bounds,
) -> bool:
    terminal = path_dir / "terminal.json"
    if not terminal.is_file():
        return False
    try:
        record = json.loads(terminal.read_text())
    except (OSError, json.JSONDecodeError, UnicodeDecodeError):
        return False
    expected = {
        "schema_version": SCHEMA_VERSION,
        "mechanism_version": MECHANISM_VERSION,
        "gate_version": GATE_VERSION,
        "key": key,
        "settings": asdict(settings),
        "bounds": asdict(bounds),
    }
    if any(record.get(field) != value for field, value in expected.items()):
        return False
    if record.get("status") not in TERMINAL_STATUSES:
        return False
    artifacts = record.get("artifacts")
    if not isinstance(artifacts, dict) or not artifacts:
        return False
    for relative, digest in artifacts.items():
        artifact = (path_dir / relative).resolve()
        if not artifact.is_relative_to(path_dir.resolve()) or not artifact.is_file():
            return False
        if sha256_path(artifact) != digest:
            return False
    return True


def _terminal(
    path_dir: Path,
    *,
    base: dict[str, Any],
    status: str,
    stage: str,
    reason: str | None,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    record = {
        **base,
        "status": status,
        "stage": stage,
        "reason": reason,
        "artifacts": _artifact_hashes(path_dir),
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    if extra:
        record.update(extra)
    atomic_json(path_dir / "terminal.json", record)
    return record


def run_path(
    run_dir: Path,
    settings: DftSettings,
    bounds: Bounds,
    *,
    model: str,
    n_water: int,
    family: str,
) -> dict[str, Any]:
    key = path_key(model, n_water, family)
    path_dir = path_directory(run_dir, model, n_water, family)
    if terminal_record_is_reusable(path_dir, key=key, settings=settings, bounds=bounds):
        return json.loads((path_dir / "terminal.json").read_text())
    if path_dir.exists():
        path_dir.rename(path_dir.with_name(f"{path_dir.name}.stale-{time.time_ns()}"))
    path_dir.mkdir(parents=True)
    base = {
        "schema_version": SCHEMA_VERSION,
        "mechanism_version": MECHANISM_VERSION,
        "gate_version": GATE_VERSION,
        "key": key,
        "model": model,
        "water_count": n_water,
        "family": family,
        "settings": asdict(settings),
        "bounds": asdict(bounds),
    }
    stage = "build-endpoints"
    try:
        ends = concerted_acid_relay_endpoints(
            SI_MODELS[model](), n_water=n_water, family=family
        )
        save_xyz(ends.reactant, path_dir / "reactant-seed.xyz")
        save_xyz(ends.product, path_dir / "product-seed.xyz")

        stage = "optimize-reactant"
        reactant_result = optimize_bounded(
            ends.reactant, settings, max_steps=bounds.endpoint_steps
        )
        save_xyz(reactant_result.cluster, path_dir / "reactant-optimized.xyz")
        if not reactant_result.converged:
            raise RuntimeError("reactant geometry optimization exhausted its bound")
        if reason := endpoint_gate_reason(
            reactant_result.cluster, ends, n_water=n_water, endpoint="reactant"
        ):
            return _terminal(
                path_dir, base=base, status="rejected", stage=stage, reason=reason
            )

        stage = "reactant-hessian"
        reactant_frequency = cached_frequency(
            path_dir / "reactant-hessian.npz", reactant_result.cluster, settings
        )
        if reactant_frequency.n_imaginary:
            return _terminal(
                path_dir,
                base=base,
                status="rejected",
                stage=stage,
                reason=(
                    "reactant is not a minimum: imaginary modes "
                    f"{reactant_frequency.imaginary_cm.tolist()}"
                ),
            )

        stage = "optimize-product"
        product_result = optimize_bounded(
            ends.product, settings, max_steps=bounds.endpoint_steps
        )
        save_xyz(product_result.cluster, path_dir / "product-optimized.xyz")
        if not product_result.converged:
            raise RuntimeError("product geometry optimization exhausted its bound")
        if reason := endpoint_gate_reason(
            product_result.cluster, ends, n_water=n_water, endpoint="product"
        ):
            return _terminal(
                path_dir, base=base, status="rejected", stage=stage, reason=reason
            )

        stage = "product-hessian"
        product_frequency = cached_frequency(
            path_dir / "product-hessian.npz", product_result.cluster, settings
        )
        if product_frequency.n_imaginary:
            return _terminal(
                path_dir,
                base=base,
                status="rejected",
                stage=stage,
                reason=(
                    "product is not a minimum: imaginary modes "
                    f"{product_frequency.imaginary_cm.tolist()}"
                ),
            )

        stage = "ci-neb"
        neb_root = path_dir / "neb-checkpoints"
        crest = neb_ts_guess(
            reactant_result.cluster,
            product_result.cluster,
            settings,
            n_images=bounds.neb_images,
            pre_relax_steps=bounds.neb_pre_steps,
            max_steps=bounds.neb_climb_steps,
            checkpoint_dir=neb_root,
            checkpoint_interval=5,
            climb_optimizer="ode",
        )
        save_xyz(crest, path_dir / "climbing-band-crest.xyz")

        stage = "saddle"
        active_indices = sorted(
            {
                0,
                1,
                ends.ow_index,
                ends.transferred_h_index,
                *ends.relay_h_indices,
                *ends.solvent_oxygen_indices,
            }
        )
        mode = reaction_path_vector(
            reactant_result.cluster,
            product_result.cluster,
            active_indices=active_indices,
        )
        saddle = find_ts(
            crest,
            settings,
            max_steps=bounds.saddle_steps,
            trajectory=str(path_dir / "saddle.traj"),
            initial_mode=mode,
            internal=False,
        )
        save_xyz(saddle, path_dir / "saddle.xyz")

        stage = "saddle-hessian"
        saddle_frequency = cached_frequency(
            path_dir / "saddle-hessian.npz", saddle, settings
        )
        mode_receipt = coupled_mode_components(saddle, saddle_frequency, ends)
        atomic_json(path_dir / "coupled-mode.json", mode_receipt)
        if not mode_receipt["accepted"]:
            return _terminal(
                path_dir,
                base=base,
                status="rejected",
                stage=stage,
                reason=str(mode_receipt["reason"]),
            )

        stage = "full-irc"
        backward, forward = full_irc(
            saddle,
            settings,
            max_steps=bounds.irc_steps,
            trajectory=path_dir / "full-irc.traj",
            logfile=path_dir / "full-irc.log",
        )
        save_xyz(backward, path_dir / "irc-backward.xyz")
        save_xyz(forward, path_dir / "irc-forward.xyz")
        if reason := irc_channel_reason(backward, forward, ends, n_water=n_water):
            return _terminal(
                path_dir, base=base, status="rejected", stage=stage, reason=reason
            )

        stage = "barrier"
        reactant_thermo = phase1.thermo_result(reactant_frequency, 298.15)
        saddle_thermo = phase1.thermo_result(saddle_frequency, 298.15)
        barrier_kj = float(saddle_thermo.gibbs - reactant_thermo.gibbs)
        if not math.isfinite(barrier_kj) or barrier_kj <= 0.0:
            return _terminal(
                path_dir,
                base=base,
                status="rejected",
                stage=stage,
                reason=f"non-positive/non-finite barrier {barrier_kj}",
            )
        return _terminal(
            path_dir,
            base=base,
            status="accepted",
            stage="completed",
            reason=None,
            extra={
                "barrier_kj_mol": barrier_kj,
                "barrier_kcal_mol": barrier_kj / 4.184,
                "imaginary_cm": float(saddle_frequency.imaginary_cm[0]),
                "coupled_mode": mode_receipt,
            },
        )
    except Exception as exc:
        return _terminal(
            path_dir,
            base=base,
            status="failed",
            stage=f"incomplete-{stage}",
            reason=f"{type(exc).__name__}: {exc}",
        )


def _expected_manifest(
    settings: DftSettings, bounds: Bounds, log_path: str
) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "mechanism_version": MECHANISM_VERSION,
        "gate_version": GATE_VERSION,
        "settings": asdict(settings),
        "bounds": asdict(bounds),
        "water_counts": list(CONCERTED_ACID_WATER_COUNTS),
        "families": list(CONCERTED_ACID_FAMILIES),
        "log": log_path,
        "paths": {},
        "summary": None,
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
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text())
        for field in (
            "schema_version",
            "mechanism_version",
            "gate_version",
            "settings",
            "bounds",
            "water_counts",
            "families",
        ):
            if manifest.get(field) != expected[field]:
                raise RuntimeError(f"existing A1e manifest has different {field}")
        manifest["log"] = log_path
    else:
        manifest = expected
        atomic_json(manifest_path, manifest)

    si_records: list[dict[str, Any]] = []
    matched: list[dict[str, Any]] = []
    for n_water in CONCERTED_ACID_WATER_COUNTS:
        for family in CONCERTED_ACID_FAMILIES:
            si_record = run_path(
                run_dir,
                settings,
                bounds,
                model="si",
                n_water=n_water,
                family=family,
            )
            manifest["paths"][si_record["key"]] = {
                "terminal": str(
                    (
                        path_directory(run_dir, "si", n_water, family) / "terminal.json"
                    ).relative_to(run_dir)
                ),
                "sha256": sha256_path(
                    path_directory(run_dir, "si", n_water, family) / "terminal.json"
                ),
                "status": si_record["status"],
            }
            atomic_json(manifest_path, manifest)
            si_records.append(si_record)
            if si_record["status"] != "accepted":
                continue
            al_record = run_path(
                run_dir,
                settings,
                bounds,
                model="al",
                n_water=n_water,
                family=family,
            )
            manifest["paths"][al_record["key"]] = {
                "terminal": str(
                    (
                        path_directory(run_dir, "al", n_water, family) / "terminal.json"
                    ).relative_to(run_dir)
                ),
                "sha256": sha256_path(
                    path_directory(run_dir, "al", n_water, family) / "terminal.json"
                ),
                "status": al_record["status"],
            }
            atomic_json(manifest_path, manifest)
            if al_record["status"] == "accepted":
                si_barrier = float(si_record["barrier_kj_mol"])
                al_barrier = float(al_record["barrier_kj_mol"])
                matched.append(
                    {
                        "water_count": n_water,
                        "family": family,
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
                )

    si_counts = {
        status: sum(record["status"] == status for record in si_records)
        for status in sorted(TERMINAL_STATUSES)
    }
    if matched:
        verdict = "matched-si-al-barrier-ready"
    elif si_counts["rejected"] == 4:
        verdict = "all-four-si-paths-conclusive-rejection"
    elif si_counts["failed"] or si_counts["blocked"]:
        verdict = "incomplete-si-campaign"
    else:
        verdict = "si-path-accepted-without-matched-al-barrier"
    manifest["summary"] = {
        "verdict": verdict,
        "si_path_count": len(si_records),
        "si_status_counts": si_counts,
        "matched": matched,
        "next_mechanism_card_required": verdict
        == "all-four-si-paths-conclusive-rejection",
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    atomic_json(manifest_path, manifest)
    print("A1E_SUMMARY " + json.dumps(manifest["summary"], sort_keys=True), flush=True)
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
        / f"acid-concerted-relay-v{MECHANISM_VERSION}-g{GATE_VERSION}-{method}"
    )
    manifest = run_campaign(
        run_dir,
        settings,
        bounds,
        log_path=str(_ETIQUETTE_SESSION.log_path),
    )
    verdict = manifest["summary"]["verdict"]
    return 1 if verdict.startswith("incomplete-") else 0


if __name__ == "__main__":
    raise SystemExit(main())
