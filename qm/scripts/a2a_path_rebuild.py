#!/usr/bin/env python3
"""A2a: rebuild the si-neutral production path from exact r2SCAN-3c minima.

The rejected banked saddle is deliberately not an input.  The driver treats the
observed associative basin as an explicit intermediate and reconstructs two
independent elementary paths (R<->I and I<->P) with persisted pre-climb/climb
CI-NEB, reaction-aligned Cartesian Sella refinement, same-surface finite-
difference Hessians at two predeclared steps, and non-zero full Gonzalez--
Schlegel IRCs into exact typed neighboring basins.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    _ETIQUETTE = bootstrap_cli(
        "a2a-path-rebuild",
        default_run_root="/mnt/data/vsletten/dissertation-data/task208-a2a-path-rebuild",
    )

import numpy as np

from quarry.clusters import Cluster
from quarry.pipeline import (
    HARTREE_TO_KJ,
    DftSettings,
    FrequencyResult,
    frequencies_finite_difference,
    frequency_geometry_fingerprint,
    frequency_settings_fingerprint,
)
from quarry.rates import rate_from_thermo
from quarry.store import Store, geometry_hash
from quarry.ts import (
    find_ts,
    full_irc,
    neb_ts_guess,
    relax_at_fixed_distances,
)
from scripts import production_energetics as a2

PATH_GATE_VERSION = "a2a-sequential-path-v1"
SELLA_MODE_STRATEGY = "conditioned-crest-local-ci-neb-tangent-v1"
REACTION_COORDINATE_PAIRS = ((1, 15), (0, 17), (0, 1))
HESSIAN_STEPS_BOHR = (1.0e-3, 2.0e-3)
SIGNIFICANT_IMAGINARY_FLOOR_CM = 30.0
MODE_COSINE_MINIMUM = 0.80
REACTION_VECTOR_OVERLAP_MINIMUM = 0.20
FREQUENCY_RELATIVE_DRIFT_LIMIT = 0.15
FREQUENCY_ABSOLUTE_DRIFT_LIMIT_CM = 10.0
IRC_MINIMUM_CARTESIAN_DISPLACEMENT_A = 0.05

REACTANT_BASIN = (False, True, False)
ASSOCIATIVE_BASIN = (True, True, True)
PRODUCT_BASIN = (True, False, True)


@dataclass(frozen=True)
class SegmentSpec:
    slug: str
    start_role: str
    end_role: str


SEGMENTS = (
    SegmentSpec("addition", "reactant", "intermediate"),
    SegmentSpec("cleavage", "intermediate", "product"),
)


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    a2.atomic_json(path, payload)


def load_reference(
    path: Path, *, role: str, template: Cluster | None = None
) -> Cluster:
    if template is None:
        result = a2.parse_xyz(path.read_text(), name=role, charge=0, spin=0)
    else:
        result = a2.load_xyz_like(path, template, name=role)
    if result.charge != 0 or result.spin != 0:
        raise ValueError(f"{role}: si-neutral reference must be a closed-shell neutral")
    return result


def typed_identity_payload(cluster: Cluster, attacker_index: int) -> dict[str, Any]:
    basin, hydrogen_owners, heavy_bonds = a2.endpoint_identity(cluster, attacker_index)
    return {
        "basin": list(basin),
        "hydrogen_owners": [list(pair) for pair in hydrogen_owners],
        "heavy_bonds": [list(pair) for pair in heavy_bonds],
    }


def require_basin(
    cluster: Cluster, attacker_index: int, expected: tuple[bool, bool, bool]
) -> None:
    actual = a2.si_neutral_signature(cluster, attacker_index)
    # endpoint_identity also enforces finite coordinates, the minimum-pair gate,
    # complete physical-H ownership, and an unambiguous heavy-atom graph.
    a2.endpoint_identity(cluster, attacker_index)
    if actual != expected:
        raise RuntimeError(f"{cluster.name}: expected basin {expected}, found {actual}")


def frequency_payload(result: FrequencyResult, *, step_bohr: float) -> dict[str, Any]:
    payload = a2.frequency_payload(result)
    payload["hessian_method"] = "finite-difference-gradient"
    payload["hessian_step_bohr"] = step_bohr
    payload["imaginary_mode"] = (
        result.imaginary_mode.tolist() if result.imaginary_mode is not None else None
    )
    payload["imaginary_modes"] = (
        result.imaginary_modes.tolist() if result.imaginary_modes is not None else None
    )
    return payload


def frequency_from_payload(payload: dict[str, Any]) -> FrequencyResult:
    base = a2.frequency_from_payload(payload)
    imaginary_mode = payload.get("imaginary_mode")
    imaginary_modes = payload.get("imaginary_modes")
    return replace(
        base,
        imaginary_mode=(
            np.asarray(imaginary_mode, dtype=float)
            if imaginary_mode is not None
            else None
        ),
        imaginary_modes=(
            np.asarray(imaginary_modes, dtype=float)
            if imaginary_modes is not None
            else None
        ),
    )


def checkpoint_fd_frequency(
    path: Path,
    cluster: Cluster,
    settings: DftSettings,
    *,
    step_bohr: float,
) -> FrequencyResult:
    expected_geometry = frequency_geometry_fingerprint(cluster)
    expected_settings = frequency_settings_fingerprint(settings)
    if path.exists():
        payload = json.loads(path.read_text())
        expected = (
            expected_geometry,
            expected_settings,
            "finite-difference-gradient",
            step_bohr,
        )
        actual = (
            payload.get("geometry_fingerprint"),
            payload.get("settings_fingerprint"),
            payload.get("hessian_method"),
            payload.get("hessian_step_bohr"),
        )
        if actual != expected:
            raise ValueError(f"{path.name}: cached finite-difference identity drift")
        return frequency_from_payload(payload)
    result = frequencies_finite_difference(cluster, settings, step_bohr=step_bohr)
    atomic_json(path, frequency_payload(result, step_bohr=step_bohr))
    return result


def _one_significant_mode(result: FrequencyResult) -> tuple[float, np.ndarray]:
    indices = np.flatnonzero(result.imaginary_cm > SIGNIFICANT_IMAGINARY_FLOOR_CM)
    if indices.size != 1:
        raise RuntimeError(
            "expected exactly one significant imaginary mode above "
            f"{SIGNIFICANT_IMAGINARY_FLOOR_CM:.0f} cm^-1, found "
            f"{result.imaginary_cm.tolist()}"
        )
    if result.imaginary_modes is None:
        raise RuntimeError("finite-difference receipt lacks imaginary mode vectors")
    index = int(indices[0])
    if index >= len(result.imaginary_modes):
        raise RuntimeError("imaginary frequency/mode cardinality drift")
    mode = np.asarray(result.imaginary_modes[index], dtype=float)
    norm = float(np.linalg.norm(mode))
    if mode.ndim != 2 or mode.shape[1] != 3 or not np.isfinite(norm) or norm < 1e-12:
        raise RuntimeError("significant imaginary mode is non-finite or zero")
    return float(result.imaginary_cm[index]), mode / norm


def require_hessian_step_stability(
    primary: FrequencyResult,
    secondary: FrequencyResult,
    reaction_vector: np.ndarray,
) -> dict[str, float]:
    primary_cm, primary_mode = _one_significant_mode(primary)
    secondary_cm, secondary_mode = _one_significant_mode(secondary)
    vector = np.asarray(reaction_vector, dtype=float)
    vector_norm = float(np.linalg.norm(vector))
    if (
        vector.shape != primary_mode.shape
        or not np.isfinite(vector_norm)
        or vector_norm < 1e-12
    ):
        raise ValueError("reaction vector is incompatible with the Hessian modes")
    vector /= vector_norm
    mode_cosine = abs(
        float(np.vdot(primary_mode.reshape(-1), secondary_mode.reshape(-1)))
    )
    primary_overlap = abs(float(np.vdot(primary_mode.reshape(-1), vector.reshape(-1))))
    secondary_overlap = abs(
        float(np.vdot(secondary_mode.reshape(-1), vector.reshape(-1)))
    )
    drift_cm = abs(primary_cm - secondary_cm)
    drift_limit_cm = max(
        FREQUENCY_ABSOLUTE_DRIFT_LIMIT_CM,
        FREQUENCY_RELATIVE_DRIFT_LIMIT * max(primary_cm, secondary_cm),
    )
    if drift_cm > drift_limit_cm:
        raise RuntimeError(
            f"imaginary frequency changed by {drift_cm:.3f} cm^-1 across the "
            f"bounded Hessian-step check (limit {drift_limit_cm:.3f})"
        )
    if mode_cosine < MODE_COSINE_MINIMUM:
        raise RuntimeError(
            f"imaginary-mode cosine {mode_cosine:.3f} is below "
            f"{MODE_COSINE_MINIMUM:.2f}"
        )
    if min(primary_overlap, secondary_overlap) < REACTION_VECTOR_OVERLAP_MINIMUM:
        raise RuntimeError(
            "imaginary mode is not reaction-aligned at both Hessian steps: "
            f"{primary_overlap:.3f}, {secondary_overlap:.3f}"
        )
    return {
        "primary_imaginary_cm": primary_cm,
        "secondary_imaginary_cm": secondary_cm,
        "frequency_drift_cm": drift_cm,
        "frequency_drift_limit_cm": drift_limit_cm,
        "mode_cosine": mode_cosine,
        "primary_reaction_overlap": primary_overlap,
        "secondary_reaction_overlap": secondary_overlap,
    }


def newest_pre_relaxed_checkpoint(root: Path) -> Path | None:
    candidates: list[Path] = []
    if not root.is_dir():
        return None
    for path in root.iterdir():
        manifest_path = path / "manifest.json"
        if not manifest_path.is_file():
            continue
        manifest = json.loads(manifest_path.read_text())
        if (
            manifest.get("stage") == "pre-relax-final"
            and manifest.get("converged") is True
        ):
            candidates.append(path)
    return max(candidates, key=lambda value: value.stat().st_mtime_ns, default=None)


def final_neb_tangent(
    checkpoint_root: Path,
    crest: Cluster,
    template: Cluster,
    active_indices: list[int],
) -> tuple[np.ndarray, dict[str, Any]]:
    """Recover the local converged CI-NEB tangent around the exact peak image."""
    latest = json.loads((checkpoint_root / "latest.json").read_text())
    checkpoint = checkpoint_root / str(latest["checkpoint"])
    manifest = json.loads((checkpoint / "manifest.json").read_text())
    if manifest.get("stage") != "climb-final" or manifest.get("converged") is not True:
        raise ValueError("Sella refinement requires a converged climb-final NEB band")
    image_paths = [checkpoint / str(name) for name in manifest.get("images", [])]
    if len(image_paths) < 3:
        raise ValueError("climb-final checkpoint does not contain a complete NEB band")
    images = [
        a2.load_xyz_like(path, template, name=f"neb-image-{index}")
        for index, path in enumerate(image_paths)
    ]
    distances = [
        float(np.sqrt(np.mean((image.coords - crest.coords) ** 2))) for image in images
    ]
    peak = int(np.argmin(distances))
    if distances[peak] > 1.0e-8:
        raise ValueError("saved NEB crest does not match any climb-final image")
    if peak in (0, len(images) - 1):
        raise ValueError("saved NEB crest is not an interior image")
    tangent = images[peak + 1].coords - images[peak - 1].coords
    masked = np.zeros_like(tangent)
    masked[np.asarray(active_indices, dtype=int)] = tangent[
        np.asarray(active_indices, dtype=int)
    ]
    if template.frozen_indices:
        masked[np.asarray(template.frozen_indices, dtype=int)] = 0.0
    norm = float(np.linalg.norm(masked))
    if not np.isfinite(norm) or norm < 1.0e-12:
        raise ValueError("local CI-NEB tangent is non-finite or zero")
    return masked / norm, {
        "strategy": SELLA_MODE_STRATEGY,
        "checkpoint": str(checkpoint),
        "peak_index": peak,
        "peak_rms_to_saved_crest_a": distances[peak],
        "left_image": image_paths[peak - 1].name,
        "right_image": image_paths[peak + 1].name,
        "active_indices": active_indices,
    }


def reaction_coordinate_targets(cluster: Cluster) -> list[tuple[int, int, float]]:
    return [
        (
            atom_i,
            atom_j,
            float(np.linalg.norm(cluster.coords[atom_i] - cluster.coords[atom_j])),
        )
        for atom_i, atom_j in REACTION_COORDINATE_PAIRS
    ]


def endpoint_displacement(ts: Cluster, endpoint: Cluster) -> float:
    if ts.symbols != endpoint.symbols:
        raise ValueError("IRC endpoint atom order differs from the transition state")
    return float(np.linalg.norm(endpoint.coords - ts.coords))


def run_segment(
    run_dir: Path,
    spec: SegmentSpec,
    start: Cluster,
    end: Cluster,
    settings: DftSettings,
    *,
    attacker_index: int,
    active_indices: list[int],
    neb_images: int,
    neb_steps: int,
    neb_pre_steps: int,
    saddle_steps: int,
    irc_steps: int,
) -> dict[str, Any]:
    segment_dir = run_dir / spec.slug
    segment_dir.mkdir(parents=True, exist_ok=True)
    settings_identity = frequency_settings_fingerprint(settings)
    neb_root = segment_dir / "neb-checkpoints"
    resume_from = newest_pre_relaxed_checkpoint(neb_root)

    def compute_crest() -> Cluster:
        kwargs: dict[str, Any] = {
            "n_images": neb_images,
            "fmax_ev_a": 0.08,
            "max_steps": neb_steps,
            "pre_relax_fmax_ev_a": 0.20,
            "pre_relax_steps": neb_pre_steps,
            "optimizer_dt": 0.02,
            "optimizer_maxstep": 0.03,
            "climb_optimizer": "ode",
            "checkpoint_dir": neb_root,
            "checkpoint_interval": 5,
        }
        if resume_from is not None:
            kwargs["resume_from"] = resume_from
        return neb_ts_guess(start, end, settings, **kwargs)

    crest = a2.checkpoint_cluster(
        segment_dir / "neb-crest.xyz",
        start,
        compute_crest,
        identity={
            "gate_version": PATH_GATE_VERSION,
            "stage": "persisted-preclimb-climb-ci-neb",
            "segment": spec.slug,
            "start_geometry": frequency_geometry_fingerprint(start),
            "end_geometry": frequency_geometry_fingerprint(end),
            "settings": settings_identity,
            "n_images": neb_images,
            "pre_relax_steps": neb_pre_steps,
            "climb_steps": neb_steps,
            "climb_optimizer": "ode",
        },
    )
    local_tangent, tangent_receipt = final_neb_tangent(
        neb_root,
        crest,
        start,
        active_indices,
    )
    atomic_json(segment_dir / "local-neb-tangent.receipt.json", tangent_receipt)
    coordinate_targets = reaction_coordinate_targets(crest)
    conditioned_crest = a2.checkpoint_cluster(
        segment_dir / "conditioned-crest.xyz",
        crest,
        lambda: relax_at_fixed_distances(
            crest,
            settings,
            fixed_distances=coordinate_targets,
            fmax_ev_a=0.02,
            max_steps=120,
            optimizer_maxstep=0.03,
            trajectory=str(segment_dir / "conditioned-crest.traj"),
            logfile=str(segment_dir / "conditioned-crest.log"),
        ),
        identity={
            "gate_version": PATH_GATE_VERSION,
            "stage": "three-coordinate-conditioned-crest",
            "segment": spec.slug,
            "settings": settings_identity,
            "fixed_distances": coordinate_targets,
            "fmax_ev_a": 0.02,
            "max_steps": 120,
            "optimizer_maxstep": 0.03,
        },
    )
    transition_state = a2.checkpoint_cluster(
        segment_dir / "transition-state.xyz",
        conditioned_crest,
        lambda: find_ts(
            conditioned_crest,
            settings,
            max_steps=saddle_steps,
            trajectory=str(segment_dir / "transition-state.sella.traj"),
            initial_mode=local_tangent,
            internal=False,
        ),
        identity={
            "gate_version": PATH_GATE_VERSION,
            "stage": "reaction-aligned-sella",
            "segment": spec.slug,
            "start_geometry": frequency_geometry_fingerprint(start),
            "end_geometry": frequency_geometry_fingerprint(end),
            "settings": settings_identity,
            "max_steps": saddle_steps,
            "active_indices": active_indices,
            "initial_mode_strategy": SELLA_MODE_STRATEGY,
            "local_tangent": tangent_receipt,
            "conditioned_crest_geometry": frequency_geometry_fingerprint(
                conditioned_crest
            ),
        },
    )
    primary_frequency = checkpoint_fd_frequency(
        segment_dir / "transition-state.fd-1e-3.frequency.json",
        transition_state,
        settings,
        step_bohr=HESSIAN_STEPS_BOHR[0],
    )
    secondary_frequency = checkpoint_fd_frequency(
        segment_dir / "transition-state.fd-2e-3.frequency.json",
        transition_state,
        settings,
        step_bohr=HESSIAN_STEPS_BOHR[1],
    )
    stability = require_hessian_step_stability(
        primary_frequency,
        secondary_frequency,
        local_tangent,
    )
    irc_endpoints = a2.checkpoint_irc_endpoints(
        segment_dir,
        transition_state,
        lambda: full_irc(
            transition_state,
            settings,
            max_steps=irc_steps,
            trajectory=segment_dir / "full-irc.traj",
            logfile=segment_dir / "full-irc.log",
        ),
        identity={
            "gate_version": PATH_GATE_VERSION,
            "stage": "full-irc",
            "segment": spec.slug,
            "algorithm": "sella-gonzalez-schlegel-full-irc-v1",
            "settings": settings_identity,
            "max_steps": irc_steps,
            "step_size_a": 0.10,
            "fmax_ev_a": 0.05,
            "fmax_inner_ev_a": 0.01,
        },
    )
    typed_irc = a2.require_si_neutral_irc(
        irc_endpoints,
        start,
        end,
        attacker_index,
    )
    displacements = [
        endpoint_displacement(transition_state, endpoint) for endpoint in irc_endpoints
    ]
    if min(displacements) <= IRC_MINIMUM_CARTESIAN_DISPLACEMENT_A:
        raise RuntimeError(
            f"{spec.slug}: full IRC contains a zero/near-zero direction: "
            f"{displacements}"
        )
    return {
        "spec": asdict(spec),
        "transition_state": transition_state,
        "frequency": primary_frequency,
        "hessian_step_stability": stability,
        "full_irc": {
            **typed_irc,
            "cartesian_displacement_a": displacements,
            "minimum_required_a": IRC_MINIMUM_CARTESIAN_DISPLACEMENT_A,
        },
        "artifacts": {
            "neb_checkpoints": str(neb_root),
            "neb_crest": str(segment_dir / "neb-crest.xyz"),
            "conditioned_crest": str(segment_dir / "conditioned-crest.xyz"),
            "local_tangent_receipt": str(
                segment_dir / "local-neb-tangent.receipt.json"
            ),
            "transition_state": str(segment_dir / "transition-state.xyz"),
            "primary_frequency": str(
                segment_dir / "transition-state.fd-1e-3.frequency.json"
            ),
            "secondary_frequency": str(
                segment_dir / "transition-state.fd-2e-3.frequency.json"
            ),
            "irc_receipt": str(segment_dir / "irc.r2scan3c.receipt.json"),
        },
    }


def record_store(
    path: Path,
    clusters: dict[str, Cluster],
    frequencies: dict[str, FrequencyResult],
    energies: dict[str, dict[str, float]],
    results: dict[str, Any],
    *,
    gpu: bool,
) -> None:
    temporary = path.with_name(f".{path.name}.{time.time_ns()}.tmp")
    temporary.unlink(missing_ok=True)
    with Store(temporary) as store:
        structure_ids: dict[str, int] = {}
        for role, cluster in clusters.items():
            structure_ids[role] = store.add_structure(
                f"si-neutral-{role}",
                cluster.formula,
                a2.exact_xyz(cluster),
                charge=cluster.charge,
                spin=cluster.spin,
            )
        for role, result in frequencies.items():
            job = store.add_job(
                structure_ids[role],
                "freq",
                a2.R2SCAN3C_METHOD,
                "gpu4pyscf" if gpu else "pyscf",
                detail=json.dumps(
                    {
                        "hessian": "central-difference-full-composite-gradient",
                        "step_bohr": HESSIAN_STEPS_BOHR[0],
                        "path_gate_version": PATH_GATE_VERSION,
                    },
                    sort_keys=True,
                ),
            )
            store.set_job_status(job, "done")
            store.add_result(job, "electronic", result.electronic_hartree, "hartree")
            store.add_result(job, "imaginary_count", result.n_imaginary, "count")
        for method, role_values in energies.items():
            for role, value in role_values.items():
                job = store.add_job(
                    structure_ids[role],
                    "sp",
                    method,
                    "gpu4pyscf" if gpu else "pyscf",
                    detail=json.dumps(
                        {"path_gate_version": PATH_GATE_VERSION}, sort_keys=True
                    ),
                )
                store.set_job_status(job, "done")
                store.add_result(job, "electronic", value, "hartree")
        for segment, payload in results["segments"].items():
            role = f"{segment}-transition-state"
            job = store.add_job(
                structure_ids[role],
                "analysis",
                a2.PRODUCTION_METHOD,
                "quarry",
                detail=json.dumps({"temperature_k": a2.TEMPERATURE_K}, sort_keys=True),
            )
            store.set_job_status(job, "done")
            store.add_result(job, "local_dg_kj", payload["local_dg_kj"], "kJ/mol")
            store.add_result(
                job,
                "reactant_referenced_dg_kj",
                payload["reactant_referenced_dg_kj"],
                "kJ/mol",
            )
    temporary.replace(path)


def run(args: argparse.Namespace) -> int:
    run_dir = args.run_dir.resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    reactant_source = load_reference(args.reactant_reference.resolve(), role="reactant")
    intermediate_source = load_reference(
        args.intermediate_reference.resolve(),
        role="intermediate",
        template=reactant_source,
    )
    product_source = load_reference(
        args.product_reference.resolve(),
        role="product",
        template=reactant_source,
    )
    for value in active_indices(args.active_indices, len(reactant_source.symbols)):
        if value in reactant_source.frozen_indices:
            raise ValueError("active-indices contains a frozen atom")
    selected_active_indices = active_indices(
        args.active_indices, len(reactant_source.symbols)
    )

    sources = {
        role: {
            "path": str(path.resolve()),
            "sha256": a2.sha256_path(path.resolve()),
        }
        for role, path in {
            "reactant": args.reactant_reference,
            "intermediate": args.intermediate_reference,
            "product": args.product_reference,
        }.items()
    }
    atomic_json(
        run_dir / "source-receipt.json",
        {
            "gate_version": PATH_GATE_VERSION,
            "sources": sources,
            "rejected_banked_transition_state_used": False,
        },
    )

    r2scan3c, production, b3lyp_d4 = a2.settings(use_gpu=args.gpu)
    atomic_json(
        run_dir / "settings.json",
        {
            "gate_version": PATH_GATE_VERSION,
            "r2scan3c": asdict(r2scan3c),
            "production": asdict(production),
            "b3lyp_d4": asdict(b3lyp_d4),
            "hessian_steps_bohr": list(HESSIAN_STEPS_BOHR),
            "stability_contract": {
                "significant_imaginary_floor_cm": SIGNIFICANT_IMAGINARY_FLOOR_CM,
                "mode_cosine_minimum": MODE_COSINE_MINIMUM,
                "reaction_vector_overlap_minimum": REACTION_VECTOR_OVERLAP_MINIMUM,
                "frequency_relative_drift_limit": FREQUENCY_RELATIVE_DRIFT_LIMIT,
                "frequency_absolute_drift_limit_cm": FREQUENCY_ABSOLUTE_DRIFT_LIMIT_CM,
                "irc_minimum_cartesian_displacement_a": (
                    IRC_MINIMUM_CARTESIAN_DISPLACEMENT_A
                ),
            },
            "bounds": {
                "minimum_steps": args.minimum_steps,
                "neb_images": args.neb_images,
                "neb_pre_steps": args.neb_pre_steps,
                "neb_steps": args.neb_steps,
                "saddle_steps": args.saddle_steps,
                "irc_steps": args.irc_steps,
            },
            "active_indices": selected_active_indices,
        },
    )

    references = {
        "reactant": reactant_source,
        "intermediate": intermediate_source,
        "product": product_source,
    }
    minima: dict[str, Cluster] = {}
    for role, source in references.items():
        minima[role] = a2.checkpoint_cluster(
            run_dir / f"{role}.r2scan3c.xyz",
            source,
            lambda source=source, role=role: a2.optimize_minimum(
                source,
                r2scan3c,
                max_steps=args.minimum_steps,
                trajectory=run_dir / f"{role}.r2scan3c.traj",
            ),
            identity={
                "gate_version": PATH_GATE_VERSION,
                "stage": "minimum",
                "role": role,
                "settings": frequency_settings_fingerprint(r2scan3c),
                "max_steps": args.minimum_steps,
            },
        )
    for role, expected in {
        "reactant": REACTANT_BASIN,
        "intermediate": ASSOCIATIVE_BASIN,
        "product": PRODUCT_BASIN,
    }.items():
        require_basin(minima[role], args.attacker_index, expected)

    frequencies: dict[str, FrequencyResult] = {}
    for role, minimum in minima.items():
        result = checkpoint_fd_frequency(
            run_dir / f"{role}.fd-1e-3.frequency.json",
            minimum,
            r2scan3c,
            step_bohr=HESSIAN_STEPS_BOHR[0],
        )
        if result.n_imaginary != 0:
            raise RuntimeError(
                f"{role}: r2SCAN-3c endpoint is not a zero-index minimum"
            )
        frequencies[role] = result

    segment_runs: dict[str, dict[str, Any]] = {}
    clusters = dict(minima)
    for spec in SEGMENTS:
        segment = run_segment(
            run_dir,
            spec,
            minima[spec.start_role],
            minima[spec.end_role],
            r2scan3c,
            attacker_index=args.attacker_index,
            active_indices=selected_active_indices,
            neb_images=args.neb_images,
            neb_steps=args.neb_steps,
            neb_pre_steps=args.neb_pre_steps,
            saddle_steps=args.saddle_steps,
            irc_steps=args.irc_steps,
        )
        transition_role = f"{spec.slug}-transition-state"
        clusters[transition_role] = segment.pop("transition_state")
        frequencies[transition_role] = segment.pop("frequency")
        segment_runs[spec.slug] = segment

    method_settings: dict[str, DftSettings] = {
        a2.PRODUCTION_METHOD: production,
        a2.B3LYP_D4_METHOD: b3lyp_d4,
    }
    energies: dict[str, dict[str, float]] = {}
    for method, method_settings_value in method_settings.items():
        slug = method.split("/")[0].replace("(", "-").replace(")", "")
        energies[method] = {}
        for role, cluster in clusters.items():
            energies[method][role] = a2.checkpoint_energy(
                run_dir / f"{role}.{slug}.energy.json",
                cluster,
                method_settings_value,
                method,
            )

    reactant_thermo = a2.thermo(
        frequencies["reactant"], energies[a2.PRODUCTION_METHOD]["reactant"]
    )
    segment_results: dict[str, dict[str, Any]] = {}
    for spec in SEGMENTS:
        transition_role = f"{spec.slug}-transition-state"
        start_thermo = a2.thermo(
            frequencies[spec.start_role],
            energies[a2.PRODUCTION_METHOD][spec.start_role],
        )
        transition_thermo = a2.thermo(
            frequencies[transition_role],
            energies[a2.PRODUCTION_METHOD][transition_role],
        )
        imaginary = _one_significant_mode(frequencies[transition_role])[0]
        rate = rate_from_thermo(
            start_thermo,
            transition_thermo,
            imag_nu_cm=imaginary,
            tunneling="wigner",
        )
        payload = segment_runs[spec.slug]
        payload.update(
            {
                "start_role": spec.start_role,
                "end_role": spec.end_role,
                "local_dg_kj": rate.dg_kj,
                "local_dh_kj": rate.dh_kj,
                "local_ds_kj_per_k": rate.ds_kj_per_k,
                "wigner_kappa": rate.kappa,
                "k_per_s": rate.k,
                "reactant_referenced_dg_kj": transition_thermo.gibbs
                - reactant_thermo.gibbs,
                "imaginary_cm": imaginary,
            }
        )
        segment_results[spec.slug] = payload

    rate_limiting = max(
        segment_results,
        key=lambda name: segment_results[name]["reactant_referenced_dg_kj"],
    )
    electronic_barriers: dict[str, dict[str, float]] = {}
    for method, values in energies.items():
        electronic_barriers[method] = {
            spec.slug: (
                values[f"{spec.slug}-transition-state"] - values[spec.start_role]
            )
            * HARTREE_TO_KJ
            for spec in SEGMENTS
        }
        electronic_barriers[method]["overall_profile"] = max(
            (values[f"{spec.slug}-transition-state"] - values["reactant"])
            * HARTREE_TO_KJ
            for spec in SEGMENTS
        )
    electronic_barriers[a2.R2SCAN3C_METHOD] = {
        spec.slug: (
            frequencies[f"{spec.slug}-transition-state"].electronic_hartree
            - frequencies[spec.start_role].electronic_hartree
        )
        * HARTREE_TO_KJ
        for spec in SEGMENTS
    }
    electronic_barriers[a2.R2SCAN3C_METHOD]["overall_profile"] = max(
        (
            frequencies[f"{spec.slug}-transition-state"].electronic_hartree
            - frequencies["reactant"].electronic_hartree
        )
        * HARTREE_TO_KJ
        for spec in SEGMENTS
    )

    results = {
        "reaction": "si-neutral",
        "mechanism": "sequential-associative",
        "path_gate_version": PATH_GATE_VERSION,
        "temperature_k": a2.TEMPERATURE_K,
        "rejected_banked_transition_state_used": False,
        "route": [
            "reactant",
            "addition-transition-state",
            "intermediate",
            "cleavage-transition-state",
            "product",
        ],
        "typed_minima": {
            role: typed_identity_payload(cluster, args.attacker_index)
            for role, cluster in minima.items()
        },
        "segments": segment_results,
        "rate_limiting_segment": rate_limiting,
        "overall_profile_dg_kj": segment_results[rate_limiting][
            "reactant_referenced_dg_kj"
        ],
        "electronic_barriers_kj": electronic_barriers,
        "geometry_hashes": {
            role: geometry_hash(a2.exact_xyz(cluster))
            for role, cluster in clusters.items()
        },
        "cc_calibration_roles": {
            "reactant": str(run_dir / "reactant.r2scan3c.xyz"),
            "ts": str(run_dir / rate_limiting / "transition-state.xyz"),
            "rate_limiting_segment": rate_limiting,
        },
    }
    atomic_json(run_dir / "results.json", results)
    record_store(
        run_dir / "store.sqlite",
        clusters,
        frequencies,
        energies,
        results,
        gpu=args.gpu,
    )
    print(json.dumps(results, indent=2, sort_keys=True))
    return 0


def active_indices(value: str, atom_count: int) -> list[int]:
    try:
        result = [int(field) for field in value.split(",") if field.strip()]
    except ValueError as exc:
        raise ValueError("active-indices must be comma-separated integers") from exc
    if not result or len(set(result)) != len(result):
        raise ValueError("active-indices must be non-empty and unique")
    if any(index < 0 or index >= atom_count for index in result):
        raise ValueError("active-indices contains an out-of-range atom")
    return result


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--reactant-reference", type=Path, required=True)
    result.add_argument("--intermediate-reference", type=Path, required=True)
    result.add_argument("--product-reference", type=Path, required=True)
    result.add_argument("--run-dir", type=Path, required=True)
    result.add_argument("--attacker-index", type=int, default=15)
    result.add_argument("--active-indices", default="0,1,2,15,16,17")
    result.add_argument("--gpu", action="store_true")
    result.add_argument("--minimum-steps", type=int, default=300)
    result.add_argument("--neb-images", type=int, default=9)
    result.add_argument("--neb-pre-steps", type=int, default=180)
    result.add_argument("--neb-steps", type=int, default=240)
    result.add_argument("--saddle-steps", type=int, default=300)
    result.add_argument("--irc-steps", type=int, default=400)
    result.add_argument("--threads", type=int, default=16)
    result.add_argument("--nice", type=int, default=10)
    result.add_argument("--log")
    return result


def execute_with_status(args: argparse.Namespace) -> int:
    status_path = args.run_dir.resolve() / "run_status.json"
    atomic_json(
        status_path,
        {
            "status": "running",
            "gate_version": PATH_GATE_VERSION,
            "started_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        },
    )
    try:
        code = run(args)
    except Exception as exc:
        atomic_json(
            status_path,
            {
                "status": "failed",
                "gate_version": PATH_GATE_VERSION,
                "error_type": type(exc).__name__,
                "error": str(exc),
                "failed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            },
        )
        raise
    atomic_json(
        status_path,
        {
            "status": "completed",
            "gate_version": PATH_GATE_VERSION,
            "completed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "results_sha256": a2.sha256_path(args.run_dir.resolve() / "results.json"),
            "store_sha256": a2.sha256_path(args.run_dir.resolve() / "store.sqlite"),
        },
    )
    return code


def main() -> int:
    args = parser().parse_args()
    for name in (
        "minimum_steps",
        "neb_images",
        "neb_pre_steps",
        "neb_steps",
        "saddle_steps",
        "irc_steps",
    ):
        if getattr(args, name) <= 0:
            raise ValueError(f"{name.replace('_', '-')} must be positive")
    if args.neb_images < 3:
        raise ValueError("neb-images must be at least 3")
    if args.attacker_index < 0:
        raise ValueError("attacker-index must be non-negative")
    return execute_with_status(args)


if __name__ == "__main__":
    raise SystemExit(main())
