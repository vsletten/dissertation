#!/usr/bin/env python3
"""Validate and close the production tier for the barrierless A2a route.

This driver consumes only the accepted R->addition-TS->I route, the classified
barrierless I->P shelf, and the fresh released product minimum.  It has no path
search, optimization, Hessian, or IRC execution route.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    _ETIQUETTE = bootstrap_cli(
        "a2a-barrierless-closeout",
        default_run_root=(
            "/mnt/data/vsletten/dissertation-data/"
            "task208-a2a-path-rebuild-20260825/production-closeout"
        ),
    )

import numpy as np

from quarry import pipeline
from quarry.clusters import Cluster
from quarry.pipeline import HARTREE_TO_KJ, DftSettings, FrequencyResult
from scripts import a2a_path_rebuild as path_driver
from scripts import production_energetics as production

CLOSEOUT_VERSION = "a2a-barrierless-production-closeout-v1"
SVP_METHOD = "b3lyp/def2-svp/df"
ROLE_PATHS = {
    "reactant": "reactant.r2scan3c.xyz",
    "intermediate": "intermediate.r2scan3c.xyz",
    "addition-transition-state": "addition/transition-state.xyz",
    "released-product": (
        "cleavage/barrierless-downhill-release-v2/released-product.xyz"
    ),
}
ACCEPTED_SHA256 = {
    "reactant": "e3b0b1fb8ab62cbf3c88df90281dcc5f4f8fb68a6cc0c34c5ca6dae75fd2f359",
    "intermediate": "bd9cab8e9198d520d1727f020a2154826513bd64fc20dbde6faec31c00d3c006",
    "addition-transition-state": (
        "d7c5ff9ae2e850f5149538ca5775e9802d3a5680b3c9e75078096fc543b5419b"
    ),
    "released-product": (
        "bf30ff37d08cd30bc523b7d286e3e5e29e2056f8e61c3e5e82b3e8732d232bab"
    ),
}


@dataclass(frozen=True)
class AcceptedRoute:
    evidence_root: Path
    clusters: dict[str, Cluster]
    frequencies: dict[str, FrequencyResult]
    validation: dict[str, Any]


def _json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path}: expected a JSON object")
    return payload


def _require_sha(path: Path, expected: str) -> None:
    actual = production.sha256_path(path)
    if actual != expected:
        raise ValueError(f"{path}: SHA-256 drift; expected {expected}, found {actual}")


def _load_frequency(path: Path) -> FrequencyResult:
    return path_driver.frequency_from_payload(_json(path))


def _classification_and_release(root: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    classification_path = (
        root / "cleavage/coupled-scan-v1/complete-grid-classification.json"
    )
    release_path = (
        root / "cleavage/barrierless-downhill-release-v2/release-receipt.json"
    )
    classification = _json(classification_path)
    release = _json(release_path)
    decision = classification.get("classification", {})
    if classification.get("cell_count") != 81:
        raise RuntimeError(
            "barrierless classification is not the complete 81-cell grid"
        )
    if decision.get("outcome") != "barrierless-shelf":
        raise RuntimeError("I->P classification is not barrierless-shelf")
    if float(decision.get("barrier_threshold_kj_mol", np.nan)) != 2.0:
        raise RuntimeError("barrierless classification threshold drift")
    if decision.get("product_cell") != [7, 7]:
        raise RuntimeError("barrierless classification product-cell drift")
    if decision.get("candidate_kind") != (
        "complete-grid-valley-no-interior-crest-above-threshold"
    ):
        raise RuntimeError("barrierless classification candidate-kind drift")
    if release.get("release_version") != "a2a-barrierless-downhill-release-v2":
        raise RuntimeError("released-product version drift")
    if release.get("status") != "completed":
        raise RuntimeError("released-product receipt is incomplete")
    if release.get("seed_cell") != decision.get("product_cell"):
        raise RuntimeError(
            "released product did not start at the classified product cell"
        )
    if release.get("zero_index_minimum") is not True:
        raise RuntimeError("released product is not a zero-index minimum")
    if release.get("typed_product_identity_matches") is not True:
        raise RuntimeError("released product does not preserve exact typed identity")
    if int(release.get("imaginary_mode_count", -1)) != 0:
        raise RuntimeError("released product has an imaginary mode")
    if float(release.get("electronic_delta_kj_mol", np.nan)) >= 0.0:
        raise RuntimeError("released product did not move downhill")
    if float(release.get("cartesian_displacement_a", 0.0)) < float(
        release.get("minimum_cartesian_displacement_a", np.inf)
    ):
        raise RuntimeError("released product did not move a nonzero accepted distance")
    return classification, release


def load_accepted_route(evidence_root: Path, attacker_index: int = 15) -> AcceptedRoute:
    root = evidence_root.resolve()
    if not root.is_dir():
        raise FileNotFoundError(root)
    for role, relative in ROLE_PATHS.items():
        _require_sha(root / relative, ACCEPTED_SHA256[role])

    reactant = path_driver.load_reference(
        root / ROLE_PATHS["reactant"], role="reactant"
    )
    clusters = {
        "reactant": reactant,
        "intermediate": path_driver.load_reference(
            root / ROLE_PATHS["intermediate"],
            role="intermediate",
            template=reactant,
        ),
        "addition-transition-state": path_driver.load_reference(
            root / ROLE_PATHS["addition-transition-state"],
            role="addition-transition-state",
            template=reactant,
        ),
        "released-product": path_driver.load_reference(
            root / ROLE_PATHS["released-product"],
            role="released-product",
            template=reactant,
        ),
    }
    path_driver.require_basin(
        clusters["reactant"], attacker_index, path_driver.REACTANT_BASIN
    )
    path_driver.require_basin(
        clusters["intermediate"], attacker_index, path_driver.ASSOCIATIVE_BASIN
    )
    path_driver.require_basin(
        clusters["released-product"], attacker_index, path_driver.PRODUCT_BASIN
    )

    frequencies = {
        "reactant": _load_frequency(root / "reactant.fd-1e-3.frequency.json"),
        "intermediate": _load_frequency(root / "intermediate.fd-1e-3.frequency.json"),
        "addition-transition-state": _load_frequency(
            root / "addition/transition-state.fd-1e-3.frequency.json"
        ),
        "released-product": _load_frequency(
            root / "cleavage/barrierless-downhill-release-v2/"
            "released-product.fd.frequency.json"
        ),
    }
    for role in ("reactant", "intermediate", "released-product"):
        if frequencies[role].n_imaginary != 0:
            raise RuntimeError(f"{role} is not a zero-index minimum")

    secondary_ts = _load_frequency(
        root / "addition/transition-state.fd-2e-3.frequency.json"
    )
    neb_crest = path_driver.load_reference(
        root / "addition/neb-crest.xyz",
        role="addition-neb-crest",
        template=reactant,
    )
    tangent, tangent_receipt = path_driver.final_neb_tangent(
        root / "addition/neb-checkpoints",
        neb_crest,
        reactant,
        [0, 1, 2, 15, 16, 17],
    )
    hessian_stability = path_driver.require_hessian_step_stability(
        frequencies["addition-transition-state"], secondary_ts, tangent
    )

    irc_receipt = _json(root / "addition/irc.r2scan3c.receipt.json")
    expected_ts_fingerprint = production.frequency_geometry_fingerprint(
        clusters["addition-transition-state"]
    )
    if irc_receipt.get("transition_state_geometry_fingerprint") != (
        expected_ts_fingerprint
    ):
        raise RuntimeError("accepted addition IRC transition-state drift")
    endpoint_basins: list[tuple[bool, bool, bool]] = []
    endpoint_displacements: list[float] = []
    for direction in ("forward", "reverse"):
        direction_receipt = irc_receipt.get("directions", {}).get(direction, {})
        endpoint_path = root / "addition" / str(direction_receipt.get("path", ""))
        _require_sha(endpoint_path, str(direction_receipt.get("sha256", "")))
        endpoint = path_driver.load_reference(
            endpoint_path,
            role=f"irc-{direction}",
            template=reactant,
        )
        basin = production.si_neutral_signature(endpoint, attacker_index)
        production.endpoint_identity(endpoint, attacker_index)
        endpoint_basins.append(basin)
        endpoint_displacements.append(
            path_driver.endpoint_displacement(
                clusters["addition-transition-state"], endpoint
            )
        )
    if set(endpoint_basins) != {
        path_driver.REACTANT_BASIN,
        path_driver.ASSOCIATIVE_BASIN,
    }:
        raise RuntimeError(f"accepted addition IRC basin drift: {endpoint_basins}")
    if min(endpoint_displacements) < path_driver.IRC_MINIMUM_CARTESIAN_DISPLACEMENT_A:
        raise RuntimeError("accepted addition IRC contains a zero-step direction")

    classification, release = _classification_and_release(root)
    validation = {
        "closeout_version": CLOSEOUT_VERSION,
        "validated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "route": [
            "reactant",
            "addition-transition-state",
            "intermediate",
            "barrierless-shelf",
            "released-product",
        ],
        "rejected_banked_transition_state_used": False,
        "cleavage": {
            "verified_saddle": False,
            "classification": classification["classification"],
            "classification_sha256": production.sha256_path(
                root / "cleavage/coupled-scan-v1/complete-grid-classification.json"
            ),
            "release_receipt_sha256": production.sha256_path(
                root / "cleavage/barrierless-downhill-release-v2/release-receipt.json"
            ),
            "released_product_delta_kj_mol": release["electronic_delta_kj_mol"],
        },
        "addition_transition_state": {
            "hessian_step_stability": hessian_stability,
            "local_tangent": tangent_receipt,
            "irc_endpoint_basins": [list(value) for value in endpoint_basins],
            "irc_displacements_a": endpoint_displacements,
        },
        "geometry_sha256": dict(ACCEPTED_SHA256),
    }
    return AcceptedRoute(root, clusters, frequencies, validation)


def method_settings(use_gpu: bool) -> dict[str, DftSettings]:
    _, wb97m_v, b3lyp_d4 = production.settings(use_gpu=use_gpu)
    return {
        SVP_METHOD: DftSettings(
            xc="b3lyp",
            basis="def2-svp",
            density_fit=True,
            use_gpu=use_gpu,
        ),
        production.PRODUCTION_METHOD: wb97m_v,
        production.B3LYP_D4_METHOD: b3lyp_d4,
    }


def checkpoint_closeout_energy(
    path: Path,
    cluster: Cluster,
    settings: DftSettings,
    method: str,
) -> float:
    """Run a receipt-bound SCF with a bounded Newton fallback.

    The GPU SMD wB97M-V direct-DIIS route exhausted its 150-cycle bound on the
    accepted reactant.  Newton is an SCF solver change, not a method change; it
    is predeclared here and persisted in the receipt rather than silently
    loosening convergence or accepting the unconverged direct value.
    """
    if path.exists():
        return production.checkpoint_energy(path, cluster, settings, method)

    expected_geometry = production.frequency_geometry_fingerprint(cluster)
    expected_settings = production.frequency_settings_fingerprint(settings)
    mf = pipeline._make_scf(pipeline.build_mol(cluster, settings), settings)
    convergence_route = "direct-diis"
    if method == production.PRODUCTION_METHOD:
        convergence_route = "newton-first"
        mf = mf.newton()
        mf.max_cycle = 100
        value = mf.kernel()
    else:
        value = mf.kernel()
        if not mf.converged:
            density = mf.make_rdm1()
            mf = mf.newton()
            mf.max_cycle = 100
            value = mf.kernel(dm0=density)
            convergence_route = "direct-diis-then-newton"
    if not mf.converged:
        raise RuntimeError(
            f"SCF did not converge for {cluster.name} via {convergence_route}"
        )
    value = float(value)
    if not np.isfinite(value):
        raise RuntimeError(f"{path.name}: computed electronic energy is non-finite")
    production.atomic_json(
        path,
        {
            "method": method,
            "electronic_hartree": value,
            "geometry_fingerprint": expected_geometry,
            "settings_fingerprint": expected_settings,
            "scf_contract": "bounded-direct-diis-or-newton-v1",
            "convergence_route": convergence_route,
            "converged": True,
        },
    )
    return value


def run_dft(route: AcceptedRoute, *, use_gpu: bool) -> dict[str, Any]:
    path_driver.preflight_gpu_contraction_engine(use_gpu)
    output_root = route.evidence_root / "production-closeout"
    energy_root = output_root / "energies"
    energy_root.mkdir(parents=True, exist_ok=True)
    production.atomic_json(output_root / "route-validation.json", route.validation)

    energies: dict[str, dict[str, float]] = {}
    settings_by_method = method_settings(use_gpu)
    for method, settings in settings_by_method.items():
        slug = method.split("/")[0].replace("(", "-").replace(")", "")
        energies[method] = {}
        for role, cluster in route.clusters.items():
            energies[method][role] = checkpoint_closeout_energy(
                energy_root / f"{role}.{slug}.energy.json",
                cluster,
                settings,
                method,
            )

    reactant = "reactant"
    transition_state = "addition-transition-state"
    r2scan_barrier = (
        route.frequencies[transition_state].electronic_hartree
        - route.frequencies[reactant].electronic_hartree
    ) * HARTREE_TO_KJ
    barriers = {
        method: (values[transition_state] - values[reactant]) * HARTREE_TO_KJ
        for method, values in energies.items()
    }
    wb97m_reactant_thermo = production.thermo(
        route.frequencies[reactant],
        energies[production.PRODUCTION_METHOD][reactant],
    )
    wb97m_ts_thermo = production.thermo(
        route.frequencies[transition_state],
        energies[production.PRODUCTION_METHOD][transition_state],
    )
    payload = {
        "closeout_version": CLOSEOUT_VERSION,
        "status": "dft-completed",
        "route_validation_sha256": production.sha256_path(
            output_root / "route-validation.json"
        ),
        "methods": {key: asdict(value) for key, value in settings_by_method.items()},
        "energies_hartree": energies,
        "electronic_barriers_kj_mol": {
            production.R2SCAN3C_METHOD: r2scan_barrier,
            **barriers,
        },
        "same_geometry_barrier_shifts_kj_mol": {
            "svp_to_r2scan3c": r2scan_barrier - barriers[SVP_METHOD],
            "r2scan3c_to_wb97m_v": (
                barriers[production.PRODUCTION_METHOD] - r2scan_barrier
            ),
            "r2scan3c_to_b3lyp_d4": (
                barriers[production.B3LYP_D4_METHOD] - r2scan_barrier
            ),
        },
        "production_dg_dagger_kj_mol": (
            wb97m_ts_thermo.gibbs - wb97m_reactant_thermo.gibbs
        ),
        "cleavage_verified_saddle": False,
        "rejected_banked_transition_state_used": False,
    }
    production.atomic_json(output_root / "dft-summary.json", payload)
    return payload


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command", required=True)
    for command in ("validate", "dft"):
        subparser = subparsers.add_parser(command)
        subparser.add_argument("--evidence-root", type=Path, required=True)
        subparser.add_argument("--attacker-index", type=int, default=15)
        if command == "dft":
            subparser.add_argument("--gpu", action="store_true")
    return result


def main() -> int:
    args = parser().parse_args()
    if args.attacker_index < 0:
        raise ValueError("attacker-index must be non-negative")
    route = load_accepted_route(args.evidence_root, args.attacker_index)
    output_root = route.evidence_root / "production-closeout"
    output_root.mkdir(parents=True, exist_ok=True)
    production.atomic_json(output_root / "route-validation.json", route.validation)
    if args.command == "validate":
        print(json.dumps(route.validation, indent=2, sort_keys=True))
        return 0
    payload = run_dft(route, use_gpu=args.gpu)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
