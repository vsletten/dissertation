from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from scripts import d2c_input_bundle as bundle
from scripts import d2c_sct_campaign as campaign

FIXED_GIT_SHA = "a" * 40
FIXED_DEPENDENCIES = {
    "ase": "3.29.0",
    "cuda_device_count": "1",
    "cuda_driver": "13000",
    "cuda_runtime": "12090",
    "cupy": "14.1.1",
    "cupy_distribution": "cupy-cuda12x",
    "geometric": "1.1",
    "gpu4pyscf": "1.8.1",
    "gpu4pyscf_distribution": "gpu4pyscf-cuda12x",
    "gpu_backend": "gpu4pyscf",
    "numpy": "2.3.2",
    "pyscf": "2.10.0",
    "sella": "2.5.1",
}
BUNDLE_ROOT = Path(__file__).parents[1] / "data" / "D2c-instanton-tier" / "d2b-inputs"


def _preflight(run_root: Path, bundle_root: Path = BUNDLE_ROOT) -> dict:
    return campaign.create_preflight_receipt(
        bundle_root,
        run_root,
        git_sha=FIXED_GIT_SHA,
        dependency_versions=FIXED_DEPENDENCIES,
        created_utc="2026-09-07T12:00:00Z",
    )


def test_dry_run_identity_and_exact_four_route_inventory(tmp_path: Path):
    first = _preflight(tmp_path / "first")
    second = _preflight(tmp_path / "second")
    changed_gpu = dict(FIXED_DEPENDENCIES)
    changed_gpu["gpu4pyscf"] = "1.8.2"
    third = campaign.create_preflight_receipt(
        BUNDLE_ROOT,
        tmp_path / "third",
        git_sha=FIXED_GIT_SHA,
        dependency_versions=changed_gpu,
        created_utc="2026-09-07T12:00:00Z",
    )

    assert first["state"] == "pending"
    assert first["identity"] == second["identity"]
    assert first["identity"] != third["identity"]
    assert first["campaign"]["git_sha"] == FIXED_GIT_SHA
    assert first["campaign"]["bundle_manifest_sha256"] == bundle.sha256_path(
        BUNDLE_ROOT / "manifest.json"
    )
    assert list(first["routes"]) == list(bundle.ROUTES)
    assert len(first["routes"]) == 4
    assert first["campaign"]["dependencies"] == FIXED_DEPENDENCIES
    assert first["campaign"]["bounds"]["sct"][
        "straightness_tolerance_per_angstrom"
    ] == pytest.approx(1.0e-12)
    assert first["campaign"]["bounds"]["hessian"] == {
        "cartesian_hessian_units": "hartree / bohr^2",
        "coverage": "every retained IRC point including TS and endpoints",
        "maximum_retained_points_per_direction": 201,
        "required_transverse_mode_count": "3N-7",
        "transverse_eigenvalue_units": "hartree / bohr^2 / amu",
        "transverse_negative_eigenvalue_tolerance": 1.0e-8,
        "transverse_negative_eigenvalue_tolerance_units": ("hartree / bohr^2 / amu"),
    }
    for route, route_receipt in first["routes"].items():
        assert route_receipt["symbols"]
        assert len(route_receipt["atom_identity_labels"]) == len(
            route_receipt["symbols"]
        )
        assert len(set(route_receipt["atom_identity_labels"])) == len(
            route_receipt["symbols"]
        )
        assert len(route_receipt["atom_mapping_sha256"]) == 64
        assert (
            route_receipt["canonical_transition_state_geometry_sha256"]
            == (campaign.TRUSTED_CANONICAL_TS_GEOMETRY_SHA256[route])
        )
        assert len(route_receipt["masses_amu"]) == len(route_receipt["symbols"])
        assert all(
            stage["state"] == "pending" for stage in route_receipt["stages"].values()
        )
    assert first["routes"]["h-co-1w-oside"]["atom_identity_labels"][:3] == [
        "carbon_monoxide_carbon",
        "carbon_monoxide_oxygen",
        "incoming_hydrogen",
    ]
    abstraction_labels = first["routes"]["h-h2co-h2-hco-1w"]["atom_identity_labels"]
    assert abstraction_labels[2] == "abstracted_formaldehyde_hydrogen"
    assert abstraction_labels[4] == "incoming_hydrogen"
    assert first["campaign_stages"]["branching_common_reference_gate"]["required"]
    assert first["campaign_stages"]["final_freeze"]["required"]
    assert first["accepted_result"] is None


def test_dry_run_refuses_bundle_drift_without_writing_receipt(tmp_path: Path):
    copied = tmp_path / "bundle"
    shutil.copytree(BUNDLE_ROOT, copied)
    target = copied / bundle.ROUTES[0] / "ts.xyz"
    target.write_text(target.read_text() + "\n")
    run_root = tmp_path / "run"

    with pytest.raises(ValueError, match="bundled input drifted"):
        _preflight(run_root, copied)

    assert not (run_root / campaign.PREFLIGHT_RECEIPT).exists()


def test_preflight_rejects_same_element_permutation_despite_rehashed_manifest(
    tmp_path: Path,
):
    copied = tmp_path / "bundle"
    shutil.copytree(BUNDLE_ROOT, copied)
    route = "h-h2co-ch3o-1w"
    target = copied / route / "ts.xyz"
    lines = target.read_text().splitlines()
    assert lines[4].split()[0] == lines[5].split()[0] == "H"
    lines[4], lines[5] = lines[5], lines[4]
    target.write_text("\n".join(lines) + "\n")

    manifest_path = copied / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    file_receipt = manifest["routes"][route]["files"]["ts.xyz"]
    file_receipt["bytes"] = target.stat().st_size
    file_receipt["sha256"] = bundle.sha256_path(target)
    manifest["routes"][route]["canonical_checkpoint_geometry_sha256"]["ts.xyz"] = (
        bundle.geometry_hash_xyz(target)
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    bundle.verify_bundle(copied)

    run_root = tmp_path / "run"
    with pytest.raises(ValueError, match="canonical transition-state geometry"):
        _preflight(run_root, copied)
    assert not (run_root / campaign.PREFLIGHT_RECEIPT).exists()


def test_preflight_receipt_is_non_overwriting(tmp_path: Path):
    run_root = tmp_path / "run"
    first = _preflight(run_root)
    receipt_path = run_root / campaign.PREFLIGHT_RECEIPT
    original_bytes = receipt_path.read_bytes()

    with pytest.raises(FileExistsError, match="already exists"):
        _preflight(run_root)

    assert receipt_path.read_bytes() == original_bytes
    assert json.loads(original_bytes)["identity"] == first["identity"]


def test_preflight_refuses_active_run_claim_then_recovers(tmp_path: Path):
    run_root = tmp_path / "run"
    with campaign._exclusive_run_claim(run_root):
        with pytest.raises(RuntimeError, match="already claimed"):
            _preflight(run_root)
        assert not (run_root / campaign.PREFLIGHT_RECEIPT).exists()

    receipt = _preflight(run_root)
    assert receipt["state"] == "pending"


def test_preflight_refuses_and_preserves_stale_accepted_result(tmp_path: Path):
    run_root = tmp_path / "run"
    run_root.mkdir()
    stale = run_root / "accepted-result.json"
    stale.write_text('{"state":"accepted","identity":"stale"}\n')
    original_bytes = stale.read_bytes()

    with pytest.raises(FileExistsError, match="stale accepted/final result"):
        _preflight(run_root)

    assert stale.read_bytes() == original_bytes
    assert not (run_root / campaign.PREFLIGHT_RECEIPT).exists()

    stale.unlink()
    recovered = _preflight(run_root)
    assert recovered["state"] == "pending"


def test_preflight_refuses_and_preserves_final_freeze_receipt(tmp_path: Path):
    run_root = tmp_path / "run"
    run_root.mkdir()
    stale = run_root / "final-freeze.json"
    stale.write_bytes(b"untrusted final freeze receipt\n")
    original_bytes = stale.read_bytes()

    with pytest.raises(FileExistsError, match="stale accepted/final result"):
        _preflight(run_root)

    assert stale.read_bytes() == original_bytes
    assert not (run_root / campaign.PREFLIGHT_RECEIPT).exists()


def test_preflight_refuses_and_preserves_declared_route_stage_receipt(
    tmp_path: Path,
):
    run_root = tmp_path / "run"
    stale = run_root / bundle.ROUTES[0] / "irc-forward" / "receipt.json"
    stale.parent.mkdir(parents=True)
    stale.write_bytes(b"not JSON and must not be parsed or replaced\n")
    original_bytes = stale.read_bytes()

    with pytest.raises(FileExistsError, match="stale stage receipt"):
        _preflight(run_root)

    assert stale.read_bytes() == original_bytes
    assert not (run_root / campaign.PREFLIGHT_RECEIPT).exists()


def test_cli_refuses_execution_without_dry_run(tmp_path: Path, capsys):
    with pytest.raises(SystemExit) as excinfo:
        campaign.main(
            ["--run-root", str(tmp_path / "run")],
            git_sha=FIXED_GIT_SHA,
            dependency_versions=FIXED_DEPENDENCIES,
        )

    assert excinfo.value.code == 2
    assert "--dry-run is required" in capsys.readouterr().err
    assert not (tmp_path / "run").exists()


def test_cli_enforces_compute_etiquette_bounds(tmp_path: Path, capsys):
    with pytest.raises(SystemExit):
        campaign.main(
            [
                "--dry-run",
                "--run-root",
                str(tmp_path / "run"),
                "--threads",
                "17",
            ],
            git_sha=FIXED_GIT_SHA,
            dependency_versions=FIXED_DEPENDENCIES,
        )
    assert "--threads must be <= 16" in capsys.readouterr().err

    with pytest.raises(SystemExit):
        campaign.main(
            [
                "--dry-run",
                "--run-root",
                str(tmp_path / "run"),
                "--nice",
                "9",
            ],
            git_sha=FIXED_GIT_SHA,
            dependency_versions=FIXED_DEPENDENCIES,
        )
    assert "--nice must be >= 10" in capsys.readouterr().err


def test_git_identity_refuses_dirty_worktree(monkeypatch, tmp_path: Path):
    calls = []

    def fake_run(command, **kwargs):
        calls.append(command)
        return subprocess.CompletedProcess(
            command, 0, stdout="?? untracked.py\n", stderr=""
        )

    monkeypatch.setattr(campaign.subprocess, "run", fake_run)
    with pytest.raises(RuntimeError, match="requires a clean tracked and untracked"):
        campaign._git_sha(tmp_path)

    assert calls == [["git", "status", "--porcelain", "--untracked-files=all"]]


def test_dependency_inventory_fails_closed_without_gpu_backend(monkeypatch):
    available = {name: "1.0" for name in campaign.DEPENDENCY_DISTRIBUTIONS}

    def fake_version(distribution: str) -> str:
        if distribution in available:
            return available[distribution]
        raise campaign.importlib.metadata.PackageNotFoundError(distribution)

    monkeypatch.setattr(campaign.importlib.metadata, "version", fake_version)
    with pytest.raises(RuntimeError, match="GPU4PySCF distribution"):
        campaign._dependency_versions()


def test_injected_dependency_inventory_must_include_exact_gpu_identity(tmp_path: Path):
    incomplete = dict(FIXED_DEPENDENCIES)
    incomplete.pop("cuda_runtime")

    with pytest.raises(ValueError, match="exact and complete"):
        campaign.create_preflight_receipt(
            BUNDLE_ROOT,
            tmp_path / "run",
            git_sha=FIXED_GIT_SHA,
            dependency_versions=incomplete,
        )


def _strict_gate_fixture(
    eigenvalues: np.ndarray,
    *,
    fmax: float = 0.001,
    route_mode_index: int = 0,
):
    from quarry import reaction_path
    from quarry.clusters import Cluster
    from quarry.native_hessian import HARTREE_TO_EV, NativeHessianResult
    from quarry.pipeline import BOHR_TO_ANGSTROM, frequency_geometry_fingerprint

    cluster = Cluster(
        name="strict-gate",
        symbols=["C", "O", "H"],
        coords=np.array([[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [0.1, 0.9, 0.4]]),
    )
    masses = np.array([campaign.ISOTOPIC_MASSES_AMU[s] for s in cluster.symbols])
    rigid = reaction_path._rigid_motion_basis(cluster.coords, masses)
    vibrational = np.linalg.svd(rigid.T, full_matrices=True)[2][6:].T
    mass_weighted = vibrational @ np.diag(eigenvalues) @ vibrational.T
    factors = np.repeat(np.sqrt(masses), 3)
    hessian = factors[:, None] * mass_weighted * factors[None, :]
    gradient = np.zeros_like(cluster.coords)
    gradient[0, 0] = fmax * BOHR_TO_ANGSTROM / HARTREE_TO_EV
    native = NativeHessianResult(
        electronic_hartree=-100.0,
        gradient_hartree_per_bohr=gradient,
        physical_fmax_ev_per_angstrom=fmax,
        cartesian_hessian_hartree_per_bohr2=hessian,
        requested_backend="gpu4pyscf",
        actual_backend="gpu4pyscf",
        gpu_fallback_used=False,
        geometry_fingerprint=frequency_geometry_fingerprint(cluster),
        settings_fingerprint="settings",
    )
    displacement = (
        vibrational[:, route_mode_index] / np.repeat(np.sqrt(masses), 3)
    ).reshape(cluster.coords.shape)
    reactant = cluster.coords - 0.01 * displacement
    product = cluster.coords + 0.01 * displacement
    mapped_path = reaction_path.build_mass_scaled_path(
        np.stack((reactant, cluster.coords, product)),
        masses,
        transition_state_index=1,
    )
    route_vector = (
        mapped_path.mass_scaled_coordinates[2] - mapped_path.mass_scaled_coordinates[0]
    )
    return cluster, masses, native, route_vector, reactant, product


def _validate_strict_gate(fixture):
    cluster, masses, native, reaction_vector, reactant, product = fixture
    return campaign.validate_transition_state_gate(
        cluster,
        masses,
        native,
        expected_settings_fingerprint="settings",
        reaction_vector_mass_scaled=reaction_vector,
        reaction_vector_source="route-fixture:mapped-reactant-product-displacement",
        mapped_reactant_coordinates_angstrom=reactant,
        mapped_product_coordinates_angstrom=product,
    )


def test_strict_transition_state_gate_accepts_one_deep_imaginary_mode():
    receipt = _validate_strict_gate(_strict_gate_fixture(np.array([-0.02, 0.01, 0.03])))

    assert receipt["accepted"] is True
    assert receipt["imaginary_mode_count"] == 1
    assert receipt["imaginary_wavenumber_cm"] >= 200.0
    assert receipt["mapped_reaction_vector_overlap"] == pytest.approx(1.0)
    assert receipt["reaction_vector_source"].startswith("route-fixture:")
    assert len(receipt["reaction_vector_sha256"]) == 64
    assert len(receipt["gradient_sha256"]) == 64
    assert len(receipt["canonical_hessian_sha256"]) == 64


@pytest.mark.parametrize(
    "eigenvalues",
    [
        np.array([0.01, 0.02, 0.03]),
        np.array([-0.02, -1.0e-5, 0.03]),
        np.array([-0.02, -1.0e-10, 0.03]),
    ],
)
def test_strict_transition_state_gate_rejects_zero_second_or_noise_floor_modes(
    eigenvalues: np.ndarray,
):
    with pytest.raises(ValueError, match="exactly one significant negative"):
        _validate_strict_gate(_strict_gate_fixture(eigenvalues))


def test_strict_transition_state_gate_rejects_nonstationary_geometry():
    with pytest.raises(ValueError, match="not stationary"):
        _validate_strict_gate(
            _strict_gate_fixture(np.array([-0.02, 0.01, 0.03]), fmax=0.02)
        )


def test_strict_transition_state_gate_rejects_stale_geometry_settings_and_masses():
    from dataclasses import replace

    cluster, masses, native, reaction_vector, reactant, product = _strict_gate_fixture(
        np.array([-0.02, 0.01, 0.03])
    )
    route_kwargs = {
        "reaction_vector_mass_scaled": reaction_vector,
        "reaction_vector_source": "route-fixture:mapped-reactant-product-displacement",
        "mapped_reactant_coordinates_angstrom": reactant,
        "mapped_product_coordinates_angstrom": product,
    }
    with pytest.raises(ValueError, match="geometry fingerprint"):
        campaign.validate_transition_state_gate(
            cluster,
            masses,
            replace(native, geometry_fingerprint="wrong"),
            expected_settings_fingerprint="settings",
            **route_kwargs,
        )
    with pytest.raises(ValueError, match="settings fingerprint"):
        campaign.validate_transition_state_gate(
            cluster,
            masses,
            native,
            expected_settings_fingerprint="different",
            **route_kwargs,
        )
    with pytest.raises(ValueError, match="isotopic standard"):
        campaign.validate_transition_state_gate(
            cluster,
            masses + 0.001,
            native,
            expected_settings_fingerprint="settings",
            **route_kwargs,
        )


def test_strict_transition_state_gate_rejects_orthogonal_deep_spectator_mode():
    cluster, masses, native, reaction_vector, reactant, product = _strict_gate_fixture(
        np.array([0.01, -0.02, 0.03]), route_mode_index=0
    )
    with pytest.raises(ValueError, match="imaginary mode overlap"):
        campaign.validate_transition_state_gate(
            cluster,
            masses,
            native,
            expected_settings_fingerprint="settings",
            reaction_vector_mass_scaled=reaction_vector,
            reaction_vector_source="route-fixture:mapped-reactant-product-displacement",
            mapped_reactant_coordinates_angstrom=reactant,
            mapped_product_coordinates_angstrom=product,
        )


def test_strict_transition_state_gate_rejects_unbound_vector():
    cluster, masses, native, reaction_vector, reactant, product = _strict_gate_fixture(
        np.array([-0.02, 0.01, 0.03])
    )
    with pytest.raises(ValueError, match="reaction vector source"):
        campaign.validate_transition_state_gate(
            cluster,
            masses,
            native,
            expected_settings_fingerprint="settings",
            reaction_vector_mass_scaled=reaction_vector,
            reaction_vector_source="",
            mapped_reactant_coordinates_angstrom=reactant,
            mapped_product_coordinates_angstrom=product,
        )


def test_strict_transition_state_gate_rejects_vector_not_derived_from_mapped_basins():
    cluster, masses, native, _, reactant, product = _strict_gate_fixture(
        np.array([-0.02, 0.01, 0.03])
    )
    from quarry import reaction_path

    rigid = reaction_path._rigid_motion_basis(cluster.coords, masses)
    unrelated = np.linalg.svd(rigid.T, full_matrices=True)[2][7]
    with pytest.raises(ValueError, match="does not match the receipt-bound"):
        campaign.validate_transition_state_gate(
            cluster,
            masses,
            native,
            expected_settings_fingerprint="settings",
            reaction_vector_mass_scaled=unrelated,
            reaction_vector_source="route-fixture:forged",
            mapped_reactant_coordinates_angstrom=reactant,
            mapped_product_coordinates_angstrom=product,
        )


def test_strict_transition_state_gate_rejects_projection_overflow():
    cluster, masses, native, _, reactant, product = _strict_gate_fixture(
        np.array([-0.02, 0.01, 0.03])
    )
    overflowing = np.array(
        [-1e308, -1e308, -1e308, 1e308, -1e308, -1e308, -1e308, 1e308, 1e308]
    )
    with pytest.raises(ValueError, match="projection must be finite"):
        campaign.validate_transition_state_gate(
            cluster,
            masses,
            native,
            expected_settings_fingerprint="settings",
            reaction_vector_mass_scaled=overflowing,
            reaction_vector_source="route-fixture:overflow",
            mapped_reactant_coordinates_angstrom=reactant,
            mapped_product_coordinates_angstrom=product,
        )
