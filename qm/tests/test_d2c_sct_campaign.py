from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from quarry.clusters import Cluster
from quarry.ts import IrcDirectionPath, IrcPoint, SellaIrcTrace
from scripts import d2c_input_bundle as bundle
from scripts import d2c_sct_campaign as campaign
from scripts.production_energetics import load_xyz_like
from scripts.surface_rate_protocol import reactions

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


def _route_template(route: str):
    return reactions(gpu=True, basis="def2-svp")[route].cluster


def _frozen_endpoint(route: str, filename: str):
    template = _route_template(route)
    return load_xyz_like(
        BUNDLE_ROOT / route / filename,
        template,
        name=f"{route}-{filename}",
    )


def _trace(
    route: str,
    forward_terminal_file: str,
    reverse_terminal_file: str,
    *,
    reverse_ts_offset: float = 0.0,
) -> tuple[Cluster, SellaIrcTrace]:
    transition_state = load_xyz_like(
        BUNDLE_ROOT / route / "ts.xyz",
        _route_template(route),
        name=f"{route}-qualified-ts",
    )
    terminals = {
        "forward": _frozen_endpoint(route, forward_terminal_file).coords,
        "reverse": _frozen_endpoint(route, reverse_terminal_file).coords,
    }
    directions = []
    for direction, sign in (("forward", 1), ("reverse", -1)):
        start = transition_state.coords.copy()
        if direction == "reverse":
            start[0, 0] += reverse_ts_offset
        terminal = terminals[direction]
        points = (
            IrcPoint(0, start, -10.0, 0.01),
            IrcPoint(1, 0.5 * (transition_state.coords + terminal), -10.1, 0.02),
            IrcPoint(2, terminal, -10.2, 0.01),
        )
        directions.append(IrcDirectionPath(direction, sign, points))
    masses = np.array(
        [campaign.ISOTOPIC_MASSES_AMU[s] for s in transition_state.symbols]
    )
    return transition_state, SellaIrcTrace(masses, tuple(directions))


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


@pytest.mark.parametrize(
    ("route", "reactant_file", "product_file"),
    [
        ("h-co-1w-oside", "irc_fwd.xyz", "irc_back.xyz"),
        ("h-co-1w-cside", "irc_back.xyz", "irc_fwd.xyz"),
        ("h-h2co-ch3o-1w", "irc_fwd.xyz", "irc_back.xyz"),
        ("h-h2co-h2-hco-1w", "irc_fwd.xyz", "irc_back.xyz"),
    ],
)
def test_frozen_endpoint_classifier_truth_table(
    route: str, reactant_file: str, product_file: str
):
    reactant = campaign.classify_endpoint_basin(
        route, _frozen_endpoint(route, reactant_file)
    )
    product = campaign.classify_endpoint_basin(
        route, _frozen_endpoint(route, product_file)
    )

    assert reactant.basin == "reactant"
    assert product.basin == "product"
    assert set(reactant.covalent_edges) == set(
        campaign.ENDPOINT_GRAPH_CONTRACTS[route]["reactant"]
    )
    assert set(product.covalent_edges) == set(
        campaign.ENDPOINT_GRAPH_CONTRACTS[route]["product"]
    )
    assert (3, 4) in reactant.covalent_edges or (5, 6) in reactant.covalent_edges
    assert (3, 5) in reactant.covalent_edges or (5, 7) in reactant.covalent_edges


@pytest.mark.parametrize("distance", [1.25, 1.40, 1.55])
def test_endpoint_classifier_rejects_unknown_collision_gray_zone_and_boundaries(
    distance: float,
):
    route = "h-h2co-ch3o-1w"
    endpoint = _frozen_endpoint(route, "irc_fwd.xyz")

    unknown_coords = endpoint.coords.copy()
    unknown_coords[6] += np.array([8.0, 0.0, 0.0])
    with pytest.raises(ValueError, match="unassigned endpoint covalent graph"):
        campaign.classify_endpoint_basin(
            route, replace(endpoint, coords=unknown_coords)
        )

    collision_coords = endpoint.coords.copy()
    collision_coords[4] = collision_coords[0]
    with pytest.raises(ValueError, match="collision floor"):
        campaign.classify_endpoint_basin(
            route, replace(endpoint, coords=collision_coords)
        )

    gray_coords = endpoint.coords.copy()
    gray_coords[4] = gray_coords[0] + np.array([0.0, 0.0, distance])
    with pytest.raises(ValueError, match="gray-zone"):
        campaign.classify_endpoint_basin(route, replace(endpoint, coords=gray_coords))


def test_endpoint_classifier_rejects_ambiguous_graph_contract(monkeypatch):
    route = "h-co-1w-cside"
    endpoint = _frozen_endpoint(route, "irc_back.xyz")
    contract = campaign.ENDPOINT_GRAPH_CONTRACTS[route]
    monkeypatch.setitem(contract, "product", contract["reactant"])

    with pytest.raises(ValueError, match="ambiguous endpoint covalent graph"):
        campaign.classify_endpoint_basin(route, endpoint)


def test_endpoint_classifier_rejects_symbols_state_frozen_and_nonfinite():
    route = "h-co-1w-cside"
    endpoint = _frozen_endpoint(route, "irc_back.xyz")
    bad_coords = endpoint.coords.copy()
    bad_coords[0, 0] = np.inf

    with pytest.raises(ValueError, match="ordered symbols"):
        campaign.classify_endpoint_basin(
            route, replace(endpoint, symbols=["O", *endpoint.symbols[1:]])
        )
    with pytest.raises(ValueError, match="electronic state"):
        campaign.classify_endpoint_basin(route, replace(endpoint, spin=0))
    with pytest.raises(ValueError, match="frozen indices"):
        campaign.classify_endpoint_basin(route, replace(endpoint, frozen_indices=[0]))
    with pytest.raises(ValueError, match="finite"):
        campaign.classify_endpoint_basin(route, replace(endpoint, coords=bad_coords))


@pytest.mark.parametrize(
    ("field", "malformed"),
    [("charge", False), ("spin", True)],
)
def test_boolean_cluster_state_cannot_publish_path_and_corrected_retry_succeeds(
    tmp_path: Path, field: str, malformed: bool
):
    route = "h-co-1w-cside"
    transition_state, trace = _trace(route, "irc_back.xyz", "irc_fwd.xyz")
    if field == "charge":
        bad_transition_state = replace(transition_state, charge=malformed)
        bad_endpoint = replace(
            _frozen_endpoint(route, "irc_back.xyz"), charge=malformed
        )
    else:
        bad_transition_state = replace(transition_state, spin=malformed)
        bad_endpoint = replace(_frozen_endpoint(route, "irc_back.xyz"), spin=malformed)
    run_root = tmp_path / field

    with pytest.raises(ValueError, match="electronic state.*integer"):
        campaign.classify_endpoint_basin(route, bad_endpoint)
    with pytest.raises(ValueError, match="electronic state.*integer"):
        campaign.publish_typed_irc_path(
            run_root,
            campaign_identity="c" * 64,
            route=route,
            atom_mapping_sha256="d" * 64,
            qualified_transition_state=bad_transition_state,
            trace=trace,
        )

    assert not (run_root / route / "path").exists()
    corrected = campaign.publish_typed_irc_path(
        run_root,
        campaign_identity="c" * 64,
        route=route,
        atom_mapping_sha256="d" * 64,
        qualified_transition_state=transition_state,
        trace=trace,
    )
    assert corrected.receipt_path.is_file()


def test_preflight_binds_trusted_typed_endpoint_evidence(tmp_path: Path):
    receipt = _preflight(tmp_path / "run")

    assert receipt["schema"] == "d2c-sct-campaign-preflight-v3"
    assert receipt["campaign"]["endpoint_classification_policy"] == (
        campaign.endpoint_classification_policy_payload()
    )
    for route, route_record in receipt["routes"].items():
        endpoints = route_record["frozen_endpoint_evidence"]
        assert {record["basin"] for record in endpoints.values()} == {
            "reactant",
            "product",
        }
        assert (
            endpoints
            == receipt["campaign"]["routes"][route]["frozen_endpoint_evidence"]
        )
        for filename, endpoint in endpoints.items():
            assert (
                endpoint["trusted_geometry_sha256"]
                == (
                    campaign.TRUSTED_FROZEN_ENDPOINT_EVIDENCE[route][filename][
                        "geometry_sha256"
                    ]
                )
            )
            assert (
                endpoint["trusted_file_sha256"]
                == (
                    campaign.TRUSTED_FROZEN_ENDPOINT_EVIDENCE[route][filename][
                        "file_sha256"
                    ]
                )
            )


@pytest.mark.parametrize(
    ("forward_file", "reverse_file", "reactant_direction", "product_direction"),
    [
        ("irc_fwd.xyz", "irc_back.xyz", "forward", "reverse"),
        ("irc_back.xyz", "irc_fwd.xyz", "reverse", "forward"),
    ],
)
def test_typed_trace_orientation_is_chemical_not_sella_direction(
    forward_file: str,
    reverse_file: str,
    reactant_direction: str,
    product_direction: str,
):
    route = "h-co-1w-oside"
    transition_state, trace = _trace(route, forward_file, reverse_file)

    oriented = campaign.orient_sella_trace(route, transition_state, trace)

    assert oriented.transition_state_index == 2
    assert oriented.coordinates_angstrom.shape == (5, 6, 3)
    assert (
        np.count_nonzero(
            np.all(
                oriented.coordinates_angstrom == transition_state.coords, axis=(1, 2)
            )
        )
        == 1
    )
    assert (
        campaign.classify_endpoint_basin(
            route, replace(transition_state, coords=oriented.coordinates_angstrom[0])
        ).basin
        == "reactant"
    )
    assert (
        campaign.classify_endpoint_basin(
            route, replace(transition_state, coords=oriented.coordinates_angstrom[-1])
        ).basin
        == "product"
    )
    assert oriented.point_provenance[0] == {
        "source_sella_direction": reactant_direction,
        "source_outer_step": 2,
    }
    assert oriented.point_provenance[-1] == {
        "source_sella_direction": product_direction,
        "source_outer_step": 2,
    }
    assert oriented.point_provenance[2] == {
        "source_sella_directions": ["forward", "reverse"],
        "source_outer_steps": [0, 0],
        "transition_state": True,
    }
    assert np.all(np.diff(oriented.mass_scaled_path.coordinate_angstrom) > 0.0)


def test_typed_trace_rejects_same_basin_unknown_terminal_and_ts_disagreement():
    route = "h-co-1w-cside"
    transition_state, same_basin = _trace(route, "irc_back.xyz", "irc_back.xyz")
    with pytest.raises(ValueError, match="one reactant and one product"):
        campaign.orient_sella_trace(route, transition_state, same_basin)

    transition_state, unknown = _trace(route, "irc_back.xyz", "irc_fwd.xyz")
    directions = list(unknown.directions)
    bad_terminal = directions[1].points[-1].coordinates_angstrom.copy()
    bad_terminal[5] += np.array([8.0, 0.0, 0.0])
    bad_points = (
        *directions[1].points[:-1],
        replace(directions[1].points[-1], coordinates_angstrom=bad_terminal),
    )
    directions[1] = replace(directions[1], points=bad_points)
    with pytest.raises(ValueError, match="unassigned endpoint covalent graph"):
        campaign.orient_sella_trace(
            route, transition_state, replace(unknown, directions=tuple(directions))
        )

    transition_state, disagreeing = _trace(
        route, "irc_back.xyz", "irc_fwd.xyz", reverse_ts_offset=1.0e-6
    )
    with pytest.raises(ValueError, match="TS copies"):
        campaign.orient_sella_trace(route, transition_state, disagreeing)


def test_atomic_path_publication_resume_and_tamper_fail_closed(tmp_path: Path):
    route = "h-co-1w-cside"
    transition_state, trace = _trace(route, "irc_back.xyz", "irc_fwd.xyz")
    run_root = tmp_path / "run"
    kwargs = {
        "campaign_identity": "c" * 64,
        "route": route,
        "atom_mapping_sha256": "d" * 64,
        "qualified_transition_state": transition_state,
    }

    published = campaign.publish_typed_irc_path(run_root, trace=trace, **kwargs)
    path_dir = run_root / route / "path"
    assert set(path.name for path in path_dir.iterdir()) == {
        "coordinates.f64",
        "receipt.json",
    }
    assert published.receipt["schema"] == "d2c-typed-irc-path-v1"
    assert published.receipt["stage"] == "typed_irc_path"
    assert len(published.receipt_sha256) == 64

    resumed = campaign.publish_typed_irc_path(run_root, trace=None, **kwargs)
    assert resumed.receipt_sha256 == published.receipt_sha256
    assert np.array_equal(resumed.coordinates_angstrom, published.coordinates_angstrom)

    coordinates_path = path_dir / "coordinates.f64"
    tampered = bytearray(coordinates_path.read_bytes())
    tampered[-1] ^= 1
    coordinates_path.write_bytes(tampered)
    tampered_bytes = coordinates_path.read_bytes()
    with pytest.raises(ValueError, match="coordinates hash"):
        campaign.publish_typed_irc_path(run_root, trace=None, **kwargs)
    assert coordinates_path.read_bytes() == tampered_bytes


def test_atomic_path_failure_injection_leaves_no_canonical_partial(tmp_path: Path):
    route = "h-co-1w-oside"
    transition_state, trace = _trace(route, "irc_fwd.xyz", "irc_back.xyz")
    run_root = tmp_path / "run"

    def fail(stage: str) -> None:
        if stage == "before_path_commit":
            raise RuntimeError("injected path crash")

    with pytest.raises(RuntimeError, match="injected path crash"):
        campaign.publish_typed_irc_path(
            run_root,
            campaign_identity="c" * 64,
            route=route,
            atom_mapping_sha256="d" * 64,
            qualified_transition_state=transition_state,
            trace=trace,
            _failure_injector=fail,
        )

    route_root = run_root / route
    assert not (route_root / "path").exists()
    assert not list(route_root.glob(".path.*.tmp"))


def test_path_precommit_rejects_boolean_outer_step_and_corrected_retry_succeeds(
    tmp_path: Path,
):
    route = "h-co-1w-cside"
    transition_state, trace = _trace(route, "irc_back.xyz", "irc_fwd.xyz")
    directions = list(trace.directions)
    points = list(directions[0].points)
    points[1] = replace(points[1], outer_step=True)
    directions[0] = replace(directions[0], points=tuple(points))
    malformed_trace = replace(trace, directions=tuple(directions))
    run_root = tmp_path / "run"
    publication_kwargs = {
        "campaign_identity": "c" * 64,
        "route": route,
        "atom_mapping_sha256": "d" * 64,
        "qualified_transition_state": transition_state,
    }

    with pytest.raises(ValueError, match="outer_step.*integer"):
        campaign.publish_typed_irc_path(
            run_root, trace=malformed_trace, **publication_kwargs
        )

    assert not (run_root / route / "path").exists()
    corrected = campaign.publish_typed_irc_path(
        run_root, trace=trace, **publication_kwargs
    )
    assert corrected.receipt_path.is_file()


def test_path_resume_rejects_unexpected_file_and_symlink(tmp_path: Path):
    route = "h-co-1w-oside"
    transition_state, trace = _trace(route, "irc_fwd.xyz", "irc_back.xyz")
    kwargs = {
        "campaign_identity": "c" * 64,
        "route": route,
        "atom_mapping_sha256": "d" * 64,
        "qualified_transition_state": transition_state,
    }
    run_root = tmp_path / "run"
    campaign.publish_typed_irc_path(run_root, trace=trace, **kwargs)
    unexpected = run_root / route / "path" / "extra"
    unexpected.write_text("untrusted")
    with pytest.raises(ValueError, match="unexpected path artifact"):
        campaign.publish_typed_irc_path(run_root, trace=None, **kwargs)

    other_root = tmp_path / "other"
    other_root.mkdir()
    symlink_root = tmp_path / "linked"
    symlink_root.symlink_to(other_root, target_is_directory=True)
    with pytest.raises(ValueError, match="symlink"):
        campaign.publish_typed_irc_path(symlink_root, trace=trace, **kwargs)


def test_path_resume_preserves_malformed_receipt(tmp_path: Path):
    route = "h-co-1w-oside"
    transition_state, trace = _trace(route, "irc_fwd.xyz", "irc_back.xyz")
    run_root = tmp_path / "run"
    kwargs = {
        "campaign_identity": "c" * 64,
        "route": route,
        "atom_mapping_sha256": "d" * 64,
        "qualified_transition_state": transition_state,
    }
    campaign.publish_typed_irc_path(run_root, trace=trace, **kwargs)
    receipt_path = run_root / route / "path" / "receipt.json"
    receipt_path.write_bytes(b"not-json\n")
    original = receipt_path.read_bytes()

    with pytest.raises(ValueError, match="malformed typed path receipt"):
        campaign.publish_typed_irc_path(run_root, trace=None, **kwargs)
    assert receipt_path.read_bytes() == original


def _published_path(run_root: Path):
    route = "h-co-1w-cside"
    transition_state, trace = _trace(route, "irc_back.xyz", "irc_fwd.xyz")
    path = campaign.publish_typed_irc_path(
        run_root,
        campaign_identity="c" * 64,
        route=route,
        atom_mapping_sha256="d" * 64,
        qualified_transition_state=transition_state,
        trace=trace,
    )
    return route, transition_state, path


def _native_result(
    cluster, *, requested="gpu4pyscf", actual="gpu4pyscf", fallback=False
):
    point_index = int(cluster.name.rsplit("-", 1)[-1])
    gradient = np.zeros_like(cluster.coords)
    gradient[0, 0] = (point_index + 1) * 1.0e-5
    fmax = float(np.max(np.linalg.norm(gradient, axis=1)))
    fmax *= campaign.HARTREE_TO_EV / campaign.BOHR_TO_ANGSTROM
    return campaign.NativeHessianResult(
        electronic_hartree=-100.0 - point_index,
        gradient_hartree_per_bohr=gradient,
        physical_fmax_ev_per_angstrom=fmax,
        cartesian_hessian_hartree_per_bohr2=np.eye(3 * len(cluster.symbols)) * 0.01,
        requested_backend=requested,
        actual_backend=actual,
        gpu_fallback_used=fallback,
        geometry_fingerprint=campaign.frequency_geometry_fingerprint(cluster),
        settings_fingerprint="settings-v1",
    )


def _publish_hessians(run_root: Path, path, transition_state, evaluator, **kwargs):
    return campaign.publish_path_hessians(
        run_root,
        campaign_identity="c" * 64,
        route="h-co-1w-cside",
        atom_mapping_sha256="d" * 64,
        qualified_transition_state=transition_state,
        path=path,
        settings_fingerprint="settings-v1",
        backend_policy={
            "requested_backend": "gpu4pyscf",
            "allow_cpu_fallback": True,
        },
        evaluator=evaluator,
        **kwargs,
    )


def _rewrite_canonical_json(path: Path, mutate) -> None:
    payload = json.loads(path.read_text())
    mutate(payload)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def test_stale_private_path_temporary_is_removed_not_promoted(tmp_path: Path):
    route = "h-co-1w-cside"
    transition_state, trace = _trace(route, "irc_back.xyz", "irc_fwd.xyz")
    run_root = tmp_path / "run"
    publication_kwargs = {
        "campaign_identity": "c" * 64,
        "route": route,
        "atom_mapping_sha256": "d" * 64,
        "qualified_transition_state": transition_state,
    }
    original = campaign.publish_typed_irc_path(
        run_root, trace=trace, **publication_kwargs
    )
    route_root = run_root / route
    temporary = route_root / ".path.123.456.tmp"
    temporary.mkdir()
    (temporary / "coordinates.f64").write_bytes(b"private coordinates")
    (temporary / "receipt.json").write_bytes(b"private receipt\n")

    resumed = campaign.publish_typed_irc_path(
        run_root, trace=None, **publication_kwargs
    )

    assert not temporary.exists()
    assert resumed.receipt_sha256 == original.receipt_sha256
    assert resumed.receipt_path.read_bytes() != b"private receipt\n"


def test_stale_private_hessian_temporaries_do_not_brick_resume(tmp_path: Path):
    run_root = tmp_path / "run"
    _, transition_state, path = _published_path(run_root)
    hessian_root = run_root / "h-co-1w-cside" / "hessians"
    points_root = hessian_root / "points"
    points_root.mkdir(parents=True)
    point_temporary = points_root / ".000000.123.456.tmp"
    point_temporary.mkdir()
    (point_temporary / "gradient.f64").write_bytes(b"private gradient")
    (point_temporary / "hessian.f64").write_bytes(b"private Hessian")
    (point_temporary / "receipt.json").write_bytes(b"private point receipt\n")
    aggregate_temporary = hessian_root / ".receipt.json.123.456.tmp"
    aggregate_temporary.write_bytes(b"private aggregate receipt\n")

    published = _publish_hessians(run_root, path, transition_state, _native_result)

    assert not point_temporary.exists()
    assert not aggregate_temporary.exists()
    assert published.receipt["accepted"] is True
    assert published.receipt_path.read_bytes() != b"private aggregate receipt\n"


def test_path_publication_destination_race_never_clobbers(monkeypatch, tmp_path: Path):
    route = "h-co-1w-cside"
    transition_state, trace = _trace(route, "irc_back.xyz", "irc_fwd.xyz")
    run_root = tmp_path / "run"
    canonical = run_root / route / "path"
    original_rename = campaign._renameat2_noreplace

    def create_destination_then_rename(source: Path, destination: Path) -> None:
        destination.mkdir()
        (destination / "racer-marker").write_bytes(b"racer owns destination")
        original_rename(source, destination)

    monkeypatch.setattr(
        campaign, "_renameat2_noreplace", create_destination_then_rename
    )
    with pytest.raises(FileExistsError):
        campaign.publish_typed_irc_path(
            run_root,
            campaign_identity="c" * 64,
            route=route,
            atom_mapping_sha256="d" * 64,
            qualified_transition_state=transition_state,
            trace=trace,
        )

    assert (canonical / "racer-marker").read_bytes() == b"racer owns destination"
    assert not list((run_root / route).glob(".path.*.tmp"))


def test_hessian_point_destination_race_never_clobbers(monkeypatch, tmp_path: Path):
    run_root = tmp_path / "run"
    _, transition_state, path = _published_path(run_root)
    point_dir = run_root / "h-co-1w-cside" / "hessians" / "points" / "000000"
    original_rename = campaign._renameat2_noreplace

    def create_destination_then_rename(source: Path, destination: Path) -> None:
        destination.mkdir()
        (destination / "racer-marker").write_bytes(b"racer owns point")
        original_rename(source, destination)

    monkeypatch.setattr(
        campaign, "_renameat2_noreplace", create_destination_then_rename
    )
    with pytest.raises(FileExistsError):
        _publish_hessians(run_root, path, transition_state, _native_result)

    assert (point_dir / "racer-marker").read_bytes() == b"racer owns point"
    assert not list(point_dir.parent.glob(".000000.*.tmp"))


def test_hessian_aggregate_destination_race_never_clobbers(monkeypatch, tmp_path: Path):
    run_root = tmp_path / "run"
    _, transition_state, path = _published_path(run_root)
    aggregate = run_root / "h-co-1w-cside" / "hessians" / "receipt.json"
    original_rename = campaign._renameat2_noreplace

    def create_destination_then_rename(source: Path, destination: Path) -> None:
        if destination == aggregate:
            destination.write_bytes(b"racer owns aggregate\n")
        original_rename(source, destination)

    monkeypatch.setattr(
        campaign, "_renameat2_noreplace", create_destination_then_rename
    )
    with pytest.raises(FileExistsError):
        _publish_hessians(run_root, path, transition_state, _native_result)

    assert aggregate.read_bytes() == b"racer owns aggregate\n"
    assert not list(aggregate.parent.glob(".receipt.json.*.tmp"))


def test_regular_file_publication_fallback_is_atomic_and_non_overwriting(
    monkeypatch, tmp_path: Path
):
    def unavailable(_source: Path, _destination: Path) -> None:
        raise NotImplementedError

    monkeypatch.setattr(campaign, "_renameat2_noreplace", unavailable)
    source = tmp_path / "temporary"
    destination = tmp_path / "receipt.json"
    source.write_bytes(b"new receipt")
    campaign._publish_noreplace(source, destination, source_is_directory=False)
    assert destination.read_bytes() == b"new receipt"
    assert not source.exists()

    source.write_bytes(b"must not overwrite")
    with pytest.raises(FileExistsError):
        campaign._publish_noreplace(source, destination, source_is_directory=False)
    assert destination.read_bytes() == b"new receipt"
    assert source.read_bytes() == b"must not overwrite"


def test_directory_publication_fails_closed_without_renameat2(
    monkeypatch, tmp_path: Path
):
    def unavailable(_source: Path, _destination: Path) -> None:
        raise NotImplementedError

    monkeypatch.setattr(campaign, "_renameat2_noreplace", unavailable)
    source = tmp_path / "temporary"
    source.mkdir()

    with pytest.raises(RuntimeError, match="requires renameat2"):
        campaign._publish_noreplace(
            source, tmp_path / "canonical", source_is_directory=True
        )
    assert source.is_dir()
    assert not (tmp_path / "canonical").exists()


@pytest.mark.parametrize("kind", ["path-symlink", "point-file"])
def test_malformed_private_temporary_is_preserved_and_rejected(
    tmp_path: Path, kind: str
):
    run_root = tmp_path / "run"
    route = "h-co-1w-cside"
    transition_state, trace = _trace(route, "irc_back.xyz", "irc_fwd.xyz")
    if kind == "path-symlink":
        route_root = run_root / route
        route_root.mkdir(parents=True)
        outside = tmp_path / "outside"
        outside.mkdir()
        temporary = route_root / ".path.123.456.tmp"
        temporary.symlink_to(outside, target_is_directory=True)

        def call():
            return campaign.publish_typed_irc_path(
                run_root,
                campaign_identity="c" * 64,
                route=route,
                atom_mapping_sha256="d" * 64,
                qualified_transition_state=transition_state,
                trace=trace,
            )
    else:
        path = campaign.publish_typed_irc_path(
            run_root,
            campaign_identity="c" * 64,
            route=route,
            atom_mapping_sha256="d" * 64,
            qualified_transition_state=transition_state,
            trace=trace,
        )
        points_root = run_root / route / "hessians" / "points"
        points_root.mkdir(parents=True)
        temporary = points_root / ".000000.123.456.tmp"
        temporary.write_bytes(b"not a directory")

        def call():
            return _publish_hessians(run_root, path, transition_state, _native_result)

    with pytest.raises(ValueError, match="owned temporary must be a real confined"):
        call()
    assert temporary.exists() or temporary.is_symlink()


@pytest.mark.parametrize(
    "mutate",
    [
        lambda receipt: receipt["points"][0].__setitem__("unexpected", "field"),
        lambda receipt: receipt.__setitem__("accepted", 1),
        lambda receipt: receipt.__setitem__("charge", False),
        lambda receipt: receipt.__setitem__("spin_2s", True),
        lambda receipt: receipt.__setitem__("point_count", 5.0),
        lambda receipt: receipt["points"][0].__setitem__("index", False),
        lambda receipt: receipt["points"][1]["source"].__setitem__(
            "source_outer_step", True
        ),
    ],
    ids=[
        "extra-point-key",
        "int-for-bool",
        "bool-for-charge-int",
        "bool-for-spin-int",
        "float-for-count-int",
        "bool-for-index-int",
        "bool-for-provenance-int",
    ],
)
def test_path_receipt_rejects_extra_keys_and_primitive_type_substitution(
    tmp_path: Path, mutate
):
    run_root = tmp_path / "run"
    route, transition_state, _ = _published_path(run_root)
    receipt_path = run_root / route / "path" / "receipt.json"
    _rewrite_canonical_json(receipt_path, mutate)
    malformed = receipt_path.read_bytes()

    with pytest.raises(ValueError):
        campaign.publish_typed_irc_path(
            run_root,
            campaign_identity="c" * 64,
            route=route,
            atom_mapping_sha256="d" * 64,
            qualified_transition_state=transition_state,
            trace=None,
        )
    assert receipt_path.read_bytes() == malformed


@pytest.mark.parametrize(
    "mutate",
    [
        lambda receipt: receipt.__setitem__("accepted", 1),
        lambda receipt: receipt.__setitem__("spin_2s", True),
        lambda receipt: receipt.__setitem__("point_index", False),
        lambda receipt: receipt.__setitem__("point_count", 5.0),
        lambda receipt: receipt.__setitem__("is_transition_state", 0),
        lambda receipt: receipt.__setitem__("gpu_fallback_used", 0),
        lambda receipt: receipt["path_point_provenance"].__setitem__(
            "source_outer_step", 2.0
        ),
        lambda receipt: receipt["path_point_provenance"].__setitem__(
            "unexpected", "field"
        ),
    ],
    ids=[
        "int-for-accepted-bool",
        "bool-for-spin-int",
        "bool-for-point-index",
        "float-for-point-count",
        "int-for-ts-bool",
        "int-for-fallback-bool",
        "float-for-provenance-index",
        "extra-provenance-key",
    ],
)
def test_hessian_point_receipt_rejects_identity_and_provenance_type_substitution(
    tmp_path: Path, mutate
):
    run_root = tmp_path / "run"
    _, transition_state, path = _published_path(run_root)
    _publish_hessians(run_root, path, transition_state, _native_result)
    receipt_path = (
        run_root / "h-co-1w-cside" / "hessians" / "points" / "000000" / "receipt.json"
    )
    _rewrite_canonical_json(receipt_path, mutate)
    malformed = receipt_path.read_bytes()

    with pytest.raises(ValueError):
        _publish_hessians(
            run_root,
            path,
            transition_state,
            lambda _cluster: pytest.fail("malformed receipt called evaluator"),
        )
    assert receipt_path.read_bytes() == malformed


def test_aggregate_receipt_rejects_child_identity_bool_for_int(tmp_path: Path):
    run_root = tmp_path / "run"
    _, transition_state, path = _published_path(run_root)
    _publish_hessians(run_root, path, transition_state, _native_result)
    receipt_path = run_root / "h-co-1w-cside" / "hessians" / "receipt.json"
    _rewrite_canonical_json(
        receipt_path,
        lambda receipt: receipt["ordered_child_hash_chain"][0].__setitem__(
            "point_index", False
        ),
    )

    with pytest.raises(ValueError):
        _publish_hessians(
            run_root,
            path,
            transition_state,
            lambda _cluster: pytest.fail("malformed aggregate called evaluator"),
        )


def test_run_root_and_route_path_escapes_are_rejected(tmp_path: Path):
    route = "h-co-1w-cside"
    transition_state, trace = _trace(route, "irc_back.xyz", "irc_fwd.xyz")
    kwargs = {
        "campaign_identity": "c" * 64,
        "atom_mapping_sha256": "d" * 64,
        "qualified_transition_state": transition_state,
        "trace": trace,
    }

    with pytest.raises(ValueError, match="path escape"):
        campaign.publish_typed_irc_path(
            tmp_path / "safe" / ".." / "escaped", route=route, **kwargs
        )
    with pytest.raises(ValueError, match="unsupported typed endpoint route"):
        campaign.publish_typed_irc_path(tmp_path / "run", route="../escaped", **kwargs)
    assert not (tmp_path / "escaped").exists()


@pytest.mark.parametrize("artifact", ["root-file", "point-symlink", "point-file"])
def test_hessian_tree_rejects_symlink_and_unexpected_artifacts(
    tmp_path: Path, artifact: str
):
    run_root = tmp_path / "run"
    _, transition_state, path = _published_path(run_root)
    hessian_root = run_root / "h-co-1w-cside" / "hessians"
    points_root = hessian_root / "points"
    points_root.mkdir(parents=True)
    if artifact == "root-file":
        (hessian_root / "unexpected").write_bytes(b"untrusted")
        match = "unexpected Hessian stage artifact"
    elif artifact == "point-symlink":
        outside = tmp_path / "outside"
        outside.mkdir()
        (points_root / "000000").symlink_to(outside, target_is_directory=True)
        match = "real directory"
    else:
        (points_root / "unexpected").write_bytes(b"untrusted")
        match = "unexpected Hessian point artifact"

    with pytest.raises(ValueError, match=match):
        _publish_hessians(
            run_root,
            path,
            transition_state,
            lambda _cluster: pytest.fail("invalid tree called evaluator"),
        )


def test_hessian_data_symlink_is_rejected_without_following_it(tmp_path: Path):
    run_root = tmp_path / "run"
    _, transition_state, path = _published_path(run_root)
    _publish_hessians(run_root, path, transition_state, _native_result)
    hessian_path = (
        run_root / "h-co-1w-cside" / "hessians" / "points" / "000000" / "hessian.f64"
    )
    outside = tmp_path / "outside-hessian.f64"
    outside.write_bytes(b"untrusted outside data")
    hessian_path.unlink()
    hessian_path.symlink_to(outside)

    with pytest.raises(ValueError, match="cartesian_hessian must be a regular file"):
        _publish_hessians(
            run_root,
            path,
            transition_state,
            lambda _cluster: pytest.fail("symlinked Hessian called evaluator"),
        )
    assert outside.read_bytes() == b"untrusted outside data"


def test_hessian_point_failure_resumes_only_remaining_then_full_resume(tmp_path: Path):
    run_root = tmp_path / "run"
    _, transition_state, path = _published_path(run_root)
    first_calls = []

    def fail_at_two(cluster):
        point_index = int(cluster.name.rsplit("-", 1)[-1])
        first_calls.append(point_index)
        if point_index == 2:
            raise RuntimeError("injected evaluator failure")
        return _native_result(cluster)

    with pytest.raises(RuntimeError, match="injected evaluator failure"):
        _publish_hessians(run_root, path, transition_state, fail_at_two)
    points_root = run_root / "h-co-1w-cside" / "hessians" / "points"
    assert first_calls == [0, 1, 2]
    assert sorted(child.name for child in points_root.iterdir()) == ["000000", "000001"]
    assert not (run_root / "h-co-1w-cside" / "hessians" / "receipt.json").exists()

    resumed_calls = []

    def resume_evaluator(cluster):
        point_index = int(cluster.name.rsplit("-", 1)[-1])
        resumed_calls.append(point_index)
        return _native_result(cluster)

    completed = _publish_hessians(run_root, path, transition_state, resume_evaluator)
    assert resumed_calls == [2, 3, 4]
    assert completed.receipt["point_count"] == 5
    assert len(completed.receipt["ordered_child_hash_chain"]) == 5
    assert len(completed.receipt_sha256) == 64

    def forbidden(_cluster):
        raise AssertionError("full aggregate resume called evaluator")

    fully_resumed = _publish_hessians(run_root, path, transition_state, forbidden)
    assert fully_resumed.receipt_sha256 == completed.receipt_sha256


def test_hessian_atomic_point_failure_has_no_canonical_partial(tmp_path: Path):
    run_root = tmp_path / "run"
    _, transition_state, path = _published_path(run_root)

    def fail(stage: str) -> None:
        if stage == "before_hessian_point_commit:2":
            raise RuntimeError("injected point commit crash")

    with pytest.raises(RuntimeError, match="injected point commit crash"):
        _publish_hessians(
            run_root,
            path,
            transition_state,
            _native_result,
            _failure_injector=fail,
        )

    points_root = run_root / "h-co-1w-cside" / "hessians" / "points"
    assert sorted(child.name for child in points_root.iterdir()) == ["000000", "000001"]
    assert not list(points_root.glob(".000002.*.tmp"))


def test_hessian_aggregate_failure_resumes_without_evaluator(tmp_path: Path):
    run_root = tmp_path / "run"
    _, transition_state, path = _published_path(run_root)

    def fail(stage: str) -> None:
        if stage == "before_hessian_aggregate_commit":
            raise RuntimeError("injected aggregate crash")

    with pytest.raises(RuntimeError, match="injected aggregate crash"):
        _publish_hessians(
            run_root,
            path,
            transition_state,
            _native_result,
            _failure_injector=fail,
        )
    hessian_root = run_root / "h-co-1w-cside" / "hessians"
    assert len(list((hessian_root / "points").iterdir())) == 5
    assert not (hessian_root / "receipt.json").exists()

    completed = _publish_hessians(
        run_root,
        path,
        transition_state,
        lambda _cluster: pytest.fail("aggregate-only resume called evaluator"),
    )
    assert completed.receipt["accepted"] is True


def test_hessian_aggregate_is_strictly_validated_before_commit(
    monkeypatch, tmp_path: Path
):
    run_root = tmp_path / "run"
    _, transition_state, path = _published_path(run_root)
    original_json_bytes = campaign._json_bytes

    def poison_aggregate(payload):
        if payload.get("schema") == "d2c-native-path-hessians-v1":
            payload = {**payload, "accepted": 1}
        return original_json_bytes(payload)

    monkeypatch.setattr(campaign, "_json_bytes", poison_aggregate)
    with pytest.raises(ValueError, match="aggregate receipt"):
        _publish_hessians(run_root, path, transition_state, _native_result)

    aggregate_path = run_root / "h-co-1w-cside" / "hessians" / "receipt.json"
    assert not aggregate_path.exists()
    monkeypatch.setattr(campaign, "_json_bytes", original_json_bytes)
    corrected = _publish_hessians(
        run_root,
        path,
        transition_state,
        lambda _cluster: pytest.fail("corrected aggregate retry called evaluator"),
    )
    assert corrected.receipt_path == aggregate_path


def test_hessian_child_and_aggregate_tampering_block_without_evaluator(
    tmp_path: Path,
):
    run_root = tmp_path / "child"
    _, transition_state, path = _published_path(run_root)
    _publish_hessians(run_root, path, transition_state, _native_result)
    gradient_path = (
        run_root / "h-co-1w-cside" / "hessians" / "points" / "000003" / "gradient.f64"
    )
    tampered = bytearray(gradient_path.read_bytes())
    tampered[0] ^= 1
    gradient_path.write_bytes(tampered)

    with pytest.raises(ValueError, match="gradient hash"):
        _publish_hessians(
            run_root,
            path,
            transition_state,
            lambda _cluster: pytest.fail("tampered child called evaluator"),
        )

    aggregate_root = tmp_path / "aggregate"
    _, aggregate_ts, aggregate_path = _published_path(aggregate_root)
    _publish_hessians(aggregate_root, aggregate_path, aggregate_ts, _native_result)
    aggregate_receipt = aggregate_root / "h-co-1w-cside" / "hessians" / "receipt.json"
    payload = json.loads(aggregate_receipt.read_text())
    payload["ordered_child_hash_chain"][-1]["chain_sha256"] = "0" * 64
    aggregate_receipt.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    original = aggregate_receipt.read_bytes()
    with pytest.raises(ValueError, match="child hash chain"):
        _publish_hessians(
            aggregate_root,
            aggregate_path,
            aggregate_ts,
            lambda _cluster: pytest.fail("tampered aggregate called evaluator"),
        )
    assert aggregate_receipt.read_bytes() == original


@pytest.mark.parametrize("failure", ["nan", "inf", "shape", "provenance", "geometry"])
def test_hessian_invalid_result_never_publishes_point(tmp_path: Path, failure: str):
    run_root = tmp_path / failure
    _, transition_state, path = _published_path(run_root)

    def invalid(cluster):
        if failure == "provenance":
            return _native_result(cluster, requested="pyscf", actual="pyscf")
        result = _native_result(cluster)
        if failure == "geometry":
            return replace(result, geometry_fingerprint="wrong")
        object.__setattr__(
            result,
            "gradient_hartree_per_bohr",
            (
                np.full(
                    (len(cluster.symbols), 3),
                    np.nan if failure == "nan" else np.inf,
                )
                if failure in {"nan", "inf"}
                else np.zeros((1, 3))
            ),
        )
        return result

    match = {
        "nan": "non-finite",
        "inf": "non-finite",
        "shape": "shape",
        "provenance": "backend policy",
        "geometry": "geometry fingerprint",
    }[failure]
    with pytest.raises(ValueError, match=match):
        _publish_hessians(run_root, path, transition_state, invalid)
    points_root = run_root / "h-co-1w-cside" / "hessians" / "points"
    assert not list(points_root.iterdir())


def test_hessian_precommit_rejects_int_fallback_and_corrected_retry_succeeds(
    tmp_path: Path,
):
    run_root = tmp_path / "run"
    _, transition_state, path = _published_path(run_root)

    def malformed_fallback(cluster):
        result = _native_result(cluster)
        object.__setattr__(result, "gpu_fallback_used", 0)
        return result

    with pytest.raises(ValueError, match="gpu_fallback_used.*boolean"):
        _publish_hessians(run_root, path, transition_state, malformed_fallback)

    point_zero = run_root / "h-co-1w-cside" / "hessians" / "points" / "000000"
    assert not point_zero.exists()
    corrected = _publish_hessians(run_root, path, transition_state, _native_result)
    assert corrected.receipt_path.is_file()


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
