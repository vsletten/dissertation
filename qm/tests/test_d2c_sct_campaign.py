from __future__ import annotations

import hashlib
import inspect
import json
import mmap
import multiprocessing
import os
import shutil
import subprocess
import sys
from contextlib import contextmanager
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from quarry.clusters import Cluster
from quarry.ts import IrcDirectionPath, IrcExecutionContract, IrcPoint, SellaIrcTrace
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
FIXED_MODULE_MANIFEST = {
    name: {
        "origin": f"/opt/d2c-test-modules/{name.replace('.', '/')}.py",
        "sha256": hashlib.sha256(name.encode()).hexdigest(),
        "byte_count": len(name),
        "trust_class": trust_class,
        "execution_identity": {
            "kind": "python-source",
            "loaded_module_code_sha256": hashlib.sha256(
                f"module:{name}".encode()
            ).hexdigest(),
            "loaded_code_sha256": hashlib.sha256(f"loaded:{name}".encode()).hexdigest(),
            "loaded_state_sha256": hashlib.sha256(f"state:{name}".encode()).hexdigest(),
            "source_files": {
                f"/opt/d2c-test-modules/{name.replace('.', '/')}.py": hashlib.sha256(
                    name.encode()
                ).hexdigest()
            },
        },
    }
    for name, trust_class in campaign.EXECUTABLE_MODULES.items()
}
FIXED_NATIVE_PAYLOAD_MANIFEST = {
    f"lib/python3.13/site-packages/pyscf/lib/{name}": {
        "origin": f"/opt/d2c-test-native/pyscf/lib/{name}",
        "sha256": hashlib.sha256(name.encode()).hexdigest(),
        "byte_count": len(name),
        "mapped_device": 1,
        "mapped_inode": index,
    }
    for index, name in enumerate(
        sorted(campaign.REQUIRED_NATIVE_PAYLOAD_BASENAMES), start=1
    )
}
BUNDLE_ROOT = Path(__file__).parents[1] / "data" / "D2c-instanton-tier" / "d2b-inputs"


def _hold_route_claim(route_root, entered, release, outcome):
    try:
        with campaign._exclusive_route_claim(Path(route_root)):
            entered.set()
            if not release.wait(10):
                raise TimeoutError("test did not release route claim")
    except Exception as exc:
        outcome.put(("failed", type(exc).__name__, str(exc)))
    else:
        outcome.put(("released", "", ""))


def _attempt_route_claim(route_root, outcome):
    try:
        with campaign._exclusive_route_claim(Path(route_root)):
            outcome.put(("entered", "", ""))
    except Exception as exc:
        outcome.put(("blocked", type(exc).__name__, str(exc)))


@pytest.fixture(autouse=True)
def _fixed_production_runtime_identity(monkeypatch):
    """Keep production-boundary tests deterministic without weakening wrappers."""

    monkeypatch.setattr(
        campaign,
        "_current_code_dependency_identity",
        lambda: {
            "git_sha": FIXED_GIT_SHA,
            "dependencies": dict(FIXED_DEPENDENCIES),
            "executable_modules": {
                name: dict(record) for name, record in FIXED_MODULE_MANIFEST.items()
            },
            "native_payloads": {
                name: dict(record)
                for name, record in FIXED_NATIVE_PAYLOAD_MANIFEST.items()
            },
            "python": campaign.platform.python_version(),
        },
    )


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
    initialization_fingerprint = hashlib.sha256(
        f"test-sella-initialization:{route}".encode()
    ).hexdigest()
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
        directions.append(
            IrcDirectionPath(
                direction,
                sign,
                points,
                initialization_fingerprint=initialization_fingerprint,
            )
        )
    masses = np.array(
        [campaign.ISOTOPIC_MASSES_AMU[s] for s in transition_state.symbols]
    )
    return transition_state, SellaIrcTrace(
        masses,
        tuple(directions),
        IrcExecutionContract(
            algorithm="sella-gonzalez-schlegel",
            step_size_angstrom=0.05,
            maximum_steps=200,
            outer_fmax_ev_per_angstrom=0.05,
            inner_fmax_ev_per_angstrom=0.01,
        ),
        initialization_fingerprint,
    )


def _sella_initialization_state(
    transition_state: Cluster,
    trace: SellaIrcTrace,
) -> campaign.quarry_ts.SellaIrcInitializationState:
    assert trace.execution_contract is not None
    dimension = transition_state.coords.size
    x0 = transition_state.coords.reshape(-1).copy()
    masses = np.asarray(
        [campaign.ISOTOPIC_MASSES_AMU[symbol] for symbol in transition_state.symbols]
    )
    settings = campaign.DftSettings(**campaign.DFT_SETTINGS)
    h0 = np.eye(dimension)
    h0[0, 0] = -1.0
    v0ts = np.zeros(dimension)
    v0ts[0] = trace.execution_contract.step_size_angstrom / np.sqrt(masses[0])
    identity = np.eye(dimension)
    pes_current = {
        "L": np.zeros(0),
        "Ucons": np.zeros((dimension, 0)),
        "Ufree": identity,
        "Unred": identity,
        "drdx": np.zeros((0, dimension)),
        "f": -10.0,
        "g": np.zeros(dimension),
        "state_hash": x0.tobytes(),
        "x": x0,
    }
    return campaign.quarry_ts.SellaIrcInitializationState(
        symbols=tuple(transition_state.symbols),
        charge=transition_state.charge,
        spin=transition_state.spin,
        frozen_indices=tuple(sorted(transition_state.frozen_indices)),
        settings_fingerprint=hashlib.sha256(
            campaign.frequency_settings_fingerprint(settings).encode()
        ).hexdigest(),
        execution_contract=trace.execution_contract,
        x0=x0,
        masses_amu=masses,
        h0=h0,
        v0ts=v0ts,
        pes_current=pes_current,
        pes_last={"x": None, "f": None, "g": None},
    )


def _bind_trace_initialization(
    trace: SellaIrcTrace,
    initialization: campaign.quarry_ts.SellaIrcInitializationState,
) -> SellaIrcTrace:
    directions = tuple(
        replace(direction, initialization_fingerprint=initialization.fingerprint)
        for direction in trace.directions
    )
    return replace(
        trace,
        directions=directions,
        initialization_fingerprint=initialization.fingerprint,
    )


def _checkpointing_trace_runner(transition_state: Cluster, trace: SellaIrcTrace):
    initialization = _sella_initialization_state(transition_state, trace)
    bound_trace = _bind_trace_initialization(trace, initialization)

    def runner(*_args, _initialization_callback, **_kwargs):
        _initialization_callback(initialization)
        return bound_trace

    return bound_trace, runner


def _preflight(run_root: Path, bundle_root: Path = BUNDLE_ROOT) -> dict:
    return campaign.create_preflight_receipt(
        bundle_root,
        run_root,
        git_sha=FIXED_GIT_SHA,
        dependency_versions=FIXED_DEPENDENCIES,
        executable_module_manifest=FIXED_MODULE_MANIFEST,
        native_payload_manifest=FIXED_NATIVE_PAYLOAD_MANIFEST,
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
        executable_module_manifest=FIXED_MODULE_MANIFEST,
        native_payload_manifest=FIXED_NATIVE_PAYLOAD_MANIFEST,
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
        "electronic_energy_reproduction_absolute_tolerance_ev": 1.0e-6,
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
        campaign._publish_typed_irc_path(
            run_root,
            campaign_identity="c" * 64,
            route=route,
            atom_mapping_sha256="d" * 64,
            qualified_transition_state=bad_transition_state,
            trace=trace,
        )

    assert not (run_root / route / "path").exists()
    corrected = campaign._publish_typed_irc_path(
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

    assert receipt["schema"] == "d2c-sct-campaign-preflight-v6"
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

    transition_state, extra_direction = _trace(route, "irc_back.xyz", "irc_fwd.xyz")
    poisoned_third = replace(
        extra_direction.directions[1],
        points=(
            replace(
                extra_direction.directions[1].points[0],
                electronic_energy_ev=(
                    extra_direction.directions[1].points[0].electronic_energy_ev + 7.0
                ),
            ),
            *extra_direction.directions[1].points[1:],
        ),
    )
    poisoned_trace = object.__new__(SellaIrcTrace)
    object.__setattr__(poisoned_trace, "masses_amu", extra_direction.masses_amu)
    object.__setattr__(
        poisoned_trace,
        "directions",
        (*extra_direction.directions, poisoned_third),
    )
    with pytest.raises(ValueError, match="exactly two directions"):
        campaign.orient_sella_trace(route, transition_state, poisoned_trace)


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

    published = campaign._publish_typed_irc_path(run_root, trace=trace, **kwargs)
    path_dir = run_root / route / "path"
    assert set(path.name for path in path_dir.iterdir()) == {
        "coordinates.f64",
        "receipt.json",
    }
    assert published.receipt["schema"] == "d2c-typed-irc-path-v1"
    assert published.receipt["stage"] == "typed_irc_path"
    assert len(published.receipt_sha256) == 64

    resumed = campaign._publish_typed_irc_path(run_root, trace=None, **kwargs)
    assert resumed.receipt_sha256 == published.receipt_sha256
    assert np.array_equal(resumed.coordinates_angstrom, published.coordinates_angstrom)

    coordinates_path = path_dir / "coordinates.f64"
    tampered = bytearray(coordinates_path.read_bytes())
    tampered[-1] ^= 1
    coordinates_path.write_bytes(tampered)
    tampered_bytes = coordinates_path.read_bytes()
    with pytest.raises(ValueError, match="coordinates hash"):
        campaign._publish_typed_irc_path(run_root, trace=None, **kwargs)
    assert coordinates_path.read_bytes() == tampered_bytes


def test_v3_typed_path_requires_exact_five_ancestor_hashes(tmp_path: Path):
    route = "h-co-1w-cside"
    transition_state, trace = _trace(route, "irc_back.xyz", "irc_fwd.xyz")
    run_root = tmp_path / "run"

    with pytest.raises(ValueError, match="ancestor receipts are incomplete"):
        campaign._publish_typed_irc_path(
            run_root,
            campaign_identity="c" * 64,
            route=route,
            atom_mapping_sha256="d" * 64,
            qualified_transition_state=transition_state,
            trace=trace,
            ancestor_receipts={
                "preflight": "a" * 64,
                "transition_state_qualification": "b" * 64,
                "irc_forward": "c" * 64,
                "irc_reverse": "d" * 64,
            },
        )

    assert not (run_root / route / "path").exists()


def test_atomic_path_failure_injection_leaves_no_canonical_partial(tmp_path: Path):
    route = "h-co-1w-oside"
    transition_state, trace = _trace(route, "irc_fwd.xyz", "irc_back.xyz")
    run_root = tmp_path / "run"

    def fail(stage: str) -> None:
        if stage == "before_path_commit":
            raise RuntimeError("injected path crash")

    with pytest.raises(RuntimeError, match="injected path crash"):
        campaign._publish_typed_irc_path(
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
        campaign._publish_typed_irc_path(
            run_root, trace=malformed_trace, **publication_kwargs
        )

    assert not (run_root / route / "path").exists()
    corrected = campaign._publish_typed_irc_path(
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
    campaign._publish_typed_irc_path(run_root, trace=trace, **kwargs)
    unexpected = run_root / route / "path" / "extra"
    unexpected.write_text("untrusted")
    with pytest.raises(ValueError, match="unexpected path artifact"):
        campaign._publish_typed_irc_path(run_root, trace=None, **kwargs)

    other_root = tmp_path / "other"
    other_root.mkdir()
    symlink_root = tmp_path / "linked"
    symlink_root.symlink_to(other_root, target_is_directory=True)
    with pytest.raises(ValueError, match="symlink"):
        campaign._publish_typed_irc_path(symlink_root, trace=trace, **kwargs)


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
    campaign._publish_typed_irc_path(run_root, trace=trace, **kwargs)
    receipt_path = run_root / route / "path" / "receipt.json"
    receipt_path.write_bytes(b"not-json\n")
    original = receipt_path.read_bytes()

    with pytest.raises(ValueError, match="malformed typed path receipt"):
        campaign._publish_typed_irc_path(run_root, trace=None, **kwargs)
    assert receipt_path.read_bytes() == original


def _published_path(run_root: Path):
    route = "h-co-1w-cside"
    transition_state, trace = _trace(route, "irc_back.xyz", "irc_fwd.xyz")
    path = campaign._publish_typed_irc_path(
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
    return campaign._publish_path_hessians(
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
    original = campaign._publish_typed_irc_path(
        run_root, trace=trace, **publication_kwargs
    )
    route_root = run_root / route
    temporary = route_root / ".path.123.456.tmp"
    temporary.mkdir()
    (temporary / "coordinates.f64").write_bytes(b"private coordinates")
    (temporary / "receipt.json").write_bytes(b"private receipt\n")

    resumed = campaign._publish_typed_irc_path(
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
        campaign._publish_typed_irc_path(
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
        if destination.name == aggregate.name and destination.parent.name == "hessians":
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
            return campaign._publish_typed_irc_path(
                run_root,
                campaign_identity="c" * 64,
                route=route,
                atom_mapping_sha256="d" * 64,
                qualified_transition_state=transition_state,
                trace=trace,
            )
    else:
        path = campaign._publish_typed_irc_path(
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
        campaign._publish_typed_irc_path(
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
        campaign._publish_typed_irc_path(
            tmp_path / "safe" / ".." / "escaped", route=route, **kwargs
        )
    with pytest.raises(ValueError, match="unsupported typed endpoint route"):
        campaign._publish_typed_irc_path(tmp_path / "run", route="../escaped", **kwargs)
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


def _qualification_native_result(
    preflight: dict,
    transition_state: Cluster,
    reactant: np.ndarray,
    product: np.ndarray,
):
    from quarry import reaction_path

    masses = np.array(
        [campaign.ISOTOPIC_MASSES_AMU[symbol] for symbol in transition_state.symbols]
    )
    rigid = reaction_path._rigid_motion_basis(transition_state.coords, masses)
    vibrational = np.linalg.svd(rigid.T, full_matrices=True)[2][6:].T
    route_vector = campaign._mapped_route_vector(
        transition_state, masses, reactant, product
    )
    coefficients = vibrational.T @ route_vector
    coefficients /= np.linalg.norm(coefficients)
    complement = np.linalg.svd(coefficients[None, :], full_matrices=True)[2][1:].T
    eigenvectors_in_vibrational_space = np.column_stack((coefficients, complement))
    eigenvalues = np.linspace(0.01, 0.03, vibrational.shape[1])
    eigenvalues[0] = -0.02
    mass_weighted = (
        vibrational
        @ eigenvectors_in_vibrational_space
        @ np.diag(eigenvalues)
        @ eigenvectors_in_vibrational_space.T
        @ vibrational.T
    )
    factors = np.repeat(np.sqrt(masses), 3)
    hessian = factors[:, None] * mass_weighted * factors[None, :]
    _, settings_fingerprint = campaign._canonical_dft_settings(preflight)
    return campaign.NativeHessianResult(
        electronic_hartree=-100.0,
        gradient_hartree_per_bohr=np.zeros_like(transition_state.coords),
        physical_fmax_ev_per_angstrom=0.0,
        cartesian_hessian_hartree_per_bohr2=hessian,
        requested_backend="gpu4pyscf",
        actual_backend="gpu4pyscf",
        gpu_fallback_used=False,
        geometry_fingerprint=campaign.frequency_geometry_fingerprint(transition_state),
        settings_fingerprint=settings_fingerprint,
    )


def _qualified_ancestry(
    run_root: Path,
    *,
    route: str = "h-co-1w-cside",
):
    preflight = _preflight(run_root)
    transition_state, trace = _trace(route, "irc_back.xyz", "irc_fwd.xyz")
    reactant = _frozen_endpoint(route, "irc_back.xyz").coords
    product = _frozen_endpoint(route, "irc_fwd.xyz").coords
    campaign._publish_transition_state_qualification(
        run_root,
        route=route,
        qualified_transition_state=transition_state,
        native_hessian=_qualification_native_result(
            preflight, transition_state, reactant, product
        ),
        mapped_reactant_coordinates_angstrom=reactant,
        mapped_product_coordinates_angstrom=product,
    )
    return preflight, transition_state, trace


def _authoritative_ancestry(
    run_root: Path,
    *,
    route: str = "h-co-1w-cside",
):
    preflight, transition_state, trace = _qualified_ancestry(run_root, route=route)
    trace, runner = _checkpointing_trace_runner(transition_state, trace)
    campaign._run_and_publish_irc(
        run_root,
        route=route,
        _runner=runner,
    )
    return preflight, transition_state, trace


def _rewrite_direction_receipt(
    run_root: Path,
    route: str,
    direction: IrcDirectionPath,
    trace: SellaIrcTrace,
) -> None:
    ancestry = campaign._load_canonical_qualification(run_root, route)
    receipt_path = (
        run_root / route / f"irc-{direction.sella_direction}" / "receipt.json"
    )
    run_identity = json.loads(receipt_path.read_text())["irc_run_identity"]
    assert trace.execution_contract is not None
    receipt = campaign._irc_direction_receipt_payload(
        preflight=ancestry.preflight,
        preflight_receipt_sha256=ancestry.preflight_receipt_sha256,
        ts_qualification_receipt_sha256=ancestry.ts_qualification_receipt_sha256,
        route=route,
        qualified_transition_state=ancestry.qualified_transition_state,
        unstable_mode_mass_scaled=ancestry.unstable_mode_mass_scaled,
        transition_state_vibrational_basis=(
            ancestry.transition_state_vibrational_basis
        ),
        direction=direction,
        masses_amu=trace.masses_amu,
        execution_contract=trace.execution_contract,
        irc_run_identity=run_identity,
    )
    receipt_path.write_bytes(campaign._json_bytes(receipt))
    restart_path = (
        run_root / route / "irc-restart" / direction.sella_direction / "receipt.json"
    )
    restart_path.write_bytes(campaign._json_bytes(receipt))
    execution_path = run_root / route / "irc-execution" / "receipt.json"
    execution = json.loads(execution_path.read_text())
    execution["direction_receipts"][direction.sella_direction] = receipt
    execution_path.write_bytes(campaign._json_bytes(execution))


def test_bound_irc_runner_receives_every_exact_campaign_argument(tmp_path: Path):
    run_root = tmp_path / "run"
    _, transition_state, trace = _qualified_ancestry(run_root)
    initialization = _sella_initialization_state(transition_state, trace)
    trace = _bind_trace_initialization(trace, initialization)
    calls = []

    def spy(ts, settings, **kwargs):
        calls.append((ts, settings, kwargs))
        kwargs["_initialization_callback"](initialization)
        return trace

    published = campaign._run_and_publish_irc(
        run_root,
        route="h-co-1w-cside",
        _runner=spy,
    )

    assert len(calls) == 1
    observed_ts, observed_settings, observed = calls[0]
    assert np.array_equal(observed_ts.coords, transition_state.coords)
    assert observed_settings == campaign.DftSettings(**campaign.DFT_SETTINGS)
    assert set(observed) == {
        "masses_amu",
        "fmax_ev_a",
        "fmax_inner_ev_a",
        "max_steps",
        "step_size_a",
        "_completed_directions",
        "_direction_callback",
        "_initialization_state",
        "_initialization_callback",
    }
    assert np.array_equal(observed["masses_amu"], trace.masses_amu)
    assert observed["step_size_a"] == 0.05
    assert observed["max_steps"] == 200
    assert observed["fmax_ev_a"] == 0.05
    assert observed["fmax_inner_ev_a"] == 0.01
    assert len(published.run_identity) == 64
    assert (
        published.execution_receipt_sha256
        == hashlib.sha256(published.execution_receipt_path.read_bytes()).hexdigest()
    )
    assert set(published.direction_receipt_sha256) == {"forward", "reverse"}
    ancestry = campaign._load_canonical_qualification(run_root, "h-co-1w-cside")
    _, canonical_hashes, canonical_run_identity = (
        campaign._validate_canonical_irc_receipts(ancestry)
    )
    assert canonical_hashes == {
        "execution": published.execution_receipt_sha256,
        "forward": published.direction_receipt_sha256["forward"],
        "reverse": published.direction_receipt_sha256["reverse"],
    }
    assert canonical_run_identity == published.run_identity
    resumed = campaign._run_and_publish_irc(
        run_root,
        route="h-co-1w-cside",
        _runner=lambda *_args, **_kwargs: pytest.fail("full resume reran IRC"),
    )
    assert resumed.run_identity == published.run_identity
    assert resumed.direction_receipt_sha256 == published.direction_receipt_sha256


def test_bound_irc_runner_rejects_wrong_observed_arguments_before_publication(
    tmp_path: Path,
):
    run_root = tmp_path / "run"
    _, _, trace = _qualified_ancestry(run_root)
    assert trace.execution_contract is not None
    wrong = replace(
        trace,
        execution_contract=replace(
            trace.execution_contract,
            step_size_angstrom=0.10,
        ),
    )

    with pytest.raises(ValueError, match="runner observed execution contract"):
        campaign._run_and_publish_irc(
            run_root,
            route="h-co-1w-cside",
            _runner=lambda *_args, **_kwargs: wrong,
        )
    route_root = run_root / "h-co-1w-cside"
    assert not (route_root / "irc-execution").exists()
    assert not (route_root / "irc-forward").exists()
    assert not (route_root / "irc-reverse").exists()


def test_irc_direction_crash_resumes_without_runner_or_overwrite(tmp_path: Path):
    run_root = tmp_path / "run"
    _, transition_state, trace = _qualified_ancestry(run_root)
    initialization = _sella_initialization_state(transition_state, trace)
    trace = _bind_trace_initialization(trace, initialization)
    runner_calls = 0

    def runner(*_args, _initialization_callback, **_kwargs):
        nonlocal runner_calls
        runner_calls += 1
        _initialization_callback(initialization)
        return trace

    def fail_after_forward(stage: str) -> None:
        if stage == "after_irc_direction_commit:forward":
            raise RuntimeError("injected crash after forward")

    with pytest.raises(RuntimeError, match="crash after forward"):
        campaign._run_and_publish_irc(
            run_root,
            route="h-co-1w-cside",
            _runner=runner,
            _failure_injector=fail_after_forward,
        )
    forward_path = run_root / "h-co-1w-cside" / "irc-forward" / "receipt.json"
    forward_before = forward_path.read_bytes()
    assert runner_calls == 1
    assert not (run_root / "h-co-1w-cside" / "irc-reverse").exists()

    resumed = campaign._run_and_publish_irc(
        run_root,
        route="h-co-1w-cside",
        _runner=lambda *_args, **_kwargs: pytest.fail("resume reran IRC"),
    )
    assert runner_calls == 1
    assert forward_path.read_bytes() == forward_before
    assert (run_root / "h-co-1w-cside" / "irc-reverse" / "receipt.json").is_file()
    assert resumed.trace.execution_contract == trace.execution_contract


@pytest.mark.parametrize("race", ["delete", "replace"])
def test_irc_final_validation_rejects_shared_execution_race(tmp_path: Path, race: str):
    route = "h-co-1w-cside"
    run_root = tmp_path / race
    _, transition_state, trace = _qualified_ancestry(run_root, route=route)
    trace, runner = _checkpointing_trace_runner(transition_state, trace)
    execution_path = run_root / route / "irc-execution" / "receipt.json"
    replacement_raw = None

    def race_after_reverse(stage: str) -> None:
        nonlocal replacement_raw
        if stage != "after_irc_direction_commit:reverse":
            return
        if race == "delete":
            shutil.rmtree(execution_path.parent)
            return
        replacement = json.loads(execution_path.read_text())
        replacement["irc_run_identity"] = "f" * 64
        for direction in replacement["direction_receipts"].values():
            direction["irc_run_identity"] = "f" * 64
        replacement_raw = campaign._json_bytes(replacement)
        execution_path.write_bytes(replacement_raw)

    with pytest.raises(ValueError, match="IRC execution receipt"):
        campaign._run_and_publish_irc(
            run_root,
            route=route,
            _runner=runner,
            _failure_injector=race_after_reverse,
        )

    if race == "delete":
        assert not execution_path.exists()
    else:
        assert execution_path.read_bytes() == replacement_raw
    assert (run_root / route / "irc-forward" / "receipt.json").is_file()
    assert (run_root / route / "irc-reverse" / "receipt.json").is_file()


def test_public_path_rejects_standalone_irc_run_identity_divergence(tmp_path: Path):
    run_root = tmp_path / "run"
    _authoritative_ancestry(run_root)
    reverse_receipt = run_root / "h-co-1w-cside" / "irc-reverse" / "receipt.json"
    _rewrite_canonical_json(
        reverse_receipt,
        lambda receipt: receipt.__setitem__("irc_run_identity", "f" * 64),
    )

    with pytest.raises(
        ValueError,
        match="(standalone/shared execution|restart/final canonical).*receipt",
    ):
        campaign.publish_typed_irc_path(run_root, route="h-co-1w-cside")
    assert not (run_root / "h-co-1w-cside" / "path").exists()


def test_canonical_ts_qualification_persists_native_evidence_and_reruns_gate(
    tmp_path: Path,
):
    route = "h-co-1w-cside"
    run_root = tmp_path / "run"
    preflight = _preflight(run_root)
    transition_state, _ = _trace(route, "irc_back.xyz", "irc_fwd.xyz")
    reactant = _frozen_endpoint(route, "irc_back.xyz").coords
    product = _frozen_endpoint(route, "irc_fwd.xyz").coords
    published = campaign._publish_transition_state_qualification(
        run_root,
        route=route,
        qualified_transition_state=transition_state,
        native_hessian=_qualification_native_result(
            preflight, transition_state, reactant, product
        ),
        mapped_reactant_coordinates_angstrom=reactant,
        mapped_product_coordinates_angstrom=product,
    )
    qualification_root = published.receipt_path.parent
    assert {path.name for path in qualification_root.iterdir()} == {
        "coordinates.f64",
        "gradient.f64",
        "hessian.f64",
        "mapped-reactant.f64",
        "mapped-product.f64",
        "receipt.json",
    }
    assert published.receipt["schema"] == "d2c-ts-qualification-v3"
    assert published.receipt["gate_evidence"]["accepted"] is True
    assert published.receipt["native_hessian"]["electronic_hartree"] == -100.0
    assert (
        campaign._publish_transition_state_qualification(
            run_root, route=route
        ).receipt_sha256
        == published.receipt_sha256
    )
    with pytest.raises(TypeError, match="unstable_mode"):
        campaign._publish_transition_state_qualification(
            run_root, route=route, unstable_mode_mass_scaled=np.ones(18)
        )

    hessian_path = qualification_root / "hessian.f64"
    replacement = np.eye(18, dtype="<f8").tobytes()
    hessian_path.write_bytes(replacement)
    receipt = json.loads(published.receipt_path.read_text())
    receipt["native_hessian"]["cartesian_hessian"]["sha256"] = hashlib.sha256(
        replacement
    ).hexdigest()
    published.receipt_path.write_bytes(campaign._json_bytes(receipt))
    with pytest.raises(ValueError, match="exactly one significant negative"):
        campaign._publish_transition_state_qualification(run_root, route=route)


def test_ts_qualification_post_commit_preflight_race_removes_owned_child(
    tmp_path: Path,
):
    route = "h-co-1w-cside"
    run_root = tmp_path / "run"
    preflight = _preflight(run_root)
    transition_state, _ = _trace(route, "irc_back.xyz", "irc_fwd.xyz")
    reactant = _frozen_endpoint(route, "irc_back.xyz").coords
    product = _frozen_endpoint(route, "irc_fwd.xyz").coords
    preflight_path = run_root / campaign.PREFLIGHT_RECEIPT

    def mutate_after_commit(stage: str) -> None:
        if stage == "after_ts_qualification_commit":
            _rewrite_canonical_json(
                preflight_path,
                lambda receipt: receipt.__setitem__(
                    "created_utc", "2026-09-07T12:00:01Z"
                ),
            )

    with pytest.raises(ValueError, match="preflight receipt SHA-256"):
        campaign._publish_transition_state_qualification(
            run_root,
            route=route,
            qualified_transition_state=transition_state,
            native_hessian=_qualification_native_result(
                preflight, transition_state, reactant, product
            ),
            mapped_reactant_coordinates_angstrom=reactant,
            mapped_product_coordinates_angstrom=product,
            _failure_injector=mutate_after_commit,
        )
    assert not (run_root / route / "ts-qualification").exists()


def test_legacy_v3_and_self_attested_v1_roots_fail_with_fresh_root_guidance(
    tmp_path: Path,
):
    v3_root = tmp_path / "v3"
    preflight = _preflight(v3_root)
    preflight["schema"] = "d2c-sct-campaign-preflight-v3"
    (v3_root / campaign.PREFLIGHT_RECEIPT).write_bytes(campaign._json_bytes(preflight))
    with pytest.raises(ValueError, match="legacy D2c run root.*fresh v6 run root"):
        campaign._publish_transition_state_qualification(v3_root, route="h-co-1w-cside")

    v1_root = tmp_path / "v1"
    _preflight(v1_root)
    qualification = v1_root / "h-co-1w-cside" / "ts-qualification"
    qualification.mkdir(parents=True)
    (qualification / "receipt.json").write_bytes(
        campaign._json_bytes({"schema": "d2c-ts-qualification-v1"})
    )
    with pytest.raises(
        ValueError, match="legacy D2c TS qualification.*fresh v6 run root"
    ):
        campaign._publish_transition_state_qualification(v1_root, route="h-co-1w-cside")


def _authoritative_path(run_root: Path):
    _, transition_state, _ = _authoritative_ancestry(run_root)
    path = campaign.publish_typed_irc_path(
        run_root,
        route="h-co-1w-cside",
    )
    return transition_state, path


def _energy_reproducing_result(cluster, path, *, offset_ev: float = 0.0):
    point_index = int(cluster.name.rsplit("-", 1)[-1])
    result = _native_result(cluster)
    path_energy_ev = path.receipt["points"][point_index]["electronic_energy_ev"]
    canonical_settings = campaign.frequency_settings_fingerprint(
        campaign.DftSettings(**campaign.DFT_SETTINGS)
    )
    return replace(
        result,
        electronic_hartree=(path_energy_ev + offset_ev) / campaign.HARTREE_TO_EV,
        settings_fingerprint=canonical_settings,
    )


def test_authoritative_path_direct_and_resume_revalidate_canonical_ancestry(
    tmp_path: Path,
):
    run_root = tmp_path / "run"
    preflight, _, _ = _authoritative_ancestry(run_root)

    published = campaign.publish_typed_irc_path(run_root, route="h-co-1w-cside")
    assert published.receipt["schema"] == "d2c-typed-irc-path-v3"
    assert set(published.receipt["ancestor_receipts"]) == {
        "preflight",
        "transition_state_qualification",
        "irc_execution",
        "irc_forward",
        "irc_reverse",
    }
    execution_raw = (
        run_root / "h-co-1w-cside" / "irc-execution" / "receipt.json"
    ).read_bytes()
    assert (
        published.receipt["ancestor_receipts"]["irc_execution"]
        == hashlib.sha256(execution_raw).hexdigest()
    )
    resumed = campaign.publish_typed_irc_path(run_root, route="h-co-1w-cside")
    assert resumed.receipt_sha256 == published.receipt_sha256

    preflight_path = run_root / campaign.PREFLIGHT_RECEIPT
    preflight["created_utc"] = "2026-09-07T12:00:01Z"
    preflight_path.write_bytes(campaign._json_bytes(preflight))
    with pytest.raises(ValueError, match="preflight receipt SHA-256"):
        campaign.publish_typed_irc_path(run_root, route="h-co-1w-cside")


def test_authoritative_path_rejects_legacy_v2_with_fresh_root_guidance(
    tmp_path: Path,
):
    route = "h-co-1w-cside"
    run_root = tmp_path / "run"
    _authoritative_ancestry(run_root, route=route)
    campaign.publish_typed_irc_path(run_root, route=route)
    path_receipt = run_root / route / "path" / "receipt.json"
    _rewrite_canonical_json(
        path_receipt,
        lambda receipt: receipt.__setitem__("schema", "d2c-typed-irc-path-v2"),
    )

    with pytest.raises(ValueError, match="legacy.*fresh run root"):
        campaign.publish_typed_irc_path(run_root, route=route)


@pytest.mark.parametrize("resume", [False, True])
def test_authoritative_path_rejects_deleted_shared_execution_receipt(
    tmp_path: Path, resume: bool
):
    route = "h-co-1w-cside"
    run_root = tmp_path / "run"
    _authoritative_ancestry(run_root, route=route)
    published = (
        campaign.publish_typed_irc_path(run_root, route=route) if resume else None
    )
    shutil.rmtree(run_root / route / "irc-execution")

    with pytest.raises(ValueError, match="IRC execution receipt"):
        campaign.publish_typed_irc_path(run_root, route=route)
    if published is None:
        assert not (run_root / route / "path").exists()
    else:
        assert (run_root / route / "path" / "receipt.json").read_bytes() == (
            campaign._json_bytes(published.receipt)
        )


def test_authoritative_path_rejects_shared_nested_direction_divergence(
    tmp_path: Path,
):
    route = "h-co-1w-cside"
    run_root = tmp_path / "run"
    _authoritative_ancestry(run_root, route=route)
    execution_path = run_root / route / "irc-execution" / "receipt.json"

    def mutate_nested_direction(receipt):
        receipt["direction_receipts"]["forward"]["points"][-1][
            "electronic_energy_ev"
        ] -= 0.25

    _rewrite_canonical_json(execution_path, mutate_nested_direction)

    with pytest.raises(ValueError, match="restart/final canonical.*receipt"):
        campaign.publish_typed_irc_path(run_root, route=route)
    assert not (run_root / route / "path").exists()


def test_authoritative_path_api_does_not_accept_caller_identity_ts_or_trace(
    tmp_path: Path,
):
    run_root = tmp_path / "run"
    _, transition_state, trace = _authoritative_ancestry(run_root)

    with pytest.raises(TypeError, match="trace"):
        campaign.publish_typed_irc_path(
            run_root,
            route="h-co-1w-cside",
            trace=trace,  # type: ignore[call-arg]
        )
    with pytest.raises(TypeError, match="campaign_identity"):
        campaign.publish_typed_irc_path(
            run_root,
            route="h-co-1w-cside",
            campaign_identity="f" * 64,
            qualified_transition_state=transition_state,
        )
    assert not (run_root / "h-co-1w-cside" / "path").exists()


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("algorithm", "claimed-sella"),
        ("step_size_angstrom", 0.0500000001),
        ("maximum_steps", 201),
        ("maximum_retained_points", 202),
        ("outer_fmax_ev_per_angstrom", 0.0500000001),
        ("inner_fmax_ev_per_angstrom", 0.0100000001),
    ],
)
def test_authoritative_path_rejects_irc_contract_tamper(
    tmp_path: Path, field: str, value
):
    run_root = tmp_path / field
    _authoritative_ancestry(run_root)
    receipt_path = run_root / "h-co-1w-cside" / "irc-forward" / "receipt.json"
    _rewrite_canonical_json(
        receipt_path, lambda receipt: receipt.__setitem__(field, value)
    )

    with pytest.raises(ValueError, match="IRC forward"):
        campaign.publish_typed_irc_path(run_root, route="h-co-1w-cside")
    assert not (run_root / "h-co-1w-cside" / "path").exists()


def test_authoritative_path_terminal_fmax_boundary_is_strict(tmp_path: Path):
    run_root = tmp_path / "run"
    _, _, trace = _authoritative_ancestry(run_root)
    directions = list(trace.directions)
    points = list(directions[0].points)
    points[-1] = replace(
        points[-1],
        projected_fmax_ev_per_angstrom=campaign.BOUNDS["irc"][
            "outer_fmax_ev_per_angstrom"
        ],
    )
    directions[0] = replace(directions[0], points=tuple(points))
    boundary_trace = replace(trace, directions=tuple(directions))
    _rewrite_direction_receipt(run_root, "h-co-1w-cside", directions[0], boundary_trace)

    with pytest.raises(ValueError, match="terminal IRC fmax must be strictly below"):
        campaign.publish_typed_irc_path(run_root, route="h-co-1w-cside")


@pytest.mark.parametrize("point_count", [201, 202])
def test_authoritative_path_enforces_retained_point_boundary(
    tmp_path: Path, point_count: int
):
    run_root = tmp_path / str(point_count)
    _, _, trace = _authoritative_ancestry(run_root)

    def expand(direction):
        start = direction.points[0]
        terminal = direction.points[-1]
        return replace(
            direction,
            points=tuple(
                IrcPoint(
                    outer_step=index,
                    coordinates_angstrom=(
                        start.coordinates_angstrom
                        + (index / (point_count - 1))
                        * (terminal.coordinates_angstrom - start.coordinates_angstrom)
                    ),
                    electronic_energy_ev=float(-10.0 - index / 1000.0),
                    projected_fmax_ev_per_angstrom=0.01,
                )
                for index in range(point_count)
            ),
        )

    expanded = replace(
        trace, directions=tuple(expand(direction) for direction in trace.directions)
    )
    if point_count == 201:
        for direction in expanded.directions:
            _rewrite_direction_receipt(run_root, "h-co-1w-cside", direction, expanded)
        assert (
            campaign.publish_typed_irc_path(run_root, route="h-co-1w-cside").receipt[
                "point_count"
            ]
            == 401
        )
    else:
        _rewrite_direction_receipt(
            run_root, "h-co-1w-cside", expanded.directions[0], expanded
        )
        with pytest.raises(ValueError, match="maximum 201 retained points"):
            campaign.publish_typed_irc_path(run_root, route="h-co-1w-cside")


def test_authoritative_path_tangent_is_rigid_invariant_and_one_bad_side_rejects(
    tmp_path: Path,
):
    route = "h-co-1w-cside"
    invariant_root = tmp_path / "invariant"
    _, _, trace = _authoritative_ancestry(invariant_root)
    rotated_directions = list(trace.directions)
    points = list(rotated_directions[0].points)
    angle = 0.73
    rotation = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    points[1] = replace(
        points[1],
        coordinates_angstrom=points[1].coordinates_angstrom @ rotation
        + np.array([8.0, -3.0, 4.0]),
    )
    rotated_directions[0] = replace(rotated_directions[0], points=tuple(points))
    rigid_trace = replace(trace, directions=tuple(rotated_directions))
    _rewrite_direction_receipt(
        invariant_root, route, rotated_directions[0], rigid_trace
    )
    assert campaign.publish_typed_irc_path(invariant_root, route=route).receipt[
        "accepted"
    ]

    bad_root = tmp_path / "bad-side"
    _, _, bad_trace = _authoritative_ancestry(bad_root)
    ancestry = campaign._load_canonical_qualification(bad_root, route)
    bad_directions = list(bad_trace.directions)
    bad_points = list(bad_directions[0].points)
    mode = ancestry.unstable_mode_mass_scaled
    basis = ancestry.transition_state_vibrational_basis
    candidate = basis[:, 0].copy()
    candidate -= mode * np.dot(candidate, mode)
    if np.linalg.norm(candidate) < 1.0e-8:
        candidate = basis[:, 1] - mode * np.dot(basis[:, 1], mode)
    candidate /= np.linalg.norm(candidate)
    masses = np.asarray(ancestry.preflight["routes"][route]["masses_amu"])
    displacement = (candidate / np.repeat(np.sqrt(masses), 3)).reshape(
        ancestry.qualified_transition_state.coords.shape
    )
    bad_points[1] = replace(
        bad_points[1],
        coordinates_angstrom=ancestry.qualified_transition_state.coords
        + 0.01 * displacement,
    )
    bad_directions[0] = replace(bad_directions[0], points=tuple(bad_points))
    bad_trace = replace(bad_trace, directions=tuple(bad_directions))
    _rewrite_direction_receipt(bad_root, route, bad_directions[0], bad_trace)
    with pytest.raises(ValueError, match="TS-adjacent tangent overlap"):
        campaign.publish_typed_irc_path(bad_root, route=route)


def test_authoritative_hessian_computes_missing_points_and_never_reruns_resume(
    tmp_path: Path,
):
    run_root = tmp_path / "run"
    _, path = _authoritative_path(run_root)
    direct_calls = []

    def direct(cluster):
        direct_calls.append(int(cluster.name.rsplit("-", 1)[-1]))
        return _energy_reproducing_result(cluster, path)

    first = campaign._publish_authoritative_path_hessians(
        run_root,
        route="h-co-1w-cside",
        evaluator=direct,
    )
    assert direct_calls == list(range(5))
    assert first.receipt["accepted"] is True

    resume_calls = []

    def resume(cluster):
        resume_calls.append(int(cluster.name.rsplit("-", 1)[-1]))
        return _energy_reproducing_result(cluster, path)

    resumed = campaign._publish_authoritative_path_hessians(
        run_root,
        route="h-co-1w-cside",
        evaluator=resume,
    )
    assert resume_calls == []
    assert resumed.receipt["schema"] == "d2c-native-path-hessians-v2"
    assert resumed.receipt_sha256 == first.receipt_sha256


def test_authoritative_hessian_resume_validates_persisted_energy_without_evaluator(
    tmp_path: Path,
):
    run_root = tmp_path / "run"
    _, path = _authoritative_path(run_root)
    campaign._publish_authoritative_path_hessians(
        run_root,
        route="h-co-1w-cside",
        evaluator=lambda cluster: _energy_reproducing_result(cluster, path),
    )
    tolerance = campaign.BOUNDS["hessian"][
        "electronic_energy_reproduction_absolute_tolerance_ev"
    ]
    point_receipt_path = (
        run_root / "h-co-1w-cside" / "hessians" / "points" / "000000" / "receipt.json"
    )

    def forge_matching_receipt_shape(receipt):
        receipt["electronic_hartree"] += 2.0 * tolerance / campaign.HARTREE_TO_EV
        native_ev = receipt["electronic_hartree"] * campaign.HARTREE_TO_EV
        path_ev = receipt["energy_reproduction"]["path_electronic_energy_ev"]
        receipt["energy_reproduction"]["native_electronic_energy_ev"] = native_ev
        receipt["energy_reproduction"]["absolute_difference_ev"] = abs(
            native_ev - path_ev
        )

    _rewrite_canonical_json(point_receipt_path, forge_matching_receipt_shape)
    with pytest.raises(ValueError, match="persisted Hessian point 0 electronic energy"):
        campaign._publish_authoritative_path_hessians(
            run_root,
            route="h-co-1w-cside",
            evaluator=lambda _cluster: pytest.fail("cached point reached evaluator"),
        )


def test_authoritative_hessian_evaluator_must_match_preflight_settings(tmp_path: Path):
    run_root = tmp_path / "run"
    _, path = _authoritative_path(run_root)

    def wrong_settings(cluster):
        return replace(
            _energy_reproducing_result(cluster, path),
            settings_fingerprint="caller-controlled",
        )

    with pytest.raises(ValueError, match="settings fingerprint mismatch"):
        campaign._publish_authoritative_path_hessians(
            run_root,
            route="h-co-1w-cside",
            evaluator=wrong_settings,
        )
    assert not (run_root / "h-co-1w-cside" / "hessians" / "points" / "000000").exists()


@pytest.mark.parametrize("multiple", [0.999, 1.001])
def test_authoritative_hessian_energy_tolerance_boundary(
    tmp_path: Path, multiple: float
):
    run_root = tmp_path / str(multiple)
    _, path = _authoritative_path(run_root)
    tolerance = campaign.BOUNDS["hessian"][
        "electronic_energy_reproduction_absolute_tolerance_ev"
    ]

    def evaluator(cluster):
        index = int(cluster.name.rsplit("-", 1)[-1])
        offset = tolerance * multiple if index == 2 else 0.0
        return _energy_reproducing_result(cluster, path, offset_ev=offset)

    kwargs = {
        "route": "h-co-1w-cside",
        "evaluator": evaluator,
    }
    if multiple < 1.0:
        assert campaign._publish_authoritative_path_hessians(
            run_root, **kwargs
        ).receipt["accepted"]
    else:
        with pytest.raises(ValueError, match="electronic energy reproduction"):
            campaign._publish_authoritative_path_hessians(run_root, **kwargs)
        assert not (
            run_root / "h-co-1w-cside" / "hessians" / "points" / "000002"
        ).exists()


def test_authoritative_hessian_resume_rejects_irc_ancestry_tamper_before_compute(
    tmp_path: Path,
):
    run_root = tmp_path / "run"
    _, path = _authoritative_path(run_root)

    def evaluator(cluster):
        return _energy_reproducing_result(cluster, path)

    kwargs = {
        "route": "h-co-1w-cside",
    }
    campaign._publish_authoritative_path_hessians(
        run_root, evaluator=evaluator, **kwargs
    )
    receipt_path = run_root / "h-co-1w-cside" / "irc-reverse" / "receipt.json"
    _rewrite_canonical_json(
        receipt_path,
        lambda receipt: receipt.__setitem__("algorithm", "tampered"),
    )

    with pytest.raises(
        ValueError,
        match="(standalone/shared execution|restart/final canonical).*receipt",
    ):
        campaign._publish_authoritative_path_hessians(
            run_root,
            evaluator=lambda _cluster: pytest.fail("tamper reached evaluator"),
            **kwargs,
        )


def test_hessian_precommit_revalidation_rejects_shared_execution_mutation(
    tmp_path: Path,
):
    route = "h-co-1w-cside"
    run_root = tmp_path / "run"
    _, path = _authoritative_path(run_root)
    execution_path = run_root / route / "irc-execution" / "receipt.json"
    evaluation_calls = []

    def mutate_during_evaluation(cluster):
        evaluation_calls.append(int(cluster.name.rsplit("-", 1)[-1]))

        def mutate_nested_direction(receipt):
            receipt["direction_receipts"]["reverse"]["points"][-1][
                "electronic_energy_ev"
            ] -= 0.25

        _rewrite_canonical_json(execution_path, mutate_nested_direction)
        return _energy_reproducing_result(cluster, path)

    with pytest.raises(
        ValueError,
        match="(standalone/shared execution|restart/final canonical).*receipt",
    ):
        campaign._publish_authoritative_path_hessians(
            run_root,
            route=route,
            evaluator=mutate_during_evaluation,
        )

    assert evaluation_calls == [0]
    assert not (run_root / route / "hessians" / "points" / "000000").exists()


def test_authoritative_hessian_rejects_low_level_unreproduced_checkpoint(
    tmp_path: Path,
):
    run_root = tmp_path / "run"
    _, path = _authoritative_path(run_root)
    ancestry, _, ancestor_receipts = campaign._validate_authoritative_published_path(
        run_root, "h-co-1w-cside"
    )
    campaign._publish_path_hessians(
        run_root,
        campaign_identity=ancestry.campaign_identity,
        route="h-co-1w-cside",
        atom_mapping_sha256=ancestry.atom_mapping_sha256,
        qualified_transition_state=ancestry.qualified_transition_state,
        path=path,
        settings_fingerprint="settings-v1",
        backend_policy={
            "requested_backend": "gpu4pyscf",
            "allow_cpu_fallback": True,
        },
        evaluator=_native_result,
        ancestor_receipts=ancestor_receipts,
    )

    with pytest.raises(ValueError, match="legacy D2c Hessian point v1"):
        campaign._publish_authoritative_path_hessians(
            run_root,
            route="h-co-1w-cside",
            evaluator=lambda cluster: _energy_reproducing_result(cluster, path),
        )


def test_authoritative_hessian_rejects_matching_energy_legacy_hessian_tree(
    tmp_path: Path,
):
    run_root = tmp_path / "run"
    _, path = _authoritative_path(run_root)
    ancestry, _, ancestor_receipts = campaign._validate_authoritative_published_path(
        run_root, "h-co-1w-cside"
    )
    _, canonical_settings = campaign._canonical_dft_settings(ancestry.preflight)

    def matching_energy_different_hessian(cluster):
        return _energy_reproducing_result(cluster, path)

    legacy = campaign._publish_path_hessians(
        run_root,
        campaign_identity=ancestry.campaign_identity,
        route="h-co-1w-cside",
        atom_mapping_sha256=ancestry.atom_mapping_sha256,
        qualified_transition_state=ancestry.qualified_transition_state,
        path=path,
        settings_fingerprint=canonical_settings,
        backend_policy=campaign._canonical_backend_policy(ancestry.preflight),
        evaluator=matching_energy_different_hessian,
        ancestor_receipts=ancestor_receipts,
    )
    assert legacy.receipt["schema"] == "d2c-native-path-hessians-v1"
    with pytest.raises(ValueError, match="legacy D2c Hessian point v1"):
        campaign._publish_authoritative_path_hessians(
            run_root,
            route="h-co-1w-cside",
            evaluator=lambda _cluster: pytest.fail("legacy tree reached evaluator"),
        )


def test_authoritative_hessian_partial_resume_evaluates_only_missing_points(
    tmp_path: Path,
):
    run_root = tmp_path / "run"
    _, path = _authoritative_path(run_root)
    first_calls: list[int] = []

    def first_evaluator(cluster):
        first_calls.append(int(cluster.name.rsplit("-", 1)[-1]))
        return _energy_reproducing_result(cluster, path)

    def fail_before_third(stage: str) -> None:
        if stage == "before_hessian_point_commit:2":
            raise RuntimeError("stop after two cached points")

    with pytest.raises(RuntimeError, match="two cached points"):
        campaign._publish_authoritative_path_hessians(
            run_root,
            route="h-co-1w-cside",
            evaluator=first_evaluator,
            _failure_injector=fail_before_third,
        )
    assert first_calls == [0, 1, 2]

    resumed_calls: list[int] = []

    def resumed_evaluator(cluster):
        resumed_calls.append(int(cluster.name.rsplit("-", 1)[-1]))
        return _energy_reproducing_result(cluster, path)

    published = campaign._publish_authoritative_path_hessians(
        run_root,
        route="h-co-1w-cside",
        evaluator=resumed_evaluator,
    )
    assert resumed_calls == [2, 3, 4]
    assert len(published.point_results) == 5
    assert published.receipt["schema"] == "d2c-native-path-hessians-v2"


@pytest.mark.parametrize("stage", ["point", "aggregate"])
def test_authoritative_hessian_post_commit_ancestry_mutation_fails_closed(
    tmp_path: Path,
    stage: str,
):
    run_root = tmp_path / stage
    _, path = _authoritative_path(run_root)
    irc_receipt = run_root / "h-co-1w-cside" / "irc-forward" / "receipt.json"

    def mutate_after_commit(hook: str) -> None:
        target = (
            hook == "after_hessian_point_commit:0"
            if stage == "point"
            else hook == "after_hessian_aggregate_commit"
        )
        if target:
            _rewrite_canonical_json(
                irc_receipt,
                lambda receipt: receipt.__setitem__("algorithm", "raced"),
            )

    with pytest.raises(
        ValueError,
        match="(standalone/shared execution|restart/final canonical).*receipt",
    ):
        campaign._publish_authoritative_path_hessians(
            run_root,
            route="h-co-1w-cside",
            evaluator=lambda cluster: _energy_reproducing_result(cluster, path),
            _failure_injector=mutate_after_commit,
        )
    hessian_root = run_root / "h-co-1w-cside" / "hessians"
    if stage == "point":
        assert not (hessian_root / "points" / "000000").exists()
    else:
        assert not (hessian_root / "receipt.json").exists()
        assert len(list((hessian_root / "points").iterdir())) == 5


def test_post_commit_cleanup_preserves_foreign_point_replacement(tmp_path: Path):
    run_root = tmp_path / "run"
    _, path = _authoritative_path(run_root)
    point_dir = run_root / "h-co-1w-cside" / "hessians" / "points" / "000000"

    def replace_owned_point(stage: str) -> None:
        if stage == "after_hessian_point_commit:0":
            shutil.rmtree(point_dir)
            point_dir.mkdir()
            (point_dir / "foreign-marker").write_text("preserve me")
            raise RuntimeError("foreign replacement raced")

    with pytest.raises(RuntimeError, match="foreign replacement raced"):
        campaign._publish_authoritative_path_hessians(
            run_root,
            route="h-co-1w-cside",
            evaluator=lambda cluster: _energy_reproducing_result(cluster, path),
            _failure_injector=replace_owned_point,
        )
    assert (point_dir / "foreign-marker").read_text() == "preserve me"


def test_owned_directory_cleanup_preserves_same_inode_mutation_before_quarantine(
    monkeypatch, tmp_path: Path
):
    publication = tmp_path / "publication"
    publication.mkdir()
    receipt = publication / "receipt.json"
    receipt.write_bytes(b"owned")
    status = publication.stat(follow_symlinks=False)
    identity = (status.st_dev, status.st_ino)
    expected_hashes = {"receipt.json": hashlib.sha256(b"owned").hexdigest()}
    real_rename = campaign._renameat2_noreplace

    def mutate_then_rename(source: Path, destination: Path) -> None:
        if source == publication:
            receipt.write_bytes(b"foreign mutation")
        real_rename(source, destination)

    monkeypatch.setattr(campaign, "_renameat2_noreplace", mutate_then_rename)
    assert not campaign._safe_remove_owned_directory(
        publication, identity, expected_hashes
    )
    assert receipt.read_bytes() == b"foreign mutation"


def test_owned_file_cleanup_preserves_same_inode_mutation_before_quarantine(
    monkeypatch, tmp_path: Path
):
    publication = tmp_path / "receipt.json"
    publication.write_bytes(b"owned")
    status = publication.stat(follow_symlinks=False)
    identity = (status.st_dev, status.st_ino)
    real_rename = campaign._renameat2_noreplace

    def mutate_then_rename(source: Path, destination: Path) -> None:
        if source == publication:
            publication.write_bytes(b"foreign mutation")
        real_rename(source, destination)

    monkeypatch.setattr(campaign, "_renameat2_noreplace", mutate_then_rename)
    assert not campaign._safe_remove_owned_file(
        publication, identity, hashlib.sha256(b"owned").hexdigest()
    )
    assert publication.read_bytes() == b"foreign mutation"


def test_path_post_commit_ancestry_mutation_removes_only_owned_publication(
    tmp_path: Path,
):
    run_root = tmp_path / "run"
    _authoritative_ancestry(run_root)
    preflight_path = run_root / campaign.PREFLIGHT_RECEIPT

    def mutate_after_commit(stage: str) -> None:
        if stage == "after_path_commit":
            _rewrite_canonical_json(
                preflight_path,
                lambda receipt: receipt.__setitem__(
                    "created_utc", "2026-09-07T12:00:01Z"
                ),
            )

    ancestry = campaign._load_canonical_qualification(run_root, "h-co-1w-cside")
    trace, hashes, _ = campaign._validate_canonical_irc_receipts(ancestry)
    ancestors = {
        "preflight": ancestry.preflight_receipt_sha256,
        "transition_state_qualification": ancestry.ts_qualification_receipt_sha256,
        "irc_execution": hashes["execution"],
        "irc_forward": hashes["forward"],
        "irc_reverse": hashes["reverse"],
    }

    def validate_ancestry() -> None:
        campaign._load_canonical_qualification(run_root, "h-co-1w-cside")

    with pytest.raises(ValueError, match="preflight receipt SHA-256"):
        campaign._publish_typed_irc_path(
            ancestry.root,
            campaign_identity=ancestry.campaign_identity,
            route=ancestry.route,
            atom_mapping_sha256=ancestry.atom_mapping_sha256,
            qualified_transition_state=ancestry.qualified_transition_state,
            trace=trace,
            ancestor_receipts=ancestors,
            _ancestry_validator=validate_ancestry,
            _failure_injector=mutate_after_commit,
        )
    assert not (run_root / "h-co-1w-cside" / "path").exists()


def test_public_production_apis_do_not_expose_failure_injection():
    for public_api in (
        campaign.publish_transition_state_qualification,
        campaign.run_and_publish_irc,
        campaign.publish_typed_irc_path,
        campaign.publish_path_hessians,
    ):
        assert "_failure_injector" not in inspect.signature(public_api).parameters


def test_public_hessian_api_rejects_caller_supplied_settings(tmp_path: Path):
    run_root = tmp_path / "run"
    _authoritative_path(run_root)
    with pytest.raises(TypeError, match="settings_fingerprint"):
        campaign.publish_path_hessians(
            run_root,
            route="h-co-1w-cside",
            settings_fingerprint="caller-attested",
            evaluator=lambda _cluster: pytest.fail("must not evaluate"),
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


def test_public_production_apis_reject_injected_scientific_evidence(tmp_path: Path):
    run_root = tmp_path / "run"
    preflight = _preflight(run_root)
    route = "h-co-1w-cside"
    transition_state, trace = _trace(route, "irc_back.xyz", "irc_fwd.xyz")
    reactant = _frozen_endpoint(route, "irc_back.xyz").coords
    product = _frozen_endpoint(route, "irc_fwd.xyz").coords
    native = _qualification_native_result(
        preflight, transition_state, reactant, product
    )

    with pytest.raises(TypeError, match="native_hessian"):
        campaign.publish_transition_state_qualification(
            run_root,
            route=route,
            native_hessian=native,
        )
    with pytest.raises(TypeError, match="_runner"):
        campaign.run_and_publish_irc(
            run_root,
            route=route,
            _runner=lambda *_args, **_kwargs: trace,
        )
    with pytest.raises(TypeError, match="evaluator"):
        campaign.publish_path_hessians(
            run_root,
            route=route,
            evaluator=lambda _cluster: native,
        )


def test_production_boundary_rejects_code_dependency_and_endpoint_drift(
    monkeypatch, tmp_path: Path
):
    run_root = tmp_path / "run"
    preflight = _preflight(run_root)
    route = "h-co-1w-cside"

    monkeypatch.setattr(
        campaign,
        "_current_code_dependency_identity",
        lambda: {
            "git_sha": "b" * 40,
            "dependencies": FIXED_DEPENDENCIES,
            "executable_modules": FIXED_MODULE_MANIFEST,
            "python": preflight["campaign"]["python"],
        },
    )
    with pytest.raises(ValueError, match="current Git SHA"):
        campaign._validate_production_boundary(run_root, route)

    copied = tmp_path / "bundle"
    shutil.copytree(BUNDLE_ROOT, copied)
    endpoint = copied / route / "irc_back.xyz"
    endpoint.write_text(endpoint.read_text() + "\n")
    monkeypatch.setattr(campaign, "DEFAULT_BUNDLE_ROOT", copied)
    monkeypatch.setattr(
        campaign,
        "_current_code_dependency_identity",
        lambda: {
            "git_sha": FIXED_GIT_SHA,
            "dependencies": FIXED_DEPENDENCIES,
            "executable_modules": FIXED_MODULE_MANIFEST,
            "native_payloads": FIXED_NATIVE_PAYLOAD_MANIFEST,
            "python": preflight["campaign"]["python"],
        },
    )
    with pytest.raises(ValueError, match="bundle|endpoint"):
        campaign._validate_production_boundary(run_root, route)


@pytest.mark.parametrize(
    "entrypoint",
    [
        "publish_transition_state_qualification",
        "run_and_publish_irc",
        "publish_typed_irc_path",
        "publish_path_hessians",
    ],
)
def test_every_public_production_entrypoint_binds_current_identity(
    monkeypatch, tmp_path: Path, entrypoint: str
):
    def reject_boundary(_run_root: Path, _route: str):
        raise ValueError("live production identity drifted")

    monkeypatch.setattr(campaign, "_validate_production_boundary", reject_boundary)

    with pytest.raises(ValueError, match="live production identity drifted"):
        getattr(campaign, entrypoint)(tmp_path / "run", route="h-co-1w-cside")


def test_mid_irc_crash_checkpoints_completed_direction_and_resume_skips_it(
    tmp_path: Path,
):
    run_root = tmp_path / "run"
    _, transition_state, trace = _qualified_ancestry(run_root)
    initialization = _sella_initialization_state(transition_state, trace)
    trace = _bind_trace_initialization(trace, initialization)
    route = "h-co-1w-cside"
    first_calls: list[tuple[str, ...]] = []

    def crashing_runner(
        *_args,
        _completed_directions,
        _direction_callback,
        _initialization_callback,
        **_kwargs,
    ):
        first_calls.append(tuple(sorted(_completed_directions)))
        _initialization_callback(initialization)
        _direction_callback(trace.directions[0])
        raise RuntimeError("in-process crash after forward IRC")

    with pytest.raises(RuntimeError, match="crash after forward IRC"):
        campaign._run_and_publish_irc(
            run_root,
            route=route,
            _runner=crashing_runner,
        )
    checkpoint = run_root / route / "irc-restart" / "forward" / "receipt.json"
    checkpoint_before = checkpoint.read_bytes()
    initialization_receipt = json.loads(
        (
            run_root / route / "irc-restart" / "initialization" / "receipt.json"
        ).read_text()
    )
    assert (
        initialization_receipt["initialization_fingerprint"]
        == initialization.fingerprint
    )
    assert set(initialization_receipt["initialization_state"]) >= {"H0", "v0ts"}
    assert json.loads(checkpoint_before)["initialization_fingerprint"] == (
        initialization.fingerprint
    )
    assert first_calls == [()]
    assert not (run_root / route / "irc-execution").exists()

    resumed_calls: list[tuple[str, ...]] = []

    def resumed_runner(
        *_args,
        _completed_directions,
        _direction_callback,
        _initialization_state,
        **_kwargs,
    ):
        resumed_calls.append(tuple(sorted(_completed_directions)))
        assert _initialization_state.fingerprint == initialization.fingerprint
        resumed_forward = _completed_directions["forward"]
        assert len(resumed_forward.points) == len(trace.directions[0].points)
        for resumed_point, expected_point in zip(
            resumed_forward.points, trace.directions[0].points, strict=True
        ):
            assert resumed_point.outer_step == expected_point.outer_step
            assert np.array_equal(
                resumed_point.coordinates_angstrom,
                expected_point.coordinates_angstrom,
            )
        _direction_callback(trace.directions[1])
        return trace

    published = campaign._run_and_publish_irc(
        run_root,
        route=route,
        _runner=resumed_runner,
    )

    assert resumed_calls == [("forward",)]
    assert checkpoint.read_bytes() == checkpoint_before
    assert published.trace.execution_contract == trace.execution_contract
    for observed, expected in zip(
        published.trace.directions, trace.directions, strict=True
    ):
        assert observed.sella_direction == expected.sella_direction
        assert len(observed.points) == len(expected.points)
    assert (run_root / route / "irc-restart" / "reverse" / "receipt.json").is_file()


@pytest.mark.parametrize(
    ("corruption", "match"),
    [
        ("deleted", "lacks shared Sella initialization"),
        ("tampered", "Sella H0 array hash mismatch"),
        ("replaced-ancestry", "initialization TS geometry does not match"),
        ("replaced-environment", "identity does not match qualified TS"),
        ("replaced-contract", "initialization execution contract drifted"),
    ],
)
def test_complete_irc_resume_rejects_invalid_shared_initialization_before_backend(
    monkeypatch, tmp_path: Path, corruption: str, match: str
):
    route = "h-co-1w-cside"
    run_root = tmp_path / corruption
    _, transition_state, trace = _qualified_ancestry(run_root, route=route)
    initialization = _sella_initialization_state(transition_state, trace)
    trace = _bind_trace_initialization(trace, initialization)

    def checkpointing_runner(*_args, _initialization_callback, **_kwargs):
        _initialization_callback(initialization)
        return trace

    published = campaign._run_and_publish_irc(
        run_root,
        route=route,
        _runner=checkpointing_runner,
    )
    route_root = run_root / route
    initialization_root = route_root / "irc-restart" / "initialization"
    initialization_path = initialization_root / "receipt.json"

    if corruption == "deleted":
        shutil.rmtree(initialization_root)
        corrupted_initialization = None
    elif corruption == "tampered":

        def tamper_hessian(receipt):
            receipt["initialization_state"]["H0"]["values"][0] += 1.0

        _rewrite_canonical_json(initialization_path, tamper_hessian)
        corrupted_initialization = initialization_path.read_bytes()
    else:
        if corruption == "replaced-ancestry":
            replacement_x0 = initialization.x0.copy()
            replacement_x0[0] += 0.01
            replacement_pes = dict(initialization.pes_current)
            replacement_pes["x"] = replacement_x0
            replacement_pes["state_hash"] = replacement_x0.tobytes()
            replacement = replace(
                initialization,
                x0=replacement_x0,
                pes_current=replacement_pes,
            )
        elif corruption == "replaced-environment":
            replacement = replace(initialization, settings_fingerprint="f" * 64)
        else:
            replacement = replace(
                initialization,
                execution_contract=replace(
                    initialization.execution_contract,
                    maximum_steps=initialization.execution_contract.maximum_steps + 1,
                ),
            )
        initialization_receipt = json.loads(initialization_path.read_text())
        initialization_receipt["initialization_fingerprint"] = replacement.fingerprint
        initialization_receipt["initialization_state"] = (
            campaign.quarry_ts._sella_irc_initialization_payload(replacement)
        )
        initialization_path.write_bytes(campaign._json_bytes(initialization_receipt))
        for name in ("forward", "reverse"):
            _rewrite_canonical_json(
                route_root / "irc-restart" / name / "receipt.json",
                lambda receipt: receipt.__setitem__(
                    "initialization_fingerprint", replacement.fingerprint
                ),
            )
            _rewrite_canonical_json(
                route_root / f"irc-{name}" / "receipt.json",
                lambda receipt: receipt.__setitem__(
                    "initialization_fingerprint", replacement.fingerprint
                ),
            )

        def replace_execution_initialization(receipt):
            receipt["initialization_fingerprint"] = replacement.fingerprint
            for direction in receipt["direction_receipts"].values():
                direction["initialization_fingerprint"] = replacement.fingerprint

        _rewrite_canonical_json(
            published.execution_receipt_path,
            replace_execution_initialization,
        )
        corrupted_initialization = initialization_path.read_bytes()

    protected = {
        path: path.read_bytes() for path in route_root.rglob("*") if path.is_file()
    }
    backend_calls = []

    def reject_runner(*_args, **_kwargs):
        backend_calls.append("runner")
        pytest.fail("invalid completed initialization reached IRC runner")

    def reject_calculator(*_args, **_kwargs):
        backend_calls.append("calculator")
        pytest.fail("invalid completed initialization reached calculator")

    monkeypatch.setattr(campaign.quarry_ts, "_trace_sella_irc_resume", reject_runner)
    monkeypatch.setattr(campaign.quarry_ts, "make_ase_calculator", reject_calculator)

    with pytest.raises(ValueError, match=match):
        campaign.run_and_publish_irc(run_root, route=route)

    assert backend_calls == []
    assert {
        path: path.read_bytes() for path in route_root.rglob("*") if path.is_file()
    } == protected
    if corrupted_initialization is None:
        assert not initialization_root.exists()
    else:
        assert initialization_path.read_bytes() == corrupted_initialization


@pytest.mark.parametrize(
    "corruption", ["mixed-run", "tampered-initialization", "unexpected-artifact"]
)
def test_irc_restart_rejects_corrupt_checkpoint_before_recompute(
    tmp_path: Path, corruption: str
):
    run_root = tmp_path / corruption
    _, transition_state, trace = _qualified_ancestry(run_root)
    initialization = _sella_initialization_state(transition_state, trace)
    trace = _bind_trace_initialization(trace, initialization)
    route = "h-co-1w-cside"

    def stop_after_forward(
        *_args,
        _completed_directions,
        _direction_callback,
        _initialization_callback,
        **_kwargs,
    ):
        assert not _completed_directions
        _initialization_callback(initialization)
        _direction_callback(trace.directions[0])
        raise RuntimeError("checkpoint only")

    with pytest.raises(RuntimeError, match="checkpoint only"):
        campaign._run_and_publish_irc(
            run_root,
            route=route,
            _runner=stop_after_forward,
        )
    restart_root = run_root / route / "irc-restart"
    if corruption == "mixed-run":
        _rewrite_canonical_json(
            restart_root / "forward" / "receipt.json",
            lambda receipt: receipt.__setitem__("irc_run_identity", "f" * 64),
        )
        match = "mixed run identities"
    elif corruption == "tampered-initialization":

        def tamper_hessian(receipt):
            receipt["initialization_state"]["H0"]["values"][0] += 1.0

        _rewrite_canonical_json(
            restart_root / "initialization" / "receipt.json", tamper_hessian
        )
        match = "Sella H0 array hash mismatch"
    else:
        (restart_root / "foreign").write_text("must be preserved")
        match = "unexpected or incomplete"

    with pytest.raises(ValueError, match=match):
        campaign._run_and_publish_irc(
            run_root,
            route=route,
            _runner=lambda *_args, **_kwargs: pytest.fail(
                "corrupt checkpoint reached IRC runner"
            ),
        )
    assert not (run_root / route / "irc-execution").exists()
    if corruption == "unexpected-artifact":
        assert (restart_root / "foreign").read_text() == "must be preserved"


def test_irc_checkpoint_destination_race_is_non_overwriting(
    monkeypatch, tmp_path: Path
):
    run_root = tmp_path / "run"
    _, transition_state, trace = _qualified_ancestry(run_root)
    initialization = _sella_initialization_state(transition_state, trace)
    trace = _bind_trace_initialization(trace, initialization)
    route = "h-co-1w-cside"
    destination = run_root / route / "irc-restart" / "forward"
    real_rename = campaign._renameat2_noreplace

    def race(source: Path, target: Path) -> None:
        if target.name == destination.name and target.parent.name == "irc-restart":
            target.mkdir()
            (target / "foreign").write_text("racer owns checkpoint")
        real_rename(source, target)

    monkeypatch.setattr(campaign, "_renameat2_noreplace", race)

    def runner(
        *_args,
        _completed_directions,
        _direction_callback,
        _initialization_callback,
        **_kwargs,
    ):
        assert not _completed_directions
        _initialization_callback(initialization)
        _direction_callback(trace.directions[0])
        return trace

    with pytest.raises(FileExistsError):
        campaign._run_and_publish_irc(run_root, route=route, _runner=runner)
    assert (destination / "foreign").read_text() == "racer owns checkpoint"
    assert not list(destination.parent.glob(".forward.*.tmp"))
    assert not (run_root / route / "irc-execution").exists()


def test_irc_restart_rejects_reverse_without_preceding_forward(
    tmp_path: Path,
):
    run_root = tmp_path / "run"
    _, _, trace = _qualified_ancestry(run_root)
    route = "h-co-1w-cside"

    def reverse_only(*_args, _completed_directions, _direction_callback, **_kwargs):
        assert not _completed_directions
        _direction_callback(trace.directions[1])
        raise RuntimeError("stop after impossible reverse-only checkpoint")

    with pytest.raises(RuntimeError, match="reverse-only"):
        campaign._run_and_publish_irc(
            run_root,
            route=route,
            _runner=reverse_only,
        )
    with pytest.raises(ValueError, match="direction completion order"):
        campaign._run_and_publish_irc(
            run_root,
            route=route,
            _runner=lambda *_args, **_kwargs: pytest.fail(
                "malformed restart reached IRC runner"
            ),
        )


def test_irc_restart_validates_completed_direction_before_checkpoint(tmp_path: Path):
    run_root = tmp_path / "run"
    _, _, trace = _qualified_ancestry(run_root)
    route = "h-co-1w-cside"
    forward = trace.directions[0]
    bad_terminal = replace(
        forward.points[-1],
        projected_fmax_ev_per_angstrom=campaign.BOUNDS["irc"][
            "outer_fmax_ev_per_angstrom"
        ],
    )
    malformed = replace(forward, points=(*forward.points[:-1], bad_terminal))

    def runner(*_args, _direction_callback, **_kwargs):
        _direction_callback(malformed)
        pytest.fail("invalid completed direction was checkpointed")

    with pytest.raises(ValueError, match="terminal IRC fmax must be strictly below"):
        campaign._run_and_publish_irc(run_root, route=route, _runner=runner)

    assert not (run_root / route / "irc-restart" / "forward").exists()


def test_irc_restart_recovers_owned_interrupted_root_publication(tmp_path: Path):
    run_root = tmp_path / "run"
    _, transition_state, trace = _qualified_ancestry(run_root)
    _, runner = _checkpointing_trace_runner(transition_state, trace)
    route_root = run_root / "h-co-1w-cside"
    stale = route_root / ".irc-restart.123.456.tmp"
    stale.mkdir()
    (stale / "receipt.json").write_text("{}\n")

    campaign._run_and_publish_irc(
        run_root,
        route="h-co-1w-cside",
        _runner=runner,
    )

    assert not stale.exists()
    assert (route_root / "irc-restart" / "forward" / "receipt.json").is_file()
    assert (route_root / "irc-restart" / "reverse" / "receipt.json").is_file()


def test_public_boundaries_build_native_backends_internally(
    monkeypatch, tmp_path: Path
):
    route = "h-co-1w-cside"
    qualification_root = tmp_path / "qualification"
    preflight = _preflight(qualification_root)
    transition_state, _ = _trace(route, "irc_back.xyz", "irc_fwd.xyz")
    reactant = _frozen_endpoint(route, "irc_back.xyz").coords
    product = _frozen_endpoint(route, "irc_fwd.xyz").coords
    native = _qualification_native_result(
        preflight, transition_state, reactant, product
    )
    qualification_calls = []

    def qualify(cluster, settings):
        qualification_calls.append((cluster, settings))
        return native

    monkeypatch.setattr(campaign, "native_cartesian_hessian", qualify)
    qualified = campaign.publish_transition_state_qualification(
        qualification_root, route=route
    )
    assert len(qualification_calls) == 1
    assert qualified.receipt["accepted"] is True

    hessian_root = tmp_path / "hessians"
    _, path = _authoritative_path(hessian_root)
    hessian_calls = []

    def hessian(cluster, settings):
        hessian_calls.append((cluster, settings))
        return _energy_reproducing_result(cluster, path)

    monkeypatch.setattr(campaign, "native_cartesian_hessian", hessian)
    published = campaign.publish_path_hessians(hessian_root, route=route)
    assert len(hessian_calls) == len(path.coordinates_angstrom)
    assert published.receipt["accepted"] is True


def test_public_resume_rejects_runtime_drift_before_backend(
    monkeypatch, tmp_path: Path
):
    run_root = tmp_path / "run"
    _qualified_ancestry(run_root)
    route = "h-co-1w-cside"
    monkeypatch.setattr(
        campaign,
        "_current_code_dependency_identity",
        lambda: {
            "git_sha": FIXED_GIT_SHA,
            "dependencies": {**FIXED_DEPENDENCIES, "sella": "drifted"},
            "executable_modules": FIXED_MODULE_MANIFEST,
            "python": campaign.platform.python_version(),
        },
    )
    monkeypatch.setattr(
        campaign.quarry_ts,
        "_trace_sella_irc_resume",
        lambda *_args, **_kwargs: pytest.fail("runtime drift reached Sella"),
    )
    with pytest.raises(ValueError, match="current dependency identity"):
        campaign.run_and_publish_irc(run_root, route=route)


def test_trusted_xyz_fails_closed_during_path_replacement(monkeypatch, tmp_path: Path):
    route = "h-co-1w-cside"
    source = BUNDLE_ROOT / route / "ts.xyz"
    candidate = tmp_path / "ts.xyz"
    candidate.write_bytes(source.read_bytes())
    replacement = tmp_path / "replacement.xyz"
    original_read = campaign.os.read
    raced = False

    def replace_path_after_first_read(descriptor, size):
        nonlocal raced
        data = original_read(descriptor, size)
        if data and not raced:
            raced = True
            replacement.write_bytes(b"attacker-controlled replacement\n")
            replacement.replace(candidate)
        return data

    monkeypatch.setattr(campaign.os, "read", replace_path_after_first_read)
    with pytest.raises(ValueError, match="changed while its snapshot was read"):
        campaign._trusted_xyz_snapshot(
            candidate,
            _route_template(route),
            name="race-safe-ts",
            expected_file_sha256=campaign.TRUSTED_TRANSITION_STATE_FILE_SHA256[route],
            expected_geometry_sha256=campaign.TRUSTED_CANONICAL_TS_GEOMETRY_SHA256[
                route
            ],
        )

    assert raced is True
    assert candidate.read_bytes() == b"attacker-controlled replacement\n"


def test_trusted_xyz_rejects_symlink_before_read(tmp_path: Path):
    route = "h-co-1w-cside"
    link = tmp_path / "ts.xyz"
    link.symlink_to(BUNDLE_ROOT / route / "ts.xyz")
    with pytest.raises(ValueError, match="trusted regular file"):
        campaign._trusted_xyz_snapshot(
            link,
            _route_template(route),
            name="symlinked-ts",
            expected_file_sha256=campaign.TRUSTED_TRANSITION_STATE_FILE_SHA256[route],
            expected_geometry_sha256=campaign.TRUSTED_CANONICAL_TS_GEOMETRY_SHA256[
                route
            ],
        )


def test_resumed_qualification_rejects_trusted_input_fingerprint_drift(tmp_path: Path):
    route = "h-co-1w-cside"
    run_root = tmp_path / "run"
    _qualified_ancestry(run_root, route=route)
    receipt_path = run_root / route / "ts-qualification" / "receipt.json"

    def tamper(receipt):
        receipt["trusted_input_fingerprints"]["reactant"]["file_sha256"] = "f" * 64

    _rewrite_canonical_json(receipt_path, tamper)
    with pytest.raises(ValueError, match="resumed TS qualification trusted input"):
        campaign._load_canonical_qualification(run_root, route)


def test_shared_irc_validator_rejects_restart_direction_divergence(tmp_path: Path):
    route = "h-co-1w-cside"
    run_root = tmp_path / "run"
    _authoritative_ancestry(run_root, route=route)
    restart_direction = run_root / route / "irc-restart" / "forward" / "receipt.json"

    def tamper(receipt):
        receipt["points"][-1]["electronic_energy_ev"] -= 0.125

    _rewrite_canonical_json(restart_direction, tamper)
    ancestry = campaign._load_canonical_qualification(run_root, route)
    with pytest.raises(ValueError, match="restart/final canonical forward receipt"):
        campaign._validate_canonical_irc_receipts(ancestry)


def test_executable_module_manifest_hashes_concrete_repository_leaf():
    qm_root = Path(campaign.__file__).resolve().parents[1]
    code = """
import json, sys
sys.path.insert(0, '.')
from scripts import d2c_sct_campaign as campaign
campaign.EXECUTABLE_MODULES = {'quarry.clusters': 'repository'}
print(json.dumps(campaign._executable_module_manifest()))
"""
    completed = subprocess.run(
        [str(qm_root / ".venv" / "bin" / "python"), "-c", code],
        cwd=qm_root,
        check=True,
        capture_output=True,
        text=True,
    )
    manifest = json.loads(completed.stdout)
    record = manifest["quarry.clusters"]
    origin = Path(record["origin"])
    assert (
        origin
        == Path(campaign.__file__).resolve().parents[1] / "quarry" / "clusters.py"
    )
    assert record["byte_count"] == origin.stat().st_size
    assert record["sha256"] == hashlib.sha256(origin.read_bytes()).hexdigest()
    assert record["execution_identity"]["kind"] == "python-source"


def test_executable_manifest_covers_direct_sella_scipy_pyscf_and_gpu_leaves():
    required = {
        "numpy.linalg._umath_linalg",
        "scipy.linalg._flapack",
        "scipy.linalg._decomp",
        "scipy.integrate._ivp.base",
        "scipy.integrate._ivp.common",
        "scipy.integrate._dop",
        "scipy.integrate._ode",
        "scipy.integrate._odepack",
        "scipy.integrate._vode",
        "scipy.sparse.linalg._interface",
        "pyscf.dft.rks",
        "pyscf.scf.uhf",
        "pyscf.grad.uhf",
        "pyscf.hessian.rhf",
        "pyscf.lib.misc",
        "sella.eigensolvers",
        "sella.hessian_update",
        "sella.internal",
        "sella.linalg",
        "gpu4pyscf.dft.rks",
        "gpu4pyscf.scf.uhf",
        "gpu4pyscf.grad.uhf",
        "gpu4pyscf.hessian.rhf",
    }

    assert required <= campaign.EXECUTABLE_MODULES.keys()


def test_actual_nested_venv_manifest_attests_source_and_loaded_native_image():
    qm_root = Path(campaign.__file__).resolve().parents[1]
    python = qm_root / ".venv" / "bin" / "python"
    assert python.is_file()
    code = """
import json, sys
sys.path.insert(0, '.')
from scripts import d2c_sct_campaign as campaign
campaign.EXECUTABLE_MODULES = {
    'scipy.linalg._decomp': 'third-party',
    'scipy.linalg._flapack': 'third-party',
    'numpy._core._multiarray_umath': 'third-party',
    'numpy.linalg._umath_linalg': 'third-party',
    'pyscf.grad.rhf': 'third-party',
    'pyscf.lib.misc': 'third-party',
    'pyscf.scf.hf': 'third-party',
}
print(json.dumps({
    'modules': campaign._executable_module_manifest(),
    'native_payloads': campaign._native_payload_manifest(),
}, sort_keys=True))
"""
    completed = subprocess.run(
        [str(python), "-c", code],
        cwd=qm_root,
        check=True,
        capture_output=True,
        text=True,
    )
    result = json.loads(completed.stdout)
    manifest = result["modules"]

    assert Path(manifest["scipy.linalg._decomp"]["origin"]).is_relative_to(
        qm_root / ".venv"
    )
    assert manifest["scipy.linalg._decomp"]["execution_identity"]["kind"] == (
        "python-source"
    )
    assert manifest["numpy._core._multiarray_umath"]["execution_identity"]["kind"] == (
        "native-extension"
    )
    assert manifest["pyscf.grad.rhf"]["execution_identity"]["kind"] == "python-source"
    assert manifest["pyscf.scf.hf"]["execution_identity"]["kind"] == "python-source"
    for module_name in ("scipy.linalg._flapack", "numpy.linalg._umath_linalg"):
        assert manifest[module_name]["execution_identity"]["kind"] == "native-extension"
        assert Path(manifest[module_name]["origin"]).is_relative_to(qm_root / ".venv")
    native_payloads = result["native_payloads"]
    assert native_payloads
    assert all(
        Path(record["origin"]).is_relative_to(qm_root / ".venv")
        for record in native_payloads.values()
    )
    present = {Path(name).name for name in native_payloads}
    assert present >= campaign.REQUIRED_NATIVE_PAYLOAD_BASENAMES
    for required_name in campaign.REQUIRED_NATIVE_PAYLOAD_BASENAMES:
        record = next(
            record
            for name, record in native_payloads.items()
            if Path(name).name == required_name
        )
        origin = Path(record["origin"])
        status = origin.stat()
        assert record["sha256"] == hashlib.sha256(origin.read_bytes()).hexdigest()
        assert (record["mapped_device"], record["mapped_inode"]) == (
            status.st_dev,
            status.st_ino,
        )


def test_relative_module_launch_attests_manifest_and_restores_profiler():
    qm_root = Path(campaign.__file__).resolve().parents[1]
    python = qm_root / ".venv" / "bin" / "python"
    wrapper = """
import json
import pathlib
import sys
import types

def prior_profile(_frame, _event, _arg):
    return None

relative_source = 'scripts/d2c_sct_campaign.py'
raw = pathlib.Path(relative_source).read_bytes()
module_name = 'd2c_relative_launch'
module = types.ModuleType(module_name)
module.__file__ = relative_source
module.__package__ = ''
sys.modules[module_name] = module
sys.setprofile(prior_profile)
exec(compile(raw, relative_source, 'exec'), module.__dict__)
assert sys.getprofile() is prior_profile
origin = pathlib.Path(relative_source).resolve()
identity = module._python_source_execution_identity(
    module,
    module_name,
    origin,
    raw,
    trusted_source_roots=(
        origin.parents[1],
        pathlib.Path(sys.prefix).resolve(),
        pathlib.Path(sys.base_prefix).resolve(),
    ),
)
print(json.dumps(identity, sort_keys=True))
print('PROFILE_RESTORED')
"""
    completed = subprocess.run(
        [str(python), "-c", wrapper],
        cwd=qm_root,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr

    output = completed.stdout.splitlines()
    identity = json.loads(output[0])
    assert identity["kind"] == "python-source"
    assert identity["loaded_module_code_sha256"]
    assert identity["loaded_state_sha256"]
    assert output[-1] == "PROFILE_RESTORED"


@pytest.mark.parametrize(
    ("module_name", "before", "after"),
    [
        (
            "d2c_loaded_global_fixture",
            "VALUE = 1\ndef value():\n    return VALUE\n",
            "VALUE = 2\ndef value():\n    return VALUE\n",
        ),
        (
            "d2c_loaded_default_fixture",
            "def value(setting=1):\n    return setting\n",
            "def value(setting=2):\n    return setting\n",
        ),
        (
            "d2c_loaded_class_fixture",
            "class Config:\n    SETTING = 1\ndef value():\n    return Config.SETTING\n",
            "class Config:\n    SETTING = 2\ndef value():\n    return Config.SETTING\n",
        ),
    ],
)
def test_manifest_rejects_module_level_state_replacement_with_restored_timestamp(
    monkeypatch, tmp_path: Path, module_name: str, before: str, after: str
):
    source = tmp_path / f"{module_name}.py"
    source.write_text(before)
    monkeypatch.syspath_prepend(str(tmp_path))
    monkeypatch.setattr(campaign.sys, "prefix", str(tmp_path))
    monkeypatch.setattr(campaign, "EXECUTABLE_MODULES", {module_name: "third-party"})
    imported = campaign._import_executable_module(module_name)
    original = source.stat()
    replacement = tmp_path / "replacement.py"
    replacement.write_text(after)
    assert replacement.stat().st_size == original.st_size
    replacement.replace(source)
    os.utime(source, ns=(original.st_atime_ns, original.st_mtime_ns))
    assert imported.value() == 1

    try:
        with pytest.raises(
            RuntimeError, match="loaded Python (?:code|module state) disagrees"
        ):
            campaign._executable_module_manifest()
    finally:
        campaign._OBSERVED_MODULE_CODE.pop(module_name, None)
        sys.modules.pop(module_name, None)


def test_manifest_rejects_mutated_campaign_dft_settings(monkeypatch):
    source = Path(campaign.__file__).resolve()
    raw = source.read_bytes()
    monkeypatch.setitem(campaign.DFT_SETTINGS, "xc", "resident-only-drift")

    with pytest.raises(
        RuntimeError,
        match=(
            r"loaded Python module state disagrees with source: "
            r".*global\.DFT_SETTINGS"
        ),
    ):
        campaign._python_source_execution_identity(
            campaign,
            "scripts.d2c_sct_campaign",
            source,
            raw,
            trusted_source_roots=(
                source.parents[1],
                Path(sys.prefix).resolve(),
                Path(sys.base_prefix).resolve(),
            ),
        )


def test_manifest_rejects_mutated_literal_list_global(monkeypatch, tmp_path: Path):
    module_name = "d2c_loaded_list_fixture"
    source = tmp_path / f"{module_name}.py"
    source.write_text("VALUE = [1]\ndef value():\n    return VALUE\n")
    monkeypatch.syspath_prepend(str(tmp_path))
    monkeypatch.setattr(campaign.sys, "prefix", str(tmp_path))
    monkeypatch.setattr(campaign, "EXECUTABLE_MODULES", {module_name: "third-party"})
    imported = campaign._import_executable_module(module_name)

    try:
        first = campaign._executable_module_manifest()
        imported.VALUE.append(2)
        with pytest.raises(
            RuntimeError,
            match=r"loaded Python module state disagrees with source: .*global\.VALUE",
        ):
            campaign._executable_module_manifest()
        assert first[module_name]["execution_identity"]["loaded_state_sha256"]
    finally:
        campaign._OBSERVED_MODULE_CODE.pop(module_name, None)
        sys.modules.pop(module_name, None)


def test_manifest_rejects_source_replacement_after_module_was_loaded(
    monkeypatch, tmp_path: Path
):
    module_name = "d2c_loaded_replacement_fixture"
    source = tmp_path / f"{module_name}.py"
    source.write_text("def value():\n    return 1\n")
    monkeypatch.syspath_prepend(str(tmp_path))
    monkeypatch.setattr(campaign.sys, "prefix", str(tmp_path))
    monkeypatch.setattr(campaign, "EXECUTABLE_MODULES", {module_name: "third-party"})
    imported = campaign._import_executable_module(module_name)
    original = source.stat()
    replacement = tmp_path / "replacement.py"
    replacement.write_text("def value():\n    return 2\n")
    replacement.replace(source)
    os.utime(source, ns=(original.st_atime_ns, original.st_mtime_ns))
    assert imported.value() == 1

    try:
        with pytest.raises(
            RuntimeError, match="loaded Python code disagrees with source"
        ):
            campaign._executable_module_manifest()
    finally:
        sys.modules.pop(module_name, None)


def test_manifest_rejects_loaded_python_code_without_source_origin(
    monkeypatch, tmp_path: Path
):
    module_name = "d2c_unproven_loaded_fixture"
    source = tmp_path / f"{module_name}.py"
    source.write_text("def value():\n    return 1\n")
    monkeypatch.syspath_prepend(str(tmp_path))
    monkeypatch.setattr(campaign.sys, "prefix", str(tmp_path))
    monkeypatch.setattr(campaign, "EXECUTABLE_MODULES", {module_name: "third-party"})
    imported = campaign._import_executable_module(module_name)
    exec("def unproven():\n    return 2\n", imported.__dict__)

    try:
        with pytest.raises(
            RuntimeError, match="loaded Python code has no provable origin"
        ):
            campaign._executable_module_manifest()
    finally:
        sys.modules.pop(module_name, None)


def test_manifest_preimports_declared_modules_before_attestation_lazy_import(
    monkeypatch, tmp_path: Path
):
    first_name = "d2c_preimport_first_fixture"
    second_name = "d2c_preimport_second_fixture"
    for module_name in (first_name, second_name):
        (tmp_path / f"{module_name}.py").write_text(
            f"VALUE = {module_name!r}\ndef value():\n    return VALUE\n"
        )
    monkeypatch.syspath_prepend(str(tmp_path))
    monkeypatch.setattr(campaign.sys, "prefix", str(tmp_path))
    monkeypatch.setattr(
        campaign,
        "EXECUTABLE_MODULES",
        {first_name: "third-party", second_name: "third-party"},
    )
    original_identity = campaign._python_source_execution_identity

    def identity_with_lazy_import(module, module_name, origin, raw, **kwargs):
        if module_name == first_name:
            campaign.importlib.import_module(second_name)
            assert second_name in campaign._OBSERVED_MODULE_CODE
        return original_identity(module, module_name, origin, raw, **kwargs)

    monkeypatch.setattr(
        campaign, "_python_source_execution_identity", identity_with_lazy_import
    )
    try:
        manifest = campaign._executable_module_manifest()
        assert set(manifest) == {first_name, second_name}
    finally:
        for module_name in (first_name, second_name):
            campaign._OBSERVED_MODULE_CODE.pop(module_name, None)
            sys.modules.pop(module_name, None)


def test_manifest_fails_closed_for_python_module_loaded_before_execution_capture(
    monkeypatch, tmp_path: Path
):
    module_name = "d2c_uncaptured_loaded_fixture"
    source = tmp_path / f"{module_name}.py"
    source.write_text("def value():\n    return 1\n")
    monkeypatch.syspath_prepend(str(tmp_path))
    monkeypatch.setattr(campaign.sys, "prefix", str(tmp_path))
    monkeypatch.setattr(campaign, "EXECUTABLE_MODULES", {module_name: "third-party"})
    imported = campaign.importlib.import_module(module_name)
    assert imported.value() == 1

    try:
        with pytest.raises(
            RuntimeError, match="module-level code has no resident proof"
        ):
            campaign._executable_module_manifest()
    finally:
        sys.modules.pop(module_name, None)


def test_manifest_rejects_loaded_code_outside_trusted_source_roots(
    monkeypatch, tmp_path: Path
):
    module_name = "d2c_external_loaded_fixture"
    environment = tmp_path / "venv"
    environment.mkdir()
    source = environment / f"{module_name}.py"
    source.write_text("def value():\n    return 1\n")
    external = tmp_path / "external.py"
    external.write_text("def injected():\n    return 2\n")
    monkeypatch.syspath_prepend(str(environment))
    monkeypatch.setattr(campaign.sys, "prefix", str(environment))
    monkeypatch.setattr(campaign, "EXECUTABLE_MODULES", {module_name: "third-party"})
    imported = campaign._import_executable_module(module_name)
    namespace = {"__name__": module_name}
    exec(compile(external.read_text(), str(external), "exec"), namespace)
    imported.__dict__["injected"] = namespace["injected"]

    try:
        with pytest.raises(RuntimeError, match="escaped trusted source roots"):
            campaign._executable_module_manifest()
    finally:
        sys.modules.pop(module_name, None)


def test_native_execution_identity_rejects_replaced_mapped_file(tmp_path: Path):
    origin = tmp_path / "fixture.so"
    origin.write_bytes(b"mapped native fixture")
    with (
        origin.open("rb") as source,
        mmap.mmap(source.fileno(), 0, access=mmap.ACCESS_READ) as mapped,
    ):
        assert mapped[:6] == b"mapped"
        identity = campaign._native_extension_execution_identity("fixture", origin)
        assert identity["kind"] == "native-extension"
        replacement = tmp_path / "replacement.so"
        replacement.write_bytes(b"replacement fixture")
        replacement.replace(origin)

        with pytest.raises(
            RuntimeError,
            match="loaded native executable file identity cannot be proven",
        ):
            campaign._native_extension_execution_identity("fixture", origin)


def test_executable_module_manifest_rejects_import_origin_escape(
    monkeypatch, tmp_path: Path
):
    escaped = tmp_path / "clusters.py"
    escaped.write_text("raise RuntimeError('not the trusted module')\n")

    class FakeSpec:
        origin = str(escaped)

    class FakeModule:
        __spec__ = FakeSpec()
        __file__ = str(escaped)

    monkeypatch.setattr(
        campaign, "EXECUTABLE_MODULES", {"quarry.clusters": "repository"}
    )
    fake_module = FakeModule()
    monkeypatch.setitem(campaign.sys.modules, "quarry.clusters", fake_module)
    with pytest.raises(
        RuntimeError, match="repository executable module origin drifted"
    ):
        campaign._executable_module_manifest()


@pytest.mark.parametrize("field", ["origin", "sha256"])
def test_production_boundary_rejects_executable_module_origin_or_hash_drift(
    monkeypatch, tmp_path: Path, field: str
):
    route = "h-co-1w-cside"
    run_root = tmp_path / "run"
    preflight = _preflight(run_root)
    modules = {name: dict(record) for name, record in FIXED_MODULE_MANIFEST.items()}
    target = "sella.peswrapper"
    modules[target][field] = (
        "/opt/attacker/sella/peswrapper.py" if field == "origin" else "f" * 64
    )
    monkeypatch.setattr(
        campaign,
        "_current_code_dependency_identity",
        lambda: {
            "git_sha": FIXED_GIT_SHA,
            "dependencies": dict(FIXED_DEPENDENCIES),
            "executable_modules": modules,
            "python": preflight["campaign"]["python"],
        },
    )

    with pytest.raises(ValueError, match="current executable module identity"):
        campaign._validate_production_boundary(run_root, route)


def test_native_payload_manifest_rejects_missing_required_numerical_library():
    incomplete = dict(FIXED_NATIVE_PAYLOAD_MANIFEST)
    target = next(name for name in incomplete if Path(name).name == "libcint.so")
    del incomplete[target]

    with pytest.raises(
        ValueError, match="omits required numerical libraries: libcint.so"
    ):
        campaign._validated_native_payload_manifest(incomplete)


def test_native_payload_hash_rejects_same_size_in_place_rewrite(
    monkeypatch, tmp_path: Path
):
    payload = tmp_path / "mapped-native.so"
    payload.write_bytes(b"a" * (2 * 1024 * 1024))
    status = payload.stat(follow_symlinks=False)
    original_read = campaign.os.read
    rewritten = False

    def rewrite_after_first_read(descriptor: int, size: int) -> bytes:
        nonlocal rewritten
        data = original_read(descriptor, size)
        if data and not rewritten:
            rewritten = True
            payload.write_bytes(b"b" * status.st_size)
            os.utime(
                payload,
                ns=(status.st_atime_ns, status.st_mtime_ns + 1_000_000_000),
            )
        return data

    monkeypatch.setattr(campaign.os, "read", rewrite_after_first_read)

    with pytest.raises(RuntimeError, match="changed while hashed"):
        campaign._hash_mapped_native_payload(
            payload,
            mapped_device=status.st_dev,
            mapped_inode=status.st_ino,
        )
    assert rewritten is True
    assert payload.stat().st_size == status.st_size


def test_production_boundary_rejects_native_payload_hash_drift(
    monkeypatch, tmp_path: Path
):
    route = "h-co-1w-cside"
    run_root = tmp_path / "run"
    preflight = _preflight(run_root)
    native_payloads = {
        name: dict(record) for name, record in FIXED_NATIVE_PAYLOAD_MANIFEST.items()
    }
    target = next(name for name in native_payloads if Path(name).name == "libxc.so")
    native_payloads[target]["sha256"] = "f" * 64
    monkeypatch.setattr(
        campaign,
        "_current_code_dependency_identity",
        lambda: {
            "git_sha": FIXED_GIT_SHA,
            "dependencies": dict(FIXED_DEPENDENCIES),
            "executable_modules": FIXED_MODULE_MANIFEST,
            "native_payloads": native_payloads,
            "python": preflight["campaign"]["python"],
        },
    )

    with pytest.raises(ValueError, match="current native payload identity"):
        campaign._validate_production_boundary(run_root, route)


def test_native_ts_hessian_evaluation_and_publication_share_one_route_claim(
    monkeypatch, tmp_path: Path
):
    route = "h-co-1w-cside"
    run_root = tmp_path / "run"
    preflight = _preflight(run_root)
    active = False
    claim_calls = 0
    evaluator_observed_claim = False

    @contextmanager
    def nonreentrant_claim(route_root):
        nonlocal active, claim_calls
        assert active is False, "nested route claim would deadlock"
        route_root.mkdir(parents=True, exist_ok=True)
        claim_calls += 1
        active = True
        try:
            yield
        finally:
            active = False

    def evaluate(transition_state, _settings):
        nonlocal evaluator_observed_claim
        evaluator_observed_claim = active
        reactant = _frozen_endpoint(route, "irc_back.xyz").coords
        product = _frozen_endpoint(route, "irc_fwd.xyz").coords
        return _qualification_native_result(
            preflight, transition_state, reactant, product
        )

    monkeypatch.setattr(campaign, "_exclusive_route_claim", nonreentrant_claim)
    monkeypatch.setattr(campaign, "native_cartesian_hessian", evaluate)
    published = campaign.publish_transition_state_qualification(run_root, route=route)

    assert evaluator_observed_claim is True
    assert claim_calls == 1
    assert published.receipt["accepted"] is True


def test_oversized_initialization_receipt_rejects_before_json_parsing(tmp_path: Path):
    route = "h-co-1w-cside"
    run_root = tmp_path / "oversized-initialization"
    _, transition_state, trace = _qualified_ancestry(run_root, route=route)
    initialization = _sella_initialization_state(transition_state, trace)
    trace = _bind_trace_initialization(trace, initialization)

    def checkpointing_runner(*_args, _initialization_callback, **_kwargs):
        _initialization_callback(initialization)
        return trace

    campaign._run_and_publish_irc(
        run_root,
        route=route,
        _runner=checkpointing_runner,
    )
    receipt_path = run_root / route / "irc-restart" / "initialization" / "receipt.json"
    maximum = campaign._irc_initialization_receipt_maximum_bytes(
        len(transition_state.symbols)
    )
    receipt_path.write_bytes(b" " * (maximum + 1))
    ancestry = campaign._load_canonical_qualification(run_root, route)

    with pytest.raises(ValueError, match=rf"no larger than {maximum} bytes"):
        campaign._load_irc_restart(ancestry)


def test_route_replacement_during_native_evaluation_fails_closed(
    monkeypatch, tmp_path: Path
):
    route = "h-co-1w-cside"
    run_root = tmp_path / "route-replacement"
    preflight = _preflight(run_root)
    route_root = run_root / route
    displaced = run_root / f"{route}.displaced"
    evaluator_calls = 0

    def evaluate(transition_state, _settings):
        nonlocal evaluator_calls
        evaluator_calls += 1
        route_root.rename(displaced)
        route_root.mkdir()
        reactant = _frozen_endpoint(route, "irc_back.xyz").coords
        product = _frozen_endpoint(route, "irc_fwd.xyz").coords
        return _qualification_native_result(
            preflight, transition_state, reactant, product
        )

    monkeypatch.setattr(campaign, "native_cartesian_hessian", evaluate)

    with pytest.raises(RuntimeError, match="route directory identity changed"):
        campaign.publish_transition_state_qualification(run_root, route=route)

    assert evaluator_calls == 1
    assert not (route_root / "ts-qualification").exists()
    assert not (route_root / "ts-qualification" / "receipt.json").exists()
    assert not (displaced / "ts-qualification").exists()
    assert not (displaced / "ts-qualification" / "receipt.json").exists()
    assert (displaced / ".route.lock").is_file()


def test_route_claim_rejects_lock_file_replacement(tmp_path: Path):
    route_root = tmp_path / "route"
    displaced_lock = tmp_path / "displaced.lock"

    with (
        pytest.raises(RuntimeError, match="route lock identity changed"),
        campaign._exclusive_route_claim(route_root),
    ):
        (route_root / ".route.lock").replace(displaced_lock)
        (route_root / ".route.lock").write_bytes(b"replacement")


@pytest.mark.parametrize("replacement", ["lock", "route"])
def test_kernel_route_claim_blocks_second_process_after_namespace_replacement(
    tmp_path: Path, replacement: str
):
    context = multiprocessing.get_context("spawn")
    route_root = tmp_path / "route"
    entered = context.Event()
    release = context.Event()
    first_outcome = context.Queue()
    second_outcome = context.Queue()
    first = context.Process(
        target=_hold_route_claim,
        args=(str(route_root), entered, release, first_outcome),
    )
    second = context.Process(
        target=_attempt_route_claim,
        args=(str(route_root), second_outcome),
    )
    first.start()
    try:
        assert entered.wait(10), "first process did not enter route claim"
        if replacement == "lock":
            (route_root / ".route.lock").replace(tmp_path / "displaced.lock")
            (route_root / ".route.lock").write_bytes(b"replacement")
        else:
            route_root.rename(tmp_path / "displaced-route")
            route_root.mkdir()
        second.start()
        second.join(10)
        assert not second.is_alive(), "second process blocked instead of failing closed"
        assert second.exitcode == 0
        blocked = second_outcome.get(timeout=2)
        assert blocked[0] == "blocked"
        assert blocked[1] == "RuntimeError"
        assert "already claimed by another process" in blocked[2]
    finally:
        release.set()
        if second.pid is not None and second.is_alive():
            second.terminate()
        if first.pid is not None:
            first.join(10)
            if first.is_alive():
                first.terminate()
                first.join(5)
    assert first.exitcode == 0
    failed = first_outcome.get(timeout=2)
    assert failed[0] == "failed"
    assert failed[1] == "RuntimeError"
    expected_drift = (
        "route lock identity changed"
        if replacement == "lock"
        else "route directory identity changed"
    )
    assert expected_drift in failed[2]
