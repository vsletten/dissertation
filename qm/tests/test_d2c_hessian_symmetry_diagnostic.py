from __future__ import annotations

import fcntl
import hashlib
import json
import os
import shutil
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from quarry.clusters import Cluster
from quarry.pipeline import DftSettings
from scripts import d2c_hessian_symmetry_diagnostic as diagnostic


def _identity():
    execution = {
        "kind": "python-source",
        "loaded_module_code_sha256": "1" * 64,
        "loaded_code_sha256": "2" * 64,
        "loaded_state_sha256": "3" * 64,
        "source_files": {"/test/diagnostic.py": "4" * 64},
    }
    return {
        "git_sha": "a" * 40,
        "script_sha256": "b" * 64,
        "diagnostic_execution_identity": execution,
        "diagnostic_execution_identity_sha256": diagnostic.campaign._canonical_hash(
            execution
        ),
    }


def test_run_publishes_failed_case_before_reraising(tmp_path, monkeypatch):
    preflight = {
        "identity": {"campaign": "test"},
        "routes": {diagnostic.ROUTE: {"masses_amu": [1.0]}},
    }
    snapshots = {
        "transition_state": object(),
        "reactant": SimpleNamespace(coords=np.zeros((1, 3))),
        "product": SimpleNamespace(coords=np.zeros((1, 3))),
    }
    monkeypatch.setattr(
        diagnostic,
        "_git_identity",
        lambda: {"git_sha": "a" * 40, "script_sha256": "b" * 64},
    )
    monkeypatch.setattr(
        diagnostic.campaign,
        "_validate_production_boundary",
        lambda *_args: (preflight, "c" * 64),
    )
    monkeypatch.setattr(
        diagnostic.campaign,
        "_canonical_dft_settings",
        lambda _preflight: (object(), {}),
    )
    monkeypatch.setattr(
        diagnostic.campaign,
        "_trusted_input_fingerprints_from_route_record",
        lambda _route: {},
    )
    monkeypatch.setattr(
        diagnostic.campaign,
        "_trusted_route_input_snapshots",
        lambda *_args: snapshots,
    )
    monkeypatch.setattr(
        diagnostic,
        "reactions",
        lambda **_kwargs: {diagnostic.ROUTE: SimpleNamespace(cluster=object())},
    )

    def fail_case(*_args):
        raise RuntimeError("synthetic Hessian failure")

    monkeypatch.setattr(diagnostic, "_evaluate_case", fail_case)
    output_root = tmp_path / "diagnostic"

    with pytest.raises(RuntimeError, match="synthetic Hessian failure"):
        diagnostic.run(tmp_path / "preflight", output_root)

    status = json.loads((output_root / "status.json").read_text())
    assert status["state"] == "failed"
    assert status["current_case"] == diagnostic.CASES[0]["name"]
    assert status["failed_case"] == diagnostic.CASES[0]["name"]
    assert status["completed_cases"] == []
    assert status["error"] == "RuntimeError: synthetic Hessian failure"
    assert status["finished_utc"].endswith("Z")
    assert "diagnostic_execution_identity" not in status
    assert not (output_root / "receipt.json").exists()
    terminal = json.loads((output_root / diagnostic.TERMINAL).read_text())
    assert terminal["state"] == "failed"
    assert terminal["current_case"] == diagnostic.CASES[0]["name"]
    assert terminal["detail"] == "RuntimeError: synthetic Hessian failure"


def test_finalize_if_running_publishes_idempotent_dead_man_receipt(tmp_path):
    output_root = tmp_path / "diagnostic"
    output_root.mkdir()
    (output_root / "status.json").write_bytes(
        diagnostic._json_bytes(
            {
                "schema": diagnostic.SCHEMA,
                "state": "running",
                "route": diagnostic.ROUTE,
                "git_sha": "a" * 40,
                "script_sha256": "b" * 64,
                "completed_cases": [diagnostic.CASES[0]["name"]],
                "current_case": diagnostic.CASES[1]["name"],
            }
        )
    )

    first = diagnostic.finalize_if_running(output_root)
    second = diagnostic.finalize_if_running(output_root)

    assert first == second
    assert first["state"] == "failed"
    assert first["current_case"] == diagnostic.CASES[1]["name"]
    assert "bounded systemd unit ended" in first["detail"]
    status = json.loads((output_root / "status.json").read_text())
    assert status["state"] == "failed"
    assert status["error"] == first["detail"]


def _point_record(*, density_sha: str, fmax: float | None = None, s2: float = 0.76):
    record = {
        "backend": "pyscf-cpu",
        "electronic_hartree": -10.0,
        "scf_converged": True,
        "grid_point_count": 199560,
        "spin_2s": 1,
        "spin_square": [s2, 2.0],
        "gradient_grid_response": True,
        "density_initial_guess_sha256": density_sha,
        "density_final_sha256": density_sha,
        "geometry_fingerprint": "f" * 64,
        "elapsed_seconds": 1.0,
    }
    if fmax is not None:
        record["physical_fmax_ev_per_angstrom"] = fmax
    return record


def test_five_point_stencil_has_exact_orientation_and_72_displacements():
    plan = diagnostic._displacement_plan(6)
    assert len(plan) == 72
    assert [item["step_multiplier"] for item in plan[:4]] == [-2, -1, 1, 2]
    assert plan[0]["coordinate_index"] == 0
    assert plan[4]["coordinate_index"] == 1
    assert plan[-1] == {
        "key": "coordinate-017-p2",
        "coordinate_index": 17,
        "atom_index": 5,
        "axis_index": 2,
        "step_multiplier": 2,
        "displacement_bohr": 0.02,
    }

    dimension = 18
    linear = np.arange(dimension * dimension, dtype=float).reshape(dimension, dimension)
    cubic = np.flipud(linear + 1.0)
    gradients = {}
    for column in range(dimension):
        for multiplier in (-2, -1, 1, 2):
            displacement = 0.01 * multiplier
            gradients[(column, multiplier)] = (
                linear[:, column] * displacement + cubic[:, column] * displacement**3
            )

    h_h, h_2h, richardson = diagnostic._finite_difference_matrices(gradients, dimension)

    assert h_h == pytest.approx(linear + cubic * 0.01**2)
    assert h_2h == pytest.approx(linear + cubic * 0.02**2)
    assert richardson == pytest.approx(linear, abs=2.0e-12)
    assert richardson[3, 7] == pytest.approx(linear[3, 7])


def test_gradient_evaluator_uses_exact_settings_grid_response_and_shared_dm0(
    monkeypatch,
):
    cluster = Cluster(
        "test",
        ["H", "H", "H"],
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        spin=1,
    )
    base = DftSettings(
        xc="pwb6k",
        basis="def2-svp",
        dispersion="d3bj",
        density_fit=True,
        use_gpu=True,
    )
    settings = diagnostic._fd_settings(base)
    assert settings.use_gpu is False
    assert settings.grid_level == 5

    kernel_dm0 = []
    gradients = []

    class FakeGradient:
        grid_response = False

        def kernel(self):
            gradients.append(self.grid_response)
            return np.zeros((3, 3))

    class FakeMeanField:
        converged = True
        grids = SimpleNamespace(level=5, coords=np.empty((199560, 0)))

        def __init__(self):
            self.conv_tol = None
            self.max_cycle = None

        def kernel(self, dm0=None):
            kernel_dm0.append(dm0)
            assert self.conv_tol == 1.0e-12
            assert self.max_cycle == 150
            return -10.0

        def spin_square(self):
            return 0.76, 2.0

        def make_rdm1(self):
            return np.arange(8, dtype=float).reshape(2, 2, 2)

    monkeypatch.setattr(diagnostic, "build_mol", lambda candidate, _: candidate)
    monkeypatch.setattr(diagnostic, "_make_scf", lambda *_: FakeMeanField())
    monkeypatch.setattr(diagnostic, "_gradient_method", lambda *_: FakeGradient())

    _, _, density = diagnostic._evaluate_fd_gradient(
        cluster, settings, central_density=None
    )
    assert density is not None
    for _ in range(72):
        diagnostic._evaluate_fd_gradient(cluster, settings, central_density=density)

    assert len(kernel_dm0) == 73
    assert kernel_dm0[0] is None
    assert all(value is density for value in kernel_dm0[1:])
    assert gradients == [True] * 73


@pytest.mark.parametrize(
    ("mutator", "match"),
    [
        (lambda record: record.update(grid_point_count=199559), "grid point count"),
        (lambda record: record.update(spin_square=[0.739999, 2.0]), "S2"),
        (lambda record: record.update(spin_square=[0.80, 2.0]), None),
        (lambda record: record.update(physical_fmax_ev_per_angstrom=0.02), "strictly"),
        (lambda record: record.update(electronic_hartree=float("nan")), "finite"),
    ],
)
def test_point_gates_fail_closed_at_exact_boundaries(mutator, match):
    record = _point_record(density_sha="a" * 64, fmax=0.01)
    mutator(record)
    if match is None:
        assert diagnostic._validate_point_record(
            record, center_spin_square=None, center=True
        ) == pytest.approx(0.80)
    else:
        with pytest.raises(ValueError, match=match):
            diagnostic._validate_point_record(
                record, center_spin_square=None, center=True
            )


def test_displaced_s2_delta_is_strict_but_s2_range_is_inclusive(monkeypatch):
    monkeypatch.setitem(
        diagnostic.FINITE_DIFFERENCE_CONTRACT["thresholds"],
        "displaced_vs_center_spin_square_exclusive_maximum_delta",
        0.125,
    )
    record = _point_record(density_sha="a" * 64, s2=0.75)
    with pytest.raises(diagnostic.ScientificRejection, match="strictly"):
        diagnostic._validate_point_record(
            record, center_spin_square=0.625, center=False
        )
    record["spin_square"][0] = 0.74
    assert diagnostic._validate_point_record(
        record, center_spin_square=0.749, center=False
    ) == pytest.approx(0.74)


@pytest.mark.parametrize(
    "eigenvalues",
    [
        [-1.0e-8, 1.0e-5],
        [-1.0e-4, 1.0e-8],
        [-1.0e-4, 5.0e-6],
    ],
)
def test_spectrum_exact_eigenvalue_boundaries_fail_closed(eigenvalues):
    modes = SimpleNamespace(eigenvalues=np.asarray(eigenvalues, dtype=float))
    with pytest.raises(diagnostic.ScientificRejection):
        diagnostic._spectrum_payload("boundary", modes)


def test_maximum_overlap_mode_permutation_is_rejected_by_order_flag():
    left = SimpleNamespace(
        eigenvalues=np.array([-0.1, 0.01, 0.02]),
        mass_weighted_eigenvectors=np.eye(3).reshape(3, 1, 3),
    )
    right = SimpleNamespace(
        eigenvalues=np.array([-0.1, 0.01, 0.02]),
        mass_weighted_eigenvectors=np.eye(3)[[1, 0, 2]].reshape(3, 1, 3),
    )
    comparison = diagnostic._mode_comparison(left, right)
    assert comparison["maximum_overlap_assignment"] == [1, 0, 2]
    assert comparison["assignment_equals_eigenvalue_order"] is False


def test_overlap_gate_rejects_value_exactly_at_exclusive_threshold(monkeypatch):
    eigenvalues = np.array([-0.02, 1.0e-5])
    fd_vectors = np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]])
    reference_vectors = np.array(
        [
            [0.995, 0.0, np.sqrt(1.0 - 0.995**2), 0.0],
            [0.0, 1.0, 0.0, 0.0],
        ]
    )
    mode_sets = iter(
        [
            SimpleNamespace(
                eigenvalues=eigenvalues,
                mass_weighted_eigenvectors=fd_vectors.reshape(2, 1, 4),
            ),
            SimpleNamespace(
                eigenvalues=eigenvalues,
                mass_weighted_eigenvectors=fd_vectors.reshape(2, 1, 4),
            ),
            SimpleNamespace(
                eigenvalues=eigenvalues,
                mass_weighted_eigenvectors=fd_vectors.reshape(2, 1, 4),
            ),
            SimpleNamespace(
                eigenvalues=eigenvalues,
                mass_weighted_eigenvectors=reference_vectors.reshape(2, 1, 4),
            ),
        ]
    )
    monkeypatch.setattr(
        diagnostic, "project_vibrational_hessian", lambda *_: next(mode_sets)
    )
    matrix = np.eye(4)
    with pytest.raises(diagnostic.ScientificRejection, match="unstable overlap"):
        diagnostic._confirmation_analysis(
            SimpleNamespace(coords=np.zeros((1, 3))),
            np.ones(1),
            np.array([1.0, 0.0, 0.0, 0.0]),
            matrix,
            matrix,
            matrix,
            matrix,
        )


def test_matrix_gate_boundary_raises_scientific_rejection(monkeypatch):
    thresholds = diagnostic.FINITE_DIFFERENCE_CONTRACT["thresholds"]
    monkeypatch.setitem(
        thresholds, "raw_fd_maximum_absolute_asymmetry_exclusive_maximum", 0.125
    )
    monkeypatch.setitem(
        thresholds, "raw_fd_spectral_relative_asymmetry_exclusive_maximum", 10.0
    )
    raw_fd = np.array([[1.0, 0.0625], [-0.0625, 1.0]])

    with pytest.raises(diagnostic.ScientificRejection, match="maximum asymmetry"):
        diagnostic._confirmation_analysis(
            SimpleNamespace(coords=np.zeros((1, 3))),
            np.ones(1),
            np.ones(2),
            np.eye(2),
            np.eye(2),
            raw_fd,
            np.eye(2),
        )


def test_mode_index_gate_raises_scientific_rejection(monkeypatch):
    eigenvalues = np.array([-0.02, 1.0e-5])
    ordered = SimpleNamespace(
        eigenvalues=eigenvalues,
        mass_weighted_eigenvectors=np.eye(2).reshape(2, 1, 2),
    )
    permuted = SimpleNamespace(
        eigenvalues=eigenvalues,
        mass_weighted_eigenvectors=np.eye(2)[[1, 0]].reshape(2, 1, 2),
    )
    mode_sets = iter([ordered, ordered, ordered, permuted])
    monkeypatch.setattr(
        diagnostic, "project_vibrational_hessian", lambda *_args: next(mode_sets)
    )

    with pytest.raises(diagnostic.ScientificRejection, match="assignment differs"):
        diagnostic._confirmation_analysis(
            SimpleNamespace(coords=np.zeros((1, 3))),
            np.ones(1),
            np.ones(2),
            np.eye(2),
            np.eye(2),
            np.eye(2),
            np.eye(2),
        )


def test_matrix_relative_delta_uses_larger_observed_or_reference_scale():
    observed = np.diag([4.0, 1.0])
    reference = np.diag([2.0, 1.0])

    metrics = diagnostic._matrix_delta_metrics(observed, reference)

    assert metrics["spectral_norm_delta"] == pytest.approx(2.0)
    assert metrics["spectral_relative_delta"] == pytest.approx(0.5)


def _real_reference_inputs():
    preflight = json.loads(
        Path(
            "/mnt/data/vsletten/dissertation-data/"
            "task300-d2c-hessian-5f74662-preflight/preflight.json"
        ).read_text()
    )
    template = diagnostic.reactions(gpu=True, basis="def2-svp")[
        diagnostic.ROUTE
    ].cluster
    fingerprints = diagnostic.campaign._trusted_input_fingerprints_from_route_record(
        preflight["routes"][diagnostic.ROUTE]
    )
    inputs = diagnostic.campaign._trusted_route_input_snapshots(
        diagnostic.campaign.DEFAULT_BUNDLE_ROOT,
        diagnostic.ROUTE,
        template,
        fingerprints,
    )
    masses = np.asarray(preflight["routes"][diagnostic.ROUTE]["masses_amu"])
    reference = Path(
        "/mnt/data/vsletten/dissertation-data/task300-d2c-hessian-5f74662-diagnostic"
    )
    return reference, preflight, inputs, masses


def test_load_reference_accepts_real_pinned_artifact_with_current_inputs(tmp_path):
    reference, preflight, inputs, masses = _real_reference_inputs()

    receipt, matrix, case, binding = diagnostic._load_reference(
        reference,
        6,
        preflight=preflight,
        transition_state=inputs["transition_state"],
        reactant=inputs["reactant"],
        product=inputs["product"],
        masses=masses,
    )

    assert receipt["git_sha"] == "5f74662a079129b54cbbbe7a37aaba47c6cf6a99"
    assert matrix.shape == (18, 18)
    assert case["name"] == "D-dense-strict-reference"
    assert (
        case["post_symmetrization_gate"]["evidence"]["geometry_fingerprint"]
        == binding["transition_state_geometry_fingerprint"]
    )
    assert binding["analytic_reference_identity"]["receipt_sha256"] == (
        diagnostic.REFERENCE_RECEIPT_SHA256
    )

    copied = tmp_path / "reference"
    shutil.copytree(reference, copied)
    matrix_path = copied / "D-dense-strict-reference/total-symmetric.f64"
    raw = matrix_path.read_bytes()
    matrix_path.write_bytes(raw[:-1] + bytes([raw[-1] ^ 1]))
    with pytest.raises(ValueError, match="content SHA-256 mismatch"):
        diagnostic._load_reference(
            copied,
            6,
            preflight=preflight,
            transition_state=inputs["transition_state"],
            reactant=inputs["reactant"],
            product=inputs["product"],
            masses=masses,
        )


def test_reference_accepts_old_source_identity_but_rejects_scientific_drift():
    reference, preflight, inputs, masses = _real_reference_inputs()
    current_preflight = json.loads(json.dumps(preflight))
    current_preflight["identity"] = "0" * 64

    receipt, _matrix, _case, binding = diagnostic._load_reference(
        reference,
        6,
        preflight=current_preflight,
        transition_state=inputs["transition_state"],
        reactant=inputs["reactant"],
        product=inputs["product"],
        masses=masses,
    )

    assert receipt["campaign_identity"] != binding["campaign_identity"]
    assert binding["analytic_reference_identity"]["git_sha"] == receipt["git_sha"]
    drifted = diagnostic.replace(
        inputs["transition_state"], coords=inputs["transition_state"].coords + 1.0e-6
    )
    with pytest.raises(ValueError, match="TS geometry fingerprint drifted"):
        diagnostic._load_reference(
            reference,
            6,
            preflight=current_preflight,
            transition_state=drifted,
            reactant=inputs["reactant"],
            product=inputs["product"],
            masses=masses,
        )


def test_artifact_bundle_hashes_every_array_and_rejects_corruption(tmp_path):
    payload, receipt_sha = diagnostic._publish_record_bundle(
        tmp_path,
        diagnostic.Path("points/test"),
        {"schema": diagnostic.FINITE_DIFFERENCE_SCHEMA, "point_key": "test"},
        {"gradient.f64": (np.arange(6).reshape(2, 3), "hartree / bohr")},
    )
    receipt, observed_sha, arrays = diagnostic._read_record_bundle(
        tmp_path, diagnostic.Path("points/test")
    )
    assert receipt == payload
    assert observed_sha == receipt_sha
    raw = (tmp_path / "points/test/gradient.f64").read_bytes()
    assert (
        payload["artifacts"]["gradient.f64"]["sha256"]
        == hashlib.sha256(raw).hexdigest()
    )
    assert arrays["gradient.f64"] == pytest.approx(np.arange(6).reshape(2, 3))

    (tmp_path / "points/test/gradient.f64").write_bytes(raw[:-1] + b"x")
    with pytest.raises(ValueError, match="hash mismatch"):
        diagnostic._read_record_bundle(tmp_path, diagnostic.Path("points/test"))


def _fake_fd_environment(tmp_path, monkeypatch):
    transition_state = Cluster(
        "frozen-ts",
        ["C", "O", "H", "O", "H", "H"],
        np.array(
            [
                [0.0, 0.0, 0.0],
                [1.2, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.2],
                [1.0, 1.0, 0.2],
                [0.2, 1.0, 1.0],
            ]
        ),
        spin=1,
    )
    reactant = diagnostic.replace(
        transition_state, name="mapped-reactant", coords=transition_state.coords - 0.01
    )
    product = diagnostic.replace(
        transition_state, name="mapped-product", coords=transition_state.coords + 0.01
    )
    snapshots = {
        "transition_state": transition_state,
        "reactant": reactant,
        "product": product,
    }
    masses = [12.0, 16.0, 1.0, 16.0, 1.0, 1.0]
    fingerprints = {
        "transition_state": {
            "file_sha256": "1" * 64,
            "geometry_sha256": "2" * 64,
            "geometry_fingerprint": diagnostic.frequency_geometry_fingerprint(
                transition_state
            ),
        },
        "reactant": {
            "file_sha256": "3" * 64,
            "geometry_sha256": "4" * 64,
            "geometry_fingerprint": diagnostic.frequency_geometry_fingerprint(reactant),
        },
        "product": {
            "file_sha256": "5" * 64,
            "geometry_sha256": "6" * 64,
            "geometry_fingerprint": diagnostic.frequency_geometry_fingerprint(product),
        },
    }
    preflight = {
        "identity": "c" * 64,
        "routes": {
            diagnostic.ROUTE: {
                "masses_amu": masses,
                "symbols": list(transition_state.symbols),
                "atom_identity_labels": [f"atom-{index}" for index in range(6)],
                "atom_mapping_sha256": "7" * 64,
            }
        },
    }
    base = DftSettings(
        xc="pwb6k",
        basis="def2-svp",
        dispersion="d3bj",
        density_fit=True,
        use_gpu=True,
    )
    identity = _identity()
    monkeypatch.setattr(diagnostic, "_fd_execution_identity", lambda: identity)
    monkeypatch.setattr(
        diagnostic.campaign,
        "_validate_production_boundary",
        lambda *_: (preflight, "d" * 64),
    )
    monkeypatch.setattr(
        diagnostic.campaign, "_canonical_dft_settings", lambda _: (base, "settings")
    )
    monkeypatch.setattr(
        diagnostic.campaign,
        "_trusted_input_fingerprints_from_route_record",
        lambda _: fingerprints,
    )
    monkeypatch.setattr(
        diagnostic.campaign,
        "_trusted_route_input_snapshots",
        lambda *_: snapshots,
    )
    monkeypatch.setattr(
        diagnostic,
        "reactions",
        lambda **_: {diagnostic.ROUTE: SimpleNamespace(cluster=transition_state)},
    )
    monkeypatch.setattr(
        diagnostic,
        "_load_reference",
        lambda *_, **__: (
            {"git_sha": "e" * 40, "script_sha256": "f" * 64},
            np.eye(18),
            {
                "electronic_hartree": -10.0,
                "spin_square": [0.76, 2.0],
            },
            binding,
        ),
    )
    monkeypatch.setattr(
        diagnostic.campaign, "_mapped_route_vector", lambda *_: np.ones(18)
    )

    def confirmation(*args):
        h_h, h_2h, raw_fd = args[3:6]
        return {
            "confirmation_passed": True,
            "symmetric_matrices": {
                "H_h": 0.5 * (h_h + h_h.T),
                "H_2h": 0.5 * (h_2h + h_2h.T),
                "FD": 0.5 * (raw_fd + raw_fd.T),
            },
        }

    monkeypatch.setattr(diagnostic, "_confirmation_analysis", confirmation)
    binding = diagnostic._current_reference_binding(
        preflight,
        transition_state,
        reactant,
        product,
        np.asarray(masses),
    )
    return transition_state, identity, binding


def test_resume_reconciles_durable_point_and_never_recomputes_it(tmp_path, monkeypatch):
    transition_state, identity, binding = _fake_fd_environment(tmp_path, monkeypatch)
    preflight_root = tmp_path / "preflight"
    reference_root = tmp_path / "reference"
    output_root = tmp_path / "output"
    output_root.mkdir()
    status = diagnostic._initial_fd_status(identity, preflight_root, reference_root)
    status.update(
        {
            "preflight_receipt_sha256": "d" * 64,
            "reference_receipt_sha256": diagnostic.REFERENCE_RECEIPT_SHA256,
            "reference_matrix_sha256": diagnostic.REFERENCE_MATRIX_SHA256,
            "reference_binding": binding,
            "reference_binding_sha256": diagnostic.campaign._canonical_hash(binding),
        }
    )
    density = np.arange(8, dtype=float).reshape(2, 2, 2)
    density_sha = hashlib.sha256(density.astype("<f8").tobytes()).hexdigest()
    center = _point_record(density_sha=density_sha, fmax=0.0)
    center.update(
        {
            "schema": diagnostic.FINITE_DIFFERENCE_SCHEMA,
            "point_key": "center",
            "kind": "center",
            **diagnostic._point_bindings(status),
            "geometry_fingerprint": diagnostic.frequency_geometry_fingerprint(
                transition_state
            ),
            **diagnostic._center_reference_metrics(
                center,
                {"electronic_hartree": -10.0, "spin_square": [0.76, 2.0]},
            ),
        }
    )
    _, center_sha = diagnostic._publish_record_bundle(
        output_root,
        diagnostic.Path("points/center"),
        center,
        {
            "gradient.f64": (np.zeros((6, 3)), "hartree / bohr"),
            "density.f64": (density, "electrons"),
        },
    )
    first_definition = diagnostic._displacement_plan(6)[0]
    first = _point_record(density_sha=density_sha)
    first.update(
        {
            "schema": diagnostic.FINITE_DIFFERENCE_SCHEMA,
            "point_key": first_definition["key"],
            "kind": "displacement",
            **diagnostic._point_bindings(status),
            **first_definition,
            **diagnostic._density_continuity_metrics(density, density),
        }
    )
    first_coords = transition_state.coords.copy()
    first_coords.reshape(-1)[first_definition["coordinate_index"]] += (
        first_definition["displacement_bohr"] * diagnostic.BOHR_TO_ANGSTROM
    )
    first["geometry_fingerprint"] = diagnostic.frequency_geometry_fingerprint(
        diagnostic.replace(transition_state, coords=first_coords)
    )
    diagnostic._publish_record_bundle(
        output_root,
        diagnostic.Path("points") / first_definition["key"],
        first,
        {
            "gradient.f64": (np.zeros((6, 3)), "hartree / bohr"),
            "density.f64": (density, "electrons"),
        },
    )
    status["completed_points"] = [{"point_key": "center", "receipt_sha256": center_sha}]
    diagnostic._atomic_write(
        output_root / "status.json", diagnostic._json_bytes(status)
    )

    calls = []

    def evaluate(cluster, _settings, *, central_density):
        calls.append((cluster.coords.copy(), central_density))
        gradient = np.zeros((6, 3))
        record = _point_record(density_sha=density_sha)
        record["geometry_fingerprint"] = diagnostic.frequency_geometry_fingerprint(
            cluster
        )
        return record, gradient, density

    monkeypatch.setattr(diagnostic, "_evaluate_fd_gradient", evaluate)
    receipt = diagnostic.run_finite_difference(
        preflight_root, reference_root, output_root
    )

    assert len(calls) == 71
    assert all(call_density is calls[0][1] for _, call_density in calls)
    assert len(receipt["point_receipts"]) == 73
    assert receipt["accepted_campaign_result"] is False
    assert receipt["confirmation_passed"] is True
    receipt_raw = (output_root / "receipt.json").read_bytes()
    terminal = json.loads((output_root / diagnostic.TERMINAL).read_text())
    assert terminal["receipt_sha256"] == hashlib.sha256(receipt_raw).hexdigest()
    assert (output_root / "points" / first_definition["key"]).is_dir()
    matrix_receipt = json.loads((output_root / "matrices/receipt.json").read_text())
    assert set(matrix_receipt["artifacts"]) == {
        "H_h-raw.f64",
        "H_2h-raw.f64",
        "richardson-raw.f64",
        "H_h-symmetric.f64",
        "H_2h-symmetric.f64",
        "richardson-symmetric.f64",
    }


def test_completed_scientific_rejection_preserves_matrix_receipt(tmp_path, monkeypatch):
    _fake_fd_environment(tmp_path, monkeypatch)
    density = np.arange(8, dtype=float).reshape(2, 2, 2)
    density_sha = hashlib.sha256(density.astype("<f8").tobytes()).hexdigest()

    def evaluate(cluster, _settings, *, central_density):
        record = _point_record(density_sha=density_sha)
        if central_density is None:
            record["physical_fmax_ev_per_angstrom"] = 0.0
        record["geometry_fingerprint"] = diagnostic.frequency_geometry_fingerprint(
            cluster
        )
        return (
            record,
            np.zeros((6, 3)),
            density,
        )

    monkeypatch.setattr(diagnostic, "_evaluate_fd_gradient", evaluate)
    monkeypatch.setattr(
        diagnostic,
        "_confirmation_analysis",
        lambda *_: (_ for _ in ()).throw(
            diagnostic.ScientificRejection("synthetic scientific rejection")
        ),
    )
    output_root = tmp_path / "rejected"

    receipt = diagnostic.run_finite_difference(
        tmp_path / "preflight", tmp_path / "reference", output_root
    )

    assert receipt["state"] == "completed"
    assert receipt["accepted_campaign_result"] is False
    assert receipt["confirmation_passed"] is False
    assert receipt["analysis"] == {
        "confirmation_passed": False,
        "rejection": "ScientificRejection: synthetic scientific rejection",
    }
    assert (output_root / "matrices/richardson-raw.f64").is_file()
    terminal = json.loads((output_root / diagnostic.TERMINAL).read_text())
    assert terminal["state"] == "completed"
    assert terminal["confirmation_passed"] is False

    (output_root / diagnostic.TERMINAL).unlink()
    reconciled = diagnostic.finalize_if_running(output_root, finite_difference=True)
    assert reconciled["state"] == "completed"
    assert reconciled["confirmation_passed"] is False


def test_resume_refuses_uncertain_unreceipted_point_to_prevent_duplicate(tmp_path):
    output_root = tmp_path / "output"
    status = diagnostic._initial_fd_status(
        _identity(),
        tmp_path / "preflight",
        tmp_path / "reference",
    )
    status["current_point"] = "center"
    with (
        diagnostic._exclusive_output_claim(output_root) as claim,
        pytest.raises(RuntimeError, match="refusing a duplicate evaluation"),
    ):
        diagnostic._resume_fd_points(
            output_root,
            status,
            diagnostic._displacement_plan(6),
            claim=claim,
        )


def test_fd_failure_and_dead_man_receipts_never_accept_campaign(tmp_path):
    output_root = tmp_path / "fd"
    output_root.mkdir()
    status = diagnostic._initial_fd_status(
        _identity(),
        tmp_path / "preflight",
        tmp_path / "reference",
    )
    status["current_point"] = "coordinate-004-p1"
    diagnostic._atomic_write(
        output_root / "status.json", diagnostic._json_bytes(status)
    )

    terminal = diagnostic.finalize_if_running(output_root)

    assert terminal["schema"] == diagnostic.FINITE_DIFFERENCE_SCHEMA
    assert terminal["state"] == "failed"
    assert terminal["accepted_campaign_result"] is False
    assert terminal["confirmation_passed"] is False
    assert terminal["current_point"] == "coordinate-004-p1"


def test_fd_run_failure_always_publishes_terminal_receipt(tmp_path, monkeypatch):
    monkeypatch.setattr(
        diagnostic,
        "_run_finite_difference_locked",
        lambda *_: (_ for _ in ()).throw(RuntimeError("synthetic FD failure")),
    )
    output_root = tmp_path / "fd-run"

    with pytest.raises(RuntimeError, match="synthetic FD failure"):
        diagnostic.run_finite_difference(
            tmp_path / "preflight", tmp_path / "reference", output_root
        )

    terminal = json.loads((output_root / diagnostic.TERMINAL).read_text())
    status = json.loads((output_root / "status.json").read_text())
    assert terminal["schema"] == diagnostic.FINITE_DIFFERENCE_SCHEMA
    assert terminal["state"] == "failed"
    assert terminal["accepted_campaign_result"] is False
    assert terminal["confirmation_passed"] is False
    assert terminal["detail"] == "RuntimeError: synthetic FD failure"
    assert status["state"] == "failed"
    assert status["contract"] == diagnostic.FINITE_DIFFERENCE_CONTRACT


def _reference_binding_inputs(monkeypatch):
    transition_state = Cluster(
        "binding-ts",
        ["H", "O"],
        np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
        spin=1,
    )
    reactant = diagnostic.replace(
        transition_state, name="binding-reactant", coords=transition_state.coords - 0.02
    )
    product = diagnostic.replace(
        transition_state, name="binding-product", coords=transition_state.coords + 0.03
    )
    masses = np.array([1.00782503223, 15.99491461957])
    fingerprints = {
        "transition_state": {
            "file_sha256": "1" * 64,
            "geometry_sha256": "2" * 64,
            "geometry_fingerprint": diagnostic.frequency_geometry_fingerprint(
                transition_state
            ),
        },
        "reactant": {
            "file_sha256": "3" * 64,
            "geometry_sha256": "4" * 64,
            "geometry_fingerprint": diagnostic.frequency_geometry_fingerprint(reactant),
        },
        "product": {
            "file_sha256": "5" * 64,
            "geometry_sha256": "6" * 64,
            "geometry_fingerprint": diagnostic.frequency_geometry_fingerprint(product),
        },
    }
    preflight = {
        "identity": "c" * 64,
        "routes": {
            diagnostic.ROUTE: {
                "symbols": ["H", "O"],
                "atom_identity_labels": ["transferring-H", "acceptor-O"],
                "atom_mapping_sha256": "7" * 64,
                "masses_amu": masses.tolist(),
            }
        },
    }
    monkeypatch.setattr(
        diagnostic.campaign,
        "_trusted_input_fingerprints_from_route_record",
        lambda _record: fingerprints,
    )
    return preflight, transition_state, reactant, product, masses, fingerprints


def test_reference_binding_records_exact_current_route_inputs(monkeypatch):
    preflight, transition_state, reactant, product, masses, fingerprints = (
        _reference_binding_inputs(monkeypatch)
    )

    binding = diagnostic._current_reference_binding(
        preflight, transition_state, reactant, product, masses
    )

    assert binding == {
        "campaign_identity": preflight["identity"],
        "route": diagnostic.ROUTE,
        "symbols": ["H", "O"],
        "atom_identity_labels": ["transferring-H", "acceptor-O"],
        "atom_mapping_sha256": "7" * 64,
        "transition_state_file_sha256": "1" * 64,
        "transition_state_geometry_sha256": "2" * 64,
        "transition_state_geometry_fingerprint": fingerprints["transition_state"][
            "geometry_fingerprint"
        ],
        "mapped_reactant_file_sha256": "3" * 64,
        "mapped_reactant_geometry_sha256": "4" * 64,
        "mapped_reactant_geometry_fingerprint": fingerprints["reactant"][
            "geometry_fingerprint"
        ],
        "mapped_reactant_coordinate_sha256": hashlib.sha256(
            np.ascontiguousarray(reactant.coords, dtype="<f8").tobytes()
        ).hexdigest(),
        "mapped_product_file_sha256": "5" * 64,
        "mapped_product_geometry_sha256": "6" * 64,
        "mapped_product_geometry_fingerprint": fingerprints["product"][
            "geometry_fingerprint"
        ],
        "mapped_product_coordinate_sha256": hashlib.sha256(
            np.ascontiguousarray(product.coords, dtype="<f8").tobytes()
        ).hexdigest(),
        "isotopic_masses_amu": masses.tolist(),
        "fd_method": diagnostic.FINITE_DIFFERENCE_CONTRACT["method"],
        "analytic_reference_settings_fingerprint": (
            diagnostic._reference_settings_fingerprint()
        ),
    }


@pytest.mark.parametrize("drift", ["ts", "reactant", "product", "masses", "route"])
def test_reference_binding_rejects_each_current_identity_drift(monkeypatch, drift):
    preflight, transition_state, reactant, product, masses, _fingerprints = (
        _reference_binding_inputs(monkeypatch)
    )
    if drift == "ts":
        transition_state = diagnostic.replace(
            transition_state, coords=transition_state.coords + 0.001
        )
    elif drift == "reactant":
        reactant = diagnostic.replace(reactant, coords=reactant.coords + 0.001)
    elif drift == "product":
        product = diagnostic.replace(product, coords=product.coords + 0.001)
    elif drift == "masses":
        masses = masses.copy()
        masses[0] += 0.001
    else:
        preflight["routes"][diagnostic.ROUTE]["atom_mapping_sha256"] = "invalid"

    with pytest.raises(ValueError):
        diagnostic._current_reference_binding(
            preflight, transition_state, reactant, product, masses
        )


def test_reference_binding_rejects_method_drift(monkeypatch):
    preflight, transition_state, reactant, product, masses, _fingerprints = (
        _reference_binding_inputs(monkeypatch)
    )
    monkeypatch.setitem(
        diagnostic.FINITE_DIFFERENCE_CONTRACT["method"], "basis", "wrong-basis"
    )
    with pytest.raises(ValueError, match="method"):
        diagnostic._current_reference_binding(
            preflight, transition_state, reactant, product, masses
        )


def test_center_reference_exclusive_gates_reject_exact_binary_boundaries(monkeypatch):
    record = _point_record(density_sha="a" * 64, s2=0.78125)
    record["electronic_hartree"] = -10.0
    reference = {"electronic_hartree": -9.875, "spin_square": [0.75, 2.0]}
    thresholds = diagnostic.FINITE_DIFFERENCE_CONTRACT["thresholds"]
    monkeypatch.setitem(
        thresholds,
        "center_energy_vs_reference_exclusive_maximum_delta_hartree",
        0.125,
    )
    with pytest.raises(diagnostic.ScientificRejection, match="energy delta"):
        diagnostic._enforce_center_reference_metrics(
            diagnostic._center_reference_metrics(record, reference)
        )

    reference["electronic_hartree"] = -10.0
    monkeypatch.setitem(
        thresholds, "center_s2_vs_reference_exclusive_maximum_delta", 0.03125
    )
    with pytest.raises(diagnostic.ScientificRejection, match="S2 delta"):
        diagnostic._enforce_center_reference_metrics(
            diagnostic._center_reference_metrics(record, reference)
        )


def test_density_continuity_is_exclusive_and_rejects_malformed_evidence(monkeypatch):
    monkeypatch.setitem(
        diagnostic.FINITE_DIFFERENCE_CONTRACT["thresholds"],
        (
            "displaced_final_density_vs_center_exclusive_maximum_"
            "normalized_frobenius_delta"
        ),
        0.5,
    )
    boundary = diagnostic._density_continuity_metrics(np.array([2.0]), np.array([1.0]))
    with pytest.raises(diagnostic.ScientificRejection, match="Frobenius"):
        diagnostic._enforce_density_continuity(boundary)
    assert diagnostic._density_continuity_metrics(np.array([1.5]), np.array([1.0]))[
        "final_vs_center_normalized_frobenius_delta"
    ] == pytest.approx(1.0 / 3.0)
    with pytest.raises(ValueError, match="one finite shape"):
        diagnostic._density_continuity_metrics(
            np.array([1.0, np.nan]), np.array([1.0, 0.0])
        )
    with pytest.raises(ValueError, match="one finite shape"):
        diagnostic._density_continuity_metrics(np.ones(2), np.ones(3))


def test_loaded_diagnostic_code_sabotage_is_detected_after_source_is_restored(
    monkeypatch,
):
    def sabotaged(value, limit, *, label):
        return 0.0

    monkeypatch.setattr(diagnostic._strict_less, "__code__", sabotaged.__code__)
    with pytest.raises(RuntimeError, match="loaded Python code disagrees with source"):
        diagnostic._diagnostic_execution_identity()


def test_loaded_diagnostic_code_rejects_foreign_callable_replacement(monkeypatch):
    monkeypatch.setattr(diagnostic, "_strict_less", lambda *_args, **_kwargs: 0.0)

    with pytest.raises(RuntimeError, match="foreign owner"):
        diagnostic._diagnostic_execution_identity()


def test_loaded_diagnostic_state_rejects_threshold_contract_mutation(monkeypatch):
    original_identity = diagnostic._diagnostic_execution_identity()
    monkeypatch.setitem(
        diagnostic.FINITE_DIFFERENCE_CONTRACT["thresholds"],
        "center_fmax_ev_per_angstrom_exclusive_maximum",
        7.0,
    )

    with pytest.raises(RuntimeError, match="scientific state disagrees with source"):
        diagnostic._diagnostic_execution_identity()
    assert len(original_identity["loaded_state_sha256"]) == 64


def test_external_claim_survives_output_directory_replacement_and_blocks_second_claim(
    tmp_path,
):
    output = tmp_path / "output"
    displaced = tmp_path / "displaced-original"

    with (
        pytest.raises(RuntimeError, match="output root was replaced"),
        diagnostic._exclusive_output_claim(output) as claim,
    ):
        os.rename(output, displaced)
        output.mkdir()
        with (
            pytest.raises(RuntimeError, match="already claimed"),
            diagnostic._exclusive_output_claim(output),
        ):
            pass
        claim.verify()


def test_external_claim_prevents_duplicate_expensive_call(tmp_path, monkeypatch):
    output = tmp_path / "output"
    calls = []
    monkeypatch.setattr(
        diagnostic,
        "_run_finite_difference_locked",
        lambda *_args: calls.append("expensive"),
    )

    with (
        diagnostic._exclusive_output_claim(output),
        pytest.raises(RuntimeError, match="already claimed"),
    ):
        diagnostic.run_finite_difference(
            tmp_path / "preflight", tmp_path / "reference", output
        )
    assert calls == []


@pytest.mark.parametrize("which", ["preflight", "reference", "output"])
def test_fd_cli_paths_reject_symlink_components_before_execution(
    tmp_path, monkeypatch, which
):
    real = tmp_path / "real"
    real.mkdir()
    linked = tmp_path / "linked"
    linked.symlink_to(real, target_is_directory=True)
    paths = {
        "preflight": tmp_path / "preflight",
        "reference": tmp_path / "reference",
        "output": tmp_path / "output",
    }
    paths[which] = linked / which
    calls = []
    monkeypatch.setattr(
        diagnostic,
        "_run_finite_difference_locked",
        lambda *_args: calls.append("ran"),
    )

    with pytest.raises(ValueError, match="contains a symlink"):
        diagnostic.run_finite_difference(
            paths["preflight"], paths["reference"], paths["output"]
        )
    assert calls == []


def test_disk_ahead_bundle_schema_rejects_receipt_and_artifact_extras(tmp_path):
    relative = diagnostic.Path("points/center")
    diagnostic._publish_record_bundle(
        tmp_path,
        relative,
        {"schema": diagnostic.FINITE_DIFFERENCE_SCHEMA, "point_key": "center"},
        {"gradient.f64": (np.zeros((2, 3)), "hartree / bohr")},
    )
    receipt_path = tmp_path / relative / "receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["unexpected"] = True
    receipt_path.write_bytes(diagnostic._json_bytes(receipt))
    with pytest.raises(ValueError, match="receipt schema is not exact"):
        diagnostic._read_record_bundle(
            tmp_path,
            relative,
            expected_receipt_keys={"schema", "point_key", "artifacts"},
            expected_artifacts={"gradient.f64": ((2, 3), "hartree / bohr")},
        )


def test_matrix_bundle_exact_shape_units_and_inventory_are_enforced(tmp_path):
    relative = diagnostic.Path("matrices")
    record = {"schema": diagnostic.FINITE_DIFFERENCE_SCHEMA, "kind": "matrix"}
    diagnostic._publish_record_bundle(
        tmp_path,
        relative,
        record,
        {"H_h-raw.f64": (np.eye(2), "hartree / bohr^2")},
    )
    with pytest.raises(ValueError, match="shape or units mismatch"):
        diagnostic._read_record_bundle(
            tmp_path,
            relative,
            expected_receipt_keys={"schema", "kind", "artifacts"},
            expected_artifacts={"H_h-raw.f64": ((3, 3), "hartree / bohr^2")},
        )


def test_programming_error_in_confirmation_analysis_produces_failed_terminal(
    tmp_path, monkeypatch
):
    _fake_fd_environment(tmp_path, monkeypatch)
    density = np.arange(8, dtype=float).reshape(2, 2, 2)
    density_sha = hashlib.sha256(density.astype("<f8").tobytes()).hexdigest()

    def evaluate(cluster, _settings, *, central_density):
        record = _point_record(density_sha=density_sha)
        if central_density is None:
            record["physical_fmax_ev_per_angstrom"] = 0.0
        record["geometry_fingerprint"] = diagnostic.frequency_geometry_fingerprint(
            cluster
        )
        return record, np.zeros((6, 3)), density

    monkeypatch.setattr(diagnostic, "_evaluate_fd_gradient", evaluate)
    monkeypatch.setattr(
        diagnostic,
        "_confirmation_analysis",
        lambda *_: (_ for _ in ()).throw(ValueError("malformed matrix shape")),
    )
    output = tmp_path / "programming-error"

    with pytest.raises(ValueError, match="malformed matrix shape"):
        diagnostic.run_finite_difference(
            tmp_path / "preflight", tmp_path / "reference", output
        )
    terminal = json.loads((output / diagnostic.TERMINAL).read_text())
    assert terminal["state"] == "failed"
    assert terminal["confirmation_passed"] is False
    assert not (output / "receipt.json").exists()


def test_fd_finalizer_malformed_status_uses_explicit_fail_closed_schema(tmp_path):
    output = tmp_path / "malformed-status"
    output.mkdir()
    (output / "status.json").write_bytes(b"{")

    terminal = diagnostic.finalize_if_running(output, finite_difference=True)

    assert terminal["schema"] == diagnostic.FINITE_DIFFERENCE_SCHEMA
    assert terminal["state"] == "failed"
    assert terminal["accepted_campaign_result"] is False
    assert terminal["confirmation_passed"] is False


def test_finalizer_never_clobbers_malformed_terminal(tmp_path):
    output = tmp_path / "malformed-terminal"
    output.mkdir()
    terminal_path = output / diagnostic.TERMINAL
    original = b'{"state":"completed"}'
    terminal_path.write_bytes(original)

    with pytest.raises(ValueError):
        diagnostic.finalize_if_running(output, finite_difference=True)

    assert terminal_path.read_bytes() == original


def test_fd_cli_stdout_includes_both_scientific_verdict_fields(
    tmp_path, monkeypatch, capsys
):
    output = tmp_path / "output"
    monkeypatch.setattr(
        diagnostic,
        "run_finite_difference",
        lambda *_args: {
            "schema": diagnostic.FINITE_DIFFERENCE_SCHEMA,
            "state": "completed",
            "confirmation_passed": False,
            "accepted_campaign_result": False,
        },
    )
    monkeypatch.setattr(
        diagnostic.sys,
        "argv",
        [
            "diagnostic",
            "--finite-difference-confirmation",
            "--preflight-root",
            str(tmp_path / "preflight"),
            "--reference-root",
            str(tmp_path / "reference"),
            "--output-root",
            str(output),
        ],
    )

    assert diagnostic.main() == 0
    stdout = json.loads(capsys.readouterr().out)
    assert stdout["confirmation_passed"] is False
    assert stdout["accepted_campaign_result"] is False


def _complete_fake_fd(tmp_path, monkeypatch, name="complete"):
    _fake_fd_environment(tmp_path, monkeypatch)
    density = np.arange(8, dtype=float).reshape(2, 2, 2)
    density_sha = hashlib.sha256(density.astype("<f8").tobytes()).hexdigest()

    def evaluate(cluster, _settings, *, central_density):
        record = _point_record(density_sha=density_sha)
        if central_density is None:
            record["physical_fmax_ev_per_angstrom"] = 0.0
        record["geometry_fingerprint"] = diagnostic.frequency_geometry_fingerprint(
            cluster
        )
        return record, np.zeros((6, 3)), density

    monkeypatch.setattr(diagnostic, "_evaluate_fd_gradient", evaluate)
    output = tmp_path / name
    receipt = diagnostic.run_finite_difference(
        tmp_path / "preflight", tmp_path / "reference", output
    )
    return output, receipt


def test_finalizer_adopts_completed_receipt_and_rebuilds_stale_status(
    tmp_path, monkeypatch
):
    output, receipt = _complete_fake_fd(tmp_path, monkeypatch)
    receipt_sha = hashlib.sha256((output / "receipt.json").read_bytes()).hexdigest()
    (output / diagnostic.TERMINAL).unlink()
    (output / "status.json").write_bytes(
        diagnostic._json_bytes(
            {"schema": diagnostic.FINITE_DIFFERENCE_SCHEMA, "state": "running"}
        )
    )

    terminal = diagnostic.finalize_if_running(output, finite_difference=True)

    assert terminal["state"] == "completed"
    assert terminal["receipt_sha256"] == receipt_sha
    assert json.loads((output / "status.json").read_text()) == (
        diagnostic._completed_fd_status(receipt, receipt_sha)
    )


def test_finalizer_reconciles_terminal_published_before_completed_status(
    tmp_path, monkeypatch
):
    output, receipt = _complete_fake_fd(tmp_path, monkeypatch, "terminal-first")
    receipt_sha = hashlib.sha256((output / "receipt.json").read_bytes()).hexdigest()
    terminal_before = (output / diagnostic.TERMINAL).read_bytes()
    (output / "status.json").write_bytes(
        diagnostic._json_bytes(
            {"schema": diagnostic.FINITE_DIFFERENCE_SCHEMA, "state": "running"}
        )
    )

    terminal = diagnostic.finalize_if_running(output, finite_difference=True)

    assert diagnostic._json_bytes(terminal) == terminal_before
    assert json.loads((output / "status.json").read_text()) == (
        diagnostic._completed_fd_status(receipt, receipt_sha)
    )


def test_finalizer_rejects_inconsistent_completed_terminal_without_clobber(
    tmp_path, monkeypatch
):
    output, _receipt = _complete_fake_fd(tmp_path, monkeypatch, "bad-terminal")
    terminal_path = output / diagnostic.TERMINAL
    terminal = json.loads(terminal_path.read_text())
    terminal["git_sha"] = "0" * 40
    original = diagnostic._json_bytes(terminal)
    terminal_path.write_bytes(original)

    with pytest.raises(ValueError, match="authoritative completed FD terminal"):
        diagnostic.finalize_if_running(output, finite_difference=True)

    assert terminal_path.read_bytes() == original


def test_finalizer_rejects_rehashed_forged_matrix_reconstruction(tmp_path, monkeypatch):
    output, _receipt = _complete_fake_fd(tmp_path, monkeypatch, "forged-matrix")
    (output / diagnostic.TERMINAL).unlink()
    matrix_path = output / "matrices/richardson-raw.f64"
    forged = np.full((18, 18), 7.0, dtype="<f8").tobytes()
    matrix_path.write_bytes(forged)
    aggregate_path = output / "matrices/receipt.json"
    aggregate = json.loads(aggregate_path.read_text())
    aggregate["artifacts"]["richardson-raw.f64"]["sha256"] = hashlib.sha256(
        forged
    ).hexdigest()
    aggregate_raw = diagnostic._json_bytes(aggregate)
    aggregate_path.write_bytes(aggregate_raw)
    receipt_path = output / "receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["aggregate_receipt"]["sha256"] = hashlib.sha256(aggregate_raw).hexdigest()
    receipt_path.write_bytes(diagnostic._json_bytes(receipt))

    terminal = diagnostic.finalize_if_running(output, finite_difference=True)

    assert terminal["state"] == "failed"
    assert "73-gradient reconstruction" in terminal["detail"]
    status = json.loads((output / "status.json").read_text())
    assert status["state"] == "failed"
    assert "receipt" not in status


def test_finalizer_rejects_malformed_reference_git_sha(tmp_path, monkeypatch):
    output, _receipt = _complete_fake_fd(tmp_path, monkeypatch, "bad-reference-git")
    (output / diagnostic.TERMINAL).unlink()
    receipt_path = output / "receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["reference_git_sha"] = "not-a-git-sha"
    receipt_path.write_bytes(diagnostic._json_bytes(receipt))

    terminal = diagnostic.finalize_if_running(output, finite_difference=True)

    assert terminal["state"] == "failed"
    assert "reference_git_sha" in terminal["detail"]
    assert "receipt" not in terminal


def test_analytic_run_holds_external_claim_and_finalizer_cannot_clobber(
    tmp_path, monkeypatch
):
    output = tmp_path / "analytic-claim"
    monkeypatch.setattr(
        diagnostic,
        "_git_identity",
        lambda: {"git_sha": "a" * 40, "script_sha256": "b" * 64},
    )

    def interrupt(_preflight, _route):
        with pytest.raises(diagnostic.ConcurrentRunError, match="already claimed"):
            diagnostic.finalize_if_running(output)
        assert not (output / diagnostic.TERMINAL).exists()
        raise RuntimeError("deterministic analytic interruption")

    monkeypatch.setattr(diagnostic.campaign, "_validate_production_boundary", interrupt)
    publication_calls = []
    original_publish = diagnostic._terminal_write_noreplace

    def publish(claim, path, payload):
        publication_calls.append(payload["state"])
        return original_publish(claim, path, payload)

    monkeypatch.setattr(diagnostic, "_terminal_write_noreplace", publish)

    with pytest.raises(RuntimeError, match="deterministic analytic interruption"):
        diagnostic.run(tmp_path / "preflight", output)

    terminal_path = output / diagnostic.TERMINAL
    original = terminal_path.read_bytes()
    assert publication_calls == ["failed"]
    assert diagnostic.finalize_if_running(output)["state"] == "failed"
    assert terminal_path.read_bytes() == original


def test_resume_rejects_corrupted_durable_final_density(tmp_path, monkeypatch):
    output, _receipt = _complete_fake_fd(tmp_path, monkeypatch, "density-corruption")
    status = json.loads((output / "status.json").read_text())
    first_key = diagnostic._displacement_plan(6)[0]["key"]
    density_path = output / "points" / first_key / "density.f64"
    raw = density_path.read_bytes()
    density_path.write_bytes(raw[:-1] + bytes([raw[-1] ^ 1]))

    with (
        diagnostic._exclusive_output_claim(output) as claim,
        pytest.raises(ValueError, match="hash mismatch"),
    ):
        diagnostic._resume_fd_points(
            output,
            status,
            diagnostic._displacement_plan(6),
            claim=claim,
        )
    with pytest.raises(ValueError, match="hash mismatch"):
        diagnostic.finalize_if_running(output, finite_difference=True)


def test_resume_rejects_exact_point_ancestry_mismatch(tmp_path, monkeypatch):
    output, _receipt = _complete_fake_fd(tmp_path, monkeypatch, "ancestry-mismatch")
    status = json.loads((output / "status.json").read_text())
    first_key = diagnostic._displacement_plan(6)[0]["key"]
    point_path = output / "points" / first_key / "receipt.json"
    point = json.loads(point_path.read_text())
    point["reference_binding_sha256"] = "0" * 64
    point_raw = diagnostic._json_bytes(point)
    point_path.write_bytes(point_raw)
    status["completed_points"][1]["receipt_sha256"] = hashlib.sha256(
        point_raw
    ).hexdigest()

    with (
        diagnostic._exclusive_output_claim(output) as claim,
        pytest.raises(ValueError, match="point ancestry mismatch"),
    ):
        diagnostic._resume_fd_points(
            output,
            status,
            diagnostic._displacement_plan(6),
            claim=claim,
        )


def test_resume_rejects_exact_matrix_ancestry_mismatch(tmp_path, monkeypatch):
    output, _receipt = _complete_fake_fd(tmp_path, monkeypatch, "matrix-ancestry")
    status_path = output / "status.json"
    status = json.loads(status_path.read_text())
    status["state"] = "running"
    status["confirmation_passed"] = False
    for field in ("receipt", "receipt_sha256", "finished_utc"):
        status.pop(field)
    status_path.write_bytes(diagnostic._json_bytes(status))
    (output / diagnostic.TERMINAL).unlink()
    (output / "receipt.json").unlink()
    matrix_receipt_path = output / "matrices" / "receipt.json"
    matrix_receipt = json.loads(matrix_receipt_path.read_text())
    matrix_receipt["point_receipts"] = []
    matrix_receipt_path.write_bytes(diagnostic._json_bytes(matrix_receipt))

    with pytest.raises(ValueError, match="resumed matrix-analysis receipt"):
        diagnostic.run_finite_difference(
            tmp_path / "preflight", tmp_path / "reference", output
        )


@pytest.mark.parametrize(
    "gate",
    ["center-fmax", "center-energy", "center-s2", "displaced-s2", "density"],
)
def test_declared_point_gate_rejection_completes_with_all_artifacts(
    tmp_path, monkeypatch, gate
):
    _fake_fd_environment(tmp_path, monkeypatch)
    center_density = np.arange(1, 9, dtype=float).reshape(2, 2, 2)
    center_sha = hashlib.sha256(center_density.astype("<f8").tobytes()).hexdigest()
    displaced_density = center_density * 2.0 if gate == "density" else center_density
    displaced_sha = hashlib.sha256(
        displaced_density.astype("<f8").tobytes()
    ).hexdigest()

    def evaluate(cluster, _settings, *, central_density):
        is_center = central_density is None
        returned = center_density if is_center else displaced_density
        record = _point_record(
            density_sha=center_sha,
            s2=(
                0.7602
                if gate == "center-s2" and is_center
                else 0.77
                if gate == "displaced-s2" and not is_center
                else 0.76
            ),
        )
        record["density_final_sha256"] = center_sha if is_center else displaced_sha
        gradient = np.zeros((6, 3))
        if gate == "center-fmax" and is_center:
            threshold = diagnostic.FINITE_DIFFERENCE_CONTRACT["thresholds"][
                "center_fmax_ev_per_angstrom_exclusive_maximum"
            ]
            gradient[0, 0] = (
                threshold * diagnostic.BOHR_TO_ANGSTROM / diagnostic.HARTREE_TO_EV
            )
        if gate == "center-energy" and is_center:
            record["electronic_hartree"] = -9.999
        record["geometry_fingerprint"] = diagnostic.frequency_geometry_fingerprint(
            cluster
        )
        return record, gradient, returned

    monkeypatch.setattr(diagnostic, "_evaluate_fd_gradient", evaluate)
    output = tmp_path / f"{gate}-rejection"

    receipt = diagnostic.run_finite_difference(
        tmp_path / "preflight", tmp_path / "reference", output
    )

    assert receipt["state"] == "completed"
    assert receipt["confirmation_passed"] is False
    assert receipt["analysis"]["rejection"].startswith("ScientificRejection:")
    assert len(receipt["point_receipts"]) == 73
    assert (output / "matrices/richardson-raw.f64").is_file()
    first_key = diagnostic._displacement_plan(6)[0]["key"]
    density_path = output / "points" / first_key / "density.f64"
    assert density_path.read_bytes() == displaced_density.astype("<f8").tobytes()
    point = json.loads((density_path.parent / "receipt.json").read_text())
    assert point["density_final_sha256"] == displaced_sha
    terminal = json.loads((output / diagnostic.TERMINAL).read_text())
    assert terminal["state"] == "completed"
    assert terminal["confirmation_passed"] is False


def test_external_claim_inode_replacement_is_detected(tmp_path):
    output = tmp_path / "output"
    claim_path = tmp_path / f".{output.name}.d2c-hessian-diagnostic.claim"

    with (
        pytest.raises(RuntimeError, match="external claim was replaced"),
        diagnostic._exclusive_output_claim(output) as claim,
    ):
        displaced = tmp_path / "displaced-claim"
        claim_path.rename(displaced)
        claim_path.write_bytes(b"replacement")
        with (
            pytest.raises(
                diagnostic.ConcurrentRunError,
                match="stable output parent is already claimed",
            ),
            diagnostic._exclusive_output_claim(output),
        ):
            pass
        claim.verify()


def test_fd_identity_rejects_source_not_equal_to_committed_bytes(monkeypatch):
    identity = {
        "git_sha": "a" * 40,
        "script_sha256": hashlib.sha256(b"committed").hexdigest(),
    }
    monkeypatch.setattr(diagnostic, "_git_identity", lambda: identity)
    monkeypatch.setattr(
        diagnostic.subprocess,
        "run",
        lambda *_args, **_kwargs: SimpleNamespace(stdout=b"different"),
    )

    with pytest.raises(RuntimeError, match="not the committed Git source"):
        diagnostic._fd_execution_identity()


def test_held_legacy_claim_blocks_before_identity_or_expensive_work(
    tmp_path, monkeypatch
):
    output = tmp_path / "legacy"
    output.mkdir()
    legacy = output / ".preflight.lock"
    legacy.write_bytes(b"legacy")
    calls = []
    monkeypatch.setattr(
        diagnostic, "_fd_execution_identity", lambda: calls.append("identity")
    )

    with legacy.open("r+b") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(
            diagnostic.ConcurrentRunError,
            match="legacy output-local claim is still held",
        ):
            diagnostic.run_finite_difference(
                tmp_path / "preflight", tmp_path / "reference", output
            )

    assert calls == []
    assert not (output / "status.json").exists()
    assert not (output / diagnostic.TERMINAL).exists()
