from __future__ import annotations

import hashlib
import inspect
import json
import math
import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest

from quarry.clusters import Cluster
from quarry.pipeline import BOHR_TO_ANGSTROM, DftSettings
from scripts import d2c_mode_projected_energy_diagnostic as diagnostic

SOURCE_SNAPSHOT = {
    "git_commit": "a" * 40,
    "git_tree": "c" * 40,
    "source_files": {
        "qm/scripts/d2c_mode_projected_energy_diagnostic.py": "b" * 64,
        "qm/scripts/d2c_hessian_symmetry_diagnostic.py": "d" * 64,
    },
}
IDENTITY = {
    "source_snapshot": SOURCE_SNAPSHOT,
    "python": "3.13.1",
    "dependencies": {
        "numpy": "2.5.2",
        "pyscf": "2.14.0",
        "pyscf-dispersion": "1.5.0",
        "scipy": "1.18.0",
    },
    "python_executable": {
        "origin": "/fake/python",
        "sha256": "9" * 64,
        "byte_count": 1,
    },
    "executable_modules": {"fake": {"sha256": "e" * 64}},
    "native_payloads": {"fake.so": {"sha256": "f" * 64}},
    "accelerator": "none",
    "cpu_only": True,
}


def _fake_evidence() -> dict:
    unstable = np.zeros(18)
    unstable[0] = 1.0
    soft = np.zeros(18)
    soft[3] = 1.0
    masses = np.array([4.0, 9.0, 16.0, 25.0, 36.0, 49.0])
    cluster = Cluster(
        name="fake",
        symbols=["C", "O", "H", "O", "H", "H"],
        coords=np.arange(18, dtype=float).reshape(6, 3),
        spin=1,
    )
    settings = DftSettings(
        xc="pwb6k",
        basis="def2-svp",
        dispersion="d3bj",
        grid_level=5,
        density_fit=True,
        use_gpu=False,
    )
    return {
        "ancestry": {
            "analytic": {
                "receipt_sha256": "1" * 64,
                "accepted_campaign_result": False,
                "internal_symmetry_policy_passed": False,
                "internal_verdicts": [
                    {
                        "case": name,
                        "electronic_symmetry_accepted": False,
                        "total_symmetry_accepted": False,
                    }
                    for name in ("A", "B", "C", "D")
                ],
            },
            "finite_difference": {
                "receipt_sha256": "2" * 64,
                "terminal_sha256": "3" * 64,
                "status_sha256": "4" * 64,
                "confirmation_passed": False,
                "accepted_campaign_result": False,
            },
            "half_step": {
                "receipt_sha256": diagnostic.SOURCE_RECEIPT_SHA256,
                "terminal_sha256": diagnostic.SOURCE_TERMINAL_SHA256,
                "status_sha256": diagnostic.SOURCE_STATUS_SHA256,
                "matrix_receipt_sha256": diagnostic.SOURCE_MATRIX_RECEIPT_SHA256,
                "richardson_symmetric_matrix_sha256": diagnostic.SOURCE_MATRIX_SHA256,
                "confirmation_passed": False,
                "accepted_campaign_result": False,
            },
        },
        "bundle": {
            "root": "/fake/bundle",
            "manifest_sha256": diagnostic.BUNDLE_MANIFEST_SHA256,
            "file_count": diagnostic.BUNDLE_FILE_COUNT,
            "total_bytes": diagnostic.BUNDLE_TOTAL_BYTES,
        },
        "reference_binding": {"route": diagnostic.ROUTE},
        "reference_binding_sha256": "5" * 64,
        "authoritative_center": {
            "receipt_sha256": diagnostic.SOURCE_CENTER_RECEIPT_SHA256,
            "electronic_hartree": diagnostic.SOURCE_CENTER_ENERGY_HARTREE,
            "spin_square": diagnostic.SOURCE_CENTER_S2,
            "density_sha256": diagnostic.SOURCE_CENTER_DENSITY_SHA256,
            "energy_exclusive_tolerance_hartree": diagnostic.CONTRACT["thresholds"][
                "center_energy_vs_source_exclusive_maximum_delta_hartree"
            ],
            "s2_exclusive_tolerance": diagnostic.CONTRACT["thresholds"][
                "center_s2_vs_source_exclusive_maximum_delta"
            ],
        },
        "transition_state_geometry_fingerprint": "6" * 64,
        "masses_amu": masses.tolist(),
        "masses_sha256": hashlib.sha256(masses.astype("<f8").tobytes()).hexdigest(),
        "modes": {
            "unstable": {
                "index": 0,
                "mass_weighted_vector": unstable.tolist(),
                "mass_weighted_vector_sha256": "7" * 64,
                "source_eigenvalue_hartree_per_bohr2_amu": -0.02,
                "source_signed_wavenumber_cm": -700.0,
            },
            "lowest_positive": {
                "index": 1,
                "mass_weighted_vector": soft.tolist(),
                "mass_weighted_vector_sha256": "8" * 64,
                "source_eigenvalue_hartree_per_bohr2_amu": 0.0001,
                "source_signed_wavenumber_cm": 50.0,
            },
        },
        "transition_state": cluster,
        "settings": settings,
    }


def _seed_preflight(root: Path, evidence: dict) -> dict:
    root.mkdir()
    payload = diagnostic._preflight_payload(IDENTITY, evidence)
    (root / "preflight.json").write_bytes(diagnostic._json_bytes(payload))
    return payload


def _seed_production_preflight(root: Path) -> dict:
    try:
        identity = diagnostic._current_identity()
    except RuntimeError as exc:
        if "clean Git worktree" in str(exc):
            pytest.skip("production-path identity requires a committed clean checkout")
        raise
    evidence = diagnostic._load_source_evidence(diagnostic.SOURCE_ROOT)
    root.mkdir()
    payload = diagnostic._preflight_payload(identity, evidence)
    (root / "preflight.json").write_bytes(diagnostic._json_bytes(payload))
    return payload


def _fake_scf_boundary(
    *, fail_at: int | None = None, center_energy_offset: float = 0.0
):
    center_density = np.array(
        [[[1.0, 0.1], [0.1, 0.7]], [[0.9, 0.05], [0.05, 0.6]]], dtype=float
    )
    plan = diagnostic._point_plan(_fake_evidence())
    instances = []
    dm0_receipts = []

    class FakeScf:
        converged = True

        def __init__(self, index: int):
            self.index = index
            self.conv_tol = None
            self.max_cycle = None
            self.grids = SimpleNamespace(level=5, coords=np.empty((199560, 0)))
            self.final_density = center_density.copy()
            if index:
                self.final_density.flat[index % self.final_density.size] += 1.0e-4

        def kernel(self, dm0=None):
            if fail_at is not None and self.index + 1 == fail_at:
                raise RuntimeError("fake scalar failure")
            if self.index == 0:
                assert dm0 is None
            else:
                assert dm0 is not None
                raw = np.ascontiguousarray(dm0, dtype="<f8").tobytes()
                dm0_receipts.append((id(self), hashlib.sha256(raw).hexdigest()))
                assert raw == center_density.astype("<f8").tobytes()
            definition = plan[self.index]
            q = definition["q_bohr_sqrt_amu"]
            mode = definition["mode"]
            if mode is None:
                return diagnostic.SOURCE_CENTER_ENERGY_HARTREE + center_energy_offset
            curvature = -0.02 if mode == "unstable" else 0.0001
            gradient = 0.0001 if mode == "unstable" else -1.0e-6
            return (
                diagnostic.SOURCE_CENTER_ENERGY_HARTREE
                + gradient * q
                + 0.5 * curvature * q * q
                + 1.0e-5 * q**4
            )

        def spin_square(self):
            return diagnostic.SOURCE_CENTER_S2 + self.index * 1.0e-5, 2.009

        def make_rdm1(self):
            return self.final_density.copy()

        def density_fit(self):
            return self

    def make_scf(_mol, _settings=None):
        instance = FakeScf(len(instances))
        instances.append(instance)
        return instance

    return make_scf, instances, dm0_receipts, center_density


def _install_fake_runtime(
    monkeypatch,
    *,
    fail_at: int | None = None,
    center_energy_offset: float = 0.0,
):
    make_scf, instances, dm0_receipts, center_density = _fake_scf_boundary(
        fail_at=fail_at, center_energy_offset=center_energy_offset
    )
    monkeypatch.setattr(diagnostic._pyscf_dft, "UKS", make_scf)
    return instances, dm0_receipts, center_density


def _run_fake(tmp_path: Path, monkeypatch, *, fail_at: int | None = None):
    root = tmp_path / "run"
    _seed_production_preflight(root)
    runtime = _install_fake_runtime(monkeypatch, fail_at=fail_at)
    receipt = diagnostic.run(root)
    return root, receipt, runtime


def _observation(
    energy: float,
    *,
    center: bool,
    s2: float = diagnostic.SOURCE_CENTER_S2,
    density_delta: float = 0.0,
    root_verified: bool | None = None,
) -> dict:
    if root_verified is None:
        root_verified = center
    return {
        "backend": "pyscf-cpu",
        "electronic_hartree": energy,
        "scf_converged": True,
        "grid_point_count": 199560,
        "spin_2s": 1,
        "spin_square": [s2, 2.009],
        "scf_tolerance": 1.0e-12,
        "scf_max_cycle": 150,
        "density_initial_guess_sha256": None if center else "a" * 64,
        "density_final_sha256": "a" * 64 if center else "b" * 64,
        "density_shape": [2, 2, 2],
        "final_frobenius_norm": 1.5,
        "center_frobenius_norm": 1.5,
        "final_vs_center_normalized_frobenius_delta": density_delta,
        "robust_root_continuity_metric": (
            "center-self" if center else "unavailable-beyond-s2-and-density-proxies"
        ),
        "robust_root_continuity_verified": root_verified,
        "gradient_computed": False,
        "hessian_computed": False,
    }


def _analysis_points() -> list[dict]:
    points = []
    for definition in diagnostic._point_plan(_fake_evidence()):
        q = definition["q_bohr_sqrt_amu"]
        mode = definition["mode"]
        curvature = 0.0 if mode is None else -0.02 if mode == "unstable" else 0.0001
        energy = diagnostic.SOURCE_CENTER_ENERGY_HARTREE + 0.5 * curvature * q * q
        points.append(
            {
                "point_key": definition["point_key"],
                "geometry_sha256": hashlib.sha256(
                    np.ascontiguousarray(definition["coords"], dtype="<f8").tobytes()
                ).hexdigest(),
                "observation": _observation(energy, center=mode is None),
            }
        )
    return points


def _restore(monkeypatch, obj, name: str, value):
    monkeypatch.setattr(obj, name, value)


def _replace_failure_for_status(
    monkeypatch, state: str, *, completed_point_count: int | None = None
):
    original = diagnostic.hessian_diagnostic.os.replace

    def injected(source, destination, *, src_dir_fd=None, dst_dir_fd=None):
        if destination == "status.json" and src_dir_fd is not None:
            descriptor = diagnostic.hessian_diagnostic.os.open(
                source,
                diagnostic.hessian_diagnostic.os.O_RDONLY,
                dir_fd=src_dir_fd,
            )
            try:
                raw = diagnostic.hessian_diagnostic.os.read(descriptor, 4 * 1024 * 1024)
            finally:
                diagnostic.hessian_diagnostic.os.close(descriptor)
            payload = json.loads(raw)
            if payload["state"] == state and (
                completed_point_count is None
                or len(payload["completed_points"]) == completed_point_count
            ):
                raise OSError(f"{state} status boundary fault")
        return original(
            source,
            destination,
            src_dir_fd=src_dir_fd,
            dst_dir_fd=dst_dir_fd,
        )

    monkeypatch.setattr(diagnostic.hessian_diagnostic.os, "replace", injected)
    return original


def _noreplace_failure_for(monkeypatch, artifact: str):
    original = diagnostic.hessian_diagnostic._renameat2_noreplace_at

    def injected(directory_fd, source, destination):
        if destination == artifact:
            raise OSError(f"{artifact} publication fault")
        return original(directory_fd, source, destination)

    monkeypatch.setattr(
        diagnostic.hessian_diagnostic, "_renameat2_noreplace_at", injected
    )
    return original


def test_real_source_chain_reconstructs_modes_and_center_without_calculator(
    monkeypatch,
):
    monkeypatch.setattr(
        diagnostic._pyscf_dft,
        "UKS",
        lambda *_args, **_kwargs: pytest.fail("source validation entered calculator"),
    )
    evidence = diagnostic._load_source_evidence(diagnostic.SOURCE_ROOT)
    assert evidence["authoritative_center"] == {
        "receipt_sha256": diagnostic.SOURCE_CENTER_RECEIPT_SHA256,
        "electronic_hartree": -190.15534901824873,
        "spin_square": 0.7593962152748404,
        "density_sha256": diagnostic.SOURCE_CENTER_DENSITY_SHA256,
        "energy_exclusive_tolerance_hartree": 1.0e-8,
        "s2_exclusive_tolerance": 1.0e-4,
    }
    assert evidence["ancestry"]["analytic"]["accepted_campaign_result"] is False
    assert evidence["ancestry"]["finite_difference"]["confirmation_passed"] is False
    assert evidence["ancestry"]["half_step"]["confirmation_passed"] is False
    assert evidence["bundle"]["file_count"] == 29
    assert evidence["bundle"]["total_bytes"] == 71934
    assert {
        name: item["mass_weighted_vector_sha256"]
        for name, item in evidence["modes"].items()
    } == diagnostic.MODE_HASHES


def test_point_plan_is_exactly_nine_unique_scalar_geometries():
    plan = diagnostic._point_plan(_fake_evidence())
    hashes = {
        hashlib.sha256(
            np.ascontiguousarray(point["coords"], dtype="<f8").tobytes()
        ).hexdigest()
        for point in plan
    }
    assert tuple(point["point_key"] for point in plan) == diagnostic.POINT_ORDER
    assert len(plan) == len(hashes) == 9
    assert (
        max(point.get("maximum_atom_displacement_bohr", 0.0) for point in plan) < 0.10
    )


def test_mass_weighted_displacement_converts_bohr_to_angstrom():
    evidence = _fake_evidence()
    plan = {point["point_key"]: point for point in diagnostic._point_plan(evidence)}
    center = evidence["transition_state"].coords
    assert plan["unstable-p0.10"]["coords"][0, 0] - center[0, 0] == pytest.approx(
        0.10 / math.sqrt(4.0) * BOHR_TO_ANGSTROM
    )


def test_bundle_inventory_rejects_entry_30_before_unbounded_traversal(tmp_path):
    root = tmp_path / "bundle"
    root.mkdir()
    for index in range(30):
        (root / f"{index:02d}").write_bytes(b"x")
    with pytest.raises(ValueError, match="exceeds 29 files"):
        diagnostic._regular_tree_inventory(root)


def test_bundle_inventory_rejects_byte_71935_before_read(tmp_path):
    root = tmp_path / "bundle"
    root.mkdir()
    (root / "oversized").write_bytes(b"x" * (diagnostic.BUNDLE_TOTAL_BYTES + 1))
    with pytest.raises(ValueError, match="byte budget"):
        diagnostic._regular_tree_inventory(root)


def test_source_bundle_path_has_no_unbounded_legacy_verifier():
    source = inspect.getsource(diagnostic._load_source_evidence)
    inventory = inspect.getsource(diagnostic._regular_tree_inventory)
    verifier = inspect.getsource(diagnostic._verify_bounded_bundle)
    assert "verify_bundle(" not in source
    assert "rglob" not in source + inventory + verifier
    assert "read_bytes" not in source + inventory + verifier


def test_descriptor_relative_bundle_read_rejects_intermediate_symlink(tmp_path):
    root = tmp_path / "bundle"
    outside = tmp_path / "outside"
    root.mkdir()
    outside.mkdir()
    (outside / "payload.json").write_text("{}\n")
    (root / "route").symlink_to(outside, target_is_directory=True)
    root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
    try:
        with pytest.raises(ValueError, match="remain below the bundle root"):
            diagnostic._read_bundle_relative_at(
                root_fd,
                Path("route/payload.json"),
                label="test bundle member",
                maximum_bytes=16,
            )
    finally:
        os.close(root_fd)


def test_scalar_evaluator_uses_fresh_scf_and_exact_center_dm0(monkeypatch):
    evidence = _fake_evidence()
    make_scf, instances, dm0_receipts, center_density = _fake_scf_boundary()
    monkeypatch.setattr(diagnostic._pyscf_dft, "UKS", make_scf)
    center_record, density = diagnostic._evaluate_scalar_energy(
        evidence["transition_state"], evidence["settings"], dm0=None
    )
    displaced = replace(
        evidence["transition_state"],
        coords=diagnostic._point_plan(evidence)[1]["coords"],
    )
    displaced_record, _ = diagnostic._evaluate_scalar_energy(
        displaced, evidence["settings"], dm0=density
    )
    center_sha = hashlib.sha256(center_density.astype("<f8").tobytes()).hexdigest()
    assert len({id(instance) for instance in instances}) == 2
    assert center_record["density_final_sha256"] == center_sha
    assert displaced_record["density_initial_guess_sha256"] == center_sha
    assert dm0_receipts == [(id(instances[1]), center_sha)]
    assert displaced_record["robust_root_continuity_verified"] is False
    assert displaced_record["gradient_computed"] is False
    assert displaced_record["hessian_computed"] is False


def test_public_run_exercises_fixed_production_boundary_with_only_scf_fake(
    tmp_path, monkeypatch
):
    root, receipt, runtime = _run_fake(tmp_path, monkeypatch)
    instances, dm0_receipts, center_density = runtime
    center_raw = (root / "center-density.f64").read_bytes()
    center_sha = hashlib.sha256(center_raw).hexdigest()
    assert len(instances) == 9
    assert len({id(instance) for instance in instances}) == 9
    assert len(dm0_receipts) == 8
    assert {digest for _, digest in dm0_receipts} == {center_sha}
    assert center_raw == center_density.astype("<f8").tobytes()
    assert len(list(root.glob("point-*.json"))) == 9
    assert len(receipt["point_receipts"]) == 9
    assert receipt["analysis"]["point_count"] == 9
    assert receipt["analysis"]["robust_root_continuity_available"] is False
    assert receipt["analysis"]["same_electronic_state_verified"] is False
    assert receipt["diagnostic_passed"] is False
    assert (
        receipt["source_snapshot"]["git_commit"]
        == subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=Path(diagnostic.__file__).resolve().parents[2],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    )
    assert (
        receipt["source_snapshot"]["git_tree"]
        == subprocess.run(
            ["git", "rev-parse", "HEAD^{tree}"],
            cwd=Path(diagnostic.__file__).resolve().parents[2],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    )


def test_center_source_state_mismatch_stops_before_any_displacement(
    tmp_path, monkeypatch
):
    root = tmp_path / "run"
    _seed_production_preflight(root)
    instances, _dm0, _density = _install_fake_runtime(
        monkeypatch, center_energy_offset=1.0e-3
    )

    with pytest.raises(
        diagnostic.ScientificRejection, match="pinned source electronic state"
    ):
        diagnostic.run(root)

    assert len(instances) == 1
    assert json.loads((root / "terminal.json").read_text())["state"] == "failed"
    assert not any(root.glob("point-*.json"))
    assert not (root / "center-density.f64").exists()


def test_run_locked_has_no_evaluator_injection_parameter():
    assert tuple(inspect.signature(diagnostic._run_locked).parameters) == ("claim",)


def test_replaced_evaluator_is_rejected_before_scf(tmp_path, monkeypatch):
    evidence = _fake_evidence()
    root = tmp_path / "run"
    _seed_preflight(root, evidence)
    calls = []
    monkeypatch.setattr(
        diagnostic,
        "_evaluate_scalar_energy",
        lambda *_args, **_kwargs: calls.append(True),
    )
    with pytest.raises(RuntimeError, match="critical local callable was replaced"):
        diagnostic.run(root)
    assert calls == []
    assert set(path.name for path in root.iterdir()) == {"preflight.json"}


def test_replaced_identity_orchestrator_is_rejected_before_scf(tmp_path, monkeypatch):
    evidence = _fake_evidence()
    root = tmp_path / "run"
    _seed_preflight(root, evidence)
    monkeypatch.setattr(diagnostic, "_current_identity", lambda: IDENTITY)
    with pytest.raises(RuntimeError, match="critical local callable was replaced"):
        diagnostic.run(root)
    assert set(path.name for path in root.iterdir()) == {"preflight.json"}


def test_replaced_local_analysis_callable_is_rejected(monkeypatch):
    monkeypatch.setattr(diagnostic, "_analyze", lambda _points: {})
    with pytest.raises(RuntimeError, match="critical local callable was replaced"):
        diagnostic._assert_resident_bindings()


@pytest.mark.parametrize("name", ["_make_scf", "build_mol"])
def test_replaced_critical_alias_is_rejected(name, monkeypatch):
    monkeypatch.setattr(diagnostic, name, lambda *_args, **_kwargs: None)
    with pytest.raises(RuntimeError, match="critical resident callable was replaced"):
        diagnostic._assert_resident_bindings()


@pytest.mark.parametrize(
    "name", ["_claimed_write", "_read_bounded_regular_snapshot_at"]
)
def test_replaced_hessian_helper_is_rejected(name, monkeypatch):
    monkeypatch.setattr(
        diagnostic.hessian_diagnostic,
        name,
        lambda *_args, **_kwargs: None,
    )
    with pytest.raises(
        RuntimeError, match="critical resident helper callable was replaced"
    ):
        diagnostic._assert_resident_bindings()


def test_replaced_critical_module_is_rejected(monkeypatch):
    module_name = diagnostic.hessian_diagnostic.__name__
    monkeypatch.setitem(sys.modules, module_name, ModuleType(module_name))
    with pytest.raises(RuntimeError, match="critical resident module was replaced"):
        diagnostic._assert_resident_bindings()


def test_executable_manifest_attests_hessian_source_and_critical_callables():
    manifest = diagnostic._module_execution_manifest()
    hessian_name = diagnostic.hessian_diagnostic.__name__
    assert manifest[hessian_name]["origin"] == str(
        Path(diagnostic.hessian_diagnostic.__file__).resolve()
    )
    assert "_claimed_write" in manifest[hessian_name]["callables"]
    assert "_make_scf" in manifest[diagnostic._make_scf.__module__]["callables"]
    assert "build_mol" in manifest[diagnostic.build_mol.__module__]["callables"]
    local = manifest[diagnostic.__name__]
    assert local["origin"] == str(Path(diagnostic.__file__).resolve())
    assert "_run_locked" in local["callables"]


def test_hessian_module_foreign_on_disk_source_is_rejected(tmp_path, monkeypatch):
    foreign = tmp_path / "foreign_hessian.py"
    foreign.write_text("VALUE = 1\n")
    monkeypatch.setattr(diagnostic.hessian_diagnostic, "__file__", str(foreign))
    with pytest.raises(RuntimeError, match="module-level code disagrees with source"):
        diagnostic._module_execution_manifest()


def test_preloaded_cpu_module_replacement_is_rejected(monkeypatch):
    module = diagnostic._PRELOADED_CPU_MODULES[0]
    monkeypatch.setitem(sys.modules, module.__name__, ModuleType(module.__name__))
    with pytest.raises(RuntimeError, match="CPU dependency module was replaced"):
        diagnostic._native_payload_manifest()


def test_cpu_identity_manifest_has_only_used_cpu_stack():
    versions = diagnostic._cpu_dependency_versions()
    assert set(versions) == {"numpy", "pyscf", "pyscf-dispersion", "scipy"}
    encoded = json.dumps(versions).lower()
    assert "cuda" not in encoded
    assert "gpu" not in encoded
    assert diagnostic.CONTRACT["method"]["backend"] == "pyscf-cpu"


def test_preflight_rejects_accelerator_identity(tmp_path):
    evidence = _fake_evidence()
    root = tmp_path / "run"
    bad_identity = {**IDENTITY, "accelerator": "cuda"}
    root.mkdir()
    payload = diagnostic._preflight_payload(bad_identity, evidence)
    (root / "preflight.json").write_bytes(diagnostic._json_bytes(payload))
    with (
        diagnostic.hessian_diagnostic._exclusive_output_claim(root) as claim,
        pytest.raises(ValueError, match="not CPU-only"),
    ):
        diagnostic._validated_preflight(claim)


def test_analysis_gates_center_source_and_every_displaced_proxy_and_root():
    analysis = diagnostic._analyze(_analysis_points())
    names = {gate["name"] for gate in analysis["gates"]}
    assert "center:source_energy_delta" in names
    assert "center:source_S2_delta" in names
    for key in diagnostic.POINT_ORDER[1:]:
        assert f"{key}:S2" in names
        assert f"{key}:S2_delta_from_center" in names
        assert f"{key}:density_delta_from_center" in names
        assert f"{key}:robust_root_continuity" in names
    assert analysis["diagnostic_passed"] is False


def test_center_source_tolerance_exact_boundary_fails():
    points = _analysis_points()
    points[0]["observation"]["electronic_hartree"] += diagnostic.CONTRACT["thresholds"][
        "center_energy_vs_source_exclusive_maximum_delta_hartree"
    ]
    analysis = diagnostic._analyze(points)
    gate = next(
        gate
        for gate in analysis["gates"]
        if gate["name"] == "center:source_energy_delta"
    )
    assert gate["passed"] is False


def test_displaced_s2_and_density_exact_boundaries_fail():
    points = _analysis_points()
    thresholds = diagnostic.CONTRACT["thresholds"]
    points[1]["observation"]["spin_square"][0] += thresholds[
        "displaced_vs_center_spin_square_exclusive_maximum_delta"
    ]
    points[1]["observation"]["final_vs_center_normalized_frobenius_delta"] = thresholds[
        "displaced_final_density_vs_center_exclusive_maximum_normalized_frobenius_delta"
    ]
    analysis = diagnostic._analyze(points)
    failed = {gate["name"] for gate in analysis["gates"] if not gate["passed"]}
    assert "unstable-m0.10:S2_delta_from_center" in failed
    assert "unstable-m0.10:density_delta_from_center" in failed


def test_failed_terminal_is_published_before_best_effort_status(tmp_path, monkeypatch):
    root = tmp_path / "run"
    _seed_production_preflight(root)
    _install_fake_runtime(monkeypatch, fail_at=2)
    _replace_failure_for_status(monkeypatch, "failed")
    with pytest.raises(RuntimeError, match="fake scalar failure"):
        diagnostic.run(root)
    assert json.loads((root / "terminal.json").read_text())["state"] == "failed"
    assert json.loads((root / "status.json").read_text())["state"] == "running"


def test_failed_terminal_publication_fault_leaves_running_status(tmp_path, monkeypatch):
    root = tmp_path / "run"
    _seed_production_preflight(root)
    _install_fake_runtime(monkeypatch, fail_at=2)
    _noreplace_failure_for(monkeypatch, "terminal.json")
    with pytest.raises(OSError, match="terminal.json publication fault"):
        diagnostic.run(root)
    assert not (root / "terminal.json").exists()
    assert json.loads((root / "status.json").read_text())["state"] == "running"


def test_point_receipt_publication_fault_commits_failed_terminal(tmp_path, monkeypatch):
    root = tmp_path / "run"
    _seed_production_preflight(root)
    _install_fake_runtime(monkeypatch)
    _noreplace_failure_for(monkeypatch, "point-center.json")

    with pytest.raises(OSError, match="point-center.json publication fault"):
        diagnostic.run(root)

    terminal = json.loads((root / "terminal.json").read_text())
    assert terminal["state"] == "failed"
    assert terminal["completed_point_count"] == 0
    assert terminal["current_point"] == "center"
    assert not (root / "point-center.json").exists()


def test_receipt_publication_fault_commits_failed_terminal(tmp_path, monkeypatch):
    root = tmp_path / "run"
    _seed_production_preflight(root)
    _install_fake_runtime(monkeypatch)
    original_noreplace = _noreplace_failure_for(monkeypatch, "receipt.json")

    with pytest.raises(OSError, match="receipt.json publication fault"):
        diagnostic.run(root)

    _restore(
        monkeypatch,
        diagnostic.hessian_diagnostic,
        "_renameat2_noreplace_at",
        original_noreplace,
    )
    terminal = diagnostic.finalize_if_running(root)
    assert terminal["state"] == "failed"
    assert terminal["completed_point_count"] == 9
    assert not (root / "receipt.json").exists()


def test_finalizer_recovers_point_commit_and_terminal_publication_fault(
    tmp_path, monkeypatch
):
    root = tmp_path / "run"
    _seed_production_preflight(root)
    _install_fake_runtime(monkeypatch)
    original_replace = _replace_failure_for_status(
        monkeypatch, "running", completed_point_count=1
    )
    original_noreplace = _noreplace_failure_for(monkeypatch, "terminal.json")
    with pytest.raises(OSError, match="terminal.json publication fault"):
        diagnostic.run(root)
    assert not (root / "terminal.json").exists()
    assert (root / "point-center.json").is_file()
    _restore(monkeypatch, diagnostic.hessian_diagnostic.os, "replace", original_replace)
    _restore(
        monkeypatch,
        diagnostic.hessian_diagnostic,
        "_renameat2_noreplace_at",
        original_noreplace,
    )
    terminal = diagnostic.finalize_if_running(root)
    assert terminal["state"] == "failed"
    assert terminal["completed_point_count"] == 1
    assert json.loads((root / "status.json").read_text())["state"] == "failed"


def test_finalizer_repairs_failed_terminal_plus_running_status(tmp_path, monkeypatch):
    root = tmp_path / "run"
    _seed_production_preflight(root)
    _install_fake_runtime(monkeypatch, fail_at=2)
    original_replace = _replace_failure_for_status(monkeypatch, "failed")
    with pytest.raises(RuntimeError, match="fake scalar failure"):
        diagnostic.run(root)
    _restore(monkeypatch, diagnostic.hessian_diagnostic.os, "replace", original_replace)
    first = diagnostic.finalize_if_running(root)
    second = diagnostic.finalize_if_running(root)
    assert first == second
    assert first["state"] == "failed"
    assert json.loads((root / "status.json").read_text())["state"] == "failed"


def test_finalizer_recovers_receipt_running_without_terminal(tmp_path, monkeypatch):
    root = tmp_path / "run"
    _seed_production_preflight(root)
    _install_fake_runtime(monkeypatch)
    original_noreplace = _noreplace_failure_for(monkeypatch, "terminal.json")
    with pytest.raises(OSError, match="terminal.json publication fault"):
        diagnostic.run(root)
    assert (root / "receipt.json").is_file()
    assert not (root / "terminal.json").exists()
    assert json.loads((root / "status.json").read_text())["state"] == "running"
    _restore(
        monkeypatch,
        diagnostic.hessian_diagnostic,
        "_renameat2_noreplace_at",
        original_noreplace,
    )
    terminal = diagnostic.finalize_if_running(root)
    assert terminal["state"] == "completed"
    assert json.loads((root / "status.json").read_text())["state"] == "completed"


def test_finalizer_repairs_completed_terminal_plus_running_status(
    tmp_path, monkeypatch
):
    root = tmp_path / "run"
    _seed_production_preflight(root)
    _install_fake_runtime(monkeypatch)
    original_replace = _replace_failure_for_status(monkeypatch, "completed")
    with pytest.raises(OSError, match="completed status boundary fault"):
        diagnostic.run(root)
    assert json.loads((root / "terminal.json").read_text())["state"] == "completed"
    assert json.loads((root / "status.json").read_text())["state"] == "running"
    _restore(monkeypatch, diagnostic.hessian_diagnostic.os, "replace", original_replace)
    terminal = diagnostic.finalize_if_running(root)
    assert terminal["state"] == "completed"
    assert json.loads((root / "status.json").read_text())["state"] == "completed"


def test_running_dead_man_is_idempotent_and_binds_preflight_status(
    tmp_path, monkeypatch
):
    root = tmp_path / "run"
    preflight = _seed_production_preflight(root)
    status = diagnostic._initial_status(
        preflight, (root / "preflight.json").read_bytes()
    )
    (root / "status.json").write_bytes(diagnostic._json_bytes(status))
    first = diagnostic.finalize_if_running(root)
    second = diagnostic.finalize_if_running(root)
    terminal = json.loads((root / "terminal.json").read_text())
    failed_status = (root / "status.json").read_bytes()
    assert first == second == terminal
    assert first["state"] == "failed"
    assert first["dead_man_finalized"] is True
    assert (
        first["preflight_sha256"]
        == hashlib.sha256((root / "preflight.json").read_bytes()).hexdigest()
    )
    assert first["status_sha256"] == hashlib.sha256(failed_status).hexdigest()


def test_finalizer_recomputes_analysis_and_rejects_tampered_receipt(
    tmp_path, monkeypatch
):
    root = tmp_path / "run"
    _seed_production_preflight(root)
    _install_fake_runtime(monkeypatch)
    original_noreplace = _noreplace_failure_for(monkeypatch, "terminal.json")
    with pytest.raises(OSError, match="terminal.json publication fault"):
        diagnostic.run(root)
    _restore(
        monkeypatch,
        diagnostic.hessian_diagnostic,
        "_renameat2_noreplace_at",
        original_noreplace,
    )
    receipt = json.loads((root / "receipt.json").read_text())
    receipt["analysis"]["center_energy_hartree"] += 1.0
    (root / "receipt.json").write_bytes(diagnostic._json_bytes(receipt))
    with pytest.raises(ValueError, match="analysis is not reproducible"):
        diagnostic.finalize_if_running(root)
    assert not (root / "terminal.json").exists()


def test_finalizer_rejects_orphan_evidence(tmp_path, monkeypatch):
    root = tmp_path / "run"
    preflight = _seed_production_preflight(root)
    status = diagnostic._initial_status(
        preflight, (root / "preflight.json").read_bytes()
    )
    (root / "status.json").write_bytes(diagnostic._json_bytes(status))
    (root / "point-center.json").write_text("orphan\n")
    with pytest.raises(ValueError, match="malformed or orphan evidence"):
        diagnostic.finalize_if_running(root)
    assert not (root / "terminal.json").exists()


def test_finalizer_rejects_terminal_timestamp_and_source_hash_tamper(
    tmp_path, monkeypatch
):
    root, _receipt, _runtime = _run_fake(tmp_path, monkeypatch)
    terminal = json.loads((root / "terminal.json").read_text())
    terminal["finished_utc"] = "not-a-time"
    (root / "terminal.json").write_bytes(diagnostic._json_bytes(terminal))
    with pytest.raises(ValueError, match="canonical UTC timestamp"):
        diagnostic.finalize_if_running(root)


def test_finalizer_rejects_terminal_source_hash_tamper(tmp_path, monkeypatch):
    root, _receipt, _runtime = _run_fake(tmp_path, monkeypatch)
    terminal = json.loads((root / "terminal.json").read_text())
    terminal["source_receipt_sha256"] = "0" * 64
    (root / "terminal.json").write_bytes(diagnostic._json_bytes(terminal))
    with pytest.raises(ValueError, match="contract/source hash drifted"):
        diagnostic.finalize_if_running(root)


def test_preflight_tamper_fails_before_scf(tmp_path, monkeypatch):
    root = tmp_path / "run"
    _seed_production_preflight(root)
    payload = json.loads((root / "preflight.json").read_text())
    payload["contract"]["method"]["grid_level"] = 4
    (root / "preflight.json").write_bytes(diagnostic._json_bytes(payload))
    factory, instances, _receipts, _density = _fake_scf_boundary()
    monkeypatch.setattr(diagnostic._pyscf_dft, "UKS", factory)
    with pytest.raises(ValueError, match="preflight contract"):
        diagnostic.run(root)
    assert instances == []


@pytest.mark.parametrize(
    ("extra", "message"),
    [
        (["--threads", "0"], "threads must be between 1 and 16"),
        (["--threads", "17"], "threads must be between 1 and 16"),
        (["--nice", "9"], "nice must be at least 10"),
        (["--gpu", "0"], "GPU/CUDA flags are forbidden"),
    ],
)
def test_cli_resource_controls_fail_before_calculator(extra, message, tmp_path):
    script = Path(diagnostic.__file__).resolve()
    result = subprocess.run(
        [
            sys.executable,
            str(script),
            *extra,
            "--run",
            "--output-root",
            str(tmp_path / "no-run"),
            "--log",
            str(tmp_path / "resource.log"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 2
    assert message in result.stderr


def test_cli_is_cpu_only_and_exposes_no_gradient_hessian_or_accelerator_flags(tmp_path):
    script = Path(diagnostic.__file__).resolve()
    log_path = tmp_path / "help.log"
    result = subprocess.run(
        [
            sys.executable,
            str(script),
            "--help",
            "--threads",
            "16",
            "--nice",
            "10",
            "--log",
            str(log_path),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0
    assert "never evaluates a gradient or hessian" in result.stdout.lower()
    for forbidden in ("--gradient", "--hessian", "--gpu", "--cuda"):
        assert forbidden not in result.stdout.lower()
    assert log_path.is_file()
