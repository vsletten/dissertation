"""Pure contracts for the resumable A2a coupled-coordinate scan."""

from __future__ import annotations

import json
from dataclasses import replace
from itertools import pairwise
from pathlib import Path

import numpy as np
import pytest

from quarry.clusters import Cluster
from quarry.pipeline import DftSettings, frequency_geometry_fingerprint
from scripts import a2a_cleavage_scan as scan


def endpoint(name: str, *, si_obr: float, hw_obr: float, h17_obr: float) -> Cluster:
    symbols = ["C"] * 18
    symbols[scan.BRIDGE_INDEX] = "O"
    symbols[scan.SI_INDEX] = "Si"
    symbols[scan.OW_INDEX] = "O"
    symbols[scan.HW_INDEX] = "H"
    symbols[scan.DEGENERATE_HW_INDEX] = "H"
    coords = np.array([[20.0 + 3.0 * index, 10.0, 0.0] for index in range(18)])
    coords[scan.BRIDGE_INDEX] = [0.0, 0.0, 0.0]
    coords[scan.SI_INDEX] = [si_obr, 0.0, 0.0]
    coords[scan.OW_INDEX] = [0.0, 4.0, 0.0]
    coords[scan.HW_INDEX] = [0.0, hw_obr, 0.0]
    coords[scan.DEGENERATE_HW_INDEX] = [0.0, 0.0, h17_obr]
    return Cluster(name, symbols, coords)


def owner_identity(cluster: Cluster, _attacker_index: int):
    owners = (
        (scan.HW_INDEX, scan.OW_INDEX),
        (scan.DEGENERATE_HW_INDEX, scan.BRIDGE_INDEX),
    )
    return ((True, True, True), owners, ((scan.BRIDGE_INDEX, scan.SI_INDEX),))


def test_default_grid_contains_exact_endpoints_extension_and_adjacent_traversal():
    grid = scan.build_grid(1.95, 3.73, 2.98, 1.85)

    assert len(grid.si_obr_targets_a) == scan.DEFAULT_AXIS_POINTS == 9
    assert grid.product_index == 7
    assert grid.si_obr_targets_a[0] == 1.95
    assert grid.si_obr_targets_a[grid.product_index] == 3.73
    assert grid.hw_obr_targets_a[0] == 2.98
    assert grid.hw_obr_targets_a[grid.product_index] == 1.85
    assert grid.si_obr_targets_a[-1] > 3.73
    assert grid.hw_obr_targets_a[-1] < 1.85
    assert len(grid.traversal) == 81
    assert grid.traversal[0] == (0, 0)
    assert grid.traversal[8] == (0, 8)
    assert grid.traversal[9] == (1, 8)
    assert all(
        abs(left[0] - right[0]) + abs(left[1] - right[1]) == 1
        for left, right in pairwise(grid.traversal)
    )


def test_coordinate_mapping_explicitly_selects_h16_and_rejects_h17_degeneracy(
    monkeypatch,
):
    intermediate = endpoint("I", si_obr=1.95, hw_obr=2.98, h17_obr=0.98)
    product = endpoint("P", si_obr=3.73, hw_obr=1.85, h17_obr=0.96)
    monkeypatch.setattr(scan.a2, "endpoint_identity", owner_identity)

    mapping = scan.validate_coordinate_mapping(intermediate, product)

    assert mapping.si_index == 1
    assert mapping.bridge_index == 0
    assert mapping.ow_index == 15
    assert mapping.hw_index == 16
    assert mapping.si_obr_pair == (1, 0)
    assert mapping.hw_obr_pair == (16, 0)
    assert mapping.hw_contraction_a == pytest.approx(1.13)

    with pytest.raises(ValueError, match="H17"):
        scan.validate_coordinate_mapping(
            intermediate,
            product,
            hw_index=scan.DEGENERATE_HW_INDEX,
        )

    degenerate = endpoint("P", si_obr=3.73, hw_obr=2.96, h17_obr=0.70)
    with pytest.raises(ValueError, match="degenerate H16/H17"):
        scan.validate_coordinate_mapping(intermediate, degenerate)


def fake_energy_checkpoint(calls: list[str]):
    def checkpoint(
        path: Path, cluster: Cluster, settings: DftSettings, method: str
    ) -> float:
        expected = {
            "method": method,
            "geometry_fingerprint": frequency_geometry_fingerprint(cluster),
            "settings_fingerprint": scan.frequency_settings_fingerprint(settings),
        }
        if path.exists():
            payload = json.loads(path.read_text())
            if any(payload.get(key) != value for key, value in expected.items()):
                raise ValueError("fake energy checkpoint drift")
            return float(payload["electronic_hartree"])
        calls.append(path.parent.name)
        value = -100.0 + float(cluster.coords[scan.SI_INDEX, 0]) / 100.0
        scan.atomic_json(path, {**expected, "electronic_hartree": value})
        return value

    return checkpoint


def test_scan_resumes_without_recompute_rejects_drift_and_chains_adjacent_seeds(
    monkeypatch, tmp_path
):
    intermediate = endpoint("I", si_obr=1.95, hw_obr=2.98, h17_obr=0.98)
    product = endpoint("P", si_obr=3.73, hw_obr=1.85, h17_obr=0.96)
    settings = DftSettings(xc="r2scan", basis="def2-mtzvpp", composite="r2scan3c")
    monkeypatch.setattr(scan.a2, "endpoint_identity", owner_identity)
    relax_calls: list[tuple[float, float]] = []
    energy_calls: list[str] = []

    def relax(cluster, _settings, *, fixed_distances, **_kwargs):
        targets = {(left, right): target for left, right, target in fixed_distances}
        si_target = targets[(scan.SI_INDEX, scan.BRIDGE_INDEX)]
        hw_target = targets[(scan.HW_INDEX, scan.BRIDGE_INDEX)]
        relax_calls.append((si_target, hw_target))
        coords = cluster.coords.copy()
        coords[scan.SI_INDEX] = [si_target, 0.0, 0.0]
        coords[scan.HW_INDEX] = [0.0, hw_target, 0.0]
        return replace(cluster, coords=coords, name=f"cell-{len(relax_calls)}")

    kwargs = {
        "axis_points": 3,
        "max_steps": 5,
        "relax_fn": relax,
        "checkpoint_energy_fn": fake_energy_checkpoint(energy_calls),
    }
    first = scan.run_scan_grid(
        tmp_path / "scan", intermediate, product, settings, **kwargs
    )
    second = scan.run_scan_grid(
        tmp_path / "scan", intermediate, product, settings, **kwargs
    )

    assert len(first) == len(second) == 9
    assert len(relax_calls) == len(energy_calls) == 9
    assert [record.seed_geometry_fingerprint for record in first[1:]] == [
        record.output_geometry_fingerprint for record in first[:-1]
    ]
    first_receipt = json.loads(Path(first[0].receipt_path).read_text())
    assert first_receipt["topology"]["valid_typed_identity"] is True
    assert first_receipt["topology"]["h16_owner"] == scan.OW_INDEX
    assert first_receipt["coarse_grid_point_is_verified_saddle"] is False
    assert all(
        abs(left.row - right.row) + abs(left.column - right.column) == 1
        for left, right in pairwise(first)
    )

    changed_product = replace(product, coords=product.coords.copy())
    changed_product.coords[scan.SI_INDEX, 0] += 0.01
    with pytest.raises(ValueError, match="scan manifest drift"):
        scan.run_scan_grid(
            tmp_path / "scan",
            intermediate,
            changed_product,
            settings,
            **kwargs,
        )


def test_classifier_prioritizes_proton_first_local_minimum():
    energies = np.zeros((5, 5))
    energies[1, 2] = -1.0
    energies[2, :] = 5.0
    topologies = {
        (1, 2): {
            "valid_typed_identity": True,
            "h16_owner": scan.BRIDGE_INDEX,
            "si_obr_bonded": True,
        }
    }

    result = scan.classify_complete_grid(
        energies,
        product_index=3,
        cell_topologies=topologies,
    )

    assert result["outcome"] == scan.PROTON_FIRST_MINIMUM
    assert result["proton_first_minimum"]["cell"] == [1, 2]
    assert result["verified_saddle"] is False


def test_classifier_does_not_infer_proton_transfer_from_grid_indices():
    energies = np.zeros((5, 5))
    energies[1, 2] = -1.0

    result = scan.classify_complete_grid(energies, product_index=3)

    assert result["outcome"] == scan.BARRIERLESS_SHELF


def test_classifier_does_not_route_through_extension_cells():
    energies = np.zeros((5, 5))
    energies[2, :4] = 4.0
    energies[2, 2] = 3.0
    energies[4, :] = -10.0
    energies[:, 4] = -10.0

    result = scan.classify_complete_grid(energies, product_index=3)

    assert result["minimax_bottleneck"]["height_kj_mol"] == pytest.approx(3.0)
    assert result["outcome"] == scan.INTERIOR_CREST


def test_classifier_finds_interior_minimax_crest_above_threshold():
    energies = np.zeros((5, 5))
    energies[2, :] = 4.0
    energies[2, 2] = 3.0

    result = scan.classify_complete_grid(energies, product_index=3)

    assert result["outcome"] == scan.INTERIOR_CREST
    assert result["crest"]["cell"] == [2, 2]
    assert result["crest"]["height_kj_mol"] == pytest.approx(3.0)
    assert result["verified_saddle"] is False
    assert "not-verified-saddle" in result["candidate_kind"]


def test_classifier_calls_smooth_valley_barrierless_shelf():
    energies = np.zeros((5, 5))
    energies[1, 1] = 1.0

    result = scan.classify_complete_grid(energies, product_index=3)

    assert result["outcome"] == scan.BARRIERLESS_SHELF
    assert result["verified_saddle"] is False


def test_classifier_treats_threshold_equality_as_barrierless():
    energies = np.zeros((5, 5))
    energies[2, :] = 3.0
    energies[2, 2] = scan.DEFAULT_BARRIER_THRESHOLD_KJ_MOL

    result = scan.classify_complete_grid(energies, product_index=3)

    assert result["minimax_bottleneck"]["height_kj_mol"] == pytest.approx(2.0)
    assert result["outcome"] == scan.BARRIERLESS_SHELF


def test_classifier_does_not_promote_boundary_bottleneck_to_interior_crest():
    energies = np.zeros((5, 5))
    energies[2, :] = 4.0
    energies[2, 0] = 3.0

    result = scan.classify_complete_grid(energies, product_index=3)

    assert result["minimax_bottleneck"]["cell"] == [2, 0]
    assert result["minimax_bottleneck"]["height_kj_mol"] == pytest.approx(3.0)
    assert result["outcome"] == scan.BARRIERLESS_SHELF
    assert result["verified_saddle"] is False


def test_parser_carries_gpu_memory_override_and_fixed_scan_location(tmp_path):
    args = scan.parser().parse_args(
        [
            "--intermediate-reference",
            str(tmp_path / "I.xyz"),
            "--product-reference",
            str(tmp_path / "P.xyz"),
            "--run-dir",
            str(tmp_path),
            "--gpu",
            "--gpu-mem-gb",
            "12",
        ]
    )

    assert args.gpu is True
    assert args.gpu_mem_gb == 12.0
    assert scan.scan_root(args.run_dir) == tmp_path / "cleavage" / "coupled-scan-v1"
