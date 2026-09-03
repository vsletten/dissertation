from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
MODULE_PATH = SCRIPTS / "e3a_campaign.py"
SPEC = importlib.util.spec_from_file_location("e3a_campaign", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
campaign = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = campaign
SPEC.loader.exec_module(campaign)


def test_potential_uses_explicit_lorentz_berthelot_pairs() -> None:
    lines = campaign.potential_lines()
    assert not any("pair_modify" in line for line in lines)
    pair_lines = [line for line in lines if line.startswith("pair_coeff")]
    # Ten types: ten self pairs plus C(10, 2) unlike pairs.
    assert len(pair_lines) == 55


def test_minimize_input_persists_relaxed_state_and_local_freeze(tmp_path: Path) -> None:
    path = tmp_path / "in.min"
    campaign.write_minimize_input(
        path,
        data_name="seed.data",
        output_prefix="relaxed",
        ftol=0.01,
        max_steps=500,
        local_mobile_ids=(675, 741, 743),
    )
    text = path.read_text(encoding="utf-8")
    assert "group           local-mobile id 675 741 743" in text
    assert "fix             local_hold local-frozen setforce 0.0 0.0 0.0" in text
    assert "unfix           local_hold" in text
    assert "minimize        0.0 0.01 500 5000" in text
    assert "write_data      relaxed.data" in text
    assert (
        "write_dump      all custom relaxed.dump id type x y z modify sort id" in text
    )


def test_parse_dump_and_write_neb_coordinates(tmp_path: Path) -> None:
    dump = tmp_path / "relaxed.dump"
    dump.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        "ITEM: ATOMS id type x y z\n2 1 4 5 6\n1 1 1 2 3\n",
        encoding="utf-8",
    )
    coords = campaign.parse_dump(dump)
    assert coords == {2: (4.0, 5.0, 6.0), 1: (1.0, 2.0, 3.0)}
    final = tmp_path / "final.coords"
    campaign.write_final_coords(final, coords)
    assert final.read_text(encoding="utf-8").splitlines() == [
        "2",
        "1 1 2 3",
        "2 4 5 6",
    ]


def test_parse_neb_receipt(tmp_path: Path) -> None:
    screen = tmp_path / "screen.log"
    row = [
        13000,
        0.009,
        0.001,
        1,
        2,
        3,
        42.5,
        41.5,
        1,
        0,
        -100,
        0.1,
        -90,
        0.2,
        -80,
        0.3,
        -70,
        0.4,
        -60,
        0.5,
        -70,
        0.6,
        -80,
        1,
        -99,
    ]
    screen.write_text(" ".join(str(value) for value in row) + "\n", encoding="utf-8")
    result = campaign.parse_neb(screen, 0.01)
    assert result["computed_forward_barrier_kcal_mol"] == 42.5
    assert result["computed_reverse_barrier_kcal_mol"] == 41.5
    assert result["converged_to_requested_ftol"] is True


def test_neb_writer_uses_compatible_every_interval(tmp_path: Path) -> None:
    path = tmp_path / "in.neb"
    campaign.write_neb_input(
        path,
        data_name="initial.data",
        final_name="final.coords",
        ftol=0.01,
        relax_steps=100,
        climb_steps=200,
    )
    assert (
        "neb             1.0e-10 0.01 100 200 100 final final.coords"
        in path.read_text(encoding="utf-8")
    )


WORKBOOK = Path(
    "/mnt/data/vsletten/dissertation-data/e3a-classical-neb-barriers-20260902/"
    "sources/nteme-2022-supplement/Results Of MD Simulations.xlsx"
)


@pytest.mark.skipif(not WORKBOOK.exists(), reason="hash-pinned workbook unavailable")
def test_real_workbook_exposes_campaign_models() -> None:
    pytest.importorskip("openpyxl")
    models, _source = campaign.build_campaign_models(WORKBOOK)
    assert tuple(model.name for model in models) == campaign.MODEL_ORDER
    replication = models[0]
    assert replication.metadata["published_barrier_kcal_mol"] == pytest.approx(
        67.644151
    )
    assert abs(campaign.contract.net_charge(replication.atoms)) <= 1.0e-8
    zone_dimensions = [
        (
            model.metadata["transformation"]["radius_a_angstrom"],
            model.metadata["transformation"]["radius_b_angstrom"],
            model.metadata["transformation"]["peak_opening_angstrom"],
        )
        for model in models
        if model.name.startswith("extended-defect-zone")
    ]
    assert zone_dimensions == [(6.0, 3.0, 0.5), (8.0, 5.0, 1.0), (10.0, 7.0, 1.5)]


def test_prepare_output_dir_refuses_nonempty_directory(tmp_path: Path) -> None:
    stale = tmp_path / "out"
    stale.mkdir()
    (stale / "old-model").mkdir()
    (stale / "old-model" / "result.json").write_text("{}\n", encoding="utf-8")
    with pytest.raises(ValueError, match="new or empty directory"):
        campaign.prepare_output_dir(stale)
    assert (stale / "old-model" / "result.json").is_file()


def test_prepare_output_dir_creates_missing_or_empty_directory(tmp_path: Path) -> None:
    missing = tmp_path / "nested" / "out"
    campaign.prepare_output_dir(missing)
    assert missing.is_dir()
    assert list(missing.iterdir()) == []
    campaign.prepare_output_dir(missing)
    assert missing.is_dir()
    assert list(missing.iterdir()) == []


def test_prepare_output_dir_refuses_existing_file(tmp_path: Path) -> None:
    path = tmp_path / "out"
    path.write_text("stale\n", encoding="utf-8")
    with pytest.raises(ValueError, match="new or empty directory"):
        campaign.prepare_output_dir(path)
    assert path.read_text(encoding="utf-8") == "stale\n"


def test_write_campaign_manifest_excludes_self_even_when_leftover_exists(
    tmp_path: Path,
) -> None:
    out = tmp_path / "campaign"
    out.mkdir()
    nested = out / "model" / "result.json"
    nested.parent.mkdir()
    nested.write_text("{}\n", encoding="utf-8")
    (out / "campaign-result.json").write_text('{"schema": "e3a"}\n', encoding="utf-8")
    leftover = out / "manifest.json"
    leftover.write_text('{"files": [{"path": "stale"}]}\n', encoding="utf-8")

    listed = campaign.write_campaign_manifest(out)

    paths = [entry["path"] for entry in listed]
    assert "manifest.json" not in paths
    assert paths == ["campaign-result.json", "model/result.json"]
    payload = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    assert [entry["path"] for entry in payload["files"]] == paths
    for entry in payload["files"]:
        path = out / entry["path"]
        assert entry["bytes"] == path.stat().st_size
        assert entry["sha256"] == campaign.sha256(path)
