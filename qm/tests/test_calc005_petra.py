"""Petra round-trip and production-emission refusal for CALC-005."""

from __future__ import annotations

import shutil
import subprocess
import tomllib
from pathlib import Path

import pytest

from scripts.calc005_si_attachment_verify import validate_petra_boundary

REPO = Path(__file__).resolve().parents[2]
PETRA = REPO / "petra"
DECK = PETRA / "examples/kaolinite.toml"


def test_live_deck_has_exact_reversible_pair_and_unpublished_calc005():
    report = validate_petra_boundary(REPO)
    assert report["adsorb_si"]["center_transition"] == "200->205"
    assert report["desorb_si"]["center_transition"] == "205->200"
    assert report["desorb_si"]["osa_transition"] == "402->404"
    assert report["legacy_si_ladder_kcal_mol"] == [0.0, 6.0, 12.0, 18.0, 24.0]
    assert report["calc005_status"] == "needed"
    assert report["production_value_emitted"] is False
    assert len(report["legacy_actions_cpp_sha256"]) == 64
    assert len(report["legacy_envrn_cpp_sha256"]) == 64


@pytest.mark.parametrize(
    ("old", "new"),
    [
        (
            'name = "adsorb-si"\ncenter = { kind = "Si", state = ["empty"] }\n'
            'guards = [{ distance = 1, state = ["@o_occ"], frozen = false, min = 1 }]',
            'name = "adsorb-si"\ncenter = { kind = "Si", state = ["empty"] }\n'
            'guards = [{ distance = 2, state = ["@o_occ"], frozen = false, min = 1 }]',
        ),
        (
            'select = { kind = "Oss", state = ["Oss.empty", "Oss.si1"], '
            "frozen = false }",
            'select = { kind = "Osa", state = ["Oss.empty", "Oss.si1"], '
            "frozen = false }",
        ),
    ],
)
def test_petra_boundary_rejects_tampered_distance_or_neighbor_kind(tmp_path, old, new):
    deck_target = tmp_path / "petra/examples/kaolinite.toml"
    calculations_target = tmp_path / "qm/CALCULATIONS.md"
    legacy_target = tmp_path / "legacy/cpp-model"
    deck_target.parent.mkdir(parents=True)
    calculations_target.parent.mkdir(parents=True)
    legacy_target.mkdir(parents=True)
    source = DECK.read_text()
    assert source.count(old) == 1
    deck_target.write_text(source.replace(old, new, 1))
    calculations_target.write_text((REPO / "qm/CALCULATIONS.md").read_text())
    for name in ("actions.cpp", "envrn.cpp"):
        (legacy_target / name).write_text(
            (REPO / "legacy/cpp-model" / name).read_text()
        )
    with pytest.raises(RuntimeError, match="exact selector"):
        validate_petra_boundary(tmp_path)


@pytest.mark.parametrize(
    ("name", "old", "new", "message"),
    [
        ("actions.cpp", "state = 205;", "state = 204;", "legacy Si transitions"),
        ("envrn.cpp", "return (x + y);", "return (x + y + 1);", "Check200 ladder"),
    ],
)
def test_petra_boundary_rejects_legacy_source_tampering(
    tmp_path, name, old, new, message
):
    shutil.copytree(REPO / "petra", tmp_path / "petra")
    calculations = tmp_path / "qm/CALCULATIONS.md"
    calculations.parent.mkdir(parents=True)
    calculations.write_text((REPO / "qm/CALCULATIONS.md").read_text())
    legacy = tmp_path / "legacy/cpp-model"
    legacy.mkdir(parents=True)
    for source_name in ("actions.cpp", "envrn.cpp"):
        source = (REPO / "legacy/cpp-model" / source_name).read_text()
        if source_name == name:
            assert old in source
            source = source.replace(old, new, 1)
        (legacy / source_name).write_text(source)
    with pytest.raises(RuntimeError, match=message):
        validate_petra_boundary(tmp_path)


@pytest.mark.petra
def test_temporary_existing_reversible_pair_round_trips_without_calc005_emission(
    tmp_path,
):
    if shutil.which("cargo") is None:
        pytest.skip("cargo not available")
    temporary = tmp_path / "kaolinite-calc005-refusal.toml"
    text = DECK.read_text()
    document = tomllib.loads(text)
    assert {reaction["name"] for reaction in document["reactions"]} >= {
        "adsorb-si",
        "desorb-si",
    }
    temporary.write_text(text.replace("steps = 20000", "steps = 1", 1))
    out = tmp_path / "out"
    result = subprocess.run(
        [
            "cargo",
            "run",
            "--quiet",
            "-p",
            "petra-cli",
            "--",
            str(temporary),
            "--steps",
            "1",
            "--seed",
            "42",
            "--out",
            str(out),
        ],
        cwd=PETRA,
        capture_output=True,
        text=True,
        timeout=600,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert DECK.read_text() == text
    assert validate_petra_boundary(REPO)["production_value_emitted"] is False
