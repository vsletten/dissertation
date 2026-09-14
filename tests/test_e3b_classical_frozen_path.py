from __future__ import annotations

import argparse
import importlib.util
import stat
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

E3B_SPEC = importlib.util.spec_from_file_location(
    "e3b_periodic_dft", SCRIPTS / "e3b_periodic_dft.py"
)
assert E3B_SPEC is not None and E3B_SPEC.loader is not None
e3b = importlib.util.module_from_spec(E3B_SPEC)
sys.modules[E3B_SPEC.name] = e3b
E3B_SPEC.loader.exec_module(e3b)

SPEC = importlib.util.spec_from_file_location(
    "e3b_classical_frozen_path", SCRIPTS / "e3b_classical_frozen_path.py"
)
assert SPEC is not None and SPEC.loader is not None
classical = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = classical
SPEC.loader.exec_module(classical)

MODEL = "dehydroxylate-lattice"


def _data() -> str:
    return (
        "fixture\n\n"
        "2 atoms\n\n"
        "1 atom types\n\n"
        "0 10 xlo xhi\n0 10 ylo yhi\n0 10 zlo zhi\n\n"
        "Masses\n\n1 1.0 # H\n\n"
        "Atoms # full\n\n"
        "11 1 1 0.0 1.0 2.0 3.0 # H\n"
        "22 1 1 0.0 4.0 5.0 6.0 # H\n\n"
    )


def _input() -> str:
    return (
        "units real\n"
        "atom_style full\n"
        "boundary p p p\n"
        "read_data initial.seed.data\n"
        "pair_style zero 10.0\n"
        "pair_coeff * *\n"
        "thermo_style custom step pe fnorm fmax press\n"
        "thermo 100\n"
        "min_style fire\n"
        "minimize 0 0.01 100 1000\n"
        "write_data initial.relaxed.data\n"
        "write_dump all custom initial.relaxed.dump id type x y z\n"
    )


def _prepared(root: Path) -> Path:
    prepared = root / "prepared"
    model = prepared / "models" / MODEL
    images = model / "images"
    matched = model / "matched-classical"
    images.mkdir(parents=True)
    matched.mkdir()
    atom_map = {
        "atoms": [
            {"cp2k_index": 1, "source_id": 11, "source_type": "H"},
            {"cp2k_index": 2, "source_id": 22, "source_type": "H"},
        ]
    }
    e3b.write_json(model / "atom-map.json", atom_map)
    for index in range(classical.IMAGE_COUNT):
        (images / f"replica-{index:02d}.xyz").write_text(
            f"2\nimage {index}\nH {index}.0 2.0 3.0\nH 4.0 5.0 6.0\n"
        )
    (matched / "initial.seed.data").write_text(_data())
    (matched / "in.min.initial").write_text(_input())
    e3b.write_json(
        prepared / "preparation.json",
        {
            "schema": "e3b-cp2k-preparation-v1",
            "models": {MODEL: {"atom_identity_sha256": "fixture-identity"}},
        },
    )
    e3b.write_manifest(prepared)
    return prepared


def _fake_lammps(path: Path) -> Path:
    path.write_text(
        "#!/usr/bin/env python3\n"
        "import os, pathlib, re, sys\n"
        "name = sys.argv[sys.argv.index('-in') + 1]\n"
        "index = int(re.search(r'(\\d\\d)', name).group(1))\n"
        "if os.environ.get('FAIL_IMAGE') == str(index):\n"
        "    print('ERROR fixture failure')\n"
        "    raise SystemExit(2)\n"
        "print('Step PotEng Fnorm Fmax Press')\n"
        "print(f'0 {-10.0 + index} 1.0 0.5 0.0')\n"
        "print('Loop time of 0.1 on 1 procs for 0 steps')\n"
    )
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def _args(tmp_path: Path, prepared: Path, lammps: Path) -> argparse.Namespace:
    runtime = tmp_path / "runtime"
    return argparse.Namespace(
        prepared_root=prepared,
        runtime_root=runtime,
        receipt=runtime / "receipt.json",
        model=MODEL,
        lammps=str(lammps),
        per_image_timeout=5.0,
    )


def test_render_data_binds_xyz_order_to_source_ids() -> None:
    atom_map = [
        {"cp2k_index": 1, "source_id": 11, "source_type": "H"},
        {"cp2k_index": 2, "source_id": 22, "source_type": "H"},
    ]
    coordinates = {11: (7.0, 8.0, 9.0), 22: (1.0, 2.0, 3.0)}
    rendered = classical.render_lammps_data(_data(), coordinates)
    assert "11 1 1 0.0 7.000000000000 8.000000000000 9.000000000000 # H" in rendered
    assert "22 1 1 0.0 1.000000000000 2.000000000000 3.000000000000 # H" in rendered
    assert len(atom_map) == len(coordinates)


def test_render_static_input_removes_relaxation() -> None:
    rendered = classical.render_static_input(_input(), "image-03.data")
    assert "read_data       image-03.data" in rendered
    assert "thermo          1" in rendered
    assert "run             0" in rendered
    assert "minimize" not in rendered
    assert "write_data" not in rendered
    assert "write_dump" not in rendered


def test_profile_receipt_emits_like_for_like_classical_rise(tmp_path: Path) -> None:
    prepared = _prepared(tmp_path)
    manifest_before = e3b.sha256(prepared / "manifest.json")
    receipt = classical.run_profile(
        _args(tmp_path, prepared, _fake_lammps(tmp_path / "lmp"))
    )
    assert receipt["status"] == "complete"
    assert receipt["profile_rise_kcal_mol"] == 7.0
    assert len(receipt["images"]) == classical.IMAGE_COUNT
    assert (
        receipt["method"]["basis"] == "matched-cell-classical-frozen-path-single-points"
    )
    assert (
        "a production activation barrier"
        in receipt["claim_boundary"]["does_not_support"]
    )
    assert e3b.sha256(prepared / "manifest.json") == manifest_before
    assert (tmp_path / "runtime" / "receipt.json").is_file()


def test_failed_image_emits_no_numeric_profile(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prepared = _prepared(tmp_path)
    lammps = _fake_lammps(tmp_path / "lmp")
    monkeypatch.setenv("FAIL_IMAGE", "3")
    receipt = classical.run_profile(_args(tmp_path, prepared, lammps))
    assert receipt["status"] == "incomplete-image"
    assert len(receipt["images"]) == 4
    assert receipt["profile_energies_kcal_mol"] is None
    assert receipt["profile_rise_kcal_mol"] is None
    assert receipt["images"][-1]["energy_kcal_mol"] is None
