from __future__ import annotations

import argparse
import importlib.util
import os
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
    "e3b_frozen_path", SCRIPTS / "e3b_frozen_path.py"
)
assert SPEC is not None and SPEC.loader is not None
frozen = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = frozen
SPEC.loader.exec_module(frozen)

MODEL = "dehydroxylate-lattice"


def _prepared(root: Path) -> Path:
    prepared = root / "prepared"
    model = prepared / "models" / MODEL
    images = model / "images"
    data = prepared / "cp2k-data"
    images.mkdir(parents=True)
    data.mkdir()
    (data / "BASIS_SET").write_text("fixture\n")
    (model / "smoke.inp").write_text(
        "&GLOBAL\n"
        f"  PROJECT e3b-{MODEL}-smoke\n"
        "  RUN_TYPE ENERGY_FORCE\n"
        "&END GLOBAL\n"
        "&DFT\n"
        "  BASIS_SET_FILE_NAME ../../cp2k-data/BASIS_SET\n"
        "  &SCF\n"
        "    SCF_GUESS ATOMIC\n"
        "  &END SCF\n"
        "&END DFT\n"
        "&TOPOLOGY\n"
        "  COORD_FILE_NAME initial.xyz\n"
        "&END TOPOLOGY\n"
    )
    for index in range(frozen.IMAGE_COUNT):
        (images / f"replica-{index:02d}.xyz").write_text(
            f"1\nimage {index}\nH {index}.0 0 0\n"
        )
    e3b.write_json(
        prepared / "preparation.json",
        {
            "schema": "e3b-cp2k-preparation-v1",
            "models": {MODEL: {"atom_identity_sha256": "fixture-identity"}},
        },
    )
    e3b.write_manifest(prepared)
    return prepared


def _fake_docker(path: Path) -> Path:
    path.write_text(
        "#!/usr/bin/env python3\n"
        "import os, pathlib, re, sys\n"
        "args = sys.argv[1:]\n"
        "if args[0] == 'inspect':\n"
        "    raise SystemExit(1)\n"
        "if args[0] == 'rm':\n"
        "    raise SystemExit(0)\n"
        "cid = pathlib.Path(args[args.index('--cidfile') + 1])\n"
        "cid.write_text('fixture-container')\n"
        "name = args[args.index('-i') + 1]\n"
        "index = int(re.search(r'(\\d\\d)', name).group(1))\n"
        "fail = os.environ.get('FAIL_IMAGE')\n"
        "if fail is not None and index == int(fail):\n"
        "    print('SCF run NOT converged')\n"
        "    raise SystemExit(2)\n"
        "print('SCF run converged in 3 steps')\n"
        "print(f' ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: {-100 + index * 0.001:.12f}')\n"
        "print(' PROGRAM ENDED AT fixture')\n"
    )
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def _args(tmp_path: Path, prepared: Path, docker: Path) -> argparse.Namespace:
    runtime = tmp_path / "runtime"
    return argparse.Namespace(
        prepared_root=prepared,
        runtime_root=runtime,
        receipt=runtime / "receipt.json",
        model=MODEL,
        image="fixture/cp2k@sha256:" + "a" * 64,
        cp2k="/fixture/cp2k.psmp",
        docker=str(docker),
        mpi_ranks=4,
        omp_threads=4,
        memory_gib=48,
        timeout_seconds=120.0,
        per_image_timeout=5.0,
        systemd_unit="fixture.service",
    )


def _controller(_unit: str) -> dict[str, str]:
    return {
        "Id": "fixture.service",
        "ActiveState": "active",
        "SubState": "running",
        "RuntimeMaxUSec": "4h",
        "MemoryMax": "56G",
        "Nice": "10",
        "CPUQuotaPerSecUSec": "16s",
        "ControlGroup": "/fixture",
        "InvocationID": "fixture-invocation",
    }


def test_render_input_binds_image_and_reuses_restart_after_first() -> None:
    template = "  PROJECT old\n  SCF_GUESS ATOMIC\n  COORD_FILE_NAME initial.xyz\n"
    first = frozen.render_input(template, MODEL, 0)
    later = frozen.render_input(template, MODEL, 3)
    assert f"PROJECT e3b-{MODEL}-frozen-profile" in first
    assert "COORD_FILE_NAME images/replica-00.xyz" in first
    assert "SCF_GUESS ATOMIC" in first
    assert "COORD_FILE_NAME images/replica-03.xyz" in later
    assert "SCF_GUESS RESTART" in later


@pytest.mark.parametrize(
    "template, message",
    [
        ("SCF_GUESS ATOMIC\nCOORD_FILE_NAME x\n", "PROJECT"),
        ("PROJECT x\nSCF_GUESS ATOMIC\n", "COORD_FILE_NAME"),
        ("PROJECT x\nCOORD_FILE_NAME x\n", "SCF_GUESS"),
    ],
)
def test_render_input_fails_closed_on_malformed_template(
    template: str, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        frozen.render_input(template, MODEL, 0)


def test_profile_receipt_emits_only_frozen_path_rise(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prepared = _prepared(tmp_path)
    manifest_before = e3b.sha256(prepared / "manifest.json")
    docker = _fake_docker(tmp_path / "docker")
    monkeypatch.setattr(frozen, "systemd_readback", _controller)
    receipt = frozen.run_profile(_args(tmp_path, prepared, docker))
    assert receipt["status"] == "complete"
    assert len(receipt["images"]) == 8
    assert receipt["profile_rise_kcal_mol"] == pytest.approx(
        0.007 * e3b.HARTREE_TO_KCAL_MOL
    )
    assert "DFT minimum-energy path" in receipt["claim_boundary"]["does_not_support"]
    assert "barrier_kcal_mol" not in receipt
    assert all(item["cleanup"]["container_absent"] for item in receipt["images"])
    assert e3b.sha256(prepared / "manifest.json") == manifest_before
    assert (tmp_path / "runtime" / "receipt.json").is_file()


def test_failed_image_emits_no_profile_number(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prepared = _prepared(tmp_path)
    docker = _fake_docker(tmp_path / "docker")
    monkeypatch.setattr(frozen, "systemd_readback", _controller)
    monkeypatch.setenv("FAIL_IMAGE", "3")
    receipt = frozen.run_profile(_args(tmp_path, prepared, docker))
    assert receipt["status"] == "incomplete-image"
    assert len(receipt["images"]) == 4
    assert receipt["profile_energies_hartree"] is None
    assert receipt["profile_rise_kcal_mol"] is None
    assert receipt["failure"] == "image 03 did not converge"
    assert not any(path.suffix == ".cid" for path in (tmp_path / "runtime").iterdir())
    assert os.path.isfile(tmp_path / "runtime" / "receipt.json")
