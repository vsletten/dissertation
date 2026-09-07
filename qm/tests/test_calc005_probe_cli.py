"""Executable reproduction gate for the CALC-005 CPU proof."""

import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
QM = REPO / "qm"
DECK = REPO / "petra" / "examples" / "kaolinite.toml"
SCRIPT = QM / "scripts" / "calc005_si_attachment.py"


def test_probe_cli_emits_complete_deterministic_json():
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "probe",
            "--deck",
            str(DECK),
            "--repeat",
            "2",
        ],
        cwd=QM,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(result.stdout[result.stdout.index("{") :])

    assert payload["schema"] == "calc005-live-pair-probe-v1"
    assert len(payload["calc005_sha256"]) == 64
    assert len(payload["crystal_sha256"]) == 64
    assert len(payload["driver_sha256"]) == 64
    assert payload["two_pass_deterministic"] is True
    assert len(payload["pairs"]) == 4
    assert all(row["atom_conservation"] for row in payload["pairs"])
