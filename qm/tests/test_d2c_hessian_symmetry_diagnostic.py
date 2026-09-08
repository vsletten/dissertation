from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pytest

from scripts import d2c_hessian_symmetry_diagnostic as diagnostic


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
    assert not (output_root / "receipt.json").exists()
