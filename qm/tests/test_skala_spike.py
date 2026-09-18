import hashlib
import json

import pytest

from quarry.store import Store, geometry_hash
from scripts import skala_spike as spike


def _xyz(name: str, offset: float) -> str:
    return (
        f"2\n{name}\n"
        f"H {offset:.6f} 0.000000 0.000000\n"
        f"H {offset + 0.740000:.6f} 0.000000 0.000000\n"
    )


def test_load_source_structures_requires_exact_hashes_and_immutable_read(
    monkeypatch, tmp_path
):
    source = tmp_path / "store.sqlite"
    expected = {}
    with Store(source) as store:
        for source_id, role in enumerate(
            (
                "reactant",
                "intermediate",
                "addition-transition-state",
                "released-product",
            ),
            start=1,
        ):
            xyz = _xyz(role, source_id * 0.1)
            name = f"fixture-{role}"
            inserted = store.add_structure(name, "H8O8Si2", xyz)
            assert inserted == source_id
            expected[source_id] = (role, name, geometry_hash(xyz))
    monkeypatch.setattr(spike, "SOURCE_STORE_SHA256", spike.sha256_path(source))
    monkeypatch.setattr(spike, "EXPECTED_STRUCTURES", expected)

    connection = spike.load_source_structures(source)

    assert [structure.source_id for structure in connection] == [1, 2, 3, 4]
    assert [structure.role for structure in connection] == [
        "reactant",
        "intermediate",
        "addition-transition-state",
        "released-product",
    ]
    assert hashlib.sha256(source.read_bytes()).hexdigest() == spike.SOURCE_STORE_SHA256


def test_load_source_structures_refuses_source_hash_drift(tmp_path):
    source = tmp_path / "store.sqlite"
    source.write_bytes(b"not the evidence store")
    with pytest.raises(RuntimeError, match="A2a store SHA-256 drift"):
        spike.load_source_structures(source)


def test_derive_rows_recomputes_barrier_and_predeclared_verdict():
    jobs = []
    energies = {
        "reactant": -10.0,
        "intermediate": -9.99,
        "addition-transition-state": -10.0
        + spike.FOCAL_BARRIER_KJ_MOL / spike.HARTREE_TO_KJ,
        "released-product": -10.02,
    }
    for method_key in spike.METHODS:
        for role, energy in energies.items():
            jobs.append(
                {
                    "method_key": method_key,
                    "role": role,
                    "status": "done",
                    "energy_hartree": energy,
                }
            )

    rows = spike.derive_rows(jobs)

    assert set(rows) == set(spike.METHODS)
    for row in rows.values():
        assert row["status"] == "complete"
        assert row["barrier_kj_mol"] == pytest.approx(
            spike.FOCAL_BARRIER_KJ_MOL, abs=1e-10
        )
        assert row["comparison_to_focal"]["verdict"] == (
            "adopt as a survey-tier single-point functional candidate"
        )


def test_derive_rows_keeps_failures_honest():
    jobs = [
        {
            "method_key": "skala-1.1-def2-tzvp",
            "role": "reactant",
            "status": "failed",
        }
    ]

    rows = spike.derive_rows(jobs)

    row = rows["skala-1.1-def2-tzvp"]
    assert row["status"] == "incomplete"
    assert row["failed_roles"] == ["reactant"]
    assert "barrier_kj_mol" not in row


def test_manifest_payload_is_strict_json(tmp_path):
    path = tmp_path / "receipt.json"
    spike.atomic_json(path, {"finite": 1.0, "nested": {"ok": True}})
    assert json.loads(path.read_text()) == {"finite": 1.0, "nested": {"ok": True}}
    with pytest.raises(ValueError):
        spike.atomic_json(path, {"bad": float("nan")})


def test_skala_cpu_tiny_molecule_runs_without_implicit_dispersion():
    pytest.importorskip("skala")
    pytest.importorskip("torch")
    structure = spike.SourceStructure(
        source_id=1,
        role="reactant",
        name="h2",
        formula="H2",
        charge=0,
        spin=0,
        xyz="2\nh2\nH 0 0 0\nH 0 0 0.74\n",
        geometry_hash="fixture",
    )

    energy, converged, cycles, settings = spike.run_skala_scf(
        structure, basis="def2-svp", use_gpu=False
    )

    assert converged is True
    assert energy < 0.0
    assert cycles is None or cycles > 0
    assert settings["density_fit"] is False
    assert settings["dftd3"] is False
    assert settings["grid"] == "SkalaKS default"
