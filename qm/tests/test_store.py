"""Gates for the SQLite provenance store."""

import pytest

from quarry.store import Store, geometry_hash

XYZ = """3
water
O 0.0 0.0 0.0
H 0.96 0.0 0.0
H -0.24 0.93 0.0"""


@pytest.fixture
def store(tmp_path):
    with Store(tmp_path / "quarry.sqlite") as s:
        yield s


class TestStructures:
    def test_roundtrip(self, store):
        sid = store.add_structure("water", "H2O", XYZ)
        row = store.get_structure(sid)
        assert row["formula"] == "H2O"
        assert row["xyz"] == XYZ
        assert row["charge"] == 0

    def test_idempotent_on_geometry(self, store):
        a = store.add_structure("water", "H2O", XYZ)
        b = store.add_structure("water-again", "H2O", XYZ)
        assert a == b

    def test_hash_canonicalizes_whitespace(self):
        assert geometry_hash("O  0.0 0.0 0.0") == geometry_hash("O 0.0  0.0\t0.0")
        assert geometry_hash("O 0.0 0.0 0.0") != geometry_hash("O 0.0 0.0 0.1")

    def test_missing_structure_raises(self, store):
        with pytest.raises(KeyError):
            store.get_structure(999)


class TestJobsAndResults:
    def test_result_roundtrip(self, store):
        sid = store.add_structure("water", "H2O", XYZ)
        jid = store.add_job(sid, "sp", "hf/sto-3g", "pyscf")
        store.add_result(jid, "energy", -74.96, "hartree")
        value, units = store.get_result(jid, "energy")
        assert value == pytest.approx(-74.96)
        assert units == "hartree"

    def test_result_upsert(self, store):
        sid = store.add_structure("water", "H2O", XYZ)
        jid = store.add_job(sid, "sp", "hf/sto-3g", "pyscf")
        store.add_result(jid, "energy", -74.0, "hartree")
        store.add_result(jid, "energy", -74.96, "hartree")
        assert store.get_result(jid, "energy")[0] == pytest.approx(-74.96)

    def test_status_transitions(self, store):
        sid = store.add_structure("water", "H2O", XYZ)
        jid = store.add_job(sid, "opt", "r2scan-3c", "pyscf")
        store.set_job_status(jid, "done")
        with pytest.raises(ValueError):
            store.set_job_status(jid, "finished")

    def test_provenance_stamp(self, store):
        sid = store.add_structure("water", "H2O", XYZ)
        jid = store.add_job(sid, "sp", "wb97m-v/def2-tzvpd", "pyscf")
        p = store.provenance(jid)
        assert p.method == "wb97m-v/def2-tzvpd"
        assert p.engine == "pyscf"
        assert p.geometry_hash == geometry_hash(XYZ)
        assert p.date  # ISO timestamp present
