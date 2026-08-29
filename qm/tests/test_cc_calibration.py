"""Cheap contracts for the open-source A2 coupled-cluster layer."""

from __future__ import annotations

import argparse
import json
import math
from types import SimpleNamespace

import pytest

from scripts import cc_calibration as cc


class ArrayPool:
    def new(self, _name, shape, dtype):
        return cc.np.empty(shape, dtype=dtype)


def receipt(engine: str, basis: str, scf: float, correlation: float) -> dict:
    return {
        "engine": engine,
        "basis": basis,
        "identity": {"engine": engine, "basis": basis},
        "input_sha256": "abc",
        "charge": 0,
        "spin_2s": 0,
        "scf_hartree": scf,
        "correlation_hartree": correlation,
        "total_hartree": scf + correlation,
    }


def test_frozen_core_orbitals_for_si_neutral_cluster():
    # H8O8Si2: 8 oxygen 1s orbitals + 2 silicon Ne cores (5 each).
    assert cc.frozen_core_orbitals(["H"] * 8 + ["O"] * 8 + ["Si"] * 2) == 18


def test_byteqc_frozen_core_df_transform_is_bounded_and_scoped():
    with_df = SimpleNamespace(blockdim=240, get_naoaux=lambda: 1337)
    coupled_cluster = SimpleNamespace(with_df=with_df)
    original = object()
    dfccsd = SimpleNamespace(nr_e2=original)

    with cc.byteqc_frozen_core_transform(
        coupled_cluster, dfccsd_module=dfccsd
    ) as blockdim:
        assert blockdim == 240
        assert dfccsd.nr_e2 is cc.byteqc_frozen_core_nr_e2
    assert with_df.blockdim == 240
    assert dfccsd.nr_e2 is original

    with_df.blockdim = 0
    with (
        pytest.raises(ValueError, match="block size must be positive"),
        cc.byteqc_frozen_core_transform(coupled_cluster, dfccsd_module=dfccsd),
    ):
        pass
    assert dfccsd.nr_e2 is original


def test_byteqc_frozen_core_nr_e2_is_numerical_across_full_and_tail_blocks():
    rng = cc.np.random.default_rng(220)
    nao = 5
    nmo = 3
    coefficients = rng.normal(size=(nmo, nao))
    blocks = [rng.normal(size=(rows, nao, nao)) for rows in (4, 4, 1)]
    blocks = [(block + block.transpose(0, 2, 1)) / 2 for block in blocks]
    lower = cc.np.tril_indices(nao)
    packed = [block[:, lower[0], lower[1]].copy() for block in blocks]
    calls = {"unpack_out": [], "gemm_buf": []}

    def unpack_tril(eri, out=None):
        calls["unpack_out"].append(out)
        matrices = cc.np.zeros((eri.shape[0], nao, nao))
        matrices[:, lower[0], lower[1]] = eri
        matrices[:, lower[1], lower[0]] = eri
        return matrices

    def contraction(_inda, a, _indb, b, _indc):
        return cc.np.einsum("ij,mjk->mik", a, b)

    def gemm(_transa, _transb, a, b, buf=None):
        calls["gemm_buf"].append(buf)
        value = a @ b.T
        if buf is None:
            return value
        view = buf.reshape(-1)[: value.size].reshape(value.shape)
        view[:] = value
        return view

    output = None
    results = []
    for ao_block, packed_block in zip(blocks, packed, strict=True):
        output = cc.byteqc_frozen_core_nr_e2(
            packed_block,
            coefficients,
            (0, nmo, 0, nmo),
            aosym="s2",
            mosym="s1",
            out=output,
            array_module=cc.np,
            lib_module=SimpleNamespace(unpack_tril=unpack_tril),
            contraction_module=SimpleNamespace(
                contraction=contraction,
                gemm=gemm,
            ),
        )
        expected = cc.np.einsum("ia,mab,jb->mij", coefficients, ao_block, coefficients)
        assert output == pytest.approx(expected, abs=1e-12)
        results.append(output)

    assert calls["unpack_out"] == [None, None, None]
    assert calls["gemm_buf"][0] is None
    assert calls["gemm_buf"][1] is results[0]
    assert cc.np.shares_memory(calls["gemm_buf"][2], results[2])


def test_ccsd_checkpoint_round_trip_is_atomic_and_fingerprint_bound(tmp_path):
    path = tmp_path / "restart.h5"
    identity = {
        "input_sha256": "reactant-geometry",
        "basis": "cc-pVQZ",
        "engine_commit": "byteqc-head",
        "driver_source_sha256": "driver-head",
        "mo_coeff_sha256": "orbitals",
    }
    t1 = cc.np.arange(6, dtype=float).reshape(2, 3) / 10
    t2 = cc.np.arange(36, dtype=float).reshape(2, 2, 3, 3) / 100
    metadata = cc.write_ccsd_checkpoint(
        path,
        identity=identity,
        completed_cycle=4,
        ccsd_correlation_hartree=-2.5,
        t1=t1,
        t2=t2,
    )
    cluster = SimpleNamespace(nocc=2, nmo=5, pool=ArrayPool())
    loaded = cc.load_ccsd_checkpoint(
        path,
        identity=identity,
        coupled_cluster=cluster,
    )
    assert loaded is not None
    loaded_metadata, loaded_t1, loaded_t2 = loaded
    assert loaded_metadata == metadata
    assert loaded_t1 == pytest.approx(t1)
    assert loaded_t2 == pytest.approx(t2)

    for key, value in (
        ("input_sha256", "transition-state-geometry"),
        ("basis", "cc-pVTZ"),
        ("engine_commit", "other-byteqc-head"),
        ("driver_source_sha256", "other-driver-head"),
        ("mo_coeff_sha256", "other-orbitals"),
    ):
        drifted = {**identity, key: value}
        with pytest.raises(ValueError, match="identity drift"):
            cc.load_ccsd_checkpoint(
                path,
                identity=drifted,
                coupled_cluster=cluster,
            )


def test_ccsd_checkpoint_refuses_corruption_nonfinite_and_shape_drift(tmp_path):
    import h5py

    identity = {"input_sha256": "geometry", "basis": "cc-pVQZ"}
    cluster = SimpleNamespace(nocc=2, nmo=4, pool=ArrayPool())
    path = tmp_path / "restart.h5"
    cc.write_ccsd_checkpoint(
        path,
        identity=identity,
        completed_cycle=2,
        ccsd_correlation_hartree=-1.0,
        t1=cc.np.zeros((2, 2)),
        t2=cc.np.zeros((2, 2, 2, 2)),
    )
    with h5py.File(path, "r+") as store:
        dataset = store["t2"]
        assert isinstance(dataset, h5py.Dataset)
        dataset[0, 0, 0, 0] = 1.0
    with pytest.raises(ValueError, match="checksum mismatch"):
        cc.load_ccsd_checkpoint(path, identity=identity, coupled_cluster=cluster)

    with pytest.raises(ValueError, match="non-finite"):
        cc.write_ccsd_checkpoint(
            tmp_path / "nonfinite.h5",
            identity=identity,
            completed_cycle=2,
            ccsd_correlation_hartree=-1.0,
            t1=cc.np.zeros((2, 2)),
            t2=cc.np.full((2, 2, 2, 2), cc.np.nan),
        )

    wrong_shape = tmp_path / "wrong-shape.h5"
    cc.write_ccsd_checkpoint(
        wrong_shape,
        identity=identity,
        completed_cycle=2,
        ccsd_correlation_hartree=-1.0,
        t1=cc.np.zeros((1, 2)),
        t2=cc.np.zeros((1, 1, 2, 2)),
    )
    with pytest.raises(ValueError, match="shape drift"):
        cc.load_ccsd_checkpoint(
            wrong_shape,
            identity=identity,
            coupled_cluster=cluster,
        )


def test_ccsd_checkpoint_binds_cycle_metadata_and_preflights_disk(
    tmp_path, monkeypatch
):
    import h5py

    identity = {"input_sha256": "geometry", "basis": "cc-pVQZ"}
    cluster = SimpleNamespace(nocc=1, nmo=2, pool=ArrayPool())
    path = tmp_path / "restart.h5"
    arrays = {
        "t1": cc.np.zeros((1, 1)),
        "t2": cc.np.zeros((1, 1, 1, 1)),
    }
    cc.write_ccsd_checkpoint(
        path,
        identity=identity,
        completed_cycle=2,
        ccsd_correlation_hartree=-1.0,
        **arrays,
    )
    with h5py.File(path, "r+") as store:
        manifest_raw = store.attrs["manifest_json"]
        assert isinstance(manifest_raw, str)
        manifest = json.loads(manifest_raw)
        manifest["completed_cycle"] = 99
        store.attrs["manifest_json"] = json.dumps(manifest, sort_keys=True)
    with pytest.raises(ValueError, match="manifest checksum mismatch"):
        cc.load_ccsd_checkpoint(path, identity=identity, coupled_cluster=cluster)

    monkeypatch.setattr(
        cc.shutil,
        "disk_usage",
        lambda _path: SimpleNamespace(free=0),
    )
    with pytest.raises(RuntimeError, match="insufficient disk"):
        cc.write_ccsd_checkpoint(
            tmp_path / "no-space.h5",
            identity=identity,
            completed_cycle=1,
            ccsd_correlation_hartree=-1.0,
            **arrays,
        )


def test_converged_cycle_100_marker_is_restartable_and_fingerprint_bound(tmp_path):
    checkpoint = tmp_path / "restart.h5"
    identity = {"input_sha256": "geometry", "basis": "cc-pVQZ"}
    metadata = cc.write_ccsd_checkpoint(
        checkpoint,
        identity=identity,
        completed_cycle=100,
        ccsd_correlation_hartree=-2.5,
        t1=cc.np.zeros((1, 1)),
        t2=cc.np.zeros((1, 1, 1, 1)),
    )
    marker = cc.write_ccsd_converged_marker(checkpoint, metadata)
    assert marker["completed_cycle"] == 100
    assert cc.load_ccsd_converged_marker(checkpoint, metadata) == marker

    drifted = {**metadata, "ccsd_correlation_hartree": -2.4}
    with pytest.raises(ValueError, match="identity drift"):
        cc.load_ccsd_converged_marker(checkpoint, drifted)


def test_ccsd_checkpoint_failed_replacement_preserves_last_good_state(
    tmp_path, monkeypatch
):
    path = tmp_path / "restart.h5"
    identity = {"input_sha256": "geometry"}
    arrays = {
        "t1": cc.np.zeros((2, 2)),
        "t2": cc.np.zeros((2, 2, 2, 2)),
    }
    cc.write_ccsd_checkpoint(
        path,
        identity=identity,
        completed_cycle=1,
        ccsd_correlation_hartree=-1.0,
        **arrays,
    )
    original = path.read_bytes()

    def fail_replace(_source, _destination):
        raise OSError("simulated interrupted atomic replace")

    monkeypatch.setattr(cc.os, "replace", fail_replace)
    with pytest.raises(OSError, match="interrupted atomic replace"):
        cc.write_ccsd_checkpoint(
            path,
            identity=identity,
            completed_cycle=2,
            ccsd_correlation_hartree=-1.1,
            **arrays,
        )
    assert path.read_bytes() == original


def test_byteqc_energy_checkpoint_persists_each_completed_cycle_before_next_update(
    tmp_path, monkeypatch
):
    calls = []
    monkeypatch.setattr(
        cc,
        "write_ccsd_checkpoint",
        lambda path, **kwargs: calls.append((path, kwargs)) or kwargs,
    )
    checkpoint = cc.ByteQCEnergyCheckpoint(
        tmp_path / "restart.h5", {"basis": "cc-pVQZ"}, base_cycle=4
    )
    energies = iter((-2.0, -2.1, -2.2))
    wrapped = checkpoint.wrap(lambda *_args: (next(energies), None))
    t1 = cc.np.zeros((1, 1))
    t2 = cc.np.zeros((1, 1, 1, 1))

    wrapped(t1, t2, None)
    assert calls == []
    wrapped(t1, t2, None)
    assert calls[0][1]["completed_cycle"] == 5
    assert calls[0][1]["ccsd_correlation_hartree"] == -2.1
    # A kill during the following update cannot roll back cycle 5: its atomic
    # state was written by the preceding post-DIIS energy evaluation.
    assert checkpoint.last_checkpoint is not None
    wrapped(t1, t2, None)
    assert calls[1][1]["completed_cycle"] == 6
    assert calls[1][1]["ccsd_correlation_hartree"] == -2.2


def test_restart_checkpoint_path_rejects_job_artifact_aliases(tmp_path):
    source = tmp_path / "input.xyz"
    output = tmp_path / "receipt.json"
    engine_log = tmp_path / "engine.log"
    for collision in (source, output, engine_log):
        args = argparse.Namespace(
            restart_checkpoint=collision,
            output=output,
            log=engine_log,
        )
        with pytest.raises(ValueError, match="must differ"):
            cc.restart_checkpoint_path(args, source)
    args = argparse.Namespace(restart_checkpoint=None, output=output, log=engine_log)
    assert (
        cc.restart_checkpoint_path(args, source) == tmp_path / "receipt.ccsd-restart.h5"
    )


def test_two_point_extrapolations_recover_synthetic_limits():
    hf_limit = -100.0
    hf_amplitude = 0.2
    tz_hf = hf_limit + hf_amplitude * math.exp(-cc.HF_CBS_ALPHA * 3)
    qz_hf = hf_limit + hf_amplitude * math.exp(-cc.HF_CBS_ALPHA * 4)
    assert cc.extrapolate_hf(tz_hf, qz_hf) == pytest.approx(hf_limit)

    corr_limit = -1.0
    corr_amplitude = 0.5
    tz_corr = corr_limit + corr_amplitude / 3**3
    qz_corr = corr_limit + corr_amplitude / 4**3
    assert cc.extrapolate_correlation(tz_corr, qz_corr) == pytest.approx(corr_limit)


def test_receipt_grid_requires_complete_role_basis_matrix():
    values = [
        "reactant=cc-pVTZ=r-tz.json",
        "reactant=cc-pVQZ=r-qz.json",
        "ts=cc-pVTZ=t-tz.json",
        "ts=cc-pVQZ=t-qz.json",
    ]
    grid = cc.receipt_grid(values)
    assert set(grid) == {"reactant", "ts"}
    assert set(grid["reactant"]) == {"tz", "qz"}
    with pytest.raises(ValueError, match="cover reactant/ts"):
        cc.receipt_grid(values[:-1])


def test_existing_receipt_fails_closed_on_identity_drift(tmp_path):
    path = tmp_path / "receipt.json"
    payload = receipt("byteqc-canonical-ccsd(t)", "cc-pVTZ", -10.0, -1.0)
    path.write_text(json.dumps(payload))
    assert cc.existing_receipt(
        path,
        engine="byteqc-canonical-ccsd(t)",
        basis="cc-pVTZ",
        input_sha256="abc",
        identity=payload["identity"],
    )
    with pytest.raises(ValueError, match="identity drift"):
        cc.existing_receipt(
            path,
            engine="byteqc-canonical-ccsd(t)",
            basis="cc-pVQZ",
            input_sha256="abc",
            identity=payload["identity"],
        )
    with pytest.raises(ValueError, match="expected basis cc-pVQZ"):
        cc.load_receipt(path, "byteqc-canonical-ccsd(t)", "cc-pVQZ")


def test_calibration_engines_refuse_open_shell_before_import():
    with pytest.raises(ValueError, match="closed-shell spin=0 only"):
        cc.byteqc_job(argparse.Namespace(spin=1))
    with pytest.raises(ValueError, match="supports spin=0 only"):
        cc.psi4_job(argparse.Namespace(spin=1))


def test_calibration_resource_contract_requires_exact_byteqc_gpu_lease():
    byteqc = argparse.Namespace(
        command="byteqc", gpu=False, gpu_mem_gb=16.0, gpu_memory_gb=16
    )
    with pytest.raises(ValueError, match="requires --gpu"):
        cc.validate_resource_contract(byteqc)

    byteqc.gpu = True
    byteqc.gpu_mem_gb = 12.0
    with pytest.raises(ValueError, match="must be at least"):
        cc.validate_resource_contract(byteqc)

    byteqc.gpu_mem_gb = 18.0
    cc.validate_resource_contract(byteqc)

    with pytest.raises(ValueError, match="only for ByteQC"):
        cc.validate_resource_contract(argparse.Namespace(command="psi4", gpu=True))


def test_summarize_reports_barrier_delta_and_gate(tmp_path):
    paths = {engine: {} for engine in ("canonical", "dlpno")}
    engine_names = {
        "canonical": "byteqc-canonical-ccsd(t)",
        "dlpno": "psi4-dlpno-ccsd(t)",
    }
    for engine, engine_name in engine_names.items():
        for role, offset in (("reactant", 0.0), ("ts", 0.04)):
            for basis, label in (("cc-pVTZ", "tz"), ("cc-pVQZ", "qz")):
                # DLPNO TS is 0.0001 Eh above canonical (~0.26 kJ/mol delta).
                local_delta = 0.0001 if engine == "dlpno" and role == "ts" else 0.0
                payload = receipt(
                    engine_name,
                    basis,
                    -100.0 + offset + local_delta,
                    -1.0,
                )
                path = tmp_path / f"{engine}-{role}-{label}.json"
                path.write_text(json.dumps(payload))
                paths[engine].setdefault(role, {})[label] = path
    output = tmp_path / "summary.json"
    result = cc.summarize(
        argparse.Namespace(
            canonical=paths["canonical"],
            dlpno=paths["dlpno"],
            output=output,
        )
    )
    assert result["canonical_barrier_kj"] == pytest.approx(0.04 * cc.HARTREE_TO_KJ)
    assert result["dlpno_minus_canonical_kj"] == pytest.approx(
        0.0001 * cc.HARTREE_TO_KJ
    )
    assert result["gate_pass"] is True
    assert json.loads(output.read_text()) == result

    drift_path = paths["dlpno"]["ts"]["qz"]
    drift = json.loads(drift_path.read_text())
    drift["input_sha256"] = "different-geometry"
    drift_path.write_text(json.dumps(drift))
    with pytest.raises(ValueError, match="do not share one exact geometry"):
        cc.summarize(
            argparse.Namespace(
                canonical=paths["canonical"],
                dlpno=paths["dlpno"],
                output=output,
            )
        )
