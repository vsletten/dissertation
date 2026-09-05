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


_WRITE_CCSD_CHECKPOINT = cc.write_ccsd_checkpoint


def checkpoint_scf(nmo: int) -> dict:
    return {
        "mo_coeff": cc.np.eye(nmo, dtype=cc.np.float64),
        "mo_energy": cc.np.linspace(-1.0, 1.0, nmo, dtype=cc.np.float64),
        "mo_occ": cc.np.array([2.0] * (nmo // 2) + [0.0] * (nmo - nmo // 2)),
        "scf_energy": -10.5,
    }


def checkpoint_identity(*, nmo: int, frozen: int = 0) -> dict:
    return {
        "input_sha256": "geometry",
        "electronic_structure": {
            "frozen_core_orbitals": frozen,
            "electron_count": 2 * (nmo // 2),
        },
    }


def checkpoint_mean_field(nmo: int):
    orbital_energies = checkpoint_scf(nmo)["mo_energy"]
    return SimpleNamespace(
        mol=SimpleNamespace(nao_nr=lambda: nmo, nelectron=2 * (nmo // 2)),
        get_ovlp=lambda: cc.np.eye(nmo),
        make_rdm1=lambda mo_coeff, mo_occ: mo_coeff @ cc.np.diag(mo_occ) @ mo_coeff.T,
        get_hcore=lambda: cc.np.diag(orbital_energies),
        get_veff=lambda _mol, density: cc.np.zeros_like(density),
        energy_tot=lambda **_kwargs: -10.5,
        get_fock=lambda **kwargs: kwargs["h1e"] + kwargs["vhf"],
        mo_coeff=None,
        mo_energy=None,
        mo_occ=None,
    )


def write_checkpoint(path, **kwargs):
    nocc, nvir = kwargs["t1"].shape
    return _WRITE_CCSD_CHECKPOINT(path, **checkpoint_scf(nocc + nvir), **kwargs)


def load_checkpoint(path, *, identity, coupled_cluster):
    mean_field = checkpoint_mean_field(coupled_cluster.nmo)
    metadata = cc.restore_ccsd_checkpoint_scf(
        path, identity=identity, mean_field=mean_field
    )
    assert metadata is not None
    expected_scf = checkpoint_scf(coupled_cluster.nmo)
    assert cc.np.array_equal(mean_field.mo_coeff, expected_scf["mo_coeff"])
    assert cc.np.array_equal(mean_field.mo_energy, expected_scf["mo_energy"])
    assert cc.np.array_equal(mean_field.mo_occ, expected_scf["mo_occ"])
    assert mean_field.e_tot == expected_scf["scf_energy"]
    t1, t2 = cc.load_ccsd_checkpoint_amplitudes(
        path,
        identity=identity,
        coupled_cluster=coupled_cluster,
        expected_metadata=metadata,
    )
    return metadata, t1, t2


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
        **checkpoint_identity(nmo=5),
        "input_sha256": "reactant-geometry",
        "basis": "cc-pVQZ",
        "engine_commit": "byteqc-head",
        "driver_source_sha256": "driver-head",
        "mo_coeff_sha256": "orbitals",
    }
    t1 = cc.np.arange(6, dtype=float).reshape(2, 3) / 10
    t2 = cc.np.arange(36, dtype=float).reshape(2, 2, 3, 3) / 100
    metadata = write_checkpoint(
        path,
        identity=identity,
        completed_cycle=4,
        ccsd_correlation_hartree=-2.5,
        t1=t1,
        t2=t2,
    )
    cluster = SimpleNamespace(nocc=2, nmo=5, pool=ArrayPool())
    loaded = load_checkpoint(
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
            load_checkpoint(
                path,
                identity=drifted,
                coupled_cluster=cluster,
            )


def test_ccsd_checkpoint_refuses_corruption_nonfinite_and_shape_drift(tmp_path):
    import h5py

    identity = {**checkpoint_identity(nmo=4), "basis": "cc-pVQZ"}
    cluster = SimpleNamespace(nocc=2, nmo=4, pool=ArrayPool())
    path = tmp_path / "restart.h5"
    write_checkpoint(
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
        load_checkpoint(path, identity=identity, coupled_cluster=cluster)

    with pytest.raises(ValueError, match="non-finite"):
        write_checkpoint(
            tmp_path / "nonfinite.h5",
            identity=identity,
            completed_cycle=2,
            ccsd_correlation_hartree=-1.0,
            t1=cc.np.zeros((2, 2)),
            t2=cc.np.full((2, 2, 2, 2), cc.np.nan),
        )

    wrong_shape = tmp_path / "wrong-shape.h5"
    write_checkpoint(
        wrong_shape,
        identity=identity,
        completed_cycle=2,
        ccsd_correlation_hartree=-1.0,
        t1=cc.np.zeros((1, 2)),
        t2=cc.np.zeros((1, 1, 2, 2)),
    )
    with pytest.raises(ValueError, match="shape drift"):
        load_checkpoint(
            wrong_shape,
            identity=identity,
            coupled_cluster=cluster,
        )


def test_ccsd_checkpoint_binds_cycle_metadata_and_preflights_disk(
    tmp_path, monkeypatch
):
    import h5py

    identity = {**checkpoint_identity(nmo=2), "basis": "cc-pVQZ"}
    cluster = SimpleNamespace(nocc=1, nmo=2, pool=ArrayPool())
    path = tmp_path / "restart.h5"
    arrays = {
        "t1": cc.np.zeros((1, 1)),
        "t2": cc.np.zeros((1, 1, 1, 1)),
    }
    write_checkpoint(
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
        load_checkpoint(path, identity=identity, coupled_cluster=cluster)

    monkeypatch.setattr(
        cc.shutil,
        "disk_usage",
        lambda _path: SimpleNamespace(free=0),
    )
    with pytest.raises(RuntimeError, match="insufficient disk"):
        write_checkpoint(
            tmp_path / "no-space.h5",
            identity=identity,
            completed_cycle=1,
            ccsd_correlation_hartree=-1.0,
            **arrays,
        )


def test_converged_cycle_100_marker_is_restartable_and_fingerprint_bound(tmp_path):
    checkpoint = tmp_path / "restart.h5"
    identity = {**checkpoint_identity(nmo=2), "basis": "cc-pVQZ"}
    metadata = write_checkpoint(
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
    identity = checkpoint_identity(nmo=4)
    arrays = {
        "t1": cc.np.zeros((2, 2)),
        "t2": cc.np.zeros((2, 2, 2, 2)),
    }
    metadata = write_checkpoint(
        path,
        identity=identity,
        completed_cycle=1,
        ccsd_correlation_hartree=-1.0,
        **arrays,
    )
    cc.write_ccsd_converged_marker(path, metadata)
    original = path.read_bytes()

    def fail_replace(_source, _destination):
        raise OSError("simulated interrupted atomic replace")

    monkeypatch.setattr(cc.os, "replace", fail_replace)
    with pytest.raises(OSError, match="interrupted atomic replace"):
        write_checkpoint(
            path,
            identity=identity,
            completed_cycle=2,
            ccsd_correlation_hartree=-1.1,
            **arrays,
        )
    assert path.read_bytes() == original
    assert not cc.ccsd_converged_marker_path(path).exists()


def test_byteqc_energy_checkpoint_stops_before_second_update(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(
        cc,
        "write_ccsd_checkpoint",
        lambda path, **kwargs: calls.append((path, kwargs)) or kwargs,
    )
    checkpoint = cc.ByteQCEnergyCheckpoint(
        tmp_path / "restart.h5",
        {"basis": "cc-pVQZ"},
        checkpoint_scf(2),
        base_cycle=4,
    )
    energies = iter((-2.0, -2.1))
    wrapped = checkpoint.wrap(lambda *_args: (next(energies), None))
    updates = []
    wrapped_update = checkpoint.wrap_update_amps(
        lambda *args: updates.append(args) or args
    )
    t1 = cc.np.zeros((1, 1))
    t2 = cc.np.zeros((1, 1, 1, 1))

    wrapped(t1, t2, None)
    assert calls == []
    wrapped_update(t1, t2, None)
    wrapped(t1, t2, None)
    assert calls[0][1]["completed_cycle"] == 5
    assert calls[0][1]["ccsd_correlation_hartree"] == -2.1
    assert len(updates) == 1
    assert checkpoint.last_checkpoint is not None
    with pytest.raises(cc.CCSDLaunchCheckpointed) as stopped:
        wrapped_update(t1, t2, None)
    assert stopped.value.checkpoint is checkpoint.last_checkpoint
    assert len(calls) == 1 and len(updates) == 1
    assert set(checkpoint_scf(2)) <= set(calls[0][1])


def test_byteqc_checkpoint_limit_allows_convergence_check(tmp_path, monkeypatch):
    writes = []
    monkeypatch.setattr(
        cc,
        "write_ccsd_checkpoint",
        lambda _path, **kwargs: writes.append(kwargs) or kwargs,
    )
    cluster = SimpleNamespace()
    cluster.energy = lambda *_args: (-2.0, None)
    cluster.update_amps = lambda t1, t2, _eris: (t1, t2)
    t1 = cc.np.zeros((1, 1))
    t2 = cc.np.zeros((1, 1, 1, 1))

    with cc.byteqc_energy_checkpoint(
        cluster,
        tmp_path / "restart.h5",
        {"basis": "cc-pVQZ"},
        checkpoint_scf(2),
        base_cycle=4,
    ) as checkpoint:
        cluster.energy(t1, t2, None)  # Initial MP2 energy.
        t1, t2 = cluster.update_amps(t1, t2, None)
        cluster.energy(t1, t2, None)  # Durable cycle 5.
        converged = True  # Models ByteQC's check after post-DIIS energy.
        if not converged:
            cluster.update_amps(t1, t2, None)

    assert checkpoint.checkpoints_written == 1
    assert writes[0]["completed_cycle"] == 5


def test_failed_checkpoint_write_is_not_counted_as_progress(tmp_path, monkeypatch):
    checkpoint = cc.ByteQCEnergyCheckpoint(
        tmp_path / "restart.h5",
        {"basis": "cc-pVQZ"},
        checkpoint_scf(2),
        base_cycle=4,
    )

    def fail_write(*_args, **_kwargs):
        raise OSError("write failed")

    monkeypatch.setattr(cc, "write_ccsd_checkpoint", fail_write)
    energies = iter((-2.0, -2.1))
    wrapped = checkpoint.wrap(lambda *_args: (next(energies), None))
    t1 = cc.np.zeros((1, 1))
    t2 = cc.np.zeros((1, 1, 1, 1))

    wrapped(t1, t2, None)
    with pytest.raises(OSError, match="write failed"):
        wrapped(t1, t2, None)
    assert checkpoint.checkpoints_written == 0
    assert checkpoint.last_checkpoint is None


def test_checkpoint_context_restores_methods_on_boundary(tmp_path, monkeypatch):
    monkeypatch.setattr(cc, "write_ccsd_checkpoint", lambda _path, **kwargs: kwargs)

    def original_energy(*_args):
        return (-2.0, None)

    def original_update(t1, t2, _eris):
        return (t1, t2)

    cluster = SimpleNamespace(energy=original_energy, update_amps=original_update)
    t1 = cc.np.zeros((1, 1))
    t2 = cc.np.zeros((1, 1, 1, 1))

    with (
        pytest.raises(cc.CCSDLaunchCheckpointed),
        cc.byteqc_energy_checkpoint(
            cluster,
            tmp_path / "restart.h5",
            {"basis": "cc-pVQZ"},
            checkpoint_scf(2),
            base_cycle=4,
        ),
    ):
        cluster.energy(t1, t2, None)
        cluster.update_amps(t1, t2, None)
        cluster.energy(t1, t2, None)
        cluster.update_amps(t1, t2, None)

    assert cluster.energy is original_energy
    assert cluster.update_amps is original_update


def test_v2_rejects_legacy_without_overwrite(tmp_path):
    import h5py

    path = tmp_path / "legacy.h5"
    metadata = {"schema": cc.CCSD_RESTART_SCHEMA_V1}
    with h5py.File(path, "w") as store:
        store.attrs["manifest_json"] = json.dumps(metadata)
        store.attrs["manifest_sha256"] = cc.canonical_json_sha256(metadata)
    original = path.read_bytes()
    with pytest.raises(ValueError, match="legacy v1.*not be migrated or overwritten"):
        write_checkpoint(
            path,
            identity={"input_sha256": "geometry"},
            completed_cycle=1,
            ccsd_correlation_hartree=-1.0,
            t1=cc.np.zeros((1, 1)),
            t2=cc.np.zeros((1, 1, 1, 1)),
        )
    assert path.read_bytes() == original


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("corrupt", "checksum mismatch"),
        ("nonfinite", "non-finite"),
        ("shape", "shape drift"),
        ("dtype", "dtype drift"),
    ],
)
def test_v2_rejects_adversarial_orbital_state(tmp_path, mutation, message):
    import h5py

    path = tmp_path / f"orbital-{mutation}.h5"
    identity = checkpoint_identity(nmo=2)
    write_checkpoint(
        path,
        identity=identity,
        completed_cycle=1,
        ccsd_correlation_hartree=-1.0,
        t1=cc.np.zeros((1, 1)),
        t2=cc.np.zeros((1, 1, 1, 1)),
    )
    with h5py.File(path, "r+") as store:
        if mutation == "corrupt":
            store["mo_coeff"][0, 0] = 9.0
        elif mutation == "nonfinite":
            store["mo_energy"][0] = cc.np.nan
        elif mutation == "shape":
            del store["mo_occ"]
            store.create_dataset("mo_occ", data=cc.np.zeros(3))
        else:
            del store["mo_energy"]
            store.create_dataset("mo_energy", data=cc.np.zeros(2, dtype=cc.np.float32))
    mean_field = checkpoint_mean_field(2)
    with pytest.raises(ValueError, match=message):
        cc.restore_ccsd_checkpoint_scf(path, identity=identity, mean_field=mean_field)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("occupations", "occupations are invalid"),
        ("energies", "energies are not ordered"),
        ("metric", "not orthonormal"),
        ("energy", "total energy is inconsistent"),
        ("fock", "do not diagonalize Fock"),
        ("amplitudes", "disagree with restored RHF occupations"),
    ],
)
def test_v2_rejects_physically_inconsistent_scf_state(tmp_path, mutation, message):
    state = checkpoint_scf(2)
    t1_shape = (1, 1)
    if mutation == "occupations":
        state["mo_occ"] = cc.np.array([1.0, 1.0])
    elif mutation == "energies":
        state["mo_energy"] = state["mo_energy"][::-1]
    elif mutation == "metric":
        state["mo_coeff"] = cc.np.diag([1.0, 2.0])
    elif mutation == "energy":
        state["scf_energy"] = 123.0
    elif mutation == "fock":
        theta = 0.3
        state["mo_coeff"] = cc.np.array(
            [
                [math.cos(theta), -math.sin(theta)],
                [math.sin(theta), math.cos(theta)],
            ]
        )
    else:
        t1_shape = (2, 1)
    with pytest.raises(ValueError, match=message):
        cc._validate_restored_scf_state(
            tmp_path / "restart.h5",
            identity=checkpoint_identity(nmo=2),
            mean_field=checkpoint_mean_field(2),
            mo_coeff=state["mo_coeff"],
            mo_energy=state["mo_energy"],
            mo_occ=state["mo_occ"],
            scf_energy=state["scf_energy"],
            t1_shape=t1_shape,
        )


def test_restart_identity_is_stable_but_contract_bound():
    kwargs = {
        "input_sha256": "geometry",
        "engine_commit": "02af2548",
        "pyscf_version": "2.5.0",
        "cupy_version": "13",
        "numpy_version": "2",
        "python_version": "3.12",
        "basis": "cc-pVQZ",
        "charge": 0,
        "spin_2s": 0,
        "frozen_core_orbitals": 18,
    }
    identity = cc.byteqc_restart_identity(**kwargs)
    encoded = json.dumps(identity)
    assert "mo_coeff_sha256" not in encoded and "scf_energy" not in encoded
    assert "threads" not in encoded and "gpu_memory" not in encoded
    for key in (
        "input_sha256",
        "engine_commit",
        "basis",
        "charge",
        "frozen_core_orbitals",
    ):
        drifted = dict(kwargs)
        drifted[key] = "different" if isinstance(kwargs[key], str) else kwargs[key] + 1
        assert cc.byteqc_restart_identity(**drifted) != identity


def test_resume_restores_scf_before_explicit_ccsd_construction(tmp_path):
    path = tmp_path / "restart.h5"
    identity = checkpoint_identity(nmo=2)
    write_checkpoint(
        path,
        identity=identity,
        completed_cycle=2,
        ccsd_correlation_hartree=-1.0,
        t1=cc.np.zeros((1, 1)),
        t2=cc.np.zeros((1, 1, 1, 1)),
    )
    events = []
    mean_field = checkpoint_mean_field(2)
    mean_field.kernel = lambda: pytest.fail("resume must not recompute SCF")

    def construct(mf, frozen=None, mo_coeff=None, mo_occ=None, gpulim=None):
        assert mf.converged and mf.e_tot == -10.5
        assert mo_coeff is mf.mo_coeff and mo_occ is mf.mo_occ
        events.append("ccsd")
        return SimpleNamespace(nocc=1, nmo=2, pool=ArrayPool())

    cluster, energy, _, metadata = cc.prepare_byteqc_coupled_cluster(
        mean_field=mean_field,
        cucc_module=SimpleNamespace(CCSD=construct),
        frozen=0,
        gpulim=1,
        restart_path=path,
        restart_identity=identity,
    )
    events.append("ao2mo")
    assert events == ["ccsd", "ao2mo"] and energy == -10.5 and metadata is not None
    t1, t2 = cc.load_ccsd_checkpoint_amplitudes(
        path,
        identity=identity,
        coupled_cluster=cluster,
        expected_metadata=metadata,
    )
    assert t1.shape == (1, 1) and t2.shape == (1, 1, 1, 1)


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


def test_focal_canonical_grid_requires_directive_adjusted_matrix():
    values = [
        "reactant=cc-pVTZ=r-tz.json",
        "reactant=cc-pVQZ=r-qz.json",
        "ts=cc-pVTZ=t-tz.json",
    ]
    grid = cc.focal_canonical_receipt_grid(values)
    assert set(grid["reactant"]) == {"tz", "qz"}
    assert set(grid["ts"]) == {"tz"}
    with pytest.raises(ValueError, match="reactant TZ/QZ and TS TZ exactly"):
        cc.focal_canonical_receipt_grid(values[:-1])
    with pytest.raises(ValueError, match="reactant TZ/QZ and TS TZ exactly"):
        cc.focal_canonical_receipt_grid(
            [*values, "ts=cc-pVQZ=forbidden-canonical-ts-qz.json"]
        )
    with pytest.raises(ValueError, match="duplicate receipt selector"):
        cc.focal_canonical_receipt_grid([*values, values[-1]])


def test_summarize_focal_uses_dlpno_basis_correction_and_tz_gate(tmp_path):
    engine_names = {
        "canonical": "byteqc-canonical-ccsd(t)",
        "dlpno": "psi4-dlpno-ccsd(t)",
    }
    energies = {
        "canonical": {
            "reactant": {"tz": (-100.0, -1.0), "qz": (-100.2, -1.2)},
            "ts": {"tz": (-99.95, -1.0)},
        },
        "dlpno": {
            "reactant": {"tz": (-99.99, -1.0), "qz": (-100.19, -1.19)},
            "ts": {"tz": (-99.9399, -1.0), "qz": (-100.135, -1.19)},
        },
    }
    paths = {engine: {} for engine in engine_names}
    for engine, roles in energies.items():
        for role, bases in roles.items():
            for basis, (scf, correlation) in bases.items():
                path = tmp_path / f"{engine}-{role}-{basis}.json"
                path.write_text(
                    json.dumps(
                        receipt(
                            engine_names[engine],
                            "cc-pVTZ" if basis == "tz" else "cc-pVQZ",
                            scf,
                            correlation,
                        )
                    )
                )
                paths[engine].setdefault(role, {})[basis] = path

    output = tmp_path / "focal-summary.json"
    result = cc.summarize_focal(
        argparse.Namespace(
            canonical=paths["canonical"],
            dlpno=paths["dlpno"],
            output=output,
        )
    )
    canonical_tz = cc.barrier_kj(
        receipt(engine_names["canonical"], "cc-pVTZ", -100.0, -1.0),
        receipt(engine_names["canonical"], "cc-pVTZ", -99.95, -1.0),
    )
    dlpno_tz = cc.barrier_kj(
        receipt(engine_names["dlpno"], "cc-pVTZ", -99.99, -1.0),
        receipt(engine_names["dlpno"], "cc-pVTZ", -99.9399, -1.0),
    )
    dlpno_cbs = (
        cc.cbs_total(
            receipt(engine_names["dlpno"], "cc-pVTZ", -99.9399, -1.0),
            receipt(engine_names["dlpno"], "cc-pVQZ", -100.135, -1.19),
        )
        - cc.cbs_total(
            receipt(engine_names["dlpno"], "cc-pVTZ", -99.99, -1.0),
            receipt(engine_names["dlpno"], "cc-pVQZ", -100.19, -1.19),
        )
    ) * cc.HARTREE_TO_KJ
    assert result["canonical_ts_qz_computed"] is False
    assert result["agreement_basis"] == "cc-pVTZ"
    assert result["dlpno_minus_canonical_kj"] == pytest.approx(dlpno_tz - canonical_tz)
    assert result["focal_point_barrier_kj"] == pytest.approx(
        canonical_tz + dlpno_cbs - dlpno_tz
    )
    assert result["gate_pass"] is True
    assert json.loads(output.read_text()) == result


def test_load_receipt_rejects_energy_arithmetic_drift(tmp_path):
    path = tmp_path / "receipt.json"
    payload = receipt("psi4-dlpno-ccsd(t)", "cc-pVTZ", -10.0, -1.0)
    payload["total_hartree"] += 0.01
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="total energy arithmetic mismatch"):
        cc.load_receipt(path, "psi4-dlpno-ccsd(t)", "cc-pVTZ")
