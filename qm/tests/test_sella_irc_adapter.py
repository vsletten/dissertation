"""Deterministic contract tests for the path-preserving Sella IRC adapter."""

import hashlib
from dataclasses import replace

import numpy as np
import pytest

from quarry.clusters import water
from quarry.pipeline import DftSettings

CHEAP = DftSettings(xc="hf", basis="sto-3g")


def _force(atom_count: int, fmax: float) -> np.ndarray:
    forces = np.zeros((atom_count, 3))
    forces[0, 0] = fmax
    return forces


def _exact_sella_pes_caches(x0, *, energy=-10.0, frozen_indices=()):
    dimension = len(x0)
    constrained = [
        3 * atom_index + axis for atom_index in frozen_indices for axis in range(3)
    ]
    free = [index for index in range(dimension) if index not in constrained]
    identity = np.eye(dimension)
    current = {
        "L": np.zeros(len(constrained)),
        "Ucons": identity[:, constrained],
        "Ufree": identity[:, free],
        "Unred": identity,
        "drdx": identity[constrained],
        "f": float(energy),
        "g": np.zeros(dimension),
        "state_hash": np.asarray(x0, dtype=float).tobytes(),
        "x": np.asarray(x0, dtype=float).copy(),
    }
    return current, {"f": None, "g": None, "x": None}


def _valid_initialization_state():
    import quarry.ts as ts_module

    transition_state = water()
    masses = np.array([15.99, 1.01, 1.02])
    dimension = transition_state.coords.size
    contract = ts_module.IrcExecutionContract(
        algorithm="sella-gonzalez-schlegel",
        step_size_angstrom=0.1,
        maximum_steps=7,
        outer_fmax_ev_per_angstrom=0.05,
        inner_fmax_ev_per_angstrom=0.01,
    )
    x0 = transition_state.coords.reshape(-1).copy()
    return ts_module.SellaIrcInitializationState(
        symbols=tuple(transition_state.symbols),
        charge=transition_state.charge,
        spin=transition_state.spin,
        frozen_indices=(),
        settings_fingerprint=hashlib.sha256(
            ts_module.frequency_settings_fingerprint(CHEAP).encode()
        ).hexdigest(),
        execution_contract=contract,
        x0=x0,
        masses_amu=masses,
        h0=np.diag([-1.0, *np.arange(2, dimension + 1, dtype=float)]),
        v0ts=(
            np.eye(1, dimension, 0).reshape(-1)
            * contract.step_size_angstrom
            / np.sqrt(masses[0])
        ),
        pes_current=_exact_sella_pes_caches(x0)[0],
        pes_last=_exact_sella_pes_caches(x0)[1],
    )


def _install_fake_irc(monkeypatch, routes):
    import sella

    import quarry.ts as ts_module

    instances = []
    calls = []
    constructor_masses = []
    diagonalizations = []
    transition_state = water().coords.copy()
    atom_count = len(transition_state)

    class FakePes:
        def __init__(self, optimizer):
            self.optimizer = optimizer
            self.curr = {"x": None, "f": None, "g": None}
            self.last = self.curr.copy()

        def get_f(self):
            return self.optimizer.energy_ev

        def get_projected_forces(self):
            return self.optimizer.projected_forces

        def set_x(self, target):
            self.optimizer.atoms.positions[:] = np.asarray(target).reshape((-1, 3))

        def set_H(self, target, *, initialized):
            assert initialized is True
            self.optimizer.restored_hessian = np.asarray(target).copy()

        def _calc_basis(self):
            current, _ = _exact_sella_pes_caches(
                self.optimizer.x0,
                energy=self.curr["f"],
                frozen_indices=self.optimizer.frozen_indices,
            )
            return tuple(current[key] for key in ("drdx", "Ucons", "Unred", "Ufree"))

        def _update_basis(self, basis):
            drdx, ucons, unred, ufree = basis
            self.curr.update(
                {
                    "L": np.zeros(drdx.shape[0]),
                    "Ucons": ucons,
                    "Ufree": ufree,
                    "Unred": unred,
                    "drdx": drdx,
                }
            )

    class FakeIRC:
        def __init__(self, atoms, **kwargs):
            assert "masses" in atoms.arrays
            constructor_masses.append(atoms.get_masses().copy())
            self.atoms = atoms
            self.kwargs = kwargs
            self.nsteps = 0
            self.energy_ev = -10.0
            self.projected_forces = _force(atom_count, 0.4)
            self.first = True
            self._terminal_converged = False
            self.pes = FakePes(self)
            self.x0 = transition_state.reshape(-1).copy()
            self.v0ts = None
            self.H0 = None
            self.pescurr = None
            self.peslast = None
            self.frozen_indices = tuple(
                index
                for constraint in atoms.constraints
                for index in getattr(constraint, "index", ())
            )
            self.restored_hessian = None
            self.terminal_convergence_calls = 0
            instances.append(self)

        def gradient_converged(self, gradient):
            return True  # ASE 3.29's ordinary force-only decision at the TS.

        def converged(self):
            self.terminal_convergence_calls += 1
            return not self.first and self._terminal_converged

        def irun(self, *, fmax, fmax_inner, steps, direction):
            if self.v0ts is None:
                diagonalizations.append(direction)
                dimension = atom_count * 3
                self.H0 = np.diag([-1.0, *np.arange(2, dimension + 1, dtype=float)])
                self.v0ts = np.zeros(dimension)
                self.v0ts[0] = self.kwargs["dx"] / np.sqrt(self.atoms.get_masses()[0])
                self.pescurr, self.peslast = _exact_sella_pes_caches(
                    self.x0, frozen_indices=self.frozen_indices
                )
            else:
                # Match Sella 2.5's new-direction branch: a restored v0ts must
                # bypass the TS kick/diagonalization and reinstall the exact
                # PES coordinates, caches, and Hessian captured at the TS.
                assert self.pescurr is not None
                assert self.peslast is not None
                assert self.H0 is not None
                self.pes.set_x(self.x0)
                self.pes.curr = self.pescurr.copy()
                self.pes.last = self.peslast.copy()
                self.pes.set_H(self.H0.copy(), initialized=True)
            calls.append(
                {
                    "fmax": fmax,
                    "fmax_inner": fmax_inner,
                    "steps": steps,
                    "direction": direction,
                }
            )
            self.atoms.positions[:] = transition_state
            self.energy_ev = -10.0
            self.projected_forces = _force(atom_count, 0.4)
            self.first = True
            self._terminal_converged = False

            def states():
                yield self.gradient_converged(np.zeros(atom_count * 3))
                for offset, energy_ev, projected_fmax, converged in routes[direction]:
                    self.atoms.positions[:] = transition_state + offset
                    self.energy_ev = energy_ev
                    self.projected_forces = _force(atom_count, projected_fmax)
                    self.first = False
                    self._terminal_converged = converged
                    self.nsteps += 1
                    yield self.gradient_converged(np.zeros(atom_count * 3))

            return states()

    monkeypatch.setattr(sella, "IRC", FakeIRC)
    monkeypatch.setattr(ts_module, "make_ase_calculator", lambda *args: object())
    return ts_module, instances, calls, constructor_masses, diagonalizations


def _converged_routes():
    return {
        "forward": [(0.1, -10.1, 0.2, False), (0.2, -10.2, 0.01, True)],
        "reverse": [(-0.1, -10.3, 0.15, False), (-0.3, -10.4, 0.02, True)],
    }


def test_path_adapter_preserves_both_directions_and_explicit_masses(monkeypatch):
    ts_module, instances, calls, constructor_masses, diagonalizations = (
        _install_fake_irc(monkeypatch, _converged_routes())
    )
    transition_state = replace(water(), frozen_indices=[0])
    frozen_masses = np.array([15.99, 1.01, 1.02])

    trace = ts_module.trace_sella_irc(
        transition_state,
        CHEAP,
        masses_amu=frozen_masses,
        fmax_ev_a=0.05,
        fmax_inner_ev_a=0.01,
        max_steps=7,
        step_size_a=0.12,
    )

    assert len(instances) == 1
    assert diagonalizations == ["forward"]
    assert [call["direction"] for call in calls] == ["forward", "reverse"]
    assert all(call["steps"] == 7 for call in calls)
    assert instances[0].kwargs["dx"] == pytest.approx(0.12)
    assert len(constructor_masses) == 1
    assert np.all(np.isfinite(constructor_masses[0]))
    assert np.all(constructor_masses[0] > 0.0)
    assert constructor_masses[0] == pytest.approx(frozen_masses)
    assert np.array_equal(trace.masses_amu, constructor_masses[0])
    assert trace.execution_contract == ts_module.IrcExecutionContract(
        algorithm="sella-gonzalez-schlegel",
        step_size_angstrom=0.12,
        maximum_steps=7,
        outer_fmax_ev_per_angstrom=0.05,
        inner_fmax_ev_per_angstrom=0.01,
    )
    with pytest.raises(ValueError):
        trace.masses_amu[0] = 1.0

    forward, reverse = trace.directions
    assert (forward.sella_direction, forward.algebraic_direction) == ("forward", 1)
    assert (reverse.sella_direction, reverse.algebraic_direction) == ("reverse", -1)
    assert [point.outer_step for point in forward.points] == [0, 1, 2]
    assert [point.outer_step for point in reverse.points] == [0, 1, 2]
    assert np.array_equal(
        forward.points[0].coordinates_angstrom, transition_state.coords
    )
    assert np.array_equal(
        reverse.points[0].coordinates_angstrom, transition_state.coords
    )
    assert forward.points[-1].electronic_energy_ev == pytest.approx(-10.2)
    assert reverse.points[-1].electronic_energy_ev == pytest.approx(-10.4)
    assert forward.points[-1].projected_fmax_ev_per_angstrom == pytest.approx(0.01)
    assert reverse.points[-1].projected_fmax_ev_per_angstrom == pytest.approx(0.02)
    assert not np.array_equal(
        forward.points[-1].coordinates_angstrom,
        reverse.points[-1].coordinates_angstrom,
    )
    with pytest.raises(ValueError):
        forward.points[0].coordinates_angstrom[0, 0] = 1.0


def test_path_adapter_checkpoints_each_new_direction_and_skips_completed(
    monkeypatch,
):
    ts_module, _, _, _, _ = _install_fake_irc(monkeypatch, _converged_routes())
    transition_state = water()
    masses = np.array([15.99, 1.01, 1.02])
    initialized = []
    original = ts_module._trace_sella_irc_resume(
        transition_state,
        CHEAP,
        masses_amu=masses,
        max_steps=7,
        _initialization_callback=initialized.append,
    )

    ts_module, _, calls, _, _ = _install_fake_irc(monkeypatch, _converged_routes())
    checkpointed = []
    with pytest.raises(TypeError, match="_completed_directions"):
        ts_module.trace_sella_irc(
            transition_state,
            CHEAP,
            masses_amu=masses,
            max_steps=7,
            _completed_directions={"forward": original.directions[0]},
        )
    resumed = ts_module._trace_sella_irc_resume(
        transition_state,
        CHEAP,
        masses_amu=masses,
        max_steps=7,
        _completed_directions={"forward": original.directions[0]},
        _direction_callback=checkpointed.append,
        _initialization_state=initialized[0],
    )

    assert [call["direction"] for call in calls] == ["reverse"]
    assert [direction.sella_direction for direction in checkpointed] == ["reverse"]
    assert resumed.directions[0] is original.directions[0]
    assert resumed.directions[1].sella_direction == "reverse"


def test_path_adapter_resume_restores_shared_initialization_without_diagonalizing(
    monkeypatch,
):
    ts_module, _, _, _, first_diagonalizations = _install_fake_irc(
        monkeypatch, _converged_routes()
    )
    transition_state = water()
    masses = np.array([15.99, 1.01, 1.02])
    initialized = []
    checkpointed = []

    def crash_after_forward(path):
        checkpointed.append(path)
        raise RuntimeError("simulated process loss after forward checkpoint")

    with pytest.raises(RuntimeError, match="process loss after forward"):
        ts_module._trace_sella_irc_resume(
            transition_state,
            CHEAP,
            masses_amu=masses,
            max_steps=7,
            _initialization_callback=initialized.append,
            _direction_callback=crash_after_forward,
        )

    assert first_diagonalizations == ["forward"]
    assert len(initialized) == 1
    assert [path.sella_direction for path in checkpointed] == ["forward"]
    fingerprint = initialized[0].fingerprint
    assert checkpointed[0].initialization_fingerprint == fingerprint

    (
        ts_module,
        resumed_instances,
        resumed_calls,
        _,
        resumed_diagonalizations,
    ) = _install_fake_irc(monkeypatch, _converged_routes())
    resumed = ts_module._trace_sella_irc_resume(
        transition_state,
        CHEAP,
        masses_amu=masses,
        max_steps=7,
        _completed_directions={"forward": checkpointed[0]},
        _initialization_state=initialized[0],
    )

    assert resumed_diagonalizations == []
    assert [call["direction"] for call in resumed_calls] == ["reverse"]
    assert len(resumed_instances) == 1
    resumed_irc = resumed_instances[0]
    assert np.array_equal(resumed_irc.restored_hessian, initialized[0].h0)
    assert np.array_equal(resumed_irc.v0ts, initialized[0].v0ts)
    assert np.array_equal(resumed_irc.pes.curr["x"], initialized[0].pes_current["x"])
    assert np.array_equal(resumed_irc.pes.curr["g"], initialized[0].pes_current["g"])
    assert resumed_irc.pes.curr["f"] == initialized[0].pes_current["f"]
    assert resumed.initialization_fingerprint == fingerprint
    assert {
        direction.initialization_fingerprint for direction in resumed.directions
    } == {fingerprint}


def test_restore_installs_checkpoint_before_fresh_real_sella_irc_has_pes_caches():
    from ase import Atoms
    from ase.calculators.lj import LennardJones
    from sella import IRC

    import quarry.ts as ts_module

    transition_state = water()
    masses = np.array([15.99, 1.01, 1.02])
    dimension = 3 * len(transition_state.symbols)
    contract = ts_module.IrcExecutionContract(
        algorithm="sella-gonzalez-schlegel",
        step_size_angstrom=0.1,
        maximum_steps=7,
        outer_fmax_ev_per_angstrom=0.05,
        inner_fmax_ev_per_angstrom=0.01,
    )
    x0 = transition_state.coords.reshape(-1).copy()
    state = ts_module.SellaIrcInitializationState(
        symbols=tuple(transition_state.symbols),
        charge=transition_state.charge,
        spin=transition_state.spin,
        frozen_indices=tuple(transition_state.frozen_indices or ()),
        settings_fingerprint=hashlib.sha256(
            ts_module.frequency_settings_fingerprint(CHEAP).encode()
        ).hexdigest(),
        execution_contract=contract,
        x0=x0,
        masses_amu=masses,
        h0=np.diag([-1.0, *np.arange(2, dimension + 1, dtype=float)]),
        v0ts=(
            np.eye(1, dimension, 0).reshape(-1)
            * contract.step_size_angstrom
            / np.sqrt(masses[0])
        ),
        pes_current=_exact_sella_pes_caches(x0, energy=-1.0)[0],
        pes_last=_exact_sella_pes_caches(x0, energy=-1.0)[1],
    )
    atoms = Atoms(symbols=transition_state.symbols, positions=transition_state.coords)
    atoms.set_masses(masses)
    atoms.calc = LennardJones()
    irc = IRC(atoms, logfile=None, dx=contract.step_size_angstrom)

    assert not hasattr(irc, "pescurr")
    ts_module._restore_sella_irc_initialization(
        irc,
        state,
        transition_state,
        CHEAP,
        masses,
        contract,
    )
    irc.pes.kick = lambda *_args, **_kwargs: pytest.fail(
        "restored real Sella IRC performed a second TS diagonalization"
    )
    states = irc.irun(
        fmax=contract.outer_fmax_ev_per_angstrom,
        fmax_inner=contract.inner_fmax_ev_per_angstrom,
        steps=contract.maximum_steps,
        direction="reverse",
    )
    states.close()

    assert np.array_equal(irc.H0, state.h0)
    assert np.array_equal(irc.v0ts, state.v0ts)
    expected_pes_x = state.pes_current["x"]
    assert isinstance(expected_pes_x, np.ndarray)
    assert np.array_equal(irc.pes.curr["x"], expected_pes_x)


@pytest.mark.parametrize(
    ("corruption", "match"),
    [
        ("extra-key", "exact Sella 2.5 cache fields"),
        ("x-drift", "x does not equal"),
        ("state-hash-drift", "state_hash does not equal"),
        ("nonpristine-last", "pristine previous"),
    ],
)
def test_sella_restart_rejects_adversarial_cache_schema(corruption, match):
    state = _valid_initialization_state()
    current = dict(state.pes_current)
    previous = dict(state.pes_last)
    if corruption == "extra-key":
        current["attacker"] = np.zeros(1)
    elif corruption == "x-drift":
        current["x"] = np.asarray(current["x"]).copy()
        current["x"][0] += 1.0e-3
    elif corruption == "state-hash-drift":
        current["state_hash"] = b"0" * len(state.x0.tobytes())
    else:
        previous["f"] = -1.0
    with pytest.raises(RuntimeError, match=match):
        replace(state, pes_current=current, pes_last=previous)


@pytest.mark.parametrize(
    ("corruption", "match"),
    [
        ("shape", "array shape"),
        ("byte-total", "byte total"),
        ("too-many-atoms", "atom inventory"),
    ],
)
def test_sella_restart_payload_is_strictly_bounded(corruption, match):
    import quarry.ts as ts_module

    payload = ts_module._sella_irc_initialization_payload(_valid_initialization_state())
    if corruption == "shape":
        payload["pes_current"]["g"]["shape"] = [10]
    elif corruption == "byte-total":
        payload["pes_current"]["g"]["byte_count"] += 8
    else:
        payload["symbols"] = ["H"] * (ts_module._MAX_SELLA_RESTART_ATOMS + 1)
    with pytest.raises(ValueError, match=match):
        ts_module._sella_irc_initialization_from_payload(payload)


@pytest.mark.parametrize(
    ("v0ts", "match"),
    [
        ("wrong-norm", "norm does not match IRC dx"),
        ("wrong-mode", "lowest mass-weighted H0 eigendirection"),
        ("wrong-sign", "sign convention"),
    ],
)
def test_sella_restart_binds_mass_weighted_kick_to_hessian_masses_and_dx(v0ts, match):
    state = _valid_initialization_state()
    kick = np.zeros_like(state.v0ts)
    if v0ts == "wrong-norm":
        kick[0] = (
            2.0
            * state.execution_contract.step_size_angstrom
            / np.sqrt(state.masses_amu[0])
        )
    elif v0ts == "wrong-mode":
        kick[3] = state.execution_contract.step_size_angstrom / np.sqrt(
            state.masses_amu[1]
        )
    else:
        kick[0] = -state.execution_contract.step_size_angstrom / np.sqrt(
            state.masses_amu[0]
        )
    with pytest.raises(ValueError, match=match):
        replace(state, v0ts=kick)


def test_sella_restart_rejects_wrong_mode_despite_huge_spectator_eigenvalue():
    state = _valid_initialization_state()
    h0 = np.eye(state.h0.shape[0])
    h0[0, 0] = -1.0
    h0[-1, -1] = 1.0e30
    wrong = np.zeros_like(state.v0ts)
    wrong[3] = state.execution_contract.step_size_angstrom / np.sqrt(
        state.masses_amu[1]
    )

    with pytest.raises(ValueError, match="lowest mass-weighted H0 eigendirection"):
        replace(state, h0=h0, v0ts=wrong)


@pytest.mark.parametrize("lowest", [0.0, 1.0])
def test_sella_restart_rejects_h0_without_resolvable_unstable_mode(lowest):
    state = _valid_initialization_state()
    h0 = np.eye(state.h0.shape[0])
    h0[0, 0] = lowest

    with pytest.raises(ValueError, match="lowest mass-weighted H0 eigendirection"):
        replace(state, h0=h0)


def test_sella_restore_reconstructs_and_rejects_derived_cache_drift():
    from ase import Atoms
    from ase.calculators.lj import LennardJones
    from sella import IRC

    import quarry.ts as ts_module

    state = _valid_initialization_state()
    corrupted = dict(state.pes_current)
    corrupted["Ufree"] = np.asarray(corrupted["Ufree"]).copy()
    corrupted["Ufree"][0, 0] = 0.5
    state = replace(state, pes_current=corrupted)
    transition_state = water()
    atoms = Atoms(symbols=transition_state.symbols, positions=transition_state.coords)
    atoms.set_masses(state.masses_amu)
    atoms.calc = LennardJones()
    irc = IRC(atoms, logfile=None, dx=state.execution_contract.step_size_angstrom)

    with pytest.raises(ValueError, match="Ufree does not match reconstruction"):
        ts_module._restore_sella_irc_initialization(
            irc,
            state,
            transition_state,
            CHEAP,
            state.masses_amu,
            state.execution_contract,
        )


def test_path_adapter_resume_fails_closed_without_exact_shared_initialization(
    monkeypatch,
):
    ts_module, _, _, _, _ = _install_fake_irc(monkeypatch, _converged_routes())
    transition_state = water()
    masses = np.array([15.99, 1.01, 1.02])
    original = ts_module.trace_sella_irc(
        transition_state,
        CHEAP,
        masses_amu=masses,
        max_steps=7,
    )

    with pytest.raises(ValueError, match="shared Sella initialization"):
        ts_module._trace_sella_irc_resume(
            transition_state,
            CHEAP,
            masses_amu=masses,
            max_steps=7,
            _completed_directions={"forward": original.directions[0]},
        )


@pytest.mark.parametrize("masses", [[1.0, 2.0], [1.0, np.nan, 2.0]])
def test_path_adapter_rejects_invalid_explicit_masses_before_sella(monkeypatch, masses):
    ts_module, instances, _, _, _ = _install_fake_irc(monkeypatch, _converged_routes())

    with pytest.raises(ValueError, match="masses|finite"):
        ts_module.trace_sella_irc(water(), CHEAP, masses_amu=masses)

    assert instances == []


def test_path_adapter_rejects_a_zero_step_direction(monkeypatch):
    routes = _converged_routes()
    routes["forward"] = []
    ts_module, _, _, _, _ = _install_fake_irc(monkeypatch, routes)

    with pytest.raises(RuntimeError, match="forward.*zero outer IRC steps"):
        ts_module.trace_sella_irc(water(), CHEAP, max_steps=3)


def test_path_adapter_fails_closed_on_exhaustion(monkeypatch):
    routes = _converged_routes()
    routes["reverse"] = [(-0.1, -10.1, 0.2, False)] * 3
    ts_module, instances, _, _, _ = _install_fake_irc(monkeypatch, routes)

    with pytest.raises(
        RuntimeError, match="reverse direction did not converge.*3 steps"
    ):
        ts_module.trace_sella_irc(water(), CHEAP, max_steps=3)

    assert instances[0].terminal_convergence_calls > 0


def test_path_adapter_rejects_more_than_the_bounded_point_count(monkeypatch):
    routes = _converged_routes()
    routes["forward"] = [(0.1, -10.1, 0.01, True)] * 3
    ts_module, _, _, _, _ = _install_fake_irc(monkeypatch, routes)

    with pytest.raises(RuntimeError, match="forward.*point bound"):
        ts_module.trace_sella_irc(water(), CHEAP, max_steps=2)


@pytest.mark.parametrize("bad_field", ["coordinates", "energy", "forces"])
def test_path_adapter_rejects_nonfinite_points(monkeypatch, bad_field):
    coords = np.nan if bad_field == "coordinates" else 0.1
    energy = np.nan if bad_field == "energy" else -10.1
    force = np.nan if bad_field == "forces" else 0.01
    routes = _converged_routes()
    routes["forward"] = [(coords, energy, force, True)]
    ts_module, _, _, _, _ = _install_fake_irc(monkeypatch, routes)

    with pytest.raises(RuntimeError, match=f"non-finite {bad_field}"):
        ts_module.trace_sella_irc(water(), CHEAP, max_steps=2)


def test_path_adapter_rejects_endpoint_above_outer_force_threshold(monkeypatch):
    routes = _converged_routes()
    routes["forward"] = [(0.1, -10.1, 0.05, True)]
    ts_module, _, _, _, _ = _install_fake_irc(monkeypatch, routes)

    with pytest.raises(RuntimeError, match="endpoint projected fmax"):
        ts_module.trace_sella_irc(water(), CHEAP, fmax_ev_a=0.05, max_steps=2)


def test_full_irc_compatibility_wrapper_returns_endpoints(monkeypatch):
    ts_module, _, _, _, _ = _install_fake_irc(monkeypatch, _converged_routes())
    transition_state = water()

    sella_forward, sella_reverse = ts_module.full_irc(
        transition_state, CHEAP, max_steps=7
    )

    assert sella_forward.name == f"{transition_state.name}-irc-forward"
    assert sella_reverse.name == f"{transition_state.name}-irc-reverse"
    assert np.array_equal(sella_forward.coords, transition_state.coords + 0.2)
    assert np.array_equal(sella_reverse.coords, transition_state.coords - 0.3)
