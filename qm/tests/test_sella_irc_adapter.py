"""Deterministic contract tests for the path-preserving Sella IRC adapter."""

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


def _install_fake_irc(monkeypatch, routes):
    import sella

    import quarry.ts as ts_module

    instances = []
    calls = []
    constructor_masses = []
    transition_state = water().coords.copy()
    atom_count = len(transition_state)

    class FakePes:
        def __init__(self, optimizer):
            self.optimizer = optimizer

        def get_f(self):
            return self.optimizer.energy_ev

        def get_projected_forces(self):
            return self.optimizer.projected_forces

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
            self.terminal_convergence_calls = 0
            instances.append(self)

        def gradient_converged(self, gradient):
            return True  # ASE 3.29's ordinary force-only decision at the TS.

        def converged(self):
            self.terminal_convergence_calls += 1
            return not self.first and self._terminal_converged

        def irun(self, *, fmax, fmax_inner, steps, direction):
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
    return ts_module, instances, calls, constructor_masses


def _converged_routes():
    return {
        "forward": [(0.1, -10.1, 0.2, False), (0.2, -10.2, 0.01, True)],
        "reverse": [(-0.1, -10.3, 0.15, False), (-0.3, -10.4, 0.02, True)],
    }


def test_path_adapter_preserves_both_directions_and_explicit_masses(monkeypatch):
    ts_module, instances, calls, constructor_masses = _install_fake_irc(
        monkeypatch, _converged_routes()
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
    assert [call["direction"] for call in calls] == ["forward", "reverse"]
    assert all(call["steps"] == 7 for call in calls)
    assert instances[0].kwargs["dx"] == pytest.approx(0.12)
    assert len(constructor_masses) == 1
    assert np.all(np.isfinite(constructor_masses[0]))
    assert np.all(constructor_masses[0] > 0.0)
    assert constructor_masses[0] == pytest.approx(frozen_masses)
    assert np.array_equal(trace.masses_amu, constructor_masses[0])
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


@pytest.mark.parametrize("masses", [[1.0, 2.0], [1.0, np.nan, 2.0]])
def test_path_adapter_rejects_invalid_explicit_masses_before_sella(monkeypatch, masses):
    ts_module, instances, _, _ = _install_fake_irc(monkeypatch, _converged_routes())

    with pytest.raises(ValueError, match="masses|finite"):
        ts_module.trace_sella_irc(water(), CHEAP, masses_amu=masses)

    assert instances == []


def test_path_adapter_rejects_a_zero_step_direction(monkeypatch):
    routes = _converged_routes()
    routes["forward"] = []
    ts_module, _, _, _ = _install_fake_irc(monkeypatch, routes)

    with pytest.raises(RuntimeError, match="forward.*zero outer IRC steps"):
        ts_module.trace_sella_irc(water(), CHEAP, max_steps=3)


def test_path_adapter_fails_closed_on_exhaustion(monkeypatch):
    routes = _converged_routes()
    routes["reverse"] = [(-0.1, -10.1, 0.2, False)] * 3
    ts_module, instances, _, _ = _install_fake_irc(monkeypatch, routes)

    with pytest.raises(
        RuntimeError, match="reverse direction did not converge.*3 steps"
    ):
        ts_module.trace_sella_irc(water(), CHEAP, max_steps=3)

    assert instances[0].terminal_convergence_calls > 0


def test_path_adapter_rejects_more_than_the_bounded_point_count(monkeypatch):
    routes = _converged_routes()
    routes["forward"] = [(0.1, -10.1, 0.01, True)] * 3
    ts_module, _, _, _ = _install_fake_irc(monkeypatch, routes)

    with pytest.raises(RuntimeError, match="forward.*point bound"):
        ts_module.trace_sella_irc(water(), CHEAP, max_steps=2)


@pytest.mark.parametrize("bad_field", ["coordinates", "energy", "forces"])
def test_path_adapter_rejects_nonfinite_points(monkeypatch, bad_field):
    coords = np.nan if bad_field == "coordinates" else 0.1
    energy = np.nan if bad_field == "energy" else -10.1
    force = np.nan if bad_field == "forces" else 0.01
    routes = _converged_routes()
    routes["forward"] = [(coords, energy, force, True)]
    ts_module, _, _, _ = _install_fake_irc(monkeypatch, routes)

    with pytest.raises(RuntimeError, match=f"non-finite {bad_field}"):
        ts_module.trace_sella_irc(water(), CHEAP, max_steps=2)


def test_path_adapter_rejects_endpoint_above_outer_force_threshold(monkeypatch):
    routes = _converged_routes()
    routes["forward"] = [(0.1, -10.1, 0.05, True)]
    ts_module, _, _, _ = _install_fake_irc(monkeypatch, routes)

    with pytest.raises(RuntimeError, match="endpoint projected fmax"):
        ts_module.trace_sella_irc(water(), CHEAP, fmax_ev_a=0.05, max_steps=2)


def test_full_irc_compatibility_wrapper_returns_endpoints(monkeypatch):
    ts_module, _, _, _ = _install_fake_irc(monkeypatch, _converged_routes())
    transition_state = water()

    sella_forward, sella_reverse = ts_module.full_irc(
        transition_state, CHEAP, max_steps=7
    )

    assert sella_forward.name == f"{transition_state.name}-irc-forward"
    assert sella_reverse.name == f"{transition_state.name}-irc-reverse"
    assert np.array_equal(sella_forward.coords, transition_state.coords + 0.2)
    assert np.array_equal(sella_reverse.coords, transition_state.coords - 0.3)
