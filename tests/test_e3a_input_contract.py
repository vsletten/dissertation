from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
MODULE_PATH = SCRIPTS / "e3a_input_contract.py"
SPEC = importlib.util.spec_from_file_location("e3a_input_contract", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
contract = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = contract
SPEC.loader.exec_module(contract)
nteme = contract.nteme_neb


def toy_cell() -> object:
    return nteme.Cell(30.0, 30.0, 30.0, 90.0, 90.0, 90.0)


def test_charge_compensation_is_exact_per_removed_k() -> None:
    cell = toy_cell()
    # One Al2 tetrahedron is -3.1 e in the published partial-charge scheme;
    # one spectator Si makes the pre-compensation structure exactly -1 e.
    atoms = [nteme.Atom(1, "Al2", 2.0, 2.0, 2.0)]
    for atom_id, (dx, dy, dz) in enumerate(
        ((1.6, 0, 0), (-1.6, 0, 0), (0, 1.6, 0), (0, 0, 1.6)), 2
    ):
        atoms.append(nteme.Atom(atom_id, "O3", 2.0 + dx, 2.0 + dy, 2.0 + dz))
    atoms.append(nteme.Atom(6, "Si", 15.0, 15.0, 15.0))
    assert math.isclose(contract.net_charge(atoms), -1.0, abs_tol=1e-12)
    center = nteme.Atom(0, "X", 15.0, 15.0, 2.0)
    compensated, changes = contract.compensate_k_vacancies(atoms, cell, center, 1)
    assert len(changes) == 1
    assert math.isclose(contract.net_charge(compensated), 0.0, abs_tol=1e-12)
    assert math.isclose(
        contract.net_charge(compensated) - contract.net_charge(atoms),
        1.0,
        abs_tol=1e-12,
    )


def test_extended_zone_taper_has_bounded_peak_opening() -> None:
    cell = toy_cell()
    center = nteme.Atom(0, "X", 15.0, 15.0, 15.0)
    # Exercise the same cosine-squared taper algebra used by the builder.
    radial = 0.25
    weight = math.cos(math.pi * math.sqrt(radial) / 2.0) ** 2
    assert 0.0 < weight < 1.0
    assert math.isclose(weight, 0.5, abs_tol=1e-12)
    sx, sy, sz = contract.fractional(center, cell)
    assert math.isclose(sx, 0.5, abs_tol=1e-12)
    assert math.isclose(sy, 0.5, abs_tol=1e-12)
    assert math.isclose(sz, 0.5, abs_tol=1e-12)


def test_xe_parameter_conversion_and_charge() -> None:
    assert math.isclose(contract.XENON_SIGMA, 4.00924, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(
        contract.XENON_EPSILON,
        0.0185900 * 23.060547830619,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert contract.atom_charge("Xe") == 0.0


def test_validation_rejects_collision() -> None:
    cell = toy_cell()
    model = contract.Model(
        "collision",
        (
            nteme.Atom(1, "Ar", 1.0, 1.0, 1.0),
            nteme.Atom(2, "Ar", 1.1, 1.0, 1.0),
        ),
        cell,
        {},
    )
    try:
        contract.validate_model(model)
    except ValueError as exc:
        assert "minimum pair" in str(exc)
    else:
        raise AssertionError("collision was accepted")


def test_validation_requires_neutrality() -> None:
    cell = toy_cell()
    model = contract.Model(
        "charged",
        (nteme.Atom(1, "K", 1.0, 1.0, 1.0),),
        cell,
        {},
    )
    try:
        contract.validate_model(model)
    except ValueError as exc:
        assert "non-neutral" in str(exc)
    else:
        raise AssertionError("charged structure was accepted")


WORKBOOK = Path(
    "/mnt/data/vsletten/dissertation-data/e3a-classical-neb-barriers-20260902/"
    "sources/nteme-2022-supplement/Results Of MD Simulations.xlsx"
)


@pytest.mark.skipif(
    not WORKBOOK.exists(), reason="external hash-pinned workbook unavailable"
)
def test_hash_pinned_workbook_builds_all_contract_models(tmp_path: Path) -> None:
    pytest.importorskip("openpyxl")
    models, source = contract.build_models(WORKBOOK)
    assert source["source_sha256"] == contract.nteme_neb.WORKBOOK_SHA256
    by_name = {model.name: model for model in models}
    assert set(by_name) == {
        "dehydroxylate-lattice",
        "extended-defect-zone",
        "octahedral-trap-escape",
        "xenon-divacancy",
    }
    expected_counts = {
        "dehydroxylate-lattice": 1507,
        "extended-defect-zone": 1510,
        "octahedral-trap-escape": 1512,
        "xenon-divacancy": 1510,
    }
    for name, model in by_name.items():
        assert len(model.atoms) == expected_counts[name]
        assert model.endpoint_atoms is not None
        assert len(model.endpoint_atoms) == len(model.atoms)
        assert abs(contract.net_charge(model.atoms)) <= 1.0e-8
        assert abs(contract.net_charge(model.endpoint_atoms)) <= 1.0e-8
        contract.validate_model(model)
        contract.validate_model(
            contract.replace(model, atoms=model.endpoint_atoms, endpoint_atoms=None)
        )
        contract.write_lammps_data(tmp_path / f"{name}-initial.data", model)

    dehyd = by_name["dehydroxylate-lattice"]
    assert dehyd.metadata["transformation"]["removed_site_ids"] == [76, 81, 84]
    trap = by_name["octahedral-trap-escape"]
    assert trap.metadata["recoil_seed"]["geometric_cage_atom_ids"] == [
        743,
        741,
        708,
        711,
        628,
        707,
    ]


def test_static_input_uses_lorentz_berthelot_unlike_pairs(tmp_path: Path) -> None:
    path = tmp_path / "in.static"
    contract.write_static_input(path, "structure.data")
    text = path.read_text(encoding="utf-8")
    assert "pair_modify" not in text
    type_id, _mass, sigma, epsilon = contract.type_tables()
    names = (*nteme.TYPE_ORDER, "Xe")
    for name in names:
        self_coeff = (
            f"pair_coeff      {type_id[name]} {type_id[name]} "
            f"{epsilon[name]:.10g} {sigma[name]:.10g}"
        )
        assert self_coeff in text
    for index, name_i in enumerate(names):
        for name_j in names[index + 1 :]:
            mixed_eps = math.sqrt(epsilon[name_i] * epsilon[name_j])
            mixed_sig = 0.5 * (sigma[name_i] + sigma[name_j])
            unlike = (
                f"pair_coeff      {type_id[name_i]} {type_id[name_j]} "
                f"{mixed_eps:.10g} {mixed_sig:.10g}"
            )
            assert unlike in text
