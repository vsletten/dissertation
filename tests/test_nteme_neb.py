from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

MODULE_PATH = Path(__file__).parents[1] / "scripts" / "nteme_neb.py"
SPEC = importlib.util.spec_from_file_location("nteme_neb", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
nteme_neb = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = nteme_neb
SPEC.loader.exec_module(nteme_neb)


def test_restricted_cell_preserves_volume() -> None:
    cell = nteme_neb.Cell(
        30.900625, 26.747522, 19.549356, 89.862673, 95.576801, 89.985803
    )
    lx, ly, lz, _xy, _xz, _yz = cell.restricted()
    expected = (
        cell.a
        * cell.b
        * cell.c
        * math.sqrt(
            1
            + 2
            * math.cos(math.radians(cell.alpha))
            * math.cos(math.radians(cell.beta))
            * math.cos(math.radians(cell.gamma))
            - math.cos(math.radians(cell.alpha)) ** 2
            - math.cos(math.radians(cell.beta)) ** 2
            - math.cos(math.radians(cell.gamma)) ** 2
        )
    )
    assert math.isclose(lx * ly * lz, expected, rel_tol=1.0e-12)


def test_prepare_divacancy_route_preserves_ids_and_sets_argon() -> None:
    atoms = [
        nteme_neb.Atom(1, "K", 0.0, 0.0, 0.0),
        nteme_neb.Atom(2, "K", 1.0, 0.0, 0.0),
        nteme_neb.Atom(3, "K", 0.0, 1.0, 0.0),
    ]
    cell = nteme_neb.Cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0)
    route = nteme_neb.Barrier("divacancy", "Ar", 1, 2, 3, 67.0, None)
    prepared, final = nteme_neb.prepare_route(atoms, cell, route)
    assert prepared == [nteme_neb.Atom(1, "Ar", 0.0, 0.0, 0.0)]
    assert final == (1.0, 0.0, 0.0)


def test_minimum_image_uses_short_periodic_hop() -> None:
    cell = nteme_neb.Cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0)
    a = nteme_neb.Atom(1, "K", 9.5, 1.0, 1.0)
    b = nteme_neb.Atom(2, "K", 0.5, 1.0, 1.0)
    dx, dy, dz = nteme_neb.minimum_image_delta(a, b, cell)
    assert math.isclose(dx, 1.0, abs_tol=1.0e-12)
    assert math.isclose(dy, 0.0, abs_tol=1.0e-12)
    assert math.isclose(dz, 0.0, abs_tol=1.0e-12)
