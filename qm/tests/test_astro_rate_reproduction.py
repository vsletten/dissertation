"""Fast gates for the D2a astrochemical-rate campaign driver."""

import numpy as np
import pytest

from scripts import astro_rate_reproduction as astro


def test_literature_fit_uses_modified_arrhenius_form():
    fit = astro.LiteratureFit(
        alpha_s=7.64e10,
        beta=1.0,
        gamma_k=1339.9,
        t0_k=153.2,
        source="fixture",
        valid_floor_k=60.0,
    )

    expected = (
        7.64e10
        * (100.0 / 300.0)
        * np.exp(-1339.9 * (100.0 + 153.2) / (100.0**2 + 153.2**2))
    )
    assert fit.rate(100.0) == pytest.approx(expected, rel=1e-14)


def test_target_set_has_four_saddles_and_barrierless_control():
    targets = astro.reactions(gpu=True, basis="def2-svp")

    assert set(targets) == {
        "h-h2co-ch3o",
        "h-h2co-ch2oh",
        "h-h2co-h2-hco",
        "oh-h2",
        "h-oh-control",
    }
    assert sum(reaction.barrierless for reaction in targets.values()) == 1
    assert all(reaction.method.use_gpu for reaction in targets.values())
    assert all(reaction.cluster.spin in (0, 1) for reaction in targets.values())


def test_literature_fits_are_positive_over_campaign_grid():
    targets = astro.reactions(gpu=False, basis="def2-svp")

    for reaction in targets.values():
        if reaction.literature_fit is None:
            continue
        rates = [reaction.literature_fit.rate(t) for t in astro.TEMPERATURES]
        assert np.all(np.isfinite(rates))
        assert np.all(np.asarray(rates) > 0.0)


def test_geometry_hash_changes_with_coordinates():
    control = astro.reactions(gpu=False, basis="sto-3g")["h-oh-control"].cluster
    moved = control.coords.copy()
    moved[-1, 0] += 0.1

    assert astro.geometry_hash(control) != astro.geometry_hash(
        type(control)(
            name=control.name,
            symbols=control.symbols,
            coords=moved,
            charge=control.charge,
            spin=control.spin,
        )
    )
