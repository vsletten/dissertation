from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts import astro_rate_reproduction as gas
from scripts import d2c_input_bundle as bundle
from scripts import surface_rate_protocol as surface


def _xyz(name: str, shift: float) -> str:
    return (
        f"3\n{name}\n"
        f"H {shift:.8f} 0.00000000 0.00000000\n"
        f"C {1.00000000 + shift:.8f} 0.00000000 0.00000000\n"
        f"O {2.20000000 + shift:.8f} 0.00000000 0.00000000\n"
    )


def _write_fixture_source(root: Path) -> str:
    (root / "results.json").write_text("{}\n")
    for route in bundle.ROUTES:
        route_root = root / route
        route_root.mkdir()
        geometries = {
            "reactant": _xyz(f"{route}-reactant", 0.0),
            "ts": _xyz(f"{route}-ts", 0.1),
            "product": _xyz(f"{route}-product", 0.2),
        }
        (route_root / "irc_back.xyz").write_text(geometries["reactant"])
        (route_root / "ts.xyz").write_text(geometries["ts"])
        (route_root / "irc_fwd.xyz").write_text(geometries["product"])
        geometry_hashes = {
            role: bundle.geometry_hash_xyz(route_root / filename)
            for role, filename in {
                "reactant": "irc_back.xyz",
                "ts": "ts.xyz",
                "product": "irc_fwd.xyz",
            }.items()
        }
        result = {
            "key": route,
            "family": route.removesuffix("-1w-oside")
            .removesuffix("-1w-cside")
            .removesuffix("-1w"),
            "site": "1w",
            "n_water": 1,
            "classification": "first-order-saddle",
            "cc_delta_source": "direct-tz",
            "method": {"xc": "pwb6k", "basis": "def2-svp", "use_gpu": True},
            "provenance": {
                "generated_utc": "2026-08-28T00:00:00Z",
                "git_sha": "1" * 40,
                "geometry_sha256": geometry_hashes,
            },
        }
        (route_root / "results.json").write_text(json.dumps(result) + "\n")
        for role in bundle.ROLES:
            cc = {
                "identity": {
                    "basis": "cc-pvtz",
                    "geometry_sha256": geometry_hashes[role],
                    "method": "uhf-uccsd(t)",
                }
            }
            (route_root / f"cc-{role}-tz.json").write_text(json.dumps(cc) + "\n")
    return bundle.sha256_path(root / "results.json")


def test_bundle_creation_is_hash_bound_and_detects_drift(tmp_path: Path):
    source = tmp_path / "source"
    source.mkdir()
    aggregate_sha = _write_fixture_source(source)
    output = tmp_path / "bundle"

    manifest = bundle.create_bundle(
        source,
        output,
        source_revision="2" * 40,
        expected_aggregate_sha256=aggregate_sha,
    )

    assert manifest["selection"]["routes"] == list(bundle.ROUTES)
    assert bundle.verify_bundle(output) == manifest
    assert set(manifest["routes"]) == set(bundle.ROUTES)
    assert manifest["source"]["revision"] == "2" * 40

    target = output / bundle.ROUTES[0] / "ts.xyz"
    target.write_text(target.read_text() + "\n")
    with pytest.raises(ValueError, match="bundled input drifted"):
        bundle.verify_bundle(output)


def test_bundle_rejects_aggregate_or_geometry_drift(tmp_path: Path):
    source = tmp_path / "source"
    source.mkdir()
    aggregate_sha = _write_fixture_source(source)

    with pytest.raises(ValueError, match="aggregate results"):
        bundle.create_bundle(
            source,
            tmp_path / "bad-aggregate",
            source_revision="2" * 40,
            expected_aggregate_sha256="0" * 64,
        )

    route = source / bundle.ROUTES[0]
    (route / "ts.xyz").write_text(_xyz("changed", 0.3))
    with pytest.raises(ValueError, match="transition-state geometry hash"):
        bundle.create_bundle(
            source,
            tmp_path / "bad-geometry",
            source_revision="2" * 40,
            expected_aggregate_sha256=aggregate_sha,
        )


def test_song_kaestner_lh_table_4_transcription():
    # Table 4: R(1) CH3O alpha=3.14e10, gamma=830 K.  The neighboring
    # R(2) CH2OH gamma is 3146 K; the old code accidentally used it as
    # R(1)'s prefactor and inflated that literature anchor by ~1000x.
    anchor_path = (
        Path(__file__).parents[1]
        / "data"
        / "D2c-instanton-tier"
        / "literature-anchors.json"
    )
    anchors = json.loads(anchor_path.read_text())["anchors"]["song_kaestner_2017_lh"]
    gas_fit = gas.reactions(gpu=False, basis="sto-3g")["h-h2co-ch3o"].literature_fit
    surface_fit = surface.reactions(gpu=False, basis="sto-3g")[
        "h-h2co-ch3o"
    ].literature_fit
    assert gas_fit is not None and surface_fit is not None
    ch3o = anchors["h_h2co_ch3o"]
    assert gas_fit.alpha_s == surface_fit.alpha_s == pytest.approx(ch3o["alpha_s^-1"])
    assert gas_fit.gamma_k == surface_fit.gamma_k == pytest.approx(ch3o["gamma_k"])
    assert gas_fit.t0_k == surface_fit.t0_k == pytest.approx(ch3o["t0_k"])
    assert (
        gas_fit.valid_floor_k
        == surface_fit.valid_floor_k
        == pytest.approx(ch3o["valid_floor_k"])
    )
