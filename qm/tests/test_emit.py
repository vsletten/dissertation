"""Emitter gates: TOML fragment shape, data.rxn format fidelity, and the
petra round-trip (emitted deck compiled by the in-repo petra CLI)."""

import shutil
import subprocess
import tomllib
from pathlib import Path

import pytest

from quarry.emit import (
    AdsorptionBlock,
    ArrheniusRate,
    ByCountModifier,
    EyringRate,
    ReactionEmit,
    RxnBlock,
    data_rxn,
    petra_fragment,
    splice_into_deck,
)
from quarry.store import Provenance

DATA = Path(__file__).parent / "data"
REPO_ROOT = Path(__file__).resolve().parents[2]
PETRA = REPO_ROOT / "petra"

PROV = Provenance(
    method="wb97m-v/def2-tzvpd/smd(water)//r2scan-3c",
    engine="pyscf",
    geometry_hash="c0ffee" * 10 + "abcd",
    date="2026-08-17T00:00:00+00:00",
)


def _sample_reactions() -> list[ReactionEmit]:
    detach = ReactionEmit(
        name="detach",
        center={"kind": "M_site", "state": ["occupied"]},
        rate=EyringRate(dh=60.0, ds=-0.02),
        produces=["M"],
        modifiers=[
            ByCountModifier(
                select={"distance": 1, "state": ["occupied"]},
                dea=[0.0, 8.4, 16.8, 25.2, 33.6, 42.0, 50.4],
                provenance=PROV,
            )
        ],
        effects=[{"target": "center", "set": "empty"}],
        provenance=PROV,
    )
    attach = ReactionEmit(
        name="attach",
        center={"kind": "M_site", "state": ["empty"]},
        guards=[{"distance": 1, "min": 1, "state": ["occupied"]}],
        rate=ArrheniusRate(prefactor=1.0e13, ea=41.84),
        effects=[{"target": "center", "set": "occupied"}],
        provenance=PROV,
    )
    # petra v1 has no auto-reverse: attachment consumes M from solution.
    attach.produces = []
    attach_consumes = attach  # explicitness for the reader
    return [detach, attach_consumes]


class TestPetraFragment:
    def test_is_valid_toml_with_expected_shape(self):
        frag = petra_fragment(_sample_reactions())
        doc = tomllib.loads(frag)
        assert len(doc["reactions"]) == 2
        detach, attach = doc["reactions"]
        assert detach["rate"]["eyring"] == {"dh": 60.0, "ds": -0.02}
        assert detach["modifiers"][0]["by_count"]["dea"][1] == 8.4
        assert attach["rate"]["arrhenius"]["prefactor"] == 1.0e13
        assert attach["guards"][0]["min"] == 1

    def test_provenance_annotations_present(self):
        frag = petra_fragment(_sample_reactions())
        assert "provenance: method=wb97m-v" in frag
        assert "geometry=sha256:c0ffeec0ffee" in frag
        assert frag.count("# provenance:") == 3  # 2 reactions + 1 modifier

    def test_splice_requires_marker_and_units(self):
        with pytest.raises(ValueError, match="marker"):
            splice_into_deck('units = "kJ/mol"', _sample_reactions())
        with pytest.raises(ValueError, match="units"):
            splice_into_deck("# {QUARRY_REACTIONS}", _sample_reactions())

    def test_spliced_template_parses(self):
        template = (DATA / "template_deck.toml").read_text()
        deck = splice_into_deck(template, _sample_reactions())
        doc = tomllib.loads(deck)
        assert doc["deck"]["units"] == "kJ/mol"
        assert [r["name"] for r in doc["reactions"]] == ["detach", "attach"]


@pytest.mark.petra
class TestPetraRoundTrip:
    def test_emitted_deck_compiles_and_runs(self, tmp_path):
        if shutil.which("cargo") is None:
            pytest.skip("cargo not available")
        template = (DATA / "template_deck.toml").read_text()
        deck_path = tmp_path / "roundtrip.toml"
        deck_path.write_text(splice_into_deck(template, _sample_reactions()))
        result = subprocess.run(
            [
                "cargo",
                "run",
                "--quiet",
                "-p",
                "petra-cli",
                "--",
                str(deck_path),
                "--steps",
                "10",
                "--seed",
                "42",
                "--out",
                str(tmp_path / "out"),
            ],
            cwd=PETRA,
            capture_output=True,
            text=True,
            timeout=600,
        )
        assert result.returncode == 0, f"petra rejected the deck:\n{result.stderr}"


class TestDataRxn:
    def _render(self) -> str:
        return data_rxn(
            temperature=8000.0,
            dmu_si=-1.0,
            dmu_al=-1.0,
            hydrolysis=[
                RxnBlock(
                    from_state=301,
                    to_state=302,
                    comment="si-o-si hydrolysis",
                    pairs=[(1.0, 2.6)] * 15,
                )
            ],
            adsorption=[
                AdsorptionBlock(100, 107, "adsorb al pre-exp, mu solid", 100.0, -14.5)
            ],
            desorption=[
                RxnBlock(
                    from_state=107,
                    comment="desorb al",
                    pairs=[(1998.0, 12.0 * i) for i in range(4)],
                )
            ],
            provenance=PROV,
        )

    def test_header_and_block_structure(self):
        text = self._render()
        lines = text.splitlines()
        assert lines[0].startswith("8000")
        assert "301\t302\t# si-o-si hydrolysis" in text
        assert "100\t107\t# adsorb al pre-exp, mu solid" in text
        assert "107\t\t# desorb al" in text

    def test_pair_counts_match_declared(self):
        text = self._render()
        lines = text.splitlines()
        i = lines.index("301\t302\t# si-o-si hydrolysis")
        n = int(lines[i + 1])
        pairs = [
            tok
            for line in lines[i + 2 : i + 2 + (n + 9) // 10]
            for tok in line.split("\t")
        ]
        assert n == 15
        assert len(pairs) == 15

    def test_legacy_parser_shape(self):
        # The legacy reader consumes whitespace-separated (k, dE) floats;
        # every pair token must split into exactly two parseable floats.
        text = self._render()
        i = text.splitlines().index("107\t\t# desorb al") + 2
        for token in text.splitlines()[i].split("\t"):
            k, de = token.split()
            float(k), float(de)

    def test_desorption_with_to_state_rejected(self):
        with pytest.raises(ValueError, match="no to_state"):
            data_rxn(
                8000.0,
                -1.0,
                -1.0,
                hydrolysis=[],
                adsorption=[],
                desorption=[RxnBlock(107, "bad", [(1.0, 0.0)], to_state=108)],
            )

    def test_hydrolysis_without_to_state_rejected(self):
        with pytest.raises(ValueError, match="needs a to_state"):
            data_rxn(
                8000.0,
                -1.0,
                -1.0,
                hydrolysis=[RxnBlock(301, "bad", [(1.0, 0.0)])],
                adsorption=[],
                desorption=[],
            )
