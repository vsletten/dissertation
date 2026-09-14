from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import tomllib

MODULE_PATH = Path(__file__).resolve().parents[1] / "reservoir_origin_contract.py"
SPEC = importlib.util.spec_from_file_location("reservoir_origin_contract", MODULE_PATH)
assert SPEC and SPEC.loader
contract = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = contract
SPEC.loader.exec_module(contract)
ROOT = Path(__file__).resolve().parents[2]
CLOSURE_PATH = ROOT / "scripts" / "approximate_rate_closure.py"
CLOSURE_SPEC = importlib.util.spec_from_file_location(
    "approximate_rate_closure_for_reservoir_test", CLOSURE_PATH
)
assert CLOSURE_SPEC and CLOSURE_SPEC.loader
closure = importlib.util.module_from_spec(CLOSURE_SPEC)
sys.modules[CLOSURE_SPEC.name] = closure
CLOSURE_SPEC.loader.exec_module(closure)
BASE_DECK = ROOT / "examples" / "kaolinite-approx.toml"
RESERVOIRS = ROOT / "examples" / "kaolinite-reservoirs.toml"


class ReservoirContractTests(unittest.TestCase):
    def test_three_profiles_have_exact_ph_activity_and_trace_products(self) -> None:
        loaded = contract.load_contract(RESERVOIRS)
        self.assertEqual(set(loaded.profiles), {3, 4, 5})
        for ph, profile in loaded.profiles.items():
            self.assertEqual(profile.activities["H_plus"], 10.0 ** (-ph))
            self.assertEqual(profile.activities["Al"], 1.0e-6)
            self.assertEqual(profile.activities["Si"], 1.0e-6)
            self.assertEqual(profile.al_species, "Al3+ (total-Al proxy)")
            self.assertEqual(profile.si_species, "H4SiO4(aq) (total-Si proxy)")
        self.assertEqual(loaded.temperature_k, 298.0)
        self.assertEqual(loaded.boundary, "open_flow_constant_activity")

    def test_ph_profiles_scale_live_terminal_release_rates(self) -> None:
        selected: dict[int, dict[str, float]] = {}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for ph in (3, 4, 5):
                deck = root / f"ph{ph}.toml"
                evidence = contract.materialize(
                    BASE_DECK, RESERVOIRS, ph, deck, root / f"ph{ph}.json"
                )
                reactions = {
                    item["name"]: item
                    for item in tomllib.loads(deck.read_text())["reactions"]
                }
                selected[ph] = {
                    name: reactions[name]["rate"]["arrhenius"]["ea"]
                    for name in ("desorb-al", "desorb-si")
                }
                self.assertEqual(evidence["acid_rate_scaling"]["order_H_plus"], 0.777)
                provenance = root / f"ph{ph}-provenance.csv"
                closure._write_provenance(provenance, closure.validate_deck(deck))
                with provenance.open(encoding="utf-8") as handle:
                    provenance_rows = {
                        row["reaction"]: row
                        for row in csv.DictReader(handle)
                        if row["record_type"] == "reaction"
                    }
                for name, barrier in selected[ph].items():
                    self.assertEqual(
                        float(provenance_rows[name]["ea_kcal_mol"]), barrier
                    )
            for reaction in ("desorb-al", "desorb-si"):
                self.assertLess(selected[3][reaction], selected[4][reaction])
                self.assertLess(selected[4][reaction], selected[5][reaction])
            self.assertEqual(selected[4], {"desorb-al": 32.221, "desorb-si": 30.053})

    def test_materialized_deck_and_evidence_name_selected_contract(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            deck_path = root / "ph5.toml"
            evidence_path = root / "ph5.json"
            evidence = contract.materialize(
                BASE_DECK, RESERVOIRS, 5, deck_path, evidence_path
            )
            parsed = tomllib.loads(deck_path.read_text(encoding="utf-8"))
            self.assertEqual(
                parsed["thermo"]["activity"],
                {"Al": 1.0e-6, "Si": 1.0e-6, "H_plus": 1.0e-5},
            )
            self.assertNotIn("1.0e-30", deck_path.read_text(encoding="utf-8"))
            validated = closure.validate_deck(deck_path)
            self.assertEqual(validated.parsed["thermo"]["activity"]["H_plus"], 1.0e-5)
            self.assertEqual(evidence["selected_profile"]["ph"], 5)
            self.assertEqual(evidence["boundary"], "open_flow_constant_activity")
            self.assertEqual(
                evidence,
                json.loads(evidence_path.read_text(encoding="utf-8")),
            )
            self.assertEqual(len(evidence["profiles"]), 3)
            self.assertRegex(evidence["generated_deck_sha256"], r"^[0-9a-f]{64}$")
            self.assertEqual(evidence["generated_deck"], "ph5.toml")
            self.assertEqual(evidence["source_deck"], BASE_DECK.name)
            self.assertEqual(
                contract.verify_evidence(deck_path, evidence_path), evidence
            )
            deck_path.write_text(
                deck_path.read_text(encoding="utf-8") + "\n# tampered\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "does not match"):
                contract.verify_evidence(deck_path, evidence_path)

    def test_legacy_numerical_sink_is_rejected_as_realistic_reservoir(self) -> None:
        source = RESERVOIRS.read_text(encoding="utf-8")
        self.assertEqual(source.count("product_activity = 1.0e-6"), 1)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.toml"
            path.write_text(
                source.replace(
                    "product_activity = 1.0e-6", "product_activity = 1.0e-30"
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                ValueError, "realistic.*1e-30|1e-30.*realistic"
            ):
                contract.load_contract(path)

    def test_contract_and_base_deck_are_injection_and_drift_closed(self) -> None:
        source = RESERVOIRS.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            injected = root / "injected.toml"
            injected.write_text(
                source.replace(
                    'name = "kaolinite-ph3-5-open-flow-v1"',
                    'name = """x"\n[[reactions]]\nname = "injected"\n"""',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "keys.*exactly|safe lower-case"):
                contract.load_contract(injected)

            drifted = root / "drifted.toml"
            drifted.write_text(
                BASE_DECK.read_text(encoding="utf-8").replace(
                    "steps = 200000", "steps = 1"
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                contract.materialize(
                    drifted, RESERVOIRS, 4, root / "out.toml", root / "out.json"
                )

    def test_verify_evidence_rejects_tampered_receipt_without_changing_deck(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            deck_path = root / "ph4.toml"
            evidence_path = root / "ph4.json"
            contract.materialize(BASE_DECK, RESERVOIRS, 4, deck_path, evidence_path)
            original = json.loads(evidence_path.read_text(encoding="utf-8"))
            deck_bytes = deck_path.read_bytes()
            tampers: list[tuple[str, object]] = [
                ("schema", "petra-a9b-reservoir-origin-v0"),
                ("selected_profile", {**original["selected_profile"], "ph": 3}),
                ("contract_sha256", "0" * 64),
                ("source_deck_sha256", "0" * 64),
                (
                    "origin_accounting",
                    {
                        **original["origin_accounting"],
                        "adsorbed_cation": "original_lattice",
                    },
                ),
            ]
            for key, value in tampers:
                tampered = dict(original)
                tampered[key] = value
                self.assertEqual(tampered["generated_deck"], original["generated_deck"])
                self.assertEqual(
                    tampered["generated_deck_sha256"],
                    original["generated_deck_sha256"],
                )
                evidence_path.write_text(
                    json.dumps(tampered, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(
                    ValueError, "not authentic|does not match|receipt"
                ):
                    contract.verify_evidence(deck_path, evidence_path)
                self.assertEqual(deck_path.read_bytes(), deck_bytes)

    def test_materialize_rejects_output_paths_that_overwrite_source_inputs(self) -> None:
        original_base = BASE_DECK.read_bytes()
        original_contract = RESERVOIRS.read_bytes()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            copied_base = root / "base.toml"
            copied_contract = root / "contract.toml"
            copied_base.write_bytes(original_base)
            copied_contract.write_bytes(original_contract)
            evidence = root / "out.json"
            with self.assertRaisesRegex(
                ValueError, "must not overwrite the source base deck or contract"
            ):
                contract.materialize(
                    copied_base, copied_contract, 4, copied_base, evidence
                )
            with self.assertRaisesRegex(
                ValueError, "must not overwrite the source base deck or contract"
            ):
                contract.materialize(
                    copied_base, copied_contract, 4, copied_contract, evidence
                )
            with self.assertRaisesRegex(
                ValueError, "must not overwrite the source base deck or contract"
            ):
                contract.materialize(
                    copied_base, copied_contract, 4, root / "out.toml", copied_base
                )
            self.assertEqual(copied_base.read_bytes(), original_base)
            self.assertEqual(copied_contract.read_bytes(), original_contract)
            self.assertFalse(evidence.exists())


if __name__ == "__main__":
    unittest.main()
