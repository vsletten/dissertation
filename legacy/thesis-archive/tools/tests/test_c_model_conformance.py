from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest

MODULE_PATH = Path(__file__).parents[1] / "run_c_model_conformance.py"
SPEC = importlib.util.spec_from_file_location("c_model_conformance", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class CompareOutputTests(unittest.TestCase):
    def write(self, root: Path, name: str, text: str) -> Path:
        path = root / name
        path.write_text(text, encoding="ascii")
        return path

    def test_identical_output_has_byte_parity(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            expected = self.write(root, "expected.dat", "0.1,2,3\n0.2,4,5\n")
            actual = self.write(root, "actual.dat", "0.1,2,3\n0.2,4,5\n")
            result = MODULE.compare_numeric_output(expected, actual, expected_columns=3)
        self.assertEqual(result["classification"], "byte_parity")
        self.assertTrue(result["structural_match"])
        self.assertEqual(result["matching_prefix_rows"], 2)

    def test_first_row_sabotage_is_detected_as_numeric_divergence(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            expected = self.write(root, "expected.dat", "0.1,2,3\n0.2,4,5\n")
            sabotaged = self.write(root, "sabotaged.dat", "0.1,2,999\n0.2,4,5\n")
            result = MODULE.compare_numeric_output(expected, sabotaged, expected_columns=3)
        self.assertEqual(result["classification"], "numeric_divergence")
        self.assertEqual(result["first_mismatch_row"], 1)
        self.assertFalse(result["byte_equal"])

    def test_late_divergence_is_not_assumed_to_be_toolchain_drift(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            expected = self.write(root, "expected.dat", "0.1,2,3\n0.2,4,5\n")
            actual = self.write(root, "actual.dat", "0.1,2,3\n0.3,4,6\n")
            result = MODULE.compare_numeric_output(expected, actual, expected_columns=3)
        self.assertEqual(result["classification"], "numeric_divergence")
        self.assertEqual(result["matching_prefix_rows"], 1)
        self.assertEqual(result["first_mismatch_row"], 2)

    def test_schema_damage_is_structural_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            expected = self.write(root, "expected.dat", "0.1,2,3\n")
            actual = self.write(root, "actual.dat", "0.1,2\n")
            result = MODULE.compare_numeric_output(expected, actual, expected_columns=3)
        self.assertEqual(result["classification"], "behavioral_mismatch")
        self.assertFalse(result["structural_match"])


class HistoricalContractTests(unittest.TestCase):
    def test_diffusion_disable_audit_accepts_historical_branch(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            (root / "envrn.c").write_text(
                "int isActive(int rxn) { if (rxn < 24) return 1; "
                "else { /* diffusion */ result = FALSE; } return result; }\n",
                encoding="ascii",
            )
            (root / "rxnlist.h").write_text("#define NDES 24\n#define NRXN 28\n", encoding="ascii")
            (root / "evtlist.c").write_text("if (isActive(s, l, i)) schedule(i);\n", encoding="ascii")
            audit = MODULE.audit_diffusion_disabled(root)
        self.assertEqual(audit["status"], "pinned_disabled")
        self.assertEqual(audit["diffusion_ids"], [24, 25, 26, 27])

    def test_diffusion_disable_audit_fails_if_branch_is_missing(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            (root / "envrn.c").write_text("int isActive(int rxn) { return 1; }\n", encoding="ascii")
            (root / "rxnlist.h").write_text("#define NDES 24\n#define NRXN 28\n", encoding="ascii")
            audit = MODULE.audit_diffusion_disabled(root)
        self.assertEqual(audit["status"], "behavioral_mismatch")

    def test_diffusion_audit_ignores_unrelated_matching_text(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            (root / "envrn.c").write_text(
                "int isActive(int rxn) { return 1; }\n"
                "int unrelated(void) { /* diffusion */ result = FALSE; return result; }\n",
                encoding="ascii",
            )
            (root / "evtlist.c").write_text("if (isActive(s, l, i)) schedule(i);\n", encoding="ascii")
            (root / "rxnlist.h").write_text("#define NDES 24\n#define NRXN 28\n", encoding="ascii")
            audit = MODULE.audit_diffusion_disabled(root)
        self.assertEqual(audit["status"], "behavioral_mismatch")

    def test_sabotage_gate_perturbs_late_plausible_row_and_passes(self) -> None:
        gate = MODULE.run_sabotage_gate()
        self.assertTrue(gate["passed"])
        self.assertEqual(gate["observed_classification"], "numeric_divergence")
        self.assertFalse(gate["conformance_gate_passed"])

    def test_drift_candidates_never_pass_the_conformance_gate(self) -> None:
        runs = [
            {
                "fixture": fixture,
                "classification": "compiler_prng_drift_candidate",
                "returncode": 0,
                "timed_out": False,
            }
            for fixture in MODULE.GOLDEN_RUNS
        ]
        passed = MODULE.conformance_gate_passes(
            runs=runs,
            diffusion={"status": "pinned_disabled"},
            duplicate={"all_byte_equal": True},
            sabotage={"passed": True},
            source_matches=True,
            fixtures_match=True,
            setup_error=None,
        )
        self.assertFalse(passed)

    def test_historical_duplicate_rejects_equal_outputs_from_different_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            left, right = MODULE.CROSS_HOST_DUPLICATE
            for fixture in (left, right):
                dest = root / fixture
                dest.mkdir(parents=True)
                for name in MODULE.INPUT_NAMES:
                    (dest / name).write_bytes(b"shared-input\n")
                for name in MODULE.OUTPUT_COLUMNS:
                    (dest / name).write_bytes(b"shared-output\n")
            (root / right / MODULE.INPUT_NAMES[0]).write_bytes(b"different-input\n")
            result = MODULE.verify_historical_duplicate(root)
        self.assertFalse(result["all_byte_equal"])


if __name__ == "__main__":
    unittest.main()
