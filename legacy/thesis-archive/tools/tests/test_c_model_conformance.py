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

    def test_sabotage_is_detected_as_behavioral_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            expected = self.write(root, "expected.dat", "0.1,2,3\n0.2,4,5\n")
            sabotaged = self.write(root, "sabotaged.dat", "0.1,2,999\n0.2,4,5\n")
            result = MODULE.compare_numeric_output(expected, sabotaged, expected_columns=3)
        self.assertEqual(result["classification"], "behavioral_mismatch")
        self.assertEqual(result["first_mismatch_row"], 1)
        self.assertFalse(result["byte_equal"])

    def test_late_divergence_is_classified_as_compiler_prng_drift(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            expected = self.write(root, "expected.dat", "0.1,2,3\n0.2,4,5\n")
            actual = self.write(root, "actual.dat", "0.1,2,3\n0.3,4,6\n")
            result = MODULE.compare_numeric_output(expected, actual, expected_columns=3)
        self.assertEqual(result["classification"], "compiler_prng_drift")
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


if __name__ == "__main__":
    unittest.main()
