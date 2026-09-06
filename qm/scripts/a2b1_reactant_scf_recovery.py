#!/usr/bin/env python3
"""One bounded, evidence-bound SCF recovery for A2b's exact reactant.

The route is intentionally different from the exhausted Newton -> direct-DIIS
sequence: a Hückel guess enters a damped, level-shifted Roothaan stabilization,
then a fresh unshifted CDIIS solve starts from that density.  Only the second
stage may publish a canonical energy.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import signal
import sys
import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

if __name__ == "__main__":
    if "--finalize-if-running" not in sys.argv:
        from quarry.etiquette import bootstrap_cli

        _ETIQUETTE = bootstrap_cli(
            "a2b1-reactant-scf-recovery",
            default_run_root=(
                "/mnt/data/vsletten/dissertation-data/"
                "task294-a2b1-wb97mv-reactant-scf-recovery-20260906"
            ),
            gpu_owner="a2b1-reactant-scf-recovery",
        )

import numpy as np

from quarry import pipeline
from quarry.clusters import Cluster
from quarry.pipeline import DftSettings
from scripts import a2b_al_neutral_production as a2b
from scripts import production_energetics as a2

RECOVERY_VERSION = "a2b1-reactant-scf-recovery-v1"
RECOVERY_CONTRACT = "a2b1-huckel-damped-level-shifted-roothaan-to-cdiis-v1"
LEGACY_RUN_ROOT = Path(
    "/mnt/data/vsletten/dissertation-data/task290-a2b-al-neutral-production-20260906"
)
DEFAULT_RUN_ROOT = Path(
    "/mnt/data/vsletten/dissertation-data/"
    "task294-a2b1-wb97mv-reactant-scf-recovery-20260906"
)
LEGACY_TERMINAL_SHA256 = (
    "c4233567cc3a88a23683514531e89af43c1fe6743a020638bc4fb20281d2fecd"
)
STABILIZATION_CYCLES = 50
FINALIZATION_CYCLES = 150
DRIVER_WALL_SECONDS = 7200
SYSTEMD_RUNTIME_MAX_SECONDS = 7800
LEVEL_SHIFT_HARTREE = 0.5
DAMPING_FACTOR = 0.5


@dataclass(frozen=True)
class InputEvidence:
    cluster: Cluster
    settings: DftSettings
    bindings: dict[str, Any]


class RecoveryTimeout(RuntimeError):
    """Raised when the driver-level hard wall expires."""


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path}: expected a JSON object")
    return payload


def _host_array(value: Any) -> np.ndarray:
    if hasattr(value, "get"):
        value = value.get()
    return np.ascontiguousarray(np.asarray(value))


def array_receipt(value: Any) -> dict[str, Any]:
    array = _host_array(value)
    digest = hashlib.sha256()
    digest.update(array.dtype.str.encode())
    digest.update(json.dumps(list(array.shape), separators=(",", ":")).encode())
    digest.update(array.tobytes(order="C"))
    return {
        "sha256": digest.hexdigest(),
        "dtype": array.dtype.str,
        "shape": list(array.shape),
    }


def _finite_float(value: Any) -> float | None:
    if value is None:
        return None
    if hasattr(value, "get"):
        value = value.get()
    result = float(np.asarray(value).item())
    return result if math.isfinite(result) else None


def _checkpoint_identity(
    source: a2b.SourceEvidence, r2scan: DftSettings
) -> dict[str, Any]:
    return {
        "a2b_version": a2b.A2B_VERSION,
        "stage": "minimum",
        "role": "reactant",
        "algorithm": "ase-bfgs-v1",
        "settings": a2.frequency_settings_fingerprint(r2scan),
        "source_sha256": source.contract.artifact_sha256,
        "max_steps": 200,
    }


def validate_inputs(legacy_root: Path, *, use_gpu: bool) -> InputEvidence:
    legacy_root = legacy_root.resolve()
    terminal_path = legacy_root / "terminal-receipt.json"
    if a2.sha256_path(terminal_path) != LEGACY_TERMINAL_SHA256:
        raise ValueError("immutable TASK-290 terminal receipt hash drift")
    terminal = _read_json(terminal_path)
    if terminal.get("outcome") != "incomplete-computational-failure":
        raise ValueError("TASK-290 terminal outcome is not the accepted finite failure")
    if terminal.get("running_record_present") is not False:
        raise ValueError("TASK-290 terminal receipt still reports a running record")
    for forbidden in ("results.json", "store.sqlite", "reactant.wb97m-v.energy.json"):
        if (legacy_root / forbidden).exists():
            raise ValueError(f"TASK-290 canonical boundary drift: {forbidden} exists")

    source = a2b.load_source_evidence(a2b.SOURCE_ROOT, attacker_index=15)
    r2scan, production, _b3lyp = a2.settings(use_gpu=use_gpu)
    contract = a2b.method_contract(use_gpu=use_gpu)
    settings_path = legacy_root / "settings.json"
    settings_payload = _read_json(settings_path)
    for key in ("a2_method_contract", "methods", "settings_fingerprints"):
        if settings_payload.get(key) != contract.get(key):
            raise ValueError(f"TASK-290 settings drift for {key}")

    geometry_path = legacy_root / "reactant.r2scan3c.xyz"
    cluster = a2.checkpoint_cluster(
        geometry_path,
        source.clusters["reactant"],
        lambda: (_ for _ in ()).throw(RuntimeError("reactant checkpoint missing")),
        identity=_checkpoint_identity(source, r2scan),
    )
    if (
        a2.si_neutral_signature(cluster, 15)
        != source.contract.expected_basins["reactant"]
    ):
        raise ValueError("reactant basin identity drift")
    if a2.endpoint_identity(cluster, 15) != a2.endpoint_identity(
        source.clusters["reactant"], 15
    ):
        raise ValueError("reactant atom/state/ownership identity drift")

    frequency_path = legacy_root / "reactant.r2scan3c.frequency.json"
    frequency = _read_json(frequency_path)
    if frequency.get("geometry_fingerprint") != a2.frequency_geometry_fingerprint(
        cluster
    ):
        raise ValueError("reactant frequency geometry drift")
    if frequency.get("settings_fingerprint") != a2.frequency_settings_fingerprint(
        r2scan
    ):
        raise ValueError("reactant frequency settings drift")
    if frequency.get("hessian_method") != "finite-difference-gradient":
        raise ValueError("reactant Hessian provenance drift")
    if frequency.get("imaginary_cm") != []:
        raise ValueError("reactant is no longer a validated r2SCAN-3c minimum")

    source_receipt_path = legacy_root / "source-receipt.json"
    source_receipt = _read_json(source_receipt_path)
    if source_receipt.get("driver_sha256") != terminal.get("driver_sha256"):
        raise ValueError("TASK-290 source/terminal driver binding drift")

    driver_path = Path(__file__).resolve()
    bindings = {
        "recovery_contract": RECOVERY_CONTRACT,
        "legacy_terminal": {
            "path": str(terminal_path),
            "sha256": LEGACY_TERMINAL_SHA256,
        },
        "legacy_source_receipt": {
            "path": str(source_receipt_path),
            "sha256": a2.sha256_path(source_receipt_path),
        },
        "legacy_settings": {
            "path": str(settings_path),
            "sha256": a2.sha256_path(settings_path),
        },
        "reactant_geometry": {
            "path": str(geometry_path),
            "sha256": a2.sha256_path(geometry_path),
            "fingerprint": a2.frequency_geometry_fingerprint(cluster),
            "charge": cluster.charge,
            "spin": cluster.spin,
        },
        "reactant_frequency": {
            "path": str(frequency_path),
            "sha256": a2.sha256_path(frequency_path),
            "imaginary_count": 0,
        },
        "production_method": a2.PRODUCTION_METHOD,
        "production_settings_fingerprint": a2.frequency_settings_fingerprint(
            production
        ),
        "driver": str(driver_path),
        "driver_sha256": a2.sha256_path(driver_path),
    }
    return InputEvidence(cluster=cluster, settings=production, bindings=bindings)


def _quarantine(paths: list[Path], run_dir: Path) -> str | None:
    existing = [path for path in paths if path.exists()]
    if not existing:
        return None
    target = run_dir / "quarantine" / f"stale-{time.time_ns()}"
    target.mkdir(parents=True, exist_ok=False)
    for index, path in enumerate(existing):
        destination = target / f"{index:02d}-{path.name}"
        shutil.move(str(path), destination)
    return str(target)


def _cycle_callback(
    receipt_path: Path, base: dict[str, Any], history: list[dict[str, Any]]
) -> Callable[[dict[str, Any]], None]:
    def callback(env: dict[str, Any]) -> None:
        record = {
            "cycle": int(env.get("cycle", len(history))),
            "electronic_hartree": _finite_float(env.get("e_tot")),
            "orbital_gradient_norm": _finite_float(env.get("norm_gorb")),
            "density_delta_norm": _finite_float(env.get("norm_ddm")),
        }
        history.append(record)
        a2.atomic_json(
            receipt_path,
            {
                **base,
                "status": "running",
                "cycles_completed": len(history),
                "cycle_diagnostics": history,
            },
        )

    return callback


def run_attempt(
    *,
    receipt_path: Path,
    evidence: InputEvidence,
    attempt_id: str,
    solver: dict[str, Any],
    dm0: Any | None,
) -> tuple[Any, dict[str, Any]]:
    mol = pipeline.build_mol(evidence.cluster, evidence.settings)
    mf = pipeline._make_scf(mol, evidence.settings)
    mf.max_cycle = int(solver["max_cycle"])
    mf.diis = bool(solver["diis"])
    mf.diis_start_cycle = int(solver["diis_start_cycle"])
    mf.damp = float(solver["damp"])
    mf.level_shift = float(solver["level_shift"])
    if dm0 is None:
        dm0 = mf.get_init_guess(key=str(solver["initial_guess"]))
        density_provenance = {
            "kind": str(solver["initial_guess"]),
            **array_receipt(dm0),
        }
    else:
        density_provenance = {
            "kind": "prior-attempt-final-density",
            **array_receipt(dm0),
        }

    base = {
        "recovery_version": RECOVERY_VERSION,
        "recovery_contract": RECOVERY_CONTRACT,
        "attempt_id": attempt_id,
        "bindings": evidence.bindings,
        "solver": solver,
        "initial_density": density_provenance,
    }
    history: list[dict[str, Any]] = []
    mf.callback = _cycle_callback(receipt_path, base, history)
    a2.atomic_json(
        receipt_path,
        {**base, "status": "starting", "cycles_completed": 0, "cycle_diagnostics": []},
    )
    raw_value: Any = None
    error: BaseException | None = None
    try:
        raw_value = mf.kernel(dm0=dm0)
    except BaseException as exc:
        error = exc
    final_density: Any | None = None
    try:
        final_density = mf.make_rdm1()
    except BaseException:
        final_density = None
    converged = bool(getattr(mf, "converged", False)) and error is None
    payload: dict[str, Any] = {
        **base,
        "status": "converged" if converged else "failed",
        "converged": converged,
        "cycles_completed": int(getattr(mf, "cycles", len(history))),
        "cycle_diagnostics": history,
        "electronic_hartree_diagnostic_only": _finite_float(raw_value),
        "final_density": array_receipt(final_density)
        if final_density is not None
        else None,
    }
    if error is not None:
        payload.update(error_type=type(error).__name__, error=str(error))
    a2.atomic_json(receipt_path, payload)
    if error is not None:
        raise error
    if final_density is None:
        raise RuntimeError(f"{attempt_id}: no final density was available")
    return final_density, payload


def _attempt_ref(path: Path) -> dict[str, str]:
    return {"path": str(path.resolve()), "sha256": a2.sha256_path(path)}


def _terminal_payload(
    run_dir: Path,
    evidence: InputEvidence | None,
    *,
    outcome: str,
    attempts: list[Path],
    error: BaseException | None = None,
    canonical_energy: Path | None = None,
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "recovery_version": RECOVERY_VERSION,
        "recovery_contract": RECOVERY_CONTRACT,
        "outcome": outcome,
        "finished_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "running_record_present": False,
        "driver_sha256": a2.sha256_path(Path(__file__).resolve()),
        "attempts": [_attempt_ref(path) for path in attempts if path.exists()],
    }
    if evidence is not None:
        payload["bindings"] = evidence.bindings
    if canonical_energy is not None and canonical_energy.exists():
        payload["canonical_energy"] = _attempt_ref(canonical_energy)
    if error is not None:
        payload.update(error_type=type(error).__name__, error=str(error))
    return payload


def _install_alarm(seconds: int) -> None:
    def expired(_signum: int, _frame: Any) -> None:
        raise RecoveryTimeout(f"recovery wall clock exhausted after {seconds} seconds")

    signal.signal(signal.SIGALRM, expired)
    signal.alarm(seconds)


def execute(args: argparse.Namespace) -> int:
    run_dir = args.run_dir.resolve()
    target_energy = args.target_energy.resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    terminal_path = run_dir / "terminal-receipt.json"
    forbidden_outputs = [run_dir / "results.json", run_dir / "store.sqlite"]
    if terminal_path.exists():
        terminal = _read_json(terminal_path)
        if terminal.get("outcome") == "success":
            a2.checkpoint_converged_energy(
                target_energy,
                validate_inputs(args.legacy_run_root, use_gpu=args.gpu).cluster,
                a2.settings(use_gpu=args.gpu)[1],
                a2.PRODUCTION_METHOD,
            )
            _quarantine(forbidden_outputs, run_dir)
            return 0
        _quarantine([target_energy, *forbidden_outputs], run_dir)
        raise RuntimeError(
            "A2b1 already has a terminal failure; no retry is authorized"
        )

    quarantine = _quarantine([target_energy, *forbidden_outputs], run_dir)
    evidence: InputEvidence | None = None
    attempt_paths = [
        run_dir / "attempt-001-huckel-roothaan.json",
        run_dir / "attempt-002-cdiis-finalization.json",
    ]
    running = {
        "recovery_version": RECOVERY_VERSION,
        "recovery_contract": RECOVERY_CONTRACT,
        "status": "running",
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "quarantine": quarantine,
        "driver_wall_seconds": int(args.wall_seconds),
    }
    a2.atomic_json(run_dir / "run_status.json", running)
    _install_alarm(int(args.wall_seconds))
    try:
        evidence = validate_inputs(args.legacy_run_root, use_gpu=args.gpu)
        stage1_solver = {
            "algorithm": "damped-level-shifted-roothaan",
            "initial_guess": "huckel",
            "max_cycle": STABILIZATION_CYCLES,
            "diis": False,
            "diis_start_cycle": STABILIZATION_CYCLES + 1,
            "damp": DAMPING_FACTOR,
            "level_shift": LEVEL_SHIFT_HARTREE,
        }
        density, stage1 = run_attempt(
            receipt_path=attempt_paths[0],
            evidence=evidence,
            attempt_id="attempt-001",
            solver=stage1_solver,
            dm0=None,
        )
        stage2_solver = {
            "algorithm": "unshifted-cdiis-finalization",
            "initial_guess": "prior-attempt-final-density",
            "max_cycle": FINALIZATION_CYCLES,
            "diis": True,
            "diis_start_cycle": 1,
            "damp": 0.0,
            "level_shift": 0.0,
        }
        _density, stage2 = run_attempt(
            receipt_path=attempt_paths[1],
            evidence=evidence,
            attempt_id="attempt-002",
            solver=stage2_solver,
            dm0=density,
        )
        if not stage2["converged"]:
            raise RuntimeError("unshifted CDIIS finalization did not converge")
        value = stage2["electronic_hartree_diagnostic_only"]
        if value is None:
            raise RuntimeError("converged finalization returned a non-finite energy")
        a2.atomic_json(
            target_energy,
            {
                "method": a2.PRODUCTION_METHOD,
                "electronic_hartree": value,
                "geometry_fingerprint": evidence.bindings["reactant_geometry"][
                    "fingerprint"
                ],
                "settings_fingerprint": evidence.bindings[
                    "production_settings_fingerprint"
                ],
                "scf_contract": RECOVERY_CONTRACT,
                "convergence_route": "huckel-roothaan-then-cdiis",
                "scf_attempts": [
                    {
                        "solver": stage1_solver["algorithm"],
                        "max_cycle": STABILIZATION_CYCLES,
                        "converged": bool(stage1["converged"]),
                    },
                    {
                        "solver": stage2_solver["algorithm"],
                        "max_cycle": FINALIZATION_CYCLES,
                        "converged": True,
                    },
                ],
                "recovery_attempt_receipts": [
                    _attempt_ref(attempt_paths[0]),
                    _attempt_ref(attempt_paths[1]),
                ],
                "recovery_driver_sha256": evidence.bindings["driver_sha256"],
                "converged": True,
            },
        )
        terminal = _terminal_payload(
            run_dir,
            evidence,
            outcome="success",
            attempts=attempt_paths,
            canonical_energy=target_energy,
        )
        completed = {
            **running,
            "status": "completed",
            "outcome": "success",
            "finished_at": terminal["finished_at"],
            "terminal_receipt_pending": True,
        }
        a2.atomic_json(run_dir / "run_status.json", completed)
        terminal["run_status_sha256"] = a2.sha256_path(run_dir / "run_status.json")
        a2.atomic_json(terminal_path, terminal)
        return 0
    except BaseException as exc:
        _quarantine(
            [target_energy, run_dir / "results.json", run_dir / "store.sqlite"], run_dir
        )
        failed = {
            **running,
            "status": "failed",
            "outcome": "incomplete-computational-failure",
            "finished_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "error_type": type(exc).__name__,
            "error": str(exc),
        }
        a2.atomic_json(run_dir / "run_status.json", failed)
        terminal = _terminal_payload(
            run_dir,
            evidence,
            outcome="incomplete-computational-failure",
            attempts=attempt_paths,
            error=exc,
        )
        terminal["run_status_sha256"] = a2.sha256_path(run_dir / "run_status.json")
        a2.atomic_json(terminal_path, terminal)
        raise
    finally:
        signal.alarm(0)


def finalize_if_running(args: argparse.Namespace) -> int:
    run_dir = args.run_dir.resolve()
    target_energy = args.target_energy.resolve()
    terminal_path = run_dir / "terminal-receipt.json"
    forbidden_outputs = [run_dir / "results.json", run_dir / "store.sqlite"]
    if terminal_path.exists():
        terminal = _read_json(terminal_path)
        cleanup = forbidden_outputs
        if terminal.get("outcome") != "success":
            cleanup = [target_energy, *cleanup]
        _quarantine(cleanup, run_dir)
        return 0
    status_path = run_dir / "run_status.json"
    prior_status = _read_json(status_path) if status_path.exists() else {}
    error = RecoveryTimeout(
        "bounded systemd unit ended before an atomic terminal receipt"
    )
    _quarantine(
        [target_energy, run_dir / "results.json", run_dir / "store.sqlite"], run_dir
    )
    failed = {
        **prior_status,
        "recovery_version": RECOVERY_VERSION,
        "recovery_contract": RECOVERY_CONTRACT,
        "status": "failed",
        "outcome": "incomplete-computational-failure",
        "finished_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "error_type": type(error).__name__,
        "error": str(error),
    }
    a2.atomic_json(status_path, failed)
    attempts = [
        run_dir / "attempt-001-huckel-roothaan.json",
        run_dir / "attempt-002-cdiis-finalization.json",
    ]
    terminal = _terminal_payload(
        run_dir,
        None,
        outcome="incomplete-computational-failure",
        attempts=attempts,
        error=error,
    )
    terminal["run_status_sha256"] = a2.sha256_path(status_path)
    a2.atomic_json(terminal_path, terminal)
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_ROOT)
    result.add_argument("--legacy-run-root", type=Path, default=LEGACY_RUN_ROOT)
    result.add_argument(
        "--target-energy",
        type=Path,
        default=LEGACY_RUN_ROOT / "reactant.wb97m-v.energy.json",
    )
    result.add_argument("--gpu", action="store_true")
    result.add_argument("--gpu-mem-gb", type=float, default=18.0)
    result.add_argument("--threads", type=int, default=16)
    result.add_argument("--nice", type=int, default=10)
    result.add_argument("--wall-seconds", type=int, default=DRIVER_WALL_SECONDS)
    result.add_argument("--finalize-if-running", action="store_true")
    result.add_argument("--log", type=Path)
    return result


def main() -> int:
    args = parser().parse_args()
    if args.threads < 1 or args.threads > 16:
        raise ValueError("threads must be in [1, 16]")
    if args.wall_seconds != DRIVER_WALL_SECONDS:
        raise ValueError(
            f"driver wall is fixed at {DRIVER_WALL_SECONDS} seconds by contract"
        )
    if args.gpu_mem_gb <= 0:
        raise ValueError("gpu-mem-gb must be positive")
    if args.finalize_if_running:
        return finalize_if_running(args)
    return execute(args)


if __name__ == "__main__":
    raise SystemExit(main())
