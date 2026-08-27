import ast
import json
import os
import re
import signal
import stat
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path
from types import SimpleNamespace

import pytest

from quarry.etiquette import (
    THREAD_ENV_VARS,
    GpuLeaseBusy,
    acquire_gpu,
    apply_cupy_memory_limit,
    bootstrap_cli,
    enforce,
    gpu_lane_available,
)


def test_enforce_sets_unset_clamps_high_and_preserves_low(monkeypatch):
    monkeypatch.delenv("OMP_NUM_THREADS", raising=False)
    monkeypatch.setenv("MKL_NUM_THREADS", "64")
    monkeypatch.setenv("OPENBLAS_NUM_THREADS", "4")
    nice_calls = []
    monkeypatch.setattr(os, "nice", nice_calls.append)

    enforce(threads=16, niceness=10)

    assert {name: os.environ[name] for name in THREAD_ENV_VARS} == {
        "OMP_NUM_THREADS": "16",
        "MKL_NUM_THREADS": "16",
        "OPENBLAS_NUM_THREADS": "4",
    }
    assert nice_calls == [10]


def test_negative_niceness_is_rejected_for_non_root(monkeypatch):
    monkeypatch.setattr(os, "geteuid", lambda: 1000)
    monkeypatch.setattr(
        os, "nice", lambda value: pytest.fail("os.nice must not be called")
    )

    with pytest.raises(ValueError, match="negative niceness requires root"):
        enforce(niceness=-1)


def test_niceness_failure_warns_and_continues(monkeypatch):
    monkeypatch.setattr(
        os, "nice", lambda value: (_ for _ in ()).throw(OSError("not permitted"))
    )

    with pytest.warns(UserWarning, match="continuing"):
        enforce(niceness=10)


def test_bootstrap_cli_applies_overrides_and_self_tees(monkeypatch, tmp_path):
    for name in THREAD_ENV_VARS:
        monkeypatch.setenv(name, "64")
    nice_calls = []
    monkeypatch.setattr(os, "nice", nice_calls.append)
    original_stdout, original_stderr = sys.stdout, sys.stderr

    session = bootstrap_cli(
        "driver",
        default_run_root=tmp_path,
        argv=["--threads", "32", "--nice", "5"],
    )
    try:
        print("durable campaign output")
        sys.stdout.flush()
    finally:
        session.close()

    assert session.threads == 32
    assert session.niceness == 5
    assert session.log_path.parent == tmp_path
    log_text = session.log_path.read_text()
    assert "durable campaign output" in log_text
    assert "OMP_NUM_THREADS=32" in log_text
    assert "MKL_NUM_THREADS=32" in log_text
    assert "OPENBLAS_NUM_THREADS=32" in log_text
    assert {os.environ[name] for name in THREAD_ENV_VARS} == {"32"}
    assert nice_calls == [5]
    assert sys.stdout is original_stdout
    assert sys.stderr is original_stderr


def test_all_campaign_entrypoints_bootstrap_before_heavy_imports():
    scripts = Path(__file__).resolve().parent.parent / "scripts"
    for path in sorted(scripts.glob("*.py")):
        source = path.read_text()
        tree = ast.parse(source)
        bootstrap_at = source.find("bootstrap_cli(")
        bootstrap_calls = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "bootstrap_cli"
        ]
        main_guards = [
            node
            for node in tree.body
            if isinstance(node, ast.If)
            and isinstance(node.test, ast.Compare)
            and isinstance(node.test.left, ast.Name)
            and node.test.left.id == "__name__"
            and len(node.test.ops) == 1
            and isinstance(node.test.ops[0], ast.Eq)
            and len(node.test.comparators) == 1
            and isinstance(node.test.comparators[0], ast.Constant)
            and node.test.comparators[0].value == "__main__"
        ]
        guarded_nodes = {
            nested
            for guard in main_guards
            for statement in guard.body
            for nested in ast.walk(statement)
        }
        heavy_imports = [
            match.start()
            for match in re.finditer(
                r"^(?:import numpy|from (?!quarry\.etiquette\b)quarry\.)",
                source,
                re.MULTILINE,
            )
        ]
        assert bootstrap_at >= 0, f"{path.name} does not bootstrap etiquette"
        assert bootstrap_calls, f"{path.name} does not call bootstrap_cli"
        assert all(call in guarded_nodes for call in bootstrap_calls), (
            f"{path.name} calls bootstrap_cli outside its __main__ guard"
        )
        assert heavy_imports, f"{path.name} has no recognized heavy import"
        assert bootstrap_at < min(heavy_imports), (
            f"{path.name} bootstraps after a heavy import"
        )
        assert '"--threads"' in source
        assert '"--nice"' in source
        assert '"--log"' in source


def test_gpu_lease_acquires_mode_0600_and_releases_exact_record(tmp_path):
    lease_path = tmp_path / "gpu" / "lease.json"
    started = datetime(2026, 8, 27, 12, 0, tzinfo=UTC)
    previous_term = signal.getsignal(signal.SIGTERM)

    lease = acquire_gpu(
        "phase1",
        16.0,
        12.0,
        lease_path=lease_path,
        pid=4242,
        now=started,
        pid_alive=lambda pid: False,
    )
    try:
        payload = json.loads(lease_path.read_text())
        assert payload == {
            "expected_gb": 16.0,
            "owner": "phase1",
            "pid": 4242,
            "started": "2026-08-27T12:00:00Z",
            "ttl": 12.0,
        }
        assert stat.S_IMODE(lease_path.stat().st_mode) == 0o600
        handler = signal.getsignal(signal.SIGTERM)
        assert getattr(handler, "__self__", None) is lease
    finally:
        lease.release()

    assert not lease_path.exists()
    assert signal.getsignal(signal.SIGTERM) == previous_term


def test_gpu_lease_refuses_a_live_owner_with_clear_message(tmp_path):
    lease_path = tmp_path / "lease.json"
    first = acquire_gpu(
        "phase1",
        16.0,
        12.0,
        lease_path=lease_path,
        pid=111,
        pid_alive=lambda pid: False,
        install_signal_handlers=False,
    )
    try:
        with pytest.raises(GpuLeaseBusy, match=r"GPU lane busy \(owner=phase1 pid=111"):
            acquire_gpu(
                "phase2",
                16.0,
                12.0,
                lease_path=lease_path,
                pid=222,
                pid_alive=lambda pid: pid == 111,
                install_signal_handlers=False,
            )
    finally:
        first.release()


def test_gpu_lease_breaks_dead_pid_and_writes_dated_note(tmp_path):
    lease_path = tmp_path / "lease.json"
    first = acquire_gpu(
        "dead-campaign",
        10.0,
        1.0,
        lease_path=lease_path,
        pid=111,
        now=datetime(2026, 8, 27, 10, 0, tzinfo=UTC),
        pid_alive=lambda pid: False,
        install_signal_handlers=False,
    )
    # Simulate a crashed process: leave its lease record behind.
    first._released = True

    second = acquire_gpu(
        "replacement",
        16.0,
        2.0,
        lease_path=lease_path,
        pid=222,
        now=datetime(2026, 8, 27, 12, 0, tzinfo=UTC),
        pid_alive=lambda pid: False,
        install_signal_handlers=False,
    )
    try:
        assert json.loads(lease_path.read_text())["owner"] == "replacement"
        note = (tmp_path / "stale-breaks.log").read_text()
        assert "2026-08-27T12:00:00Z" in note
        assert "owner=dead-campaign pid=111 reason=pid-dead" in note
    finally:
        second.release()


def test_gpu_lease_recovers_orphaned_pid_temp_and_acquires(tmp_path):
    lease_path = tmp_path / "lease.json"
    leftover = tmp_path / f".lease.json.{os.getpid()}.tmp"
    leftover.write_text("{}\n")
    leftover.chmod(0o600)

    lease = acquire_gpu(
        "replacement",
        16.0,
        1.0,
        lease_path=lease_path,
        install_signal_handlers=False,
    )
    try:
        assert json.loads(lease_path.read_text())["owner"] == "replacement"
        assert not leftover.exists()
    finally:
        lease.release()
    assert not lease_path.exists()


def test_gpu_lane_preflight_clears_only_a_dead_lease(tmp_path):
    lease_path = tmp_path / "lease.json"
    lease = acquire_gpu(
        "dead-campaign",
        16.0,
        1.0,
        lease_path=lease_path,
        pid=111,
        pid_alive=lambda pid: False,
        install_signal_handlers=False,
    )
    lease._released = True

    available, message = gpu_lane_available(
        lease_path=lease_path,
        pid_alive=lambda pid: False,
        now=datetime(2026, 8, 27, 12, 0, tzinfo=UTC),
    )

    assert available is True
    assert "broke stale owner=dead-campaign pid=111" in message
    assert not lease_path.exists()


def test_cupy_pool_limit_is_applied_and_preserves_ollama_floor(monkeypatch):
    limits = []
    pool = SimpleNamespace(set_limit=lambda *, size: limits.append(size))
    monkeypatch.setitem(
        sys.modules,
        "cupy",
        SimpleNamespace(get_default_memory_pool=lambda: pool),
    )

    assert apply_cupy_memory_limit(16.0) == 16 * 1024**3
    assert limits == [16 * 1024**3]
    with pytest.raises(ValueError, match="preserve the 6 GB ollama floor"):
        apply_cupy_memory_limit(18.1)


def test_two_gpu_driver_bootstraps_serialize_and_release(monkeypatch, tmp_path, capsys):
    lease_path = tmp_path / "lease.json"
    monkeypatch.setenv("GPU_LEASE_PATH", str(lease_path))
    pool = SimpleNamespace(set_limit=lambda *, size: None)
    monkeypatch.setitem(
        sys.modules,
        "cupy",
        SimpleNamespace(get_default_memory_pool=lambda: pool),
    )
    first = bootstrap_cli(
        "phase1",
        default_run_root=tmp_path,
        argv=["--gpu", "--gpu-mem-gb", "16", "--nice", "0"],
        gpu_owner="phase1_xiao_lasaga",
    )
    try:
        with pytest.raises(SystemExit) as exit_info:
            bootstrap_cli(
                "phase2",
                default_run_root=tmp_path,
                argv=["--gpu", "--gpu-mem-gb", "16", "--nice", "0"],
                gpu_owner="phase2_ladder",
            )
        assert exit_info.value.code == 1
        assert "GPU lane busy (owner=phase1_xiao_lasaga" in capsys.readouterr().err
        assert lease_path.exists()
    finally:
        first.close()
    assert not lease_path.exists()


def test_forked_child_cannot_release_parent_lease(tmp_path):
    lease_path = tmp_path / "lease.json"
    lease = acquire_gpu(
        "parent",
        16.0,
        1.0,
        lease_path=lease_path,
        install_signal_handlers=False,
    )
    child = os.fork()
    if child == 0:
        lease.release()
        os._exit(0)
    try:
        _, status = os.waitpid(child, 0)
        assert os.waitstatus_to_exitcode(status) == 0
        assert lease_path.exists()
        assert json.loads(lease_path.read_text())["pid"] == os.getpid()
    finally:
        lease.release()
    assert not lease_path.exists()


def test_actual_phase2_cli_refuses_live_lease_without_traceback(tmp_path):
    lease_path = tmp_path / "lease.json"
    lease = acquire_gpu(
        "parent",
        16.0,
        1.0,
        lease_path=lease_path,
        install_signal_handlers=False,
    )
    driver = Path(__file__).resolve().parent.parent / "scripts" / "phase2_ladder.py"
    env = os.environ | {"GPU_LEASE_PATH": str(lease_path)}
    try:
        completed = subprocess.run(
            [
                sys.executable,
                str(driver),
                "--family",
                "oss",
                "--dry-run",
                "--gpu",
                "--nice",
                "0",
                "--log",
                str(tmp_path / "driver.log"),
            ],
            env=env,
            text=True,
            capture_output=True,
            timeout=10,
            check=False,
        )
    finally:
        lease.release()

    assert completed.returncode == 1
    assert "GPU lane busy (owner=parent" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_sigterm_during_release_does_not_deadlock(tmp_path):
    lease_path = tmp_path / "lease.json"
    code = """
import os
import signal
import quarry.etiquette as etiquette

path = os.environ["GPU_LEASE_PATH"]
lease = etiquette.acquire_gpu("signal-test", 16.0, 1.0)
original = etiquette._read_lease

def signal_while_locked(current_path):
    os.kill(os.getpid(), signal.SIGTERM)
    return original(current_path)

etiquette._read_lease = signal_while_locked
lease.release()
raise RuntimeError("pending SIGTERM did not terminate the process")
"""
    completed = subprocess.run(
        [sys.executable, "-c", code],
        env=os.environ | {"GPU_LEASE_PATH": str(lease_path)},
        text=True,
        capture_output=True,
        timeout=5,
        check=False,
    )

    assert completed.returncode == -signal.SIGTERM
    assert not lease_path.exists()


def test_state_symlinks_fail_closed_without_touching_targets(tmp_path):
    lease_path = tmp_path / "lease.json"
    victim = tmp_path / "victim.txt"
    victim.write_text("KEEP\n")
    victim.chmod(0o644)
    lock_path = tmp_path / "lease.json.lock"
    lock_path.symlink_to(victim)

    with pytest.raises(GpuLeaseBusy, match="unsafe lock path"):
        gpu_lane_available(lease_path=lease_path)
    assert victim.read_text() == "KEEP\n"
    assert stat.S_IMODE(victim.stat().st_mode) == 0o644

    lock_path.unlink()
    lease_path.write_text(
        json.dumps(
            {
                "owner": "dead",
                "pid": 999_999_999,
                "started": "2026-08-27T12:00:00Z",
                "ttl": 1.0,
                "expected_gb": 16.0,
            }
        )
    )
    (tmp_path / "stale-breaks.log").symlink_to(victim)
    with pytest.raises(GpuLeaseBusy, match="unsafe stale-break path"):
        gpu_lane_available(lease_path=lease_path, pid_alive=lambda pid: False)
    assert victim.read_text() == "KEEP\n"
    assert stat.S_IMODE(victim.stat().st_mode) == 0o644
    assert lease_path.exists()


def test_phase_drivers_wire_gpu_lease_and_memory_override():
    scripts = Path(__file__).resolve().parent.parent / "scripts"
    owners = {
        "phase1_xiao_lasaga.py": "phase1_xiao_lasaga",
        "phase2_ladder.py": "phase2_ladder",
    }
    for filename, owner in owners.items():
        source = (scripts / filename).read_text()
        assert f'gpu_owner="{owner}"' in source
        assert '"--gpu-mem-gb"' in source
