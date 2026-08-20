import os
import re
import sys
from pathlib import Path

import pytest

from quarry.etiquette import THREAD_ENV_VARS, bootstrap_cli, enforce


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
        bootstrap_at = source.find("bootstrap_cli(")
        heavy_imports = [
            match.start()
            for match in re.finditer(
                r"^(?:import numpy|from (?!quarry\.etiquette\b)quarry\.)",
                source,
                re.MULTILINE,
            )
        ]
        assert bootstrap_at >= 0, f"{path.name} does not bootstrap etiquette"
        assert heavy_imports, f"{path.name} has no recognized heavy import"
        assert bootstrap_at < min(heavy_imports), (
            f"{path.name} bootstraps after a heavy import"
        )
        assert '"--threads"' in source
        assert '"--nice"' in source
        assert '"--log"' in source
