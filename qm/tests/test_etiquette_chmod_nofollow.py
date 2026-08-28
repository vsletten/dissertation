import os
import stat
from datetime import UTC, datetime
from pathlib import Path

from quarry import etiquette


def _payload() -> dict[str, object]:
    return {
        "owner": "phase1",
        "pid": 4242,
        "started": "2026-08-27T12:00:00Z",
        "ttl": 1.0,
        "expected_gb": 16.0,
    }


def test_write_lease_chmod_does_not_follow_post_replace_symlink(tmp_path, monkeypatch):
    path = tmp_path / "lease.json"
    victim = tmp_path / "victim.txt"
    victim.write_text("KEEP\n")
    victim.chmod(0o644)
    real_replace = os.replace

    def replace_then_swap(src, dst, *args, **kwargs):
        real_replace(src, dst, *args, **kwargs)
        dest = Path(dst)
        dest.unlink()
        dest.symlink_to(victim)

    monkeypatch.setattr(os, "replace", replace_then_swap)
    etiquette._write_lease(path, _payload())

    assert victim.read_text() == "KEEP\n"
    assert stat.S_IMODE(victim.stat().st_mode) == 0o644


def test_stale_break_chmod_does_not_follow_post_write_symlink(tmp_path, monkeypatch):
    lease_path = tmp_path / "lease.json"
    note_path = tmp_path / "stale-breaks.log"
    victim = tmp_path / "victim.txt"
    victim.write_text("KEEP\n")
    victim.chmod(0o644)
    original_fdopen = os.fdopen

    def fdopen_swap_on_close(fd, *args, **kwargs):
        stream = original_fdopen(fd, *args, **kwargs)
        closer = stream.close

        def close_and_swap():
            closer()
            if note_path.exists() and not note_path.is_symlink():
                note_path.unlink()
                note_path.symlink_to(victim)

        stream.close = close_and_swap
        return stream

    monkeypatch.setattr(os, "fdopen", fdopen_swap_on_close)
    etiquette._record_stale_break(
        lease_path,
        _payload(),
        now=datetime(2026, 8, 27, 12, 0, tzinfo=UTC),
    )

    assert victim.read_text() == "KEEP\n"
    assert stat.S_IMODE(victim.stat().st_mode) == 0o644


def test_lease_and_stale_note_modes_are_0600_via_fchmod(tmp_path, monkeypatch):
    lease_path = tmp_path / "lease.json"
    original_chmod = Path.chmod

    def refuse_path_chmod(self, *args, **kwargs):
        if self.name in {"lease.json", "stale-breaks.log"}:
            raise AssertionError("state files must not chmod by path")
        return original_chmod(self, *args, **kwargs)

    monkeypatch.setattr(Path, "chmod", refuse_path_chmod)
    etiquette._write_lease(lease_path, _payload())
    etiquette._record_stale_break(
        lease_path,
        _payload(),
        now=datetime(2026, 8, 27, 12, 0, tzinfo=UTC),
    )

    assert stat.S_IMODE(lease_path.stat().st_mode) == 0o600
    assert stat.S_IMODE((tmp_path / "stale-breaks.log").stat().st_mode) == 0o600
