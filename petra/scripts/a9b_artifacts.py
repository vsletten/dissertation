"""Generic deterministic file primitives for A9b artifacts."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from collections.abc import Iterable
from pathlib import Path


def canonical_json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_bytes_atomic(path: Path, payload: bytes, *, overwrite: bool = True) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        if not overwrite and path.exists():
            raise FileExistsError(path)
        if overwrite:
            os.replace(temporary, path)
        else:
            os.link(temporary, path)
            os.unlink(temporary)
        directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def write_json_atomic(path: Path, value: object) -> None:
    write_bytes_atomic(path, canonical_json_bytes(value))


def write_json_atomic_new(path: Path, value: object) -> None:
    write_bytes_atomic(path, canonical_json_bytes(value), overwrite=False)


def remove_file_durable(path: Path) -> bool:
    """Atomically invalidate an artifact and durably record its removal."""
    try:
        path.unlink()
    except FileNotFoundError:
        return False
    directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)
    return True


def read_jsonl(path: Path) -> list[object]:
    result: list[object] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            try:
                result.append(json.loads(line))
            except json.JSONDecodeError as exc:
                raise ValueError(f"{path}:{line_number}: invalid JSON") from exc
    return result


def exact_files(root: Path, expected: Iterable[str]) -> None:
    expected_set = set(expected)
    actual = {path.name for path in root.iterdir() if path.is_file()}
    if actual != expected_set:
        raise ValueError(
            f"{root}: file set mismatch missing={sorted(expected_set - actual)} "
            f"extra={sorted(actual - expected_set)}"
        )
