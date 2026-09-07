"""Generic crash-durable filesystem operations (no scientific policy)."""

import os
from pathlib import Path


def fsync_directory(path: Path) -> None:
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def durable_mkdir(path: Path) -> None:
    if path.is_dir():
        return
    durable_mkdir(path.parent)
    path.mkdir(exist_ok=True)
    fsync_directory(path)
    fsync_directory(path.parent)


def durable_replace(source: Path, target: Path) -> None:
    source.replace(target)
    fsync_directory(target.parent)
    if source.parent != target.parent:
        fsync_directory(source.parent)
