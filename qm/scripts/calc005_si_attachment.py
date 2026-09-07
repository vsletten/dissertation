#!/usr/bin/env python3
"""CPU-only CALC-005 state-expanded pair probe.

The production minima driver will extend this command; this slice intentionally
contains no calculator import or electronic-structure execution.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    bootstrap_cli(
        "calc005_si_attachment",
        default_run_root=Path.home() / ".local" / "state" / "quarry" / "calc005",
    )

from quarry.calc005 import probe_calc005_pairs  # noqa: E402


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    probe = subparsers.add_parser("probe")
    probe.add_argument("--deck", type=Path, required=True)
    probe.add_argument("--repeat", type=int, default=2)
    probe.add_argument("--output", type=Path)
    # Parsed first by bootstrap_cli; repeated here so the full parser accepts them.
    probe.add_argument("--threads", type=int, default=16)
    probe.add_argument("--nice", type=int, default=10)
    probe.add_argument("--log")
    return parser


def _atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        Path(temporary).replace(path)
    finally:
        temporary_path = Path(temporary)
        if temporary_path.exists():
            temporary_path.unlink()


def main() -> int:
    args = _parser().parse_args()
    if args.repeat != 2:
        raise SystemExit("--repeat must be exactly 2 for the fixed proof")
    payload = probe_calc005_pairs(args.deck)
    payload["driver_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    content = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        _atomic_write(args.output, content)
    print(content, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
