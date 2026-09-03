#!/usr/bin/env python3
"""Rebuild the archived 1999 thesis without modifying its source tree."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

GENERATED_SUFFIXES = {
    ".aux", ".toc", ".lof", ".lot", ".bbl", ".blg", ".dvi", ".pdf",
    ".ps", ".log", ".fls", ".fdb_latexmk",
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True, help="read-only thesis archive root")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = args.source.resolve()
    script_dir = Path(__file__).resolve().parent
    with tempfile.TemporaryDirectory(prefix="thesis-a8-build-") as temp:
        build = Path(temp) / "text"
        shutil.copytree(source / "text", build)
        for path in build.iterdir():
            if path.is_file() and path.suffix in GENERATED_SUFFIXES:
                path.unlink()
        for path in build.glob("*.tex"):
            text = path.read_text(encoding="latin-1")
            text = text.replace("/home/accts/vsletten/thesis", str(source))
            path.write_text(text, encoding="latin-1")
        env = os.environ.copy()
        env["TEXINPUTS"] = f"{script_dir / 'texmf'}{os.pathsep}" + env.get("TEXINPUTS", "")
        subprocess.run(
            ["latexmk", "-pdfdvi", "-interaction=nonstopmode", "-halt-on-error", "thesis.tex"],
            cwd=build,
            env=env,
            check=True,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(build / "thesis.pdf", args.output)
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
