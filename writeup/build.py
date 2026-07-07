#!/usr/bin/env python3
"""Convert the simulation write-up to Word/ODF for merging into the OB-PA manuscript.

Runs pandoc on ``sim_for_mauricio.md`` with a resource path that embeds the figures
from ``../mc``, ``../zeta`` and ``../itc``, writing ``.docx`` and ``.odt`` beside the
source (math becomes native Word/ODF equations).

    uv run writeup/build.py          # -> sim_for_mauricio.docx, .odt
    uv run writeup/build.py --pdf    # ... and a PDF proof via LibreOffice

Requires the ``pandoc`` binary on PATH; ``--pdf`` also needs ``soffice`` (LibreOffice).
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "sim_for_mauricio.md"


def _run(cmd: list[str]) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdf", action="store_true", help="also render a PDF proof via soffice")
    args = ap.parse_args()

    if shutil.which("pandoc") is None:
        sys.exit("pandoc not found on PATH (e.g. `conda install pandoc` or `brew install pandoc`)")

    resource_path = f"{HERE.parent}:{HERE}"  # so ../mc, ../zeta, ../itc figures embed
    for fmt in ("docx", "odt"):
        out = SRC.with_suffix(f".{fmt}")
        _run(["pandoc", str(SRC), "-o", str(out), f"--resource-path={resource_path}"])

    if args.pdf:
        if shutil.which("soffice") is None:
            sys.exit("soffice (LibreOffice) not found; cannot render --pdf")
        _run(
            [
                "soffice",
                "--headless",
                "--convert-to",
                "pdf",
                "--outdir",
                str(HERE),
                str(SRC.with_suffix(".docx")),
            ]
        )


if __name__ == "__main__":
    main()
