#!/usr/bin/env python3
"""Remove reproducible build products while preserving source and data files.

Run with ``uv run clean.py``. The tracked CSV analysis layer, raw simulation
directories, workbooks, and archives are deliberately outside the cleanup set.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent

GENERATED_FILES = (
    "mc/composite_I0.005.png",
    "mc/composite_I0.03.png",
    "mc/composite_I0.15.png",
    "zeta/fig_lc_arpb_schematic.png",
    "zeta/fig_lateral_arpb.png",
    "itc/fig_itc_lateral_arpb.png",
    "writeup/sim_for_mauricio.docx",
    "writeup/sim_for_mauricio.odt",
    "writeup/sim_for_mauricio.pdf",
)

SOURCE_DIRECTORIES = ("itc", "mc", "nmr", "tests", "writeup", "zeta")


def generated_paths() -> list[Path]:
    """Return existing build products and transient analysis caches."""
    paths = [ROOT / name for name in GENERATED_FILES]
    paths.extend((ROOT / "__pycache__", ROOT / ".pytest_cache", ROOT / ".ruff_cache"))
    for name in SOURCE_DIRECTORIES:
        source = ROOT / name
        paths.extend(source.rglob("__pycache__"))
        paths.extend(source.rglob(".ipynb_checkpoints"))
    return sorted({path for path in paths if path.exists()})


def main() -> None:
    """Delete known generated paths, or list them without deleting."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true", help="list paths without removing them")
    args = parser.parse_args()

    paths = generated_paths()
    for path in paths:
        print(path.relative_to(ROOT))
        if args.dry_run:
            continue
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()

    action = "Would remove" if args.dry_run else "Removed"
    print(f"{action} {len(paths)} generated path(s); data and CSV files were preserved.")


if __name__ == "__main__":
    main()
