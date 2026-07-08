#!/usr/bin/env python3
"""Stage 1 -- extract the analysis CSV layer from the raw MC output in ``mc/runs/``.

Runs the two Monte-Carlo extractors over every run directory, writing the committed CSV
layer that every figure and table is built from:

  mc/pmf_profiles.csv        PMF(r) per (pH, ionic)                 (extract_pmf.py)
  mc/window_profiles.csv     phytate/peptide charge, rotational     (extract_profiles.py)
                             mobility and peptide shape vs r
  mc/potential_profiles.csv  electric potential phi(r)              (extract_profiles.py)
  mc/simulated_zeta.csv      shear-plane zeta per condition         (extract_profiles.py)

This is the only step that needs ``mc/runs/`` (the large raw Faunus output, obtained
separately). Once the CSVs are committed, ``make_figures.py`` reproduces every figure and
table without it. Run: ``uv run extract_csvs.py``.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

MC = Path(__file__).resolve().parent / "mc"


def _run(script: str) -> None:
    print(f"\n>>> {script}")
    subprocess.run([sys.executable, str(MC / script)], cwd=MC, check=True)


if __name__ == "__main__":
    if not any((MC / "runs").glob("*/pmf.csv")):
        sys.exit("mc/runs/ not found -- obtain the raw MC output to regenerate the CSV layer.")
    _run("extract_pmf.py")  # -> mc/pmf_profiles.csv (+ pmf.png)
    _run("extract_profiles.py")  # -> window/potential profiles + simulated_zeta
    print("\nCSV layer regenerated in mc/.")
