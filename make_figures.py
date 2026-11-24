#!/usr/bin/env python3
"""Stage 2 -- regenerate every paper figure and table from the committed CSV layer.

Needs only the committed CSVs and experimental workbooks, not ``mc/runs/``. Produces:

  zeta/measured_zeta.csv           measured zeta, parsed from the workbook
  mc/composite_I{0.005,0.03,0.15}.png   MC composites (main-text + SI figures)
  zeta/fig_lc_arpb_schematic.png LC--AR--PB coupling schematic     (Fig S2)
  zeta/fig_lateral_arpb.png        zeta vs [PA]/[NaCl]           (Fig, Table S-3a)
  zeta/model_zeta.csv              the same predicted zeta scans, tabulated
  itc/itc_injections.csv          all ITC heats, one row per injection
  itc/fig_itc_lateral_arpb.png     pH-4 thermogram fit           (Fig, Table S-6)

and prints the model parameters for Tables S-3a and S-6. Run: ``uv run make_figures.py``.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def _run(rel: str) -> None:
    script = ROOT / rel
    print(f"\n>>> {rel}")
    subprocess.run([sys.executable, str(script)], cwd=script.parent, check=True)


if __name__ == "__main__":
    _run("zeta/parse_zeta_workbook.py")  # -> zeta/measured_zeta.csv
    _run("mc/plot_composites.py")  # -> mc/composite_I*.png
    _run("zeta/plot_lc_arpb_schematic.py")  # -> zeta/fig_lc_arpb_schematic.png
    _run("zeta/lateral_arpb.py")  # -> zeta/fig_lateral_arpb.png + model_zeta.csv + Table S-3a
    _run("itc/extract_thermograms.py")  # -> itc/itc_injections.csv
    _run("itc/fit_lateral_arpb.py")  # -> itc/fig_itc_lateral_arpb.png + Table S-6
    print("\nAll figures and tables regenerated.")
