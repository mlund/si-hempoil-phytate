#!/usr/bin/env python3
"""Stage 2 -- regenerate every paper figure and table from the committed CSV layer.

Needs only the committed CSVs plus the zeta workbook and the hydrodynamics libraries --
NOT ``mc/runs/``. Produces:

  zeta/measured_zeta.csv           measured zeta, parsed from the workbook
  nmr/free_phytate_diffusion.csv   free-phytate hydrodynamics (Table S-5)
  mc/composite_I{0.005,0.03,0.15}.png   MC composites (main-text + SI figures)
  zeta/fig_crpb.png                zeta vs [PA]/[NaCl]           (Fig, Table S-3a)
  itc/fig_itc_crpb.png             pH-4 thermogram fit           (Fig, Table S-6)

and prints the paper tables (N, S-3a, S-6) to stdout. Run: ``uv run make_figures.py``.
Then ``uv run writeup/build.py`` embeds the figures into the Word/ODF write-up.
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
    _run("nmr/rotational_diffusion.py")  # -> nmr/free_phytate_diffusion.csv (Table S-5)
    _run("mc/plot_composites.py")  # -> mc/composite_I*.png
    _run("zeta/charge_regulation.py")  # -> fig_crpb.png, itc/fig_itc_crpb.png + Tables N/S-3a/S-6
    print("\nAll figures and tables regenerated.")
