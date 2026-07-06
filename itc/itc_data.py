#!/usr/bin/env python3
"""Load TA NanoAnalyze .xls exports of the phytate-into-OB ITC titrations.

Each export holds NanoAnalyze's own independent-model fit (Kd, n, dH, Ka) plus the
per-injection integrated heats (``Area Data``, uJ) against mole ratio and their
fitted model. We read the heats and the titration geometry (injection volumes,
cell volume) and combine them with the concentrations from the manuscript methods
(0.30 mM PA in the syringe, 5% OB in the cell) into an absolute-concentration
titration ready for binding fits.

The oleosin molarity of "5% OB" is not tabulated, so we recover the moles of
oleosin in the cell from NanoAnalyze's mole-ratio column and the known syringe
concentration (mole_ratio = cumulative injected PA / oleosin in cell). Absolute
Ka and dH therefore inherit the ~few-% spread of that recovery; the curve shape,
n, and relative thermodynamics do not.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent

# Syringe phytate concentration, identical at every pH (manuscript Section 2.2.8).
PA_SYRINGE_M = 0.30e-3

# The .xls export whose NanoAnalyze fit the manuscript reports, per pH.
REPORTED = {4.0: "pH4/pH4_rep3.xls"}

# NanoAnalyze fit-parameter labels -> short keys.
_FIT_LABELS = {"Kd (M)": "Kd", "n": "n", "ΔH (kJ/mol)": "dH", "Ka (Mˉ¹)": "Ka"}


@dataclass
class Titration:
    """One ITC titration: per-injection heats + absolute-concentration geometry."""

    pH: float
    mole_ratio: np.ndarray  # X_t / M_t per injection (NanoAnalyze)
    heat_uJ: np.ndarray  # integrated heat per injection (Area Data)
    vinj_uL: np.ndarray  # injection volumes
    vcell_uL: float  # active cell volume
    csyr_M: float  # syringe phytate concentration
    moles_cell: float  # oleosin in the cell (recovered)
    na_fit: dict  # NanoAnalyze independent-model parameters
    na_model_uJ: np.ndarray  # NanoAnalyze per-injection model heat
    source: str

    @property
    def mt0(self) -> float:
        """Initial oleosin concentration in the cell [M]."""
        return self.moles_cell / (self.vcell_uL * 1e-6)


def _header_row(df: pd.DataFrame) -> int:
    """Row index of the per-injection table header (the one naming ``Area Data``)."""
    for i in range(len(df)):
        if df.iloc[i].astype(str).str.contains("Area Data", na=False).any():
            return i
    raise ValueError("no 'Area Data' header row found")


def _column(header: pd.Series, needle: str) -> int:
    """Index of the first header cell containing ``needle``."""
    hits = [j for j, v in enumerate(header) if needle in str(v)]
    if not hits:
        raise ValueError(f"column {needle!r} not found in {list(header)}")
    return hits[0]


def _nanoanalyze_fit(df: pd.DataFrame) -> dict:
    """Pull NanoAnalyze's independent-model parameters from the summary block."""
    flat = df.astype(str)
    fit = {}
    for label, key in _FIT_LABELS.items():
        loc = np.argwhere(flat.values == label)
        if loc.size:
            r, c = loc[0]
            fit[key] = float(df.iloc[r, c + 1])
    return fit


def load_xls(path: str | Path, pH: float) -> Titration:
    """Read one NanoAnalyze .xls export into a :class:`Titration`."""
    path = Path(path)
    df = pd.read_excel(path, header=None)

    h = _header_row(df)
    header = df.iloc[h]
    col = {
        name: _column(header, name)
        for name in (
            "Area Data",
            "Injection Volume",
            "Mole Ratio",
            "Cell Volume",
            "Independent Model",
        )
    }
    body = df.iloc[h + 1 :].copy()
    # Data rows are those whose Area Data cell is numeric (drops trailing blanks).
    heat = pd.to_numeric(body[col["Area Data"]], errors="coerce")
    body = body[heat.notna()]

    def num(name):
        return pd.to_numeric(body[col[name]], errors="coerce").to_numpy()

    heat_uJ = num("Area Data")
    vinj_uL = num("Injection Volume")
    mole_ratio = num("Mole Ratio")
    vcell_uL = float(num("Cell Volume")[0])
    na_model_uJ = num("Independent Model")

    # Recover oleosin in cell: mole_ratio = cumulative injected PA / oleosin_cell.
    cum_pa = PA_SYRINGE_M * np.cumsum(vinj_uL) * 1e-6  # mol
    moles_cell = float(np.median(cum_pa / mole_ratio))

    return Titration(
        pH=pH,
        mole_ratio=mole_ratio,
        heat_uJ=heat_uJ,
        vinj_uL=vinj_uL,
        vcell_uL=vcell_uL,
        csyr_M=PA_SYRINGE_M,
        moles_cell=moles_cell,
        na_fit=_nanoanalyze_fit(df),
        na_model_uJ=na_model_uJ,
        source=path.name,
    )


def load(pH: float = 4.0) -> Titration:
    """Load the manuscript-reported titration for ``pH``."""
    return load_xls(HERE / REPORTED[pH], pH)


if __name__ == "__main__":
    t = load(4.0)
    print(f"{t.source}: pH {t.pH}, {len(t.heat_uJ)} injections")
    print(f"  syringe PA = {t.csyr_M * 1e3:.2f} mM, cell = {t.vcell_uL:.0f} uL")
    print(
        f"  oleosin in cell = {t.moles_cell * 1e9:.1f} nmol  ->  [oleosin] = {t.mt0 * 1e6:.1f} uM"
    )
    print(f"  mole ratio {t.mole_ratio.min():.3f}..{t.mole_ratio.max():.3f}")
    print(f"  NanoAnalyze fit: {t.na_fit}")
