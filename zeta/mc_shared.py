#!/usr/bin/env python3
"""Shared MC-grounded utilities for the adsorption analyses.

This module holds constants, profile loaders, Debye/Grahame relations, and measured
zeta-potential selectors used by the LC--AR--PB and ITC analyses.
"""

from __future__ import annotations

from functools import cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import constants

# kT/e at the simulation temperature, in mV (zeta and MC potentials are in mV).
_T = 298.15
KT_MV = constants.k * _T / constants.e * 1e3
N_A = constants.N_A
# Water permittivity Faunus uses at 298.15 K; back-computed from the bjerrum_length/Å
# (7.152762) in its run output, so theory and MC share one electrostatic scale.
EPS_R = 78.3557
# Bjerrum length [A]
LB = constants.e**2 / (4 * constants.pi * constants.epsilon_0 * EPS_R * constants.k * _T) * 1e10
# Phytate hard-disk diameter [A] ~ 2 a_H, sets the 2D excluded area and the a2 range.
D_PA = 10.6

R_CONTACT = 26.0  # A, brush-periphery binding position (PMF-minimum region)

DILUTION = 1000.0  # emulsion is diluted 1000x before the zeta measurement
IONIC_MC = (0.005, 0.03, 0.15)  # MC ionic strengths available
PHS = (4.0, 6.0, 10.0)  # pH conditions carried into the manuscript

_MC = Path(__file__).resolve().parent.parent / "mc"


@cache
def _table(name: str) -> pd.DataFrame:
    return pd.read_csv(_MC / name)


def kappa(ionic: float) -> float:
    """Inverse Debye length [1/A] for a 1:1 electrolyte at ionic strength `ionic` (M)."""
    return np.sqrt(8 * np.pi * LB * ionic * N_A / 1e27)


def sigma_from_zeta(zeta_mV: float, ionic: float) -> float:
    """Grahame surface charge [e/A^2] giving potential `zeta_mV` at ionic strength `ionic`."""
    return (kappa(ionic) / (2 * np.pi * LB)) * np.sinh((zeta_mV / KT_MV) / 2)


def _ionic_bracket(ionic: float):
    """(lo, hi, t) for log-interpolating MC profiles between the bracketing simulated ionics."""
    imc = IONIC_MC
    if ionic <= imc[0]:
        return imc[0], imc[0], 0.0
    if ionic >= imc[-1]:
        return imc[-1], imc[-1], 0.0
    j = int(np.searchsorted(np.asarray(imc), ionic))
    lo, hi = imc[j - 1], imc[j]
    return lo, hi, (np.log(ionic) - np.log(lo)) / (np.log(hi) - np.log(lo))


@cache
def _exp() -> pd.DataFrame:
    """Measured zeta potentials (mean, SD, n) parsed from the tidied workbook."""
    return pd.read_csv(Path(__file__).resolve().parent / "measured_zeta.csv")


def _fig6a(pH: float):
    """(PA_mM, mean_mV, sd_mV) for PA>0 at low salt -- the zeta vs [PA] series."""
    d = _exp()
    d = d[(d.series == "PA") & (d.pH == pH) & (d.PA_mM > 0)].sort_values("PA_mM")
    return d.PA_mM.to_numpy(), d.mean_mV.to_numpy(), d.sd_mV.to_numpy()


def _fig6c(pH: float, pa_mM: float):
    """(NaCl_mM, mean_mV, sd_mV) at fixed PA -- the zeta vs ionic-strength series."""
    d = _exp()
    d = d[(d.series == "NaCl") & (d.pH == pH) & (d.PA_mM == pa_mM)].sort_values("NaCl_mM")
    return d.NaCl_mM.to_numpy(), d.mean_mV.to_numpy(), d.sd_mV.to_numpy()


# Bare-OB zeta at low salt (PA=0) anchors the far-field offset; read from the data.
ZETA_BARE = {
    pH: float(_exp().query("series=='PA' and pH==@pH and PA_mM==0").mean_mV.iloc[0]) for pH in PHS
}
