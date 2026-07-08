#!/usr/bin/env python3
"""Shared MC-grounded utilities for the adsorption analyses.

This module holds the pieces that the charge-regulated zeta model
(``charge_regulation.py``) and the NMR bound-fraction bridge (``nmr/binding_fraction.py``)
both build on: physical constants matched to the Faunus MC water model, loaders for the
extracted MC profile CSVs, the Debye/Grahame double-layer relations, the simulated
phytate valence, and the single-phytate PMF (with the shared non-electrostatic contact
well) interpolated in ionic strength.

The non-electrostatic well ``EPS_NONEL`` is a specific phosphate--Lys/Arg contact
attraction the coarse-grained MC under-represents; being salt-robust, one value carries
the binding across the 1000x-diluted zeta (~5 mM), the 15 mM ITC and the 150 mM NMR
regimes. Overridable via ``$EPS_NONEL``.
"""

from __future__ import annotations

import os
from functools import cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import constants

# kT/e at the simulation temperature, in mV (zeta and MC potentials are in mV).
_T = 298.15
KT_MV = constants.k * _T / constants.e * 1e3
N_A = constants.N_A
# Bjerrum length of the MC water model [A]; must match Faunus (298.15 K, eps_r 78.4)
# and mc/extract_pmf.py::BJERRUM_A.
LB = 7.1528
# Phytate hard-disk diameter [A] ~ 2 a_H, sets the 2D excluded area and the a2 range.
D_PA = 10.6

# Non-electrostatic contact attraction (screened Coulomb + excluded volume) modelled as a
# salt-robust Gaussian well of depth EPS_NONEL [kT] at the brush-contact position, added to
# every PMF and shared by the zeta and NMR analyses.
EPS_NONEL = float(os.environ.get("EPS_NONEL", "2.5"))  # kT (reconciles NMR/ITC/zeta binding)
R_CONTACT = 26.0  # A, brush-periphery binding position (PMF-minimum region)
W_CONTACT = 6.0  # A, contact-well half-width

DILUTION = 1000.0  # emulsion is diluted 1000x before the zeta measurement
IONIC_MC = (0.005, 0.03, 0.15)  # MC ionic strengths available
PHS = (4.0, 6.0, 10.0)  # pH conditions carried into the manuscript

_MC = Path(__file__).resolve().parent.parent / "mc"


def _nonel_well(r: np.ndarray) -> np.ndarray:
    """Non-electrostatic contact attraction added to the PMF [kT] (<= 0, deepens well)."""
    return -EPS_NONEL * np.exp(-0.5 * ((np.asarray(r) - R_CONTACT) / W_CONTACT) ** 2)


@cache
def _table(name: str) -> pd.DataFrame:
    return pd.read_csv(_MC / name)


def kappa(ionic: float) -> float:
    """Inverse Debye length [1/A] for a 1:1 electrolyte at ionic strength `ionic` (M)."""
    return np.sqrt(8 * np.pi * LB * ionic * N_A / 1e27)


def _pmf(pH: float, ionic: float) -> pd.DataFrame:
    """Bulk-referenced single-phytate PMF [kT] with the non-electrostatic well added."""
    w = _table("pmf_profiles.csv")
    w = w[(w.ionic_strength == ionic) & (w.pH == pH)].sort_values("r").copy()
    w["pmf_kT"] = w.pmf_kT + _nonel_well(w.r.to_numpy())  # add non-electrostatic contact well
    return w[w.r < 70]


def phytate_valence(pH: float, ionic: float) -> float:
    """Adsorption-weighted phytate valence |z| from the MC (window_profiles.csv, INO charge).

    The charge is read where the phytate actually binds (Boltzmann-weighted over the
    PMF well), so it is the peptide-field-regulated value, not the bulk one.
    """
    w = _pmf(pH, ionic)
    weight = np.expm1(-w.pmf_kT.to_numpy()).clip(0)
    if weight.sum() == 0:
        return 6.0
    c = _table("window_profiles.csv")
    c = c[(c.ionic_strength == ionic) & (c.pH == pH)].sort_values("r")
    z_of_r = np.interp(w.r.to_numpy(), c.r.to_numpy(), -c.q_ino.to_numpy())
    return float(np.average(z_of_r, weights=weight))


def sigma_from_zeta(zeta_mV: float, ionic: float) -> float:
    """Grahame surface charge [e/A^2] giving potential `zeta_mV` at ionic strength `ionic`."""
    return (kappa(ionic) / (2 * np.pi * LB)) * np.sinh((zeta_mV / KT_MV) / 2)


def lateral_alpha(ionic: float, z: float, gamma: float, sigma_net: float | None = None) -> float:
    """Forsman-corrected lateral coefficient alpha = 2 a2 / (kT d^2).

    The screened-Coulomb pair repulsion is replaced by its tangent below the
    correlation-hole radius R_c, following Forsman's full prescription
    R_c = min(Eq 8, Eq 12): the *density* estimate 1/sqrt(pi*gamma) (Eq 12) for the
    nearest-neighbour spacing of adsorbed phytates, and -- when a net surface charge
    `sigma_net` is supplied -- the *surface-charge* estimate sqrt(z/pi|sigma_net|)
    (Eq 8), which dominates for strongly charged surfaces (the pH-4 wall). R_c is
    capped at 2/kappa (correlations screened beyond ~two Debye lengths) and floored at
    the hard core d. The tangent piece is integrated analytically.
    """
    k = kappa(ionic)
    candidates = [2.0 / k, 1.0 / np.sqrt(np.pi * max(gamma, 1e-12))]
    if sigma_net is not None:
        candidates.append(np.sqrt(z / (np.pi * max(abs(sigma_net), 1e-12))))
    r_c = max(min(candidates), D_PA)
    e = np.exp(-k * r_c)
    u_c = z**2 * LB / r_c * e  # pair energy at R_c
    slope = z**2 * LB * (-1 / r_c**2 - k / r_c) * e  # its derivative at R_c
    outer = 2 * np.pi * z**2 * LB / k * e  # int_{R_c}^inf 2 pi r u(r) dr
    inner = (
        2
        * np.pi
        * (  # int_d^{R_c} 2 pi r [slope (r-R_c) + u_c] dr
            slope * (r_c**3 - D_PA**3) / 3
            - slope * r_c * (r_c**2 - D_PA**2) / 2
            + u_c * (r_c**2 - D_PA**2) / 2
        )
    )
    return 2 * (-0.5 * (inner + outer)) / D_PA**2


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


def _pmf_interp(pH: float, ionic: float):
    """(r_w, w0) PMF at arbitrary ionic, log-interpolated between MC ionics (nonel well in)."""
    lo, hi, t = _ionic_bracket(ionic)
    wl = _pmf(pH, lo)
    r_w, w_lo = wl.r.to_numpy(), wl.pmf_kT.to_numpy()
    if t == 0.0:
        return r_w, w_lo
    wh = _pmf(pH, hi)
    return r_w, (1 - t) * w_lo + t * np.interp(r_w, wh.r.to_numpy(), wh.pmf_kT.to_numpy())


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
