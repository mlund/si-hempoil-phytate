#!/usr/bin/env python3
"""Fit the phytate-into-OB ITC titrations (NanoAnalyze heats from extract_thermograms).

The one-site (independent / Wiseman) model reproduces NanoAnalyze's Ka, n, dH as a
validation of the pipeline. ``fit_lateral_arpb.py`` reuses this titration engine for
the LC-AR-PB surface fit. ``_bound_from_isotherm`` remains as a generic bridge for
callers that supply a tabulated coverage isotherm.

The dilution-corrected titration injects dv into the fixed cell volume V0, diluting the
cell content and overflowing the excess (MicroCal/Origin displacement heat).
Concentrations are absolute (syringe 0.30 mM PA, oleosin from extract_thermograms), so
Ka and dH come out in physical units.
"""

from __future__ import annotations

import extract_thermograms as itc_data
import numpy as np
from scipy.optimize import brentq, curve_fit


def _dilute(titr: itc_data.Titration) -> tuple[np.ndarray, np.ndarray]:
    """Total oleosin and phytate concentrations in the cell after each injection.

    Overflow model: an injection of dv into cell V0 mixes then expels dv of the
    (now diluted) content, so a cell-only species is scaled by V0/(V0+dv) and the
    injectant accumulates from the syringe concentration.
    """
    v0 = titr.vcell_uL * 1e-6
    m_t, l_t = np.empty(len(titr.vinj_uL)), np.empty(len(titr.vinj_uL))
    m, ligand = titr.mt0, 0.0
    for i, dv_uL in enumerate(titr.vinj_uL):
        dv = dv_uL * 1e-6
        m *= v0 / (v0 + dv)
        ligand = (ligand * v0 + titr.csyr_M * dv) / (v0 + dv)
        m_t[i], l_t[i] = m, ligand
    return m_t, l_t


def _heat_uJ(titr: itc_data.Titration, bound: np.ndarray, dH_kJ: float) -> np.ndarray:
    """Per-injection heat [uJ] from the bound-complex profile (Origin displacement)."""
    v0 = titr.vcell_uL * 1e-6
    dv = titr.vinj_uL * 1e-6
    prev = np.concatenate([[0.0], bound[:-1]])
    dmx = bound - prev + (dv / v0) * (bound + prev) / 2  # displacement-corrected
    return dH_kJ * 1e3 * v0 * dmx * 1e6  # kJ/mol * L * mol/L -> J -> uJ


# --- one-site (independent / Wiseman) -------------------------------------------
def onesite_bound(m_t: np.ndarray, l_t: np.ndarray, Ka: float, n: float) -> np.ndarray:
    """Bound-complex concentration [M] for the single-site model."""
    s = n * m_t  # binding sites
    b = s + l_t + 1.0 / Ka
    return 0.5 * (b - np.sqrt(np.clip(b * b - 4 * s * l_t, 0, None)))


def onesite_heat(titr, Ka, dH_kJ, n):
    """Per-injection heat [uJ] for the one-site model (Ka [1/M], dH [kJ/mol])."""
    m_t, l_t = _dilute(titr)
    return _heat_uJ(titr, onesite_bound(m_t, l_t, Ka, n), dH_kJ)


def _fit_onesite_pcov(titr: itc_data.Titration):
    """Least-squares core of the one-site fit; returns (popt=(Ka,dH,n), pcov).

    Skipped injections are removed from the residual only. The heat model still sees the
    whole titration, so their dilution and bound-complex history stay in the mass balance.
    """
    mask = titr.fit_mask
    return curve_fit(
        lambda x, Ka, dH, n: onesite_heat(titr, Ka, dH, n)[mask],
        titr.mole_ratio[mask],
        titr.heat_uJ[mask],
        p0=[2e5, -50.0, 0.33],
        maxfev=200000,
    )


def fit_onesite(titr: itc_data.Titration):
    """Independent-model fit; returns (Ka, dH_kJ, n, rms_uJ)."""
    popt, _pcov = _fit_onesite_pcov(titr)
    mask = titr.fit_mask
    residual = (titr.heat_uJ - onesite_heat(titr, *popt))[mask]
    return (*popt, float(np.sqrt(np.mean(residual**2))))


def fit_onesite_pcov(titr: itc_data.Titration):
    """One-site fit returning (popt=(Ka,dH,n), pcov) for parameter-correlation analysis."""
    return _fit_onesite_pcov(titr)


def corr_from_cov(pcov: np.ndarray) -> np.ndarray:
    """Correlation matrix from a covariance matrix."""
    d = np.sqrt(np.diag(pcov))
    return pcov / np.outer(d, d)


# --- adsorption-isotherm heat (coverage grid supplied by the caller) -------------
def _bound_from_isotherm(titr, area, grid, cov):
    """Bound phytate [M] per injection from a tabulated coverage isotherm.

    ``grid`` (free-PA concentrations) and ``cov`` (areal coverage) define the isotherm;
    the result turns that coverage into a per-injection heat after ligand mass balance.
    """
    log_grid = np.log(grid)

    def gamma_of(cfree):
        return np.interp(np.log(cfree), log_grid, cov)

    m_t, l_t = _dilute(titr)
    out = np.empty_like(l_t)
    for i, (mt, ltot) in enumerate(zip(m_t, l_t, strict=True)):
        balance = lambda cf, mt=mt, ltot=ltot: ltot - cf - mt * area * gamma_of(cf)  # noqa: E731
        out[i] = ltot if balance(1e-15) <= 0 else mt * area * gamma_of(brentq(balance, 1e-15, ltot))
    return out


if __name__ == "__main__":
    # pH-4 one-site sanity check (the only pH the manuscript fitted).
    t4 = itc_data.load(4.0)
    Ka, dH, n, rms = fit_onesite(t4)
    print(f"pH 4 one-site (independent) reproduction ({t4.source}):")
    print(f"  mine        : Ka={Ka:.2e}/M  n={n:.3f}  dH={dH:.1f} kJ/mol  RMS={rms:.1f} uJ")
    f = t4.na_fit
    print(f"  NanoAnalyze : Ka={f['Ka']:.2e}/M  n={f['n']:.3f}  dH={f['dH']:.1f} kJ/mol")
    print("  manuscript  : Ka=5e5/M  n=2.5(oleosin/PA => 0.40 PA/oleosin)  dH=-47 kJ/mol")
