#!/usr/bin/env python3
"""Fit the phytate-into-OB ITC titrations (real NanoAnalyze heats from itc_data).

Two binding models share one titration engine:

  * one-site  -- the independent (Wiseman) model of the manuscript; reproduces
    NanoAnalyze's Ka, n, dH as a validation of the pipeline.
  * gvdW      -- the analytic Gouy-Chapman surface-adsorption isotherm
    (gc_gvdw.py): K_H and the lateral coefficient a2 follow from one surface
    charge sigma_s, and a coverage->bound-PA capacity closes the mass balance.

Both use the same dilution-corrected titration: injecting dv into the fixed cell
volume V0 dilutes the cell content and overflows the excess (MicroCal/Origin
displacement heat). Concentrations are absolute (syringe 0.30 mM PA, oleosin
recovered by itc_data), so Ka and dH come out in physical units.
"""

from __future__ import annotations

import gc_gvdw as gc
import itc_data
import numpy as np
from scipy.optimize import brentq, curve_fit

NA, D_A = 6.02214076e23, 10.6


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
    m_t, l_t = _dilute(titr)
    return _heat_uJ(titr, onesite_bound(m_t, l_t, Ka, n), dH_kJ)


def fit_onesite(titr: itc_data.Titration):
    """Independent-model fit; returns (Ka, dH_kJ, n, rms_uJ)."""
    popt, _ = curve_fit(
        lambda x, Ka, dH, n: onesite_heat(titr, Ka, dH, n),
        titr.mole_ratio,
        titr.heat_uJ,
        p0=[2e5, -50.0, 0.33],
        maxfev=200000,
    )
    rms = float(np.sqrt(np.mean((titr.heat_uJ - onesite_heat(titr, *popt)) ** 2)))
    return (*popt, rms)


# --- gvdW surface adsorption -----------------------------------------------------
def _coverage(cfree_M: float, K_H: float, alpha: float, d: float = D_A) -> float:
    """gvdW 2D-isotherm coverage theta at free phytate concentration ``cfree_M``."""
    phi0 = K_H * (cfree_M * NA / 1e27) * d**2
    if phi0 <= 0:
        return 0.0

    def f(t):
        return t - phi0 * (1 - t) * np.exp(alpha * t - t / (1 - t))

    return brentq(f, 1e-13, 1 - 1e-12) if f(1e-13) < 0 else min(phi0, 1 - 1e-12)


def gvdw_bound(m_t, l_t, K_H, alpha, cap):
    """Bound phytate [M]; ``cap`` scales surface coverage to bound PA per oleosin."""
    out = np.empty_like(l_t)
    for i, (m, ltot) in enumerate(zip(m_t, l_t, strict=True)):
        balance = lambda cf, m=m, ltot=ltot: ltot - cf - cap * m * _coverage(cf, K_H, alpha)  # noqa: E731
        out[i] = ltot - brentq(balance, 0.0, ltot)
    return out


def gvdw_heat(titr, sigma_s, dH_kJ, cap, ionic, z=gc.Z_PA):
    m_t, l_t = _dilute(titr)
    K_H = gc.K_H_gc(sigma_s, -z, ionic)
    alpha = gc.alpha_gc(sigma_s, z, ionic)
    return _heat_uJ(titr, gvdw_bound(m_t, l_t, K_H, alpha, cap), dH_kJ)


def fit_gvdw(titr, mode="fitS", fixed=gc.I_BUFFER):
    """Fit the gvdW model. mode='fitS': fit sigma_s at fixed ionic strength `fixed`;
    mode='fitI': fit ionic strength at fixed sigma_s `fixed`. Returns (charge_or_I,
    dH_kJ, cap, rms_uJ)."""

    def model(x, p, dH, cap):
        sigma_s, ionic = (10**p, fixed) if mode == "fitS" else (fixed, 10**p)
        return gvdw_heat(titr, sigma_s, dH, cap, ionic)

    p0 = np.log10(1e-3) if mode == "fitS" else np.log10(0.02)
    lo, hi = (np.log10(2e-4), np.log10(3e-3)) if mode == "fitS" else (np.log10(3e-3), np.log10(0.5))
    popt, _ = curve_fit(
        model,
        titr.mole_ratio,
        titr.heat_uJ,
        p0=[p0, -50.0, 5.0],
        bounds=([lo, -1e4, 0], [hi, 0, 1e4]),
        maxfev=100000,
    )
    rms = float(np.sqrt(np.mean((titr.heat_uJ - model(titr.mole_ratio, *popt)) ** 2)))
    return (10 ** popt[0], popt[1], popt[2], rms)


if __name__ == "__main__":
    t = itc_data.load(4.0)
    print(
        f"pH {t.pH} ITC ({t.source}), {len(t.heat_uJ)} injections, [oleosin]={t.mt0 * 1e6:.1f} uM\n"
    )

    Ka, dH, n, rms = fit_onesite(t)
    print("one-site (independent):")
    print(f"  mine        : Ka={Ka:.3e}/M  n={n:.3f}  dH={dH:.1f} kJ/mol  RMS={rms:.1f} uJ")
    f = t.na_fit
    print(f"  NanoAnalyze : Ka={f['Ka']:.3e}/M  n={f['n']:.3f}  dH={f['dH']:.1f} kJ/mol\n")

    s, dHg, cap, rmsg = fit_gvdw(t, mode="fitS", fixed=gc.I_BUFFER)
    print(f"gvdW (fit sigma_s at I={gc.I_BUFFER * 1e3:.0f} mM):")
    print(
        f"  sigma_s={s:.2e} e/A^2  cap={cap:.2f} PA/oleosin  dH={dHg:.1f} kJ/mol  RMS={rmsg:.1f} uJ"
    )
