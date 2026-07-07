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

import sys
from pathlib import Path

import gc_gvdw as gc
import itc_data
import numpy as np
from scipy.optimize import brentq, curve_fit, minimize_scalar

# The MC-grounded correlation-corrected gvdW binding lives in ../zeta.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "zeta"))
import zeta_ccgvdw as cc  # noqa: E402

NA, D_A = 6.02214076e23, 10.6

# 5% OB in the ITC cell carries oleosin/phospholipid counterions plus the 3 mM
# buffer and titrated phytate, so the cell ionic strength sits well above the
# buffer; the MC-grounded thermogram is best fit near here (see fit_cpb_ionic).
I_ITC = 0.015
_IMC = (0.005, 0.03, 0.15)  # MC ionic strengths available for PMF interpolation


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


# --- MC-grounded correlation-corrected gvdW (cPB) --------------------------------
def _cpb_coverage(cfree_M: float, ionic: float, pH: float = 4.0, nit: int = 200) -> float:
    """Self-consistent cPB areal coverage at free phytate ``cfree_M``.

    Binding comes from the correlation-corrected model (K_H from the MC PMF,
    Forsman-corrected a2); the PMF and potential are interpolated in log-ionic-
    strength between the simulated values, so the ITC cell ionic strength need not
    coincide with a simulated one.
    """
    imc = np.array(_IMC)
    j = max(min(int(np.searchsorted(imc, ionic)), len(imc) - 1), 1)
    lo, hi = imc[j - 1], imc[j]
    t = (np.log(ionic) - np.log(lo)) / (np.log(hi) - np.log(lo))
    w_lo, w_hi = cc._pmf(pH, lo), cc._pmf(pH, hi)
    r_w = w_lo.r.to_numpy()
    w0 = (1 - t) * w_lo.pmf_kT.to_numpy() + t * w_hi.pmf_kT.to_numpy()
    z = (1 - t) * cc.phytate_valence(pH, lo) + t * cc.phytate_valence(pH, hi)
    r = np.linspace(0, 180, 400)
    y0 = (1 - t) * cc._potential_grid(pH, lo, r) + t * cc._potential_grid(pH, hi, r)
    u = np.zeros_like(r)
    gamma = 0.0
    for _ in range(nit):
        alpha = cc.lateral_alpha(ionic, z, gamma)
        w_eff = w0 - z * np.interp(r_w, r, u)
        gamma += 0.15 * (cc._coverage(np.trapezoid(np.expm1(-w_eff), r_w), cfree_M, alpha) - gamma)
        excess = np.expm1(-w_eff).clip(0)
        if excess.sum() == 0:
            return 0.0
        excess /= np.trapezoid(excess, r_w)
        rho = -z * gamma * np.interp(r, r_w, excess, left=0, right=0)
        u += 0.15 * (cc._linear_shift(rho, y0, ionic, r) - u)
    return gamma


def _cpb_bound(titr, area, grid, cov):
    """Bound phytate [M] per injection from the interpolated cPB isotherm."""
    log_grid = np.log(grid)

    def gamma_of(cfree):
        return np.interp(np.log(cfree), log_grid, cov)

    m_t, l_t = _dilute(titr)
    out = np.empty_like(l_t)
    for i, (mt, ltot) in enumerate(zip(m_t, l_t, strict=True)):
        balance = lambda cf, mt=mt, ltot=ltot: ltot - cf - mt * area * gamma_of(cf)  # noqa: E731
        out[i] = ltot if balance(1e-15) <= 0 else mt * area * gamma_of(brentq(balance, 1e-15, ltot))
    return out


def fit_cpb(titr, ionic: float = I_ITC):
    """Fit the MC-grounded cPB isotherm to the raw thermogram at fixed ``ionic``.

    Only the area per oleosin and dH are fit; the binding strength is read from the
    MC PMF, so this transfers the simulated adsorption to the calorimetry rather
    than fitting it. Returns (area_A2, dH_kJ, rms_uJ).
    """
    grid = np.logspace(-9, -2, 24)
    cov = np.array([_cpb_coverage(c, ionic) for c in grid])

    def model(x, area, dH):
        return _heat_uJ(titr, _cpb_bound(titr, area, grid, cov), dH)

    popt, _ = curve_fit(
        model,
        titr.mole_ratio,
        titr.heat_uJ,
        p0=[3000.0, -40.0],
        bounds=([50, -1e4], [1e5, 0]),
        maxfev=20000,
    )
    rms = float(np.sqrt(np.mean((titr.heat_uJ - model(titr.mole_ratio, *popt)) ** 2)))
    return (*popt, rms)


def fit_cpb_ionic(titr):
    """Fit the cPB isotherm with the cell ionic strength free; returns (ionic, area, dH, rms)."""
    sol = minimize_scalar(
        lambda log_I: fit_cpb(titr, ionic=10**log_I)[-1],
        bounds=(np.log10(0.006), np.log10(0.05)),
        method="bounded",
        options={"xatol": 0.03},
    )
    area, dH, rms = fit_cpb(titr, ionic=10**sol.x)
    return (10**sol.x, area, dH, rms)


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

    s, dHg, cap, rmsg = fit_gvdw(t, mode="fitS", fixed=I_ITC)
    print(
        f"gvdW (fit sigma_s at I={I_ITC * 1e3:.0f} mM): sigma_s={s:.2e} e/A^2  "
        f"cap={cap:.2f}  dH={dHg:.1f} kJ/mol  RMS={rmsg:.1f} uJ\n"
    )

    area, dHc, rmsc = fit_cpb(t, ionic=I_ITC)
    print(
        f"cPB (MC binding, I={I_ITC * 1e3:.0f} mM): area={area:.0f} A^2/oleosin "
        f"(~{area * 3e-3:.0f} peptide-equiv)  dH={dHc:.1f} kJ/mol  RMS={rmsc:.1f} uJ"
    )
    ieff, area2, dH2, rms2 = fit_cpb_ionic(t)
    print(
        f"cPB (MC binding, fit I): I={ieff * 1e3:.0f} mM  area={area2:.0f} A^2/oleosin  "
        f"dH={dH2:.1f} kJ/mol  RMS={rms2:.1f} uJ"
    )
