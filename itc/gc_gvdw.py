#!/usr/bin/env python3
"""Fully analytical gvdW binding model from Gouy-Chapman, one free parameter.

The peptide layer is a charged plane of surface charge density sigma_s (the single
free parameter, = peptide coverage x charge). Gouy-Chapman then fixes everything:

  surface potential   y0 = 2 arcsinh(2 pi lB sigma_s / kappa)        (Grahame)
  potential profile    y(x) = 4 arctanh( tanh(y0/4) e^{-kappa x} )
  K_H  = int_xmin^inf ( e^{-z_pa y(x)} - 1 ) dx                       (adsorption of the -z anion)
  a2   at local screening  I_local = I_bulk cosh(y0)  (surface counterion enhancement)

So K_H(sigma_s) and a2(sigma_s) are both analytic and self-consistent (same y0),
with the phytate valence z_pa known. No MC PMF required; sigma_s is fit to the ITC
(or read off from the simulated brush for a parameter-free prediction).
"""

import numpy as np
from scipy.integrate import quad

LB, NA, D_A = 7.1528, 6.02214076e23, 10.6
KT_MV = 25.693
X_MIN = 5.3  # closest phytate-centre approach to the plane ~ a_H [A]


def kappa(ionic):
    return np.sqrt(8 * np.pi * LB * ionic * NA / 1e27)  # 1/A


def y0_grahame(sigma_s, ionic):
    return 2 * np.arcsinh(2 * np.pi * LB * sigma_s / kappa(ionic))


def y_profile(x, y0, ionic):
    return 4 * np.arctanh(np.tanh(y0 / 4) * np.exp(-kappa(ionic) * x))


def K_H_gc(sigma_s, z_pa, ionic, x_min=X_MIN):
    """Analytic Henry constant [A]; z_pa signed (negative for the anion)."""
    y0 = y0_grahame(sigma_s, ionic)

    def integrand(x):
        return np.expm1(-z_pa * y_profile(x, y0, ionic))

    val, _ = quad(integrand, x_min, x_min + 40 / kappa(ionic), limit=200)
    return val


def alpha_gc(sigma_s, z, ionic, d=D_A):
    """Lateral alpha = 2 a2/(kT d^2), a2 screened at the surface-enhanced I_local."""
    y0 = y0_grahame(sigma_s, ionic)
    I_loc = ionic * np.cosh(y0)
    k = kappa(I_loc)
    a2 = -(np.pi * z**2 * LB / k) * np.exp(-k * d)
    return 2 * a2 / d**2


Z_PA = 5.94  # phytate valence at pH 4
SIGMA_BRUSH = 1.35e-3  # simulated brush: 3e-3 peptide/A^2 x ~0.45 e  [e/A^2]
I_BUFFER = 0.005  # ITC buffer (low salt)


def fit_itc(mode, fixed):
    """Fit the pH-4 thermogram (Fig S2C) with the analytic GC model.

    mode='fitI'   : sigma_s fixed at `fixed`, effective ionic strength fitted.
    mode='fitS'   : ionic strength fixed at `fixed`, surface charge fitted.
    Reuses the ITC titration machinery (digitised data + isotherm) from itc_fit.
    """
    import itc_fit as m
    from scipy.optimize import curve_fit

    def model(X, dH, cap, p):
        ionic, s = (10**p, fixed) if mode == "fitI" else (fixed, 10**p)
        m.K_H = K_H_gc(s, -Z_PA, ionic)
        m.ALPHA = alpha_gc(s, Z_PA, ionic)
        return dH * m.per_injection(m.gvdw_bound(X, cap))

    if mode == "fitI":
        p0, lo, hi = np.log10(0.02), np.log10(0.003), np.log10(0.5)
    else:
        p0, lo, hi = np.log10(1e-3), np.log10(2e-4), np.log10(3e-3)
    popt, _ = curve_fit(
        model, m.X, m.Q, p0=[-3000, 5, p0], bounds=([-1e5, 0, lo], [0, 1e4, hi]), maxfev=100000
    )
    rms = np.sqrt(np.mean((m.Q - model(m.X, *popt)) ** 2))
    return 10 ** popt[2], popt[1], rms


if __name__ == "__main__":
    from scipy.optimize import brentq

    ionic = I_BUFFER
    print(f"Gouy-Chapman gvdW binding, pH 4, ionic = {ionic * 1e3:.0f} mM, z = {Z_PA}\n")
    print(f"{'sigma_s(e/A^2)':>14} {'psi0(mV)':>9} {'K_H(A)':>11} {'I_local(mM)':>11} {'alpha':>7}")
    for s in (5e-4, 1e-3, 1.5e-3, 2e-3, 3e-3):
        y0 = y0_grahame(s, ionic)
        print(
            f"{s:>14.1e} {y0 * KT_MV:>8.0f}m {K_H_gc(s, -Z_PA, ionic):>11.2e} "
            f"{ionic * np.cosh(y0) * 1e3:>10.1f}m {alpha_gc(s, Z_PA, ionic):>7.0f}"
        )

    s_match = brentq(lambda s: np.log(K_H_gc(s, -Z_PA, ionic)) - np.log(4.46e8), 1e-4, 5e-3)
    print(
        f"\nvalidation: sigma_s reproducing MC K_H(5mM)=4.46e8 -> {s_match:.2e} e/A^2 "
        f"(brush {SIGMA_BRUSH:.2e})"
    )

    print("\ntwo equivalent one-parameter ITC fits (Fig S2C):")
    I_eff, cap1, r1 = fit_itc("fitI", SIGMA_BRUSH)
    print(
        f"  fix sigma_s = brush {SIGMA_BRUSH:.2e}  -> I_eff   = {I_eff * 1e3:5.1f} mM   "
        f"cap = {cap1:4.1f}  RMS = {r1:.1f}"
    )
    s_eff, cap2, r2 = fit_itc("fitS", I_BUFFER)
    print(
        f"  fix ionic = buffer {I_BUFFER * 1e3:.0f} mM      -> sigma_s = {s_eff:.2e}   "
        f"cap = {cap2:4.1f}  RMS = {r2:.1f}"
    )
