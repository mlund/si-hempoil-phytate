#!/usr/bin/env python3
"""ITC titration engine and one-site reference fit for the pH-4 thermogram (Fig S2C).

Provides the digitised data and the machinery to turn a binding curve into per-injection
heats, plus the one-site (independent) model of the manuscript as a reference. The gvdW
binding physics is supplied externally by setting the module globals K_H and ALPHA (the
analytic Gouy-Chapman model in gc_gvdw.py does this); this module stays agnostic to where
they come from.

Data are the per-injection integrated heats digitised from Fig S2C (Area Data, uJ, vs mole
ratio X = PA:oleosin). Per-injection heat Q_i = dH * V * ([bound]_i - [bound]_{i-1}).

Caveats: points are read off a figure, and the cell/oleosin concentration are estimates,
so absolute affinities are indicative; the shape and dH are the robust outputs.
"""

import numpy as np
from scipy.optimize import brentq, curve_fit

# --- digitised Fig S2C (pH 4): mole ratio, integrated heat per injection (uJ) ---
X = np.array(
    [
        0.05,
        0.09,
        0.13,
        0.17,
        0.21,
        0.25,
        0.29,
        0.33,
        0.37,
        0.41,
        0.45,
        0.50,
        0.54,
        0.58,
        0.62,
        0.67,
        0.71,
        0.75,
        0.79,
        0.84,
        0.88,
        0.92,
        0.96,
        1.00,
        1.05,
    ]
)
Q = np.array(
    [
        -105,
        -92,
        -82,
        -73,
        -68,
        -63,
        -58,
        -53,
        -50,
        -46,
        -44,
        -42,
        -37,
        -30,
        -25,
        -22,
        -18,
        -16,
        -14,
        -12,
        -11,
        -10,
        -10,
        -9,
        -10.0,
    ]
)  # uJ per injection

NA, D_A = 6.02214076e23, 10.6
M_OLEO = 0.03e-3  # oleosin in the 5% OB cell [M] (estimate; sets the mole-ratio axis)

# gvdW binding parameters -- set by the binding model (e.g. gc_gvdw.fit_itc).
K_H = 0.0  # Henry constant [A]
ALPHA = 0.0  # lateral coefficient 2 a2/(kT d^2)


def onesite_bound(X, K, n):
    """Bound ligand per macromolecule for the independent one-site model."""
    Mt = M_OLEO
    r = 1.0 / (n * K * Mt)
    disc = (1 + X / n + r) ** 2 - 4 * X / n
    return (n / 2) * (1 + X / n + r - np.sqrt(np.clip(disc, 0, None)))


def coverage_at(cfree):
    """gvdW coverage theta at free phytate concentration cfree [M] (robust bisection)."""
    if cfree <= 0:
        return 0.0
    phi0 = K_H * (cfree * NA / 1e27) * D_A**2
    if phi0 <= 0:
        return 0.0

    def f(t):
        return t - phi0 * (1 - t) * np.exp(ALPHA * t - t / (1 - t))

    if f(1e-13) >= 0:  # root below the bracket -> Henry limit theta ~ phi0
        return phi0
    return brentq(f, 1e-13, 1 - 1e-12)


def gvdw_bound(X, cap):
    """Bound phytate per oleosin from the gvdW isotherm; `cap` scales coverage to the
    bound capacity. Bound PA [M] = cap*M_oleo*theta. Uses the module globals K_H, ALPHA."""
    out = np.empty_like(X)
    for i, x in enumerate(X):
        Lt = x * M_OLEO
        if Lt <= 0:
            out[i] = 0.0
            continue

        def resid(cfree, Lt=Lt):  # mass balance: total = free + bound
            return Lt - cfree - cap * M_OLEO * coverage_at(cfree)

        cf = brentq(resid, 0.0, Lt)
        out[i] = (Lt - cf) / M_OLEO
    return out


def per_injection(bound_cumulative):
    """Convert a cumulative bound curve to per-injection increments (prepend 0)."""
    return np.diff(np.concatenate([[0.0], bound_cumulative]))


def fit_onesite():
    """Reference: independent one-site fit of the thermogram (manuscript model)."""

    def model(X, dH_scale, K, n):
        return dH_scale * per_injection(onesite_bound(X, K, n))

    popt, _ = curve_fit(model, X, Q, p0=[-3000, 1.9e5, 0.32], maxfev=200000)
    return popt, np.sqrt(np.mean((Q - model(X, *popt)) ** 2)), model(X, *popt)


if __name__ == "__main__":
    (dH, Ka, n), rms, _ = fit_onesite()
    print("one-site (independent) reference fit of Fig S2C:")
    print(f"  Ka = {Ka:.2e} /M   n = {n:.3f}   RMS = {rms:.1f} uJ")
    print("  (published: Ka=1.87e5, n=0.322, dH=-56.4 kJ/mol)")
    print("\ngvdW binding: use the analytic Gouy-Chapman model in gc_gvdw.py")
    print("  (it sets itc_fit.K_H / itc_fit.ALPHA and reuses this titration engine).")
