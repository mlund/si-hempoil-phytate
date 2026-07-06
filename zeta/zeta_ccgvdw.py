#!/usr/bin/env python3
"""Correlation-corrected gvdW model of the oil-body zeta-potential (Fig 6).

One MC-grounded model connects the calorimetry (ITC), the electrophoretic
zeta-potential (Fig 6A/C) and the NMR by extrapolating the infinite-dilution
single-phytate PMF to finite phytate concentration:

  * Henry limit          K_H = int (e^{-w(r)/kT} - 1) dr  from the MC PMF w(r).
  * finite concentration a 2D generalised-van-der-Waals isotherm with a lateral
    coefficient a2 that caps the coverage.
  * correlations         a2 is corrected for multivalent ion-ion correlations
    following Forsman (J. Phys. Chem. B 2004, 108, 9236): below the neutralising
    nearest-neighbour spacing R_c = sqrt(z / (pi sigma)) the divergent Coulomb
    repulsion is replaced by its tangent. This lets the surface over-charge --
    the pH-4 zeta reversal that mean-field Poisson-Boltzmann structurally cannot
    produce. R_c is capped at 2/kappa because correlations are screened beyond
    roughly two Debye lengths, which also removes the sigma -> 0 singularity.
  * self-consistency     the adsorbed anion layer shallows its own well and
    lowers the net surface charge (both fed back to convergence).
  * electrostatics       zeta is read from a nonlinear Poisson-Boltzmann solve of
    the adsorbed charge on the MC potential profile, so the double layer can
    invert.
  * regulated charge     the phytate valence z(r) is read from the MC (charge.csv),
    already regulated by the peptide field; no separate charge-regulation term is
    added (that would double-count -- the single-PA PMF already contains it).

The measured zeta follows a 1000-fold dilution of the emulsion, so the
zeta-relevant coverage is the isotherm re-equilibrated at [PA]/1000: the same
adsorption the ITC sees at the preparation concentration, read lower on the
isotherm. The only fitted quantity is the shear-plane position r_shear.

Caveats: pH-4 reversal is captured in trend and magnitude but not fully (the
coverage feedback uses the linear double layer while zeta uses the nonlinear one,
and the R_c = 2/kappa cap deliberately limits the softening); the pH-10 bare
anchor carries a ~10 mV offset; and Fig 6A and 6C disagree by ~12 mV at pH 6.
"""

from __future__ import annotations

from functools import cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import constants
from scipy.integrate import solve_bvp
from scipy.linalg import solve_banded
from scipy.optimize import brentq, minimize_scalar

# kT/e at the simulation temperature, in mV (zeta and MC potentials are in mV).
_T = 298.15
KT_MV = constants.k * _T / constants.e * 1e3
N_A = constants.N_A
# Bjerrum length of the MC water model [A]; must match Faunus (298.15 K, eps_r 78.4).
LB = 7.1528
# Phytate hard-disk diameter [A] ~ 2 a_H, sets the 2D excluded area and a2 range.
D_PA = 10.6

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


def phytate_valence(pH: float, ionic: float) -> float:
    """Adsorption-weighted phytate valence |z| from the MC (charge.csv, INO selection).

    The charge is read where the phytate actually binds (Boltzmann-weighted over the
    PMF well), so it is the peptide-field-regulated value, not the bulk one.
    """
    w = _pmf(pH, ionic)
    weight = np.expm1(-w.pmf_kT.to_numpy()).clip(0)
    if weight.sum() == 0:
        return 6.0
    c = _table("charge.csv")
    c = c[(c.ionic_strength == ionic) & (c.pH == pH) & (c.selection == "atomtype INO")]
    z_of_r = np.interp(w.r.to_numpy(), c.r.to_numpy(), -c.charge.to_numpy())
    return float(np.average(z_of_r, weights=weight))


def lateral_alpha(ionic: float, z: float, sigma_net: float) -> float:
    """Forsman-corrected lateral coefficient alpha = 2 a2 / (kT d^2).

    The screened-Coulomb pair repulsion is replaced by its tangent below the
    neutralising nearest-neighbour spacing R_c = sqrt(z / (pi |sigma_net|)),
    capped at 2/kappa. The tangent piece is integrated analytically.
    """
    k = kappa(ionic)
    r_c = min(np.sqrt(z / (np.pi * max(abs(sigma_net), 1e-9))), 2.0 / k)
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


def _coverage(k_henry: float, pa_M: float, alpha: float) -> float:
    """Surface excess Gamma [PA/A^2] from the 2D gvdW isotherm."""
    if k_henry <= 0 or pa_M <= 0:
        return 0.0
    phi0 = k_henry * pa_M * N_A / 1e27 * D_PA**2

    def residual(theta):
        return theta - phi0 * (1 - theta) * np.exp(alpha * theta - theta / (1 - theta))

    theta = brentq(residual, 1e-13, 1 - 1e-12) if residual(1e-13) < 0 else min(phi0, 1 - 1e-12)
    return theta / D_PA**2


def _pmf(pH: float, ionic: float) -> pd.DataFrame:
    w = _table("all_pmfs.csv")
    w = w[(w.ionic_strength == ionic) & (w.pH == pH)].sort_values("r")
    return w[w.r < 70]


def _potential_grid(pH: float, ionic: float, r: np.ndarray) -> np.ndarray:
    p = _table("potential.csv")
    p = p[(p.ionic_strength == ionic) & (p.pH == pH)].sort_values("r")
    return np.interp(r, p.r.to_numpy(), p.potential_all_mV.to_numpy() / KT_MV, right=0)


def _linear_shift(rho: np.ndarray, y0: np.ndarray, ionic: float, r: np.ndarray) -> np.ndarray:
    """Debye-Huckel response u to a fixed charge `rho`, linearised about the MC potential y0."""
    h = r[1] - r[0]
    screen = kappa(ionic) ** 2 * np.cosh(np.clip(y0, -20, 20))
    ab = np.zeros((3, r.size))
    ab[1] = -2 / h**2 - screen
    ab[0, 1:] = 1 / h**2
    ab[2, :-1] = 1 / h**2
    ab[0, 1] = 2 / h**2  # u'(0) = 0
    ab[1, -1], ab[0, -1], ab[2, -2] = 1, 0, 0  # u(L) = 0
    b = -4 * np.pi * LB * rho
    b[-1] = 0
    return solve_banded((1, 1), ab, b)


def adsorb(pH: float, ionic: float, pa_M: float, zeta_bare: float, gamma0: float = 0.0):
    """Self-consistently adsorb phytate and return (Gamma, r_grid, potential_shift u(r)).

    Coverage and the double-layer shift are relaxed together: the adsorbed anions
    shallow the well (w_eff = w - z u) and neutralise the surface (sigma_net drives
    R_c). The feedback loop uses the cheap linear double layer; the returned u(r) is
    then recomputed from the nonlinear Poisson-Boltzmann equation so that zeta can
    invert.
    """
    w = _pmf(pH, ionic)
    r_w, w0 = w.r.to_numpy(), w.pmf_kT.to_numpy()
    z = phytate_valence(pH, ionic)
    sigma_bare = sigma_from_zeta(zeta_bare, ionic)
    r = np.linspace(0, 180, 400)
    y0 = _potential_grid(pH, ionic, r)
    u = np.zeros_like(r)
    gamma, well_shape = gamma0, np.zeros_like(r_w)
    for _ in range(300):
        alpha = lateral_alpha(ionic, z, sigma_bare - z * gamma)
        w_eff = w0 - z * np.interp(r_w, r, u)
        k_henry = np.trapezoid(np.expm1(-w_eff), r_w)
        gamma += 0.15 * (_coverage(k_henry, pa_M, alpha) - gamma)
        excess = np.expm1(-w_eff).clip(0)
        if excess.sum() == 0:
            return 0.0, r, np.zeros_like(r)
        well_shape = excess / np.trapezoid(excess, r_w)
        rho = -z * gamma * np.interp(r, r_w, well_shape, left=0, right=0)
        u += 0.15 * (_linear_shift(rho, y0, ionic, r) - u)

    return gamma, *_nonlinear_shift(pH, ionic, z, gamma, well_shape, r_w, r, y0, u)


def _nonlinear_shift(pH, ionic, z, gamma, well_shape, r_w, r, y0, u_guess):
    """Nonlinear Poisson-Boltzmann shift u = y - y0 for the converged adsorbed charge."""
    rg = np.linspace(0, 150, 300)
    y0_of = lambda rr: np.interp(rr, r, y0)  # noqa: E731 -- tiny closure for solve_bvp
    source = 4 * np.pi * LB * z * gamma

    def rhs(rr, U):
        screen = kappa(ionic) ** 2 * (
            np.sinh(np.clip(y0_of(rr) + U[0], -30, 30)) - np.sinh(np.clip(y0_of(rr), -30, 30))
        )
        return np.vstack([U[1], screen + source * np.interp(rr, r_w, well_shape, left=0, right=0)])

    guess = np.vstack([np.interp(rg, r, u_guess), np.gradient(np.interp(rg, r, u_guess), rg)])
    sol = solve_bvp(rhs, lambda a, b: np.array([a[1], b[0]]), rg, guess, max_nodes=6000, tol=1e-4)
    return rg, sol.sol(rg)[0] if sol.success else np.interp(rg, r, u_guess)


# Bare-OB zeta at low salt (Fig 6A, ~0 PA) used to anchor the far-field offset [mV].
ZETA_BARE = {4.0: 38.0, 6.0: -12.0, 10.0: -56.0}
DILUTION = 1000.0  # emulsion is diluted 1000x before the zeta measurement
IONIC_MC = (0.005, 0.03, 0.15)  # MC ionic strengths available

# Digitised Fig 6A: pH -> [([PA] mM, zeta mV)] at low salt.
FIG6A = {
    4.0: [(0.15, 14), (0.75, 8), (1.5, 4), (3, 1), (4.5, 1), (7.6, -4)],
    6.0: [(0.15, -20), (0.75, -24), (1.5, -25), (3, -26), (4.5, -25), (7.6, -22)],
    10.0: [(0.15, -68), (0.75, -63), (1.5, -61), (3, -63), (4.5, -63), (7.6, -66)],
}
# Digitised Fig 6C: pH -> {"bare": [(NaCl mM, zeta)], "pa": [...] at 1.5 mM PA}.
FIG6C = {
    4.0: {
        "bare": [(5, 35), (10, 25), (20, 15), (40, 10), (60, 8), (80, 6), (100, 5)],
        "pa": [(5, 8), (10, 10), (20, 9), (40, 7), (60, 6), (80, 5), (100, 5)],
    },
    6.0: {
        "bare": [(5, -13), (20, -8), (40, -7), (60, -6), (100, -5)],
        "pa": [(5, -13), (20, -9), (40, -7), (100, -5)],
    },
    10.0: {
        "bare": [(5, -50), (10, -42), (20, -35), (40, -28), (60, -26), (80, -24), (100, -23)],
        "pa": [(5, -50), (20, -35), (60, -26), (100, -23)],
    },
}


def _profiles_6a():
    """Converged u(r) profiles for every Fig 6A point (continuation in [PA])."""
    out = {}
    for pH, points in FIG6A.items():
        gamma = 0.0
        for pa_mM, _ in points:
            gamma, rg, u = adsorb(pH, 0.005, pa_mM * 1e-3 / DILUTION, ZETA_BARE[pH], gamma0=gamma)
            out[(pH, pa_mM)] = (rg, u)
    return out


def _profiles_6c():
    """Converged u(r) profiles at 1.5 mM PA for each pH and MC ionic strength."""
    return {
        (pH, ionic): adsorb(pH, ionic, 1.5e-3 / DILUTION, ZETA_BARE[pH])[1:]
        for pH in FIG6C
        for ionic in IONIC_MC
    }


def _bare_6c(pH, nacl_mM):
    x, y = zip(*FIG6C[pH]["bare"], strict=True)
    return np.interp(nacl_mM, x, y)


def _zeta_from_profile(rg, u, r_shear):
    return np.interp(r_shear, rg, u) * KT_MV


def fit_global():
    """Fit the single shear-plane parameter to Fig 6A + 6C; return (r_shear, rms, curves)."""
    pa = _profiles_6a()
    salt = _profiles_6c()

    def zeta_6a(pH, pa_mM, r_shear):
        return ZETA_BARE[pH] + _zeta_from_profile(*pa[(pH, pa_mM)], r_shear)

    def zeta_6c(pH, nacl_mM, r_shear):
        ionic = nacl_mM / 1000 + 0.003
        shifts = [_zeta_from_profile(*salt[(pH, i)], r_shear) for i in IONIC_MC]
        return _bare_6c(pH, nacl_mM) + np.interp(np.log(ionic), np.log(IONIC_MC), shifts)

    def residuals(r_shear):
        e = [zeta_6a(pH, c, r_shear) - z for pH, pts in FIG6A.items() for c, z in pts]
        e += [zeta_6c(pH, n, r_shear) - z for pH in FIG6C for n, z in FIG6C[pH]["pa"]]
        return e

    sol = minimize_scalar(
        lambda rs: np.sqrt(np.mean(np.square(residuals(rs)))), bounds=(20, 40), method="bounded"
    )
    return sol.x, sol.fun, (zeta_6a, zeta_6c)


def plot_fig6(r_shear, rms, zeta_6a, zeta_6c, path):
    """Save the global Fig 6 comparison (model curves over the digitised data) to `path`."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {4.0: "#1f4fd8", 6.0: "#e08a00", 10.0: "#2ca02c"}
    fig = plt.figure(figsize=(13, 7))
    grid = fig.add_gridspec(2, 3)
    ax_a = fig.add_subplot(grid[0, :])
    ax_c = [fig.add_subplot(grid[1, i]) for i in range(3)]

    for pH, points in FIG6A.items():
        pa, z = zip(*points, strict=True)
        ax_a.plot(pa, z, "o", color=colors[pH], ms=8, label=f"pH {pH:g}")
        ax_a.plot(pa, [zeta_6a(pH, c, r_shear) for c in pa], "-", color=colors[pH], lw=2)
    ax_a.axhline(0, color="gray", lw=0.5)
    ax_a.set_xlabel("[PA] (mM)")
    ax_a.set_ylabel(r"$\zeta$ (mV)")
    ax_a.set_title(r"A) $\zeta$ vs [PA] (low salt)")
    ax_a.legend(ncol=3, fontsize=9)

    nacl = np.linspace(5, 100, 40)
    for ax, pH in zip(ax_c, (4.0, 6.0, 10.0), strict=True):
        xb, zb = zip(*FIG6C[pH]["bare"], strict=True)
        xr, zr = zip(*FIG6C[pH]["pa"], strict=True)
        ax.plot(xb, zb, "ks", ms=5, label="0 PA")
        ax.plot(xr, zr, "r^", ms=5, label="1.5 mM PA")
        ax.plot(nacl, [_bare_6c(pH, n) for n in nacl], "k-", lw=1)
        ax.plot(nacl, [zeta_6c(pH, n, r_shear) for n in nacl], "r-", lw=2)
        ax.axhline(0, color="gray", lw=0.5)
        ax.set_title(f"C) pH {pH:g}")
        ax.set_xlabel("[NaCl] (mM)")
    ax_c[0].set_ylabel(r"$\zeta$ (mV)")
    ax_c[0].legend(fontsize=8)

    fig.suptitle(f"Fig 6 -- correlation-corrected gvdW; r_shear={r_shear:.0f} A, RMS={rms:.1f} mV")
    fig.tight_layout()
    fig.savefig(path, dpi=140)


if __name__ == "__main__":
    r_shear, rms, (zeta_6a, zeta_6c) = fit_global()
    print(
        f"correlation-corrected gvdW, global Fig 6 fit: "
        f"r_shear = {r_shear:.1f} A, RMS = {rms:.1f} mV (33 points)\n"
    )
    print(f"Fig 6A -- zeta vs [PA], low salt:\n  {'pH':>4} {'[PA]/mM':>8} {'model':>7} {'data':>6}")
    for pH, points in FIG6A.items():
        for pa_mM, z in points:
            print(f"  {pH:>4} {pa_mM:>8} {zeta_6a(pH, pa_mM, r_shear):>+7.1f} {z:>+6}")
    print("\nFig 6C -- zeta at 1.5 mM PA vs [NaCl]:")
    print(f"  {'pH':>4} {'NaCl/mM':>8} {'model':>7} {'data':>6}")
    for pH, series in FIG6C.items():
        for nacl_mM, z in series["pa"]:
            print(f"  {pH:>4} {nacl_mM:>8} {zeta_6c(pH, nacl_mM, r_shear):>+7.1f} {z:>+6}")

    figure = Path(__file__).parent / "fig6_ccgvdw.png"
    plot_fig6(r_shear, rms, zeta_6a, zeta_6c, figure)
    print(f"\nsaved {figure.name}")
