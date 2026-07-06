#!/usr/bin/env python3
"""Global gvdW model of the zeta-potential vs (pH, [PA], [NaCl]) surface (Fig 6).

Two-potential, self-consistent model. Phytate adsorption and the measured zeta are
governed by *different* surface potentials, which is what lets the model work near
the OB isoelectric point:

  * DRIVING (adsorption)  <- the local outer cationic PATCH potential phi_bind, i.e.
    the dipolar effective charge sigma_eff = sigma_c cosh(kappa d) - sigma_a
    (positive even when the NET charge is ~0), renormalised to the simulated
    phi_pep(binding site). The adsorbed anion neutralises this local patch.
  * ZETA (far-field)      <- the NET shear-plane charge sigma_net (from the bare-OB
    zeta) plus the adsorbed PA charge z_PA * Gamma.

Both respond to the same adsorbed population Gamma, solved self-consistently, with
K_H(ionic) screening the binding (giving the NaCl dependence of Fig 6C).

sigma_net(pH) and phi_bind(pH) are fixed inputs; the only fitted quantity is the
adsorption capacity. Caveat: decoupling driving from zeta is a mean-field
approximation (the dipolar plate ties them to one sigma_eff); phi_bind is taken
from the simulated brush.
"""

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

# The analytic Gouy-Chapman binding model is shared with the ITC fit (../itc).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "itc"))
import gc_gvdw as gc  # noqa: E402

LB, NA, KT_MV, D_A = 7.1528, 6.02214076e23, 25.693, 10.6
I_BUF = 0.005
Z_PA = {4.0: 5.94, 6.0: 7.08, 10.0: 10.54}
ZETA_BARE = {4.0: 38.0, 6.0: -12.0, 10.0: -56.0}  # net far-field, bare OB [mV]
DIPOLE = {4.0: 31.6, 6.0: 34.7, 10.0: 21.8}  # peptide dipole moment (MC multipole) [e.A]
GRAFT = 3e-3  # chains/A^2


def kappa(ionic):
    return np.sqrt(8 * np.pi * LB * ionic * NA / 1e27)


def sigma_from_zeta(z_mV, ionic):
    return (kappa(ionic) / (2 * np.pi * LB)) * np.sinh((z_mV / KT_MV) / 2)


def zeta_from_sigma(sig, ionic):
    return 2 * np.arcsinh(2 * np.pi * LB * sig / kappa(ionic)) * KT_MV


SIG_NET = {pH: sigma_from_zeta(z, I_BUF) for pH, z in ZETA_BARE.items()}


def _theta(K_H, alpha, nB):
    if K_H <= 0:
        return 0.0
    phi0 = K_H * nB * D_A**2

    def f(t):
        return t - phi0 * (1 - t) * np.exp(alpha * t - t / (1 - t))

    return brentq(f, 1e-13, 1 - 1e-12) if f(1e-13) < 0 else min(phi0, 1 - 1e-12)


def gamma_sc(PA_mM, pH, ionic, a_patch, a_cap):
    """Self-consistent adsorbed PA. Both the DRIVING patch and the CAPACITY scale with
    the MC dipole p(pH) (its cationic-outer character), so both are pH-dependent:
      patch charge = a_patch * p(pH) * GRAFT   (positive; added to sigma_net)
      capacity     = a_cap   * p(pH) * GRAFT   (max Gamma = bindable sites)
    The adsorbed anion neutralises the patch, so binding self-limits."""
    if PA_mM <= 0:
        return 0.0
    nB = PA_mM * 1e-3 * NA / 1e27
    z = Z_PA[pH]
    drive0 = SIG_NET[pH] + a_patch * DIPOLE[pH] * GRAFT
    cap = a_cap * DIPOLE[pH] * GRAFT
    G = 0.0
    for _ in range(400):
        sig_drive = drive0 - z * G  # patch neutralised by adsorbed anion
        K_H = gc.K_H_gc(sig_drive, -z, ionic) if sig_drive > 0 else 0.0
        alpha = gc.alpha_gc(abs(sig_drive) + 1e-9, z, ionic)
        Gnew = min(_theta(K_H, alpha, nB) / D_A**2, cap)
        G += 0.3 * (Gnew - G)
    return G


def zeta_model(pH, PA_mM, NaCl_mM, a_patch, a_cap):
    ionic = I_BUF + NaCl_mM * 1e-3
    G = gamma_sc(PA_mM, pH, ionic, a_patch, a_cap)
    return zeta_from_sigma(SIG_NET[pH] - Z_PA[pH] * G, ionic)  # zeta from NET + adsorbed PA


# --- digitised Fig 6A (zeta vs [PA]) and 6C (zeta vs [NaCl]) ---  (pH, PA_mM, NaCl_mM, zeta_mV)
DATA = [
    *[(4.0, pa, 0, z) for pa, z in [(0.15, 14), (0.75, 8), (1.5, 4), (3, 1), (4.5, 1), (7.6, -4)]],
    *[
        (6.0, pa, 0, z)
        for pa, z in [(0.15, -20), (0.75, -24), (1.5, -25), (3, -26), (4.5, -25), (7.6, -22)]
    ],
    *[
        (10.0, pa, 0, z)
        for pa, z in [(0.15, -68), (0.75, -63), (1.5, -61), (3, -63), (4.5, -63), (7.6, -66)]
    ],
    *[
        (4.0, 0, n, z)
        for n, z in [(5, 35), (10, 25), (20, 15), (40, 10), (60, 8), (80, 6), (100, 5)]
    ],
    *[
        (4.0, 1.5, n, z)
        for n, z in [(5, 8), (10, 10), (20, 9), (40, 7), (60, 6), (80, 5), (100, 5)]
    ],
    *[(6.0, 0, n, z) for n, z in [(5, -13), (20, -8), (40, -7), (60, -6), (100, -5)]],
    *[(6.0, 1.5, n, z) for n, z in [(5, -13), (20, -9), (40, -7), (100, -5)]],
    *[
        (10.0, 0, n, z)
        for n, z in [(5, -50), (10, -42), (20, -35), (40, -28), (60, -26), (80, -24), (100, -23)]
    ],
    *[(10.0, 1.5, n, z) for n, z in [(5, -50), (20, -35), (60, -26), (100, -23)]],
]


def fit_global():
    """Global fit over ALL pH: two scale factors (a_patch, a_cap); the pH-dependence of
    both the patch driving and the capacity comes from the MC dipole p(pH)."""

    def resid(p):
        a_patch, a_cap = 10 ** p[0], 10 ** p[1]
        return [zeta_model(pH, pa, n, a_patch, a_cap) - z for pH, pa, n, z in DATA]

    sol = least_squares(resid, [np.log10(1e-2), np.log10(1e-2)])
    a_patch, a_cap = 10 ** sol.x[0], 10 ** sol.x[1]
    return a_patch, a_cap, np.sqrt(np.mean(np.array(resid(sol.x)) ** 2))


if __name__ == "__main__":
    a_patch, a_cap, rms = fit_global()
    print("two-potential zeta model; patch & capacity BOTH pH-dependent via MC dipole p(pH)")
    print(f"a_patch = {a_patch:.2e}   a_cap = {a_cap:.2e}   RMS(all pH) = {rms:.1f} mV")
    print("pH-dependent patch/capacity (a * p(pH) * graft):")
    for pH in (4.0, 6.0, 10.0):
        patch = a_patch * DIPOLE[pH] * GRAFT
        cap = a_cap * DIPOLE[pH] * GRAFT
        print(f"  pH{pH:>4}: p={DIPOLE[pH]:.1f} eA -> patch={patch:.2e}  cap={cap:.2e} e/A^2")
    print(f"\n{'pH':>4} {'PA':>4} {'NaCl':>5} {'zeta data':>9} {'zeta model':>10}")
    for pH, pa, n, z in DATA:
        print(f"{pH:>4} {pa:>4} {n:>5} {z:>9} {zeta_model(pH, pa, n, a_patch, a_cap):>10.1f}")
