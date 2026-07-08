#!/usr/bin/env python3
"""NMR bound fraction of phytate from the simulated PMF.

The MC gives the single-phytate potential of mean force w(r) against a static oleosin
brush. NMR sees a population average over all phytate, so the two meet through a bound
fraction f_bound that carries the concentration dependence the PMF cannot:

  K_H     = int (e^{-w} - 1) dr             surface partition length [A]
  f_bound = Gamma a / n_B                   a = 6 phi_OB / d43 (droplet area per volume)

The areal coverage Gamma comes from a two-dimensional generalized-van-der-Waals isotherm
(Nordholm et al. 2018) whose lateral coefficient uses the Forsman correlation correction
(mc_shared.lateral_alpha): in the dilute limit Gamma -> K_H n_B, but at finite phytate the
like-charge lateral repulsion caps it. The bound fraction then scales with the oil-body
loading through the specific area a = 6 phi_OB / d43.
"""

import sys
from pathlib import Path

import numpy as np

# The Forsman lateral term (lateral_alpha) lives in the zeta shared module in ../zeta.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "zeta"))
import mc_shared as gv  # noqa: E402


def specific_area(phi_ob: float, d43_um: float) -> float:
    """a = 6 phi / d43 [1/A] for a dispersion of spheres (surface area per volume)."""
    return 6.0 * phi_ob / (d43_um * 1e4)  # d43 um -> A


def coverage_isotherm(K_H: float, n_B: float, d2: float, z: float, ionic: float) -> float:
    """Adsorbed-layer coverage theta = Gamma*d^2 from Nordholm's gvdW 2D isotherm.

    Nordholm et al. (2018), Generalized van der Waals Theory of Molecular Fluids in Bulk
    and at Surfaces (Elsevier). Solve the implicit hard-disk (Eq 7.43) / Lennard-Jones
    (Eq 7.54) isotherm

        Gamma = Gamma_HA (1 - Gamma d^2) exp( 2 a2 Gamma/kT - Gamma d^2/(1 - Gamma d^2) )

    with Gamma_HA = K_H n_B the Henry (zero-coverage) excess. In coverage form,
    theta = phi0 (1-theta) exp(alpha*theta - theta/(1-theta)), phi0 = K_H n_B d^2,
    alpha = 2 a2/(kT d^2) (alpha < 0 for the repulsive phytate-phytate laterals).
    The Forsman-corrected alpha depends on the coverage itself (R_c = 1/sqrt(pi Gamma),
    Gamma = theta/d^2), so it is re-evaluated each fixed-point step.
    """
    phi0 = K_H * n_B * d2  # Henry-limit coverage
    if phi0 <= 0:
        return 0.0
    theta = 0.0
    for _ in range(500):  # damped fixed point; theta stays in [0, 1)
        t = min(max(theta, 0.0), 0.999)
        alpha = gv.lateral_alpha(ionic, z, t / d2)  # coverage-dependent correlation hole
        new = phi0 * (1 - t) * np.exp(alpha * t - t / (1 - t))
        theta += 0.5 * (min(new, 0.999) - theta)
        if abs(new - t) < 1e-10:
            break
    return theta


def bound_fraction(a: float, K_H: float, pa_mM: float, d_A: float, z: float, ionic: float) -> float:
    """Bound fraction from the Nordholm gvdW isotherm: crowding + lateral repulsion cap it.

    f_bound = Gamma a / n_B = (theta/d^2) a / n_B; reduces to Henry (a K_H) as theta->0.
    """
    n_B = pa_mM * 1e-3 * 6.02214076e23 / 1e27  # mM -> 1/A^3
    d2 = d_A**2
    theta = coverage_isotherm(K_H, n_B, d2, z, ionic)
    Gamma = theta / d2
    return Gamma * a / n_B


# Measured 31P DeltaR = R1rho - R1 (s^-1) at OB 50%: (free PA_aq, OB-PA). Figure 9G.
NMR_DELTA_R = {4.0: (162.0, 559.0), 6.0: (128.0, 430.0), 10.0: (5.1, 52.0)}
