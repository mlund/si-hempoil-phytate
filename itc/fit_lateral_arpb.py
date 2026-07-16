#!/usr/bin/env python3
"""Fit an ITC thermogram with the LC-AR-PB adsorption isotherm.

The calculation uses the same experimental and calorimetric assumptions as the
independent one-site fit:

* the titration dilution/overflow model comes from ``fit_thermogram.py``;
* free ligand is obtained from the per-injection mass balance;
* intrinsic affinity K, area per oleosin, and a constant molar binding enthalpy are
  fitted simultaneously;
* LC-AR-PB uses the peptide-only MC potential, MC phytate charge at the binding plane,
  nonlinear PB feedback, and the PB-matched residual lateral-correlation free energy.

The fit is diagnostic rather than evidence that all measured heat is monolayer
adsorption.  A constant dH assigns the heat of each net adsorption event to one fitted
enthalpy; coverage-dependent electrostatic/correlation enthalpies and aggregation heat
are not resolved separately.

Run from the repository root:

    uv run python itc/fit_lateral_arpb.py

The default is the pH-4 experiment at 15 mM ionic strength. ``--no-correlations``
gives the matched AR-PB baseline.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "zeta"))

import extract_thermograms as itc_data  # noqa: E402
import fit_thermogram as itc_fit  # noqa: E402
import lateral_arpb as lc  # noqa: E402

I_ITC = 0.015
AREA_OLEOSIN = 130.0 * 130.0 / 50.0  # physical grafting area, A^2 / oleosin

# Equilibrium depends on activity K*c_free.  This range covers the fit bounds
# (K=1e-6..1e10 /M) and the experimental free concentrations with saturated margins.
LOG_ACTIVITY_MIN = -18.0
LOG_ACTIVITY_MAX = 10.0
N_ACTIVITY = 701


@dataclass(frozen=True)
class CoverageIsotherm:
    """Tabulated LC-AR-PB coverage versus log10(K*c_free)."""

    pH: float
    ionic: float
    include_correlations: bool
    log_activity: np.ndarray
    gamma: np.ndarray

    def coverage(self, c_free_M: float | np.ndarray, K: float) -> float | np.ndarray:
        concentration = np.asarray(c_free_M)
        activity = K * concentration
        values = np.zeros_like(activity, dtype=float)
        positive = activity > 0.0
        if np.any(positive):
            values[positive] = np.interp(
                np.log10(activity[positive]),
                self.log_activity,
                self.gamma,
                left=self.gamma[0],
                right=self.gamma[-1],
            )
        return float(values) if values.ndim == 0 else values


@cache
def build_isotherm(
    pH: float = 4.0,
    ionic: float = I_ITC,
    include_correlations: bool = True,
) -> CoverageIsotherm:
    """Precompute the LC-AR-PB isotherm at fixed pH and ionic strength."""
    log_activity = np.linspace(LOG_ACTIVITY_MIN, LOG_ACTIVITY_MAX, N_ACTIVITY)
    activity = 10.0**log_activity
    gamma = np.array(
        [
            lc.lateral_state(
                pH,
                ionic,
                value,
                K=1.0,
                include_correlations=include_correlations,
            ).gamma
            for value in activity
        ]
    )
    # Numerical global minimization may produce sub-grid fluctuations at phase
    # coexistence.  Thermodynamic coverage cannot decrease with chemical activity;
    # retain the stable monotone envelope without smoothing a genuine jump.
    gamma = np.maximum.accumulate(gamma)
    return CoverageIsotherm(
        pH=pH,
        ionic=ionic,
        include_correlations=include_correlations,
        log_activity=log_activity,
        gamma=gamma,
    )


def bound_profile(
    titration: itc_data.Titration,
    isotherm: CoverageIsotherm,
    K: float,
    area_A2: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (bound PA, free PA) concentrations [M] for every injection.

    The balance at injection i is

        L_total = c_free + [oleosin] * area * Gamma(K*c_free).

    The oleosin and total-PA concentrations include cell overflow/dilution.
    """
    if K <= 0.0 or area_A2 <= 0.0:
        raise ValueError("K and area must be positive")

    oleosin_total, ligand_total = itc_fit._dilute(titration)
    bound = np.zeros_like(ligand_total)
    free = np.zeros_like(ligand_total)

    for index, (oleosin, total) in enumerate(zip(oleosin_total, ligand_total, strict=True)):
        if total <= 0.0:
            continue

        def balance(
            c_free: float,
            oleosin_M: float = oleosin,
            ligand_total_M: float = total,
        ) -> float:
            gamma = isotherm.coverage(c_free, K)
            return c_free + oleosin_M * area_A2 * gamma - ligand_total_M

        c_free = brentq(balance, 0.0, total, xtol=1e-15, rtol=1e-12)
        free[index] = c_free
        bound[index] = total - c_free

    return bound, free


def model_heat(
    titration: itc_data.Titration,
    isotherm: CoverageIsotherm,
    K: float,
    area_A2: float,
    dH_kJ_mol: float,
) -> np.ndarray:
    """LC-AR-PB heat per injection [uJ]."""
    bound, _free = bound_profile(titration, isotherm, K, area_A2)
    return itc_fit._heat_uJ(titration, bound, dH_kJ_mol)


@dataclass(frozen=True)
class FitResult:
    K: float
    area_A2: float
    dH_kJ_mol: float
    rms_uJ: float
    covariance: np.ndarray
    correlation: np.ndarray
    standard_error: np.ndarray
    predicted_heat_uJ: np.ndarray
    bound_M: np.ndarray
    free_M: np.ndarray
    success: bool
    message: str


def _covariance_from_jacobian(jacobian: np.ndarray, residual: np.ndarray) -> np.ndarray:
    degrees_of_freedom = max(1, residual.size - jacobian.shape[1])
    variance = float(residual @ residual) / degrees_of_freedom
    return np.linalg.pinv(jacobian.T @ jacobian) * variance


def fit_lateral_arpb(
    titration: itc_data.Titration,
    *,
    ionic: float = I_ITC,
    include_correlations: bool = True,
) -> FitResult:
    """Fit free (K, area, dH) to one ITC titration.

    Optimization variables are (log10 K, log10 area, dH) for positivity and numerical
    scaling.  Reported covariance and standard errors are transformed to physical
    (K, area, dH) coordinates.
    """
    isotherm = build_isotherm(titration.pH, ionic, include_correlations)

    def residual(parameters: np.ndarray) -> np.ndarray:
        log_K, log_area, dH = parameters
        predicted = model_heat(titration, isotherm, 10.0**log_K, 10.0**log_area, dH)
        return predicted - titration.heat_uJ

    solution = least_squares(
        residual,
        x0=np.array([2.0, np.log10(5000.0), -35.0]),
        bounds=(
            np.array([-6.0, np.log10(25.0), -1.0e4]),
            np.array([10.0, np.log10(1.0e6), -1.0e-3]),
        ),
        x_scale="jac",
        max_nfev=3000,
        xtol=1e-10,
        ftol=1e-10,
        gtol=1e-10,
    )

    log_K, log_area, dH = solution.x
    K, area = 10.0**log_K, 10.0**log_area
    predicted = model_heat(titration, isotherm, K, area, dH)
    bound, free = bound_profile(titration, isotherm, K, area)
    residual_at_solution = predicted - titration.heat_uJ
    covariance_log = _covariance_from_jacobian(solution.jac, residual_at_solution)

    # Local delta-method transformation from (log10 K, log10 A, dH) to physical units.
    transform = np.diag([np.log(10.0) * K, np.log(10.0) * area, 1.0])
    covariance = transform @ covariance_log @ transform
    standard_error = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
    scale = np.sqrt(np.clip(np.diag(covariance), 1e-300, None))
    correlation = covariance / np.outer(scale, scale)

    return FitResult(
        K=K,
        area_A2=area,
        dH_kJ_mol=float(dH),
        rms_uJ=float(np.sqrt(np.mean(residual_at_solution**2))),
        covariance=covariance,
        correlation=correlation,
        standard_error=standard_error,
        predicted_heat_uJ=predicted,
        bound_M=bound,
        free_M=free,
        success=bool(solution.success),
        message=solution.message,
    )


def plot_fit(
    titration: itc_data.Titration,
    fit: FitResult,
    output: Path,
    *,
    include_correlations: bool,
) -> None:
    """Plot experimental heats, LC-AR-PB fit, and the independent one-site fit."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    one_site_K, one_site_dH, one_site_n, _rms = itc_fit.fit_onesite(titration)
    one_site_heat = itc_fit.onesite_heat(
        titration,
        one_site_K,
        one_site_dH,
        one_site_n,
    )

    label = "LC–AR-PB" if include_correlations else "AR-PB (no residual correlations)"
    fig, axis = plt.subplots(figsize=(7, 5))
    axis.plot(
        titration.mole_ratio,
        titration.heat_uJ,
        "o",
        color="#333333",
        ms=5,
        label=f"experiment ({titration.source})",
    )
    axis.plot(
        titration.mole_ratio,
        fit.predicted_heat_uJ,
        "-",
        color="#1f4fd8",
        lw=2,
        label=label,
    )
    axis.plot(
        titration.mole_ratio,
        one_site_heat,
        "--",
        color="gray",
        lw=1.5,
        label="one-site (independent)",
    )
    axis.axhline(0.0, color="gray", lw=0.5)
    axis.set_xlabel("mole ratio PA : oleosin")
    axis.set_ylabel(r"heat per injection ($\mu$J)")
    axis.legend(fontsize=9, frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=140)
    plt.close(fig)


def _print_fit(
    titration: itc_data.Titration,
    fit: FitResult,
    include_correlations: bool,
    ionic: float,
) -> None:
    model = "LC-AR-PB" if include_correlations else "AR-PB baseline"
    errors = fit.standard_error
    print(f"{model} ITC fit: pH {titration.pH:g}, I={ionic:g} M, {titration.source}")
    print(f"  K    = {fit.K:.5g} +/- {errors[0]:.2g} /M")
    print(f"  area = {fit.area_A2:.5g} +/- {errors[1]:.2g} A^2/oleosin")
    print(
        f"         {fit.area_A2 / AREA_OLEOSIN:.1f}x physical grafting area "
        f"({AREA_OLEOSIN:.0f} A^2/oleosin)"
    )
    print(f"  dH   = {fit.dH_kJ_mol:.5g} +/- {errors[2]:.2g} kJ/mol")
    print(f"  RMS  = {fit.rms_uJ:.3g} uJ")
    print(f"  optimizer: {'converged' if fit.success else 'failed'} ({fit.message})")
    print("  parameter correlation matrix (K, area, dH):")
    for row in fit.correlation:
        print("   ", " ".join(f"{value:+.3f}" for value in row))
    max_bound_per_oleosin = float(np.max(fit.bound_M / itc_fit._dilute(titration)[0]))
    print(f"  maximum fitted bound PA/oleosin = {max_bound_per_oleosin:.3f}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pH", type=float, choices=lc.PHS, default=4.0)
    parser.add_argument("--ionic", type=float, default=I_ITC, help="ionic strength in M")
    parser.add_argument(
        "--no-correlations",
        action="store_true",
        help="fit the matched AR-PB baseline without residual lateral correlations",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=HERE / "fig_itc_lateral_arpb.png",
        help="output comparison figure",
    )
    args = parser.parse_args()

    include_correlations = not args.no_correlations
    titration = itc_data.load(args.pH)
    fit = fit_lateral_arpb(
        titration,
        ionic=args.ionic,
        include_correlations=include_correlations,
    )
    _print_fit(titration, fit, include_correlations, args.ionic)
    plot_fit(
        titration,
        fit,
        args.output,
        include_correlations=include_correlations,
    )
    print(f"  saved {args.output}")


if __name__ == "__main__":
    main()
