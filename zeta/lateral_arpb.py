#!/usr/bin/env python3
"""LC--AR--PB zeta model with peptide-only MC input and lateral correlations.

The model uses the following experimental constraints:

* the measured zero-PA zeta potential fixes the bare diffuse-layer charge;
* the intrinsic affinity K is fitted once per pH to the low-salt PA series and then
  held fixed across salt;
* the diluted analytical PA concentration is used as the free concentration.

The MC handoff and lateral thermodynamics are changed:

* the one-particle binding field is ``potential_peptide_mV`` at the binding plane,
  excluding the probe phytate potential;
* the phytate valence is the MC mean charge at the same distance Z, interpolated
  between simulated ionic strengths;
* coverage minimizes one surface grand potential containing Langmuir mixing,
  nonlinear planar PB feedback, and a residual lateral-correlation free energy;
* the residual is the Forsman-softened screened pair free energy minus the unsoftened
  Yukawa mean-field reference already represented by PB.  It is therefore not a
  second full electrostatic pair interaction.

The peptide field and the one-particle phytate charge are held fixed as coverage grows.
This is the intended dilute-probe approximation: it neglects restructuring of the brush
and charge regulation caused by neighbouring adsorbed phytates.  Langmuir site exclusion
is retained, so no additional hard-disk free energy is added.

Run directly to fit K, print the predicted Figure S3 values, and write
``fig_lateral_arpb.png`` plus ``model_zeta.csv``.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.special import xlogy

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mc_shared as gv  # noqa: E402

KT_MV = gv.KT_MV
LB = gv.LB
R_BETA = gv.R_CONTACT
I_ANCHOR = 0.005
DILUTION = gv.DILUTION
PHS = gv.PHS

# Geometric site ceiling for the Figure S3 LC--AR--PB calculation.
GAMMA_MAX = 1.0 / 11.0**2  # PA / A^2

# Scan definitions shared by the figure and the exported table, so the plotted curve and
# the tabulated numbers cannot drift apart.
PA_CURVE_MM = np.linspace(0.0, 7.7, 200)
NACL_CURVE_MM = np.linspace(1.0, 100.0, 200)
SALT_SCAN_PA_MM = 1.5  # fixed PA of the measured salt series
I_BACKGROUND = 0.003  # buffer/counterion background on top of added NaCl [M]


def _y_d(sigma: float | np.ndarray, ionic: float) -> float | np.ndarray:
    """Reduced nonlinear Gouy-Chapman potential e*psi_d/kT."""
    return 2.0 * np.arcsinh(2.0 * np.pi * LB * np.asarray(sigma) / gv.kappa(ionic))


def _psi_d(sigma: float, ionic: float) -> float:
    """Diffuse-layer potential in mV."""
    return float(_y_d(sigma, ionic) * KT_MV)


@cache
def sigma0(pH: float) -> float:
    """Bare surface charge [e/A^2], anchored to zero-PA zeta at 5 mM salt."""
    return gv.sigma_from_zeta(gv.ZETA_BARE[pH], I_ANCHOR)


def _profile_at_binding(
    table_name: str,
    column: str,
    pH: float,
    ionic: float,
) -> float:
    """Read a profile at R_BETA and log-interpolate between simulated salts."""
    lo, hi, weight = gv._ionic_bracket(ionic)

    def at(simulated_ionic: float) -> float:
        table = gv._table(table_name)
        rows = table[(table.pH == pH) & (table.ionic_strength == simulated_ionic)].sort_values("r")
        if rows.empty:
            raise ValueError(f"no MC profile for pH={pH:g}, I={simulated_ionic:g} M")
        return float(np.interp(R_BETA, rows.r.to_numpy(), rows[column].to_numpy()))

    lower = at(lo)
    return lower if weight == 0.0 else (1.0 - weight) * lower + weight * at(hi)


@cache
def peptide_potential(pH: float, ionic: float) -> float:
    """Peptide-only MC potential [mV] at the phytate binding plane."""
    return _profile_at_binding(
        "potential_profiles.csv",
        "potential_peptide_mV",
        pH,
        ionic,
    )


@cache
def phytate_valence(pH: float, ionic: float) -> float:
    """Magnitude of the MC phytate charge at the same binding-plane coordinate Z."""
    q_ino = _profile_at_binding("window_profiles.csv", "q_ino", pH, ionic)
    return -q_ino


def _yukawa_pair_integral(ionic: float, z: float) -> float:
    """Integral of the hard-core Yukawa pair energy, in kT*A^2.

    J_MF = integral_d^infinity 2*pi*r*z^2*l_B*exp(-kappa*r)/r dr.
    This is the linear mean-field reference and is not added separately to PB.
    """
    kappa = gv.kappa(ionic)
    return 2.0 * np.pi * z**2 * LB * np.exp(-kappa * gv.D_PA) / kappa


def _forsman_pair_integral(ionic: float, z: float, gamma: float) -> float:
    """Forsman-tangent screened pair integral, in kT*A^2.

    The correlation-hole scale uses the actual adsorbate density.  It is capped at two
    Debye lengths and floored at the phytate core.  No net-surface-charge estimate is
    mixed into the closure because nonlinear PB already treats that collective charge.
    """
    kappa = gv.kappa(ionic)
    spacing = 1.0 / np.sqrt(np.pi * max(gamma, 1e-30))
    r_c = max(gv.D_PA, min(2.0 / kappa, spacing))
    exponential = np.exp(-kappa * r_c)
    u_c = z**2 * LB * exponential / r_c
    slope = z**2 * LB * (-1.0 / r_c**2 - kappa / r_c) * exponential

    outer = 2.0 * np.pi * z**2 * LB * exponential / kappa
    inner = (
        2.0
        * np.pi
        * (
            slope * (r_c**3 - gv.D_PA**3) / 3.0
            - slope * r_c * (r_c**2 - gv.D_PA**2) / 2.0
            + u_c * (r_c**2 - gv.D_PA**2) / 2.0
        )
    )
    return float(inner + outer)


def residual_correlation_free_energy(
    gamma: float,
    ionic: float,
    z: float,
) -> float:
    """Residual lateral free energy beta*f_corr,res [1/A^2].

    The complete density-dependent expression is minimized, so derivatives of the
    correlation-hole radius are included automatically in the equilibrium condition.
    """
    if gamma <= 0.0:
        return 0.0
    delta_pair_integral = _forsman_pair_integral(ionic, z, gamma) - _yukawa_pair_integral(ionic, z)
    return 0.5 * gamma**2 * delta_pair_integral


def _pb_primitive(sigma: float | np.ndarray, ionic: float) -> float | np.ndarray:
    """Primitive of the reduced Grahame potential: dF/dsigma = y_d(sigma)."""
    sigma_array = np.asarray(sigma)
    scale = 2.0 * np.pi * LB / gv.kappa(ionic)
    return 2.0 * (
        sigma_array * np.arcsinh(scale * sigma_array)
        - np.sqrt(1.0 + (scale * sigma_array) ** 2) / scale
    )


def _incremental_pb_free_energy(gamma: float, pH: float, ionic: float, z: float) -> float:
    """PB free energy above the tangent at the measured bare-charge state [1/A^2].

    Subtracting the tangent makes its coverage derivative
    -z*[y_d(sigma_eff)-y_d(sigma0)], matching the incremental PB binding field while
    leaving the MC peptide-only field as the one-particle baseline.
    """
    bare = sigma0(pH)
    effective = bare - z * gamma
    return float(
        _pb_primitive(effective, ionic)
        - _pb_primitive(bare, ionic)
        - _y_d(bare, ionic) * (effective - bare)
    )


def surface_grand_potential(
    gamma: float,
    pH: float,
    ionic: float,
    c_L_M: float,
    K: float,
    *,
    include_correlations: bool = True,
    gamma_max: float = GAMMA_MAX,
) -> float:
    """Dimensionless surface grand potential beta*Omega/A [1/A^2]."""
    if not 0.0 <= gamma <= gamma_max:
        return np.inf
    if c_L_M <= 0.0 or K <= 0.0:
        return 0.0 if gamma == 0.0 else np.inf

    theta = gamma / gamma_max
    mixing = gamma_max * (xlogy(theta, theta) + xlogy(1.0 - theta, 1.0 - theta))
    activity = K * c_L_M
    z = phytate_valence(pH, ionic)
    y_peptide = peptide_potential(pH, ionic) / KT_MV

    omega = (
        mixing
        - gamma * np.log(activity)
        - z * y_peptide * gamma
        + _incremental_pb_free_energy(gamma, pH, ionic, z)
    )
    if include_correlations:
        omega += residual_correlation_free_energy(gamma, ionic, z)
    return float(omega)


def _global_minimum(objective, upper: float) -> float:
    """Locate the global minimum on [0, upper], allowing multiple local minima."""
    # Ideal mixing gives an interior minimum for every finite activity, even when its
    # coverage is many orders of magnitude below the close-packing scale.  Resolve both
    # boundary layers logarithmically so they cannot be mistaken for endpoint minima.
    tail = np.geomspace(1e-12, 1e-2, 161)
    fractions = np.unique(
        np.concatenate(
            (
                [0.0],
                tail,
                np.linspace(1e-2, 1.0 - 1e-2, 393),
                1.0 - tail,
                [1.0],
            )
        )
    )
    grid = upper * fractions
    values = np.array([objective(value) for value in grid])
    candidates: list[tuple[float, float]] = [(values[0], grid[0]), (values[-1], grid[-1])]

    for index in range(1, len(grid) - 1):
        if values[index] <= values[index - 1] and values[index] <= values[index + 1]:
            solution = minimize_scalar(
                objective,
                bounds=(grid[index - 1], grid[index + 1]),
                method="bounded",
                options={"xatol": 1e-13},
            )
            candidates.append((float(solution.fun), float(solution.x)))

    return min(candidates, key=lambda item: item[0])[1]


def _correlation_chemical_potential(gamma: float, ionic: float, z: float) -> float:
    """Numerical beta*mu_corr,res for diagnostics; equilibrium uses the full free energy."""
    if gamma <= 0.0:
        return 0.0
    step = max(1e-9, gamma * 1e-4)
    lower = max(0.0, gamma - step)
    upper = min(GAMMA_MAX, gamma + step)
    if upper == lower:
        return 0.0
    return (
        residual_correlation_free_energy(upper, ionic, z)
        - residual_correlation_free_energy(lower, ionic, z)
    ) / (upper - lower)


@dataclass(frozen=True)
class State:
    """Equilibrium surface state at one (pH, I, [PA]).

    Coverage ``gamma`` [PA/A^2] and its fraction ``theta`` of the Langmuir maximum;
    ``sigma_eff`` [e/A^2] is the bare charge less bound phytate charge; ``zeta`` and
    ``psi_binding_effective`` are in mV; ``valence`` is the phytate charge [e];
    ``beta_mu_corr`` is the residual correlation chemical potential [kT].
    """

    gamma: float
    theta: float
    sigma_eff: float
    zeta: float
    psi_binding_effective: float
    valence: float
    beta_mu_corr: float


def lateral_state(
    pH: float,
    ionic: float,
    c_L_M: float,
    K: float,
    *,
    include_correlations: bool = True,
    gamma_max: float = GAMMA_MAX,
) -> State:
    """Equilibrium coverage and zeta from global grand-potential minimization."""
    if c_L_M <= 0.0 or K <= 0.0:
        gamma = 0.0
    else:
        gamma = _global_minimum(
            lambda value: surface_grand_potential(
                value,
                pH,
                ionic,
                c_L_M,
                K,
                include_correlations=include_correlations,
                gamma_max=gamma_max,
            ),
            gamma_max,
        )

    z = phytate_valence(pH, ionic)
    bare = sigma0(pH)
    sigma_eff = bare - z * gamma
    delta_psi = _psi_d(sigma_eff, ionic) - _psi_d(bare, ionic)
    psi_binding_effective = peptide_potential(pH, ionic) + delta_psi
    beta_mu_corr = _correlation_chemical_potential(gamma, ionic, z) if include_correlations else 0.0
    return State(
        gamma=gamma,
        theta=gamma / gamma_max,
        sigma_eff=sigma_eff,
        zeta=_psi_d(sigma_eff, ionic),
        psi_binding_effective=psi_binding_effective,
        valence=z,
        beta_mu_corr=beta_mu_corr,
    )


def zeta_pa(
    pH: float,
    pa_mM: float,
    K: float,
    *,
    ionic: float = I_ANCHOR,
    include_correlations: bool = True,
) -> float:
    """Predicted zeta [mV], using the 1000-fold diluted PA concentration."""
    c_L_M = pa_mM * 1e-3 / DILUTION
    return lateral_state(
        pH,
        ionic,
        c_L_M,
        K,
        include_correlations=include_correlations,
    ).zeta


def _zeta_chi2(pH: float, K: float, include_correlations: bool = True) -> float:
    x, y, error = gv._fig6a(pH)
    residual = np.array(
        [
            zeta_pa(pH, value, K, include_correlations=include_correlations) - measured
            for value, measured in zip(x, y, strict=True)
        ]
    )
    return float(np.sum((residual / np.clip(error, 0.3, None)) ** 2))


def fit_K(pH: float, include_correlations: bool = True) -> float:
    """Fit one intrinsic K at a pH; return zero when adsorption is not identifiable."""
    solution = minimize_scalar(
        lambda log_K: _zeta_chi2(pH, 10.0**log_K, include_correlations),
        bounds=(-6.0, 8.0),
        method="bounded",
        options={"xatol": 1e-3},
    )
    K = 10.0**solution.x
    if _zeta_chi2(pH, K, include_correlations) > _zeta_chi2(pH, 0.0, include_correlations) - 1.0:
        return 0.0
    return K


def fit_K_by_pH(include_correlations: bool = True) -> dict[float, float]:
    """Return the fitted intrinsic K [1/M] for every pH in ``PHS``."""
    return {pH: fit_K(pH, include_correlations) for pH in PHS}


def scan_table(
    K_by_pH: dict[float, float],
    *,
    include_correlations: bool = True,
) -> pd.DataFrame:
    """Predicted zeta for both scans, on the plotted grid and at the measured points.

    ``sampling`` separates the dense curve used for the figure from the rows that line up
    with measurements, so one file serves both the plot and the SI table.
    """
    rows = []

    def add(scan: str, sampling: str, pH: float, pa_mM: float, nacl_mM: float) -> None:
        ionic = I_ANCHOR if scan == "PA" else nacl_mM / 1000.0 + I_BACKGROUND
        state = lateral_state(
            pH,
            ionic,
            pa_mM * 1e-3 / DILUTION,
            K_by_pH[pH],
            include_correlations=include_correlations,
        )
        rows.append(
            {
                "scan": scan,
                "sampling": sampling,
                "pH": pH,
                "PA_mM": pa_mM,
                "NaCl_mM": nacl_mM,
                "ionic_M": ionic,
                "K_per_M": K_by_pH[pH],
                "z_MC": state.valence,
                "gamma_per_A2": state.gamma,
                "theta": state.theta,
                "sigma_eff_e_per_A2": state.sigma_eff,
                "beta_mu_corr": state.beta_mu_corr,
                "zeta_mV": state.zeta,
            }
        )

    for pH in PHS:
        for pa_mM in PA_CURVE_MM:
            add("PA", "curve", pH, float(pa_mM), 0.0)
        for pa_mM in gv._fig6a(pH)[0]:
            add("PA", "measured", pH, float(pa_mM), 0.0)
        for nacl_mM in NACL_CURVE_MM:
            add("NaCl", "curve", pH, SALT_SCAN_PA_MM, float(nacl_mM))
        for nacl_mM in gv._fig6c(pH, SALT_SCAN_PA_MM)[0]:
            add("NaCl", "measured", pH, SALT_SCAN_PA_MM, float(nacl_mM))

    return pd.DataFrame(rows)


def plot_figure_s3(
    table: pd.DataFrame,
    output: Path,
) -> None:
    """Write the two-panel zeta comparison used for Figure S3."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {4.0: "#1f4fd8", 6.0: "#e08a00", 10.0: "#2ca02c"}
    fig, (axis_pa, axis_salt) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    curve = table[table.sampling == "curve"]
    for pH in PHS:
        measured_pa, measured, error = gv._fig6a(pH)
        color = colors[pH]
        axis_pa.errorbar(
            measured_pa,
            measured,
            yerr=error,
            fmt="x",
            color=color,
            ms=6,
            capsize=2,
            label=f"pH {pH:g}",
        )
        model = curve[(curve.scan == "PA") & (curve.pH == pH)].sort_values("PA_mM")
        axis_pa.plot(model.PA_mM, model.zeta_mV, "-", color=color, lw=1.8)

    axis_pa.axhline(0.0, color="gray", lw=0.5)
    axis_pa.set_ylim(-75.0, 40.0)
    axis_pa.set_xlabel("[PA] (mM)")
    axis_pa.set_ylabel(r"$\zeta$ (mV)")

    for pH in PHS:
        measured_salt, measured, error = gv._fig6c(pH, SALT_SCAN_PA_MM)
        color = colors[pH]
        axis_salt.errorbar(
            measured_salt,
            measured,
            yerr=error,
            fmt="x",
            color=color,
            ms=6,
            capsize=2,
            label=f"pH {pH:g}",
        )
        model = curve[(curve.scan == "NaCl") & (curve.pH == pH)].sort_values("NaCl_mM")
        axis_salt.plot(model.NaCl_mM, model.zeta_mV, "-", color=color, lw=1.8)

    axis_salt.axhline(0.0, color="gray", lw=0.5)
    axis_salt.set_xlabel("[NaCl] (mM)")
    axis_salt.legend(fontsize=9, frameon=False, ncol=3)
    fig.tight_layout(pad=1.0)
    fig.savefig(output, dpi=140)
    plt.close(fig)


def _print_predictions(
    K_by_pH: dict[float, float],
    table: pd.DataFrame,
    include_correlations: bool,
) -> None:
    label = "LC-AR-PB" if include_correlations else "AR-PB baseline"
    print(label)
    print(f"  binding plane Z = {R_BETA:g} A; peptide-only MC potential; MC q_PA(Z)")
    print(f"  Gamma_max = {GAMMA_MAX:.5f} /A^2; dilution = {DILUTION:g}x")
    for pH in PHS:
        print(
            f"  pH {pH:g}: K = {K_by_pH[pH]:.4g} /M, "
            f"chi2 = {_zeta_chi2(pH, K_by_pH[pH], include_correlations):.1f}"
        )

    measured = table[table.sampling == "measured"]

    print("\nFigure S3A predictions")
    print(f"  {'pH':>4} {'PA/mM':>8} {'z_MC':>7} {'theta':>8} {'mu_corr/kT':>11} {'zeta/mV':>9}")
    for row in measured[measured.scan == "PA"].itertuples():
        print(
            f"  {row.pH:>4g} {row.PA_mM:>8.2f} {row.z_MC:>7.3f} {row.theta:>8.4f} "
            f"{row.beta_mu_corr:>11.3f} {row.zeta_mV:>9.2f}"
        )

    print(f"\nFigure S3B predictions at {SALT_SCAN_PA_MM:g} mM PA")
    print(f"  {'pH':>4} {'NaCl/mM':>9} {'z_MC':>7} {'theta':>8} {'zeta/mV':>9}")
    for row in measured[measured.scan == "NaCl"].itertuples():
        print(
            f"  {row.pH:>4g} {row.NaCl_mM:>9.0f} {row.z_MC:>7.3f} {row.theta:>8.4f} "
            f"{row.zeta_mV:>9.2f}"
        )


def main() -> None:
    """Fit K per pH, print predicted zeta values and write the figure and CSV."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "fig_lateral_arpb.png",
        help="output Figure S3-style PNG",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path(__file__).resolve().parent / "model_zeta.csv",
        help="output CSV of the predicted [PA] and [NaCl] zeta scans",
    )
    parser.add_argument(
        "--no-correlations",
        action="store_true",
        help="omit the residual correlation free energy for a matched AR-PB baseline",
    )
    args = parser.parse_args()
    include_correlations = not args.no_correlations
    K_by_pH = fit_K_by_pH(include_correlations)
    table = scan_table(K_by_pH, include_correlations=include_correlations)
    _print_predictions(K_by_pH, table, include_correlations)
    table.to_csv(args.csv, index=False)
    plot_figure_s3(table, args.output)
    print(f"\nsaved {args.output}\nsaved {args.csv}")


if __name__ == "__main__":
    main()
