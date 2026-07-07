#!/usr/bin/env python3
"""Bridge the single-phytate PMF to the concentration-dependent NMR observable.

The MC gives per-molecule, position-resolved quantities against a *static* oleosin
brush: the PMF w(r), the wall alignment S_wall(r), and the relative rotational
mobility D(r)/D0. NMR instead sees a *population* average over all phytate,

    <A>_obs = f_bound * A_bound + (1 - f_bound) * A_free,

so the simulation and the experiment meet only once we (i) collapse each r-profile
into an intensive per-bound-molecule value A_bound, and (ii) supply a bound
fraction f_bound that carries the concentration dependence the PMF cannot.

Normalisation (why K, not Z): with the fitted PMF referenced to bulk (w -> 0),
exp(-w) = rho(r)/rho_bulk is the local enrichment. Averaging A over the *adsorbed
excess* (Gibbs surface excess) is box-independent, unlike dividing by the full
Z = int exp(-w) dr (which dilutes toward bulk as the box grows):

    K_H       = int (exp(-w) - 1) dr                     surface partition length [A]
    A_bound   = A_bulk + int (A - A_bulk) exp(-w) dr / K_H

Concentration enters only through f_bound, in two distinct ways:
  * OB loading, via the specific area a = 6 phi / d43  (Henry limit, PA-independent):
        f_bound = a K_H / (1 + a K_H)
  * phytate loading, via saturation once coverage is not negligible. The PMF is the
    zero-coverage limit, so this needs an isotherm; we fold in the ITC binding
    constant (Ka, n oleosin per PA) as an optional Langmuir correction.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Reuse the fitted-PMF machinery so we never depend on a possibly-stale all_pmfs.csv.
# plot_pmf lives with the Monte-Carlo scripts in ../mc; the shared gvdW lateral term
# (coverage-based Forsman correction) lives with the zeta model in ../zeta.
_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / "mc"))
sys.path.insert(0, str(_ROOT / "zeta"))
import zeta_gvdw as gv  # noqa: E402
from plot_pmf import (  # noqa: E402
    RUNS_DIR,
    box_length_z,
    debye_length,
    fit_tail,
    ionic_strength_from_name,
    ph_from_name,
)

# The MC-output CSVs (shape, rotational mobility) live with the simulation in ../mc.
MC = Path(__file__).resolve().parent.parent / "mc"


def fitted_pmf(run_dir: Path, window: float = 10.0):
    """Return (r, w) for one run, with w the bulk-referenced fitted PMF (w -> 0).

    The same salt-robust non-electrostatic contact well as the zeta/ITC models is
    added, so all three probes share one binding description (gv.EPS_NONEL).
    """
    data = np.genfromtxt(run_dir / "pmf.csv", delimiter=",", names=True)
    z, pmf, err = data["cv"], data["pmf_kT"], data["stderr_kT"]
    r = z + box_length_z(run_dir / "output.yaml") / 2  # distance from decorated wall
    debye = debye_length(ionic_strength_from_name(run_dir.name))
    _, b, _ = fit_tail(r, pmf, err, debye, window)
    return r, pmf - b + gv._nonel_well(r)


def load_profiles(ionic: float):
    """Collect {pH: dict of aligned r-profiles} for one ionic strength."""
    shape = pd.read_csv(MC / "polymer_shape.csv")
    msd = pd.read_csv(MC / "rotational_msd.csv")
    shape["Swall"] = (3 * shape.Szz / (shape.Sxx + shape.Syy + shape.Szz) - 1) / 2

    out = {}
    for run in sorted(RUNS_DIR.glob("*/pmf.csv")):
        run = run.parent
        if abs(ionic_strength_from_name(run.name) - ionic) > 1e-9:
            continue
        pH = ph_from_name(run.name)
        r, w = fitted_pmf(run)
        s = shape[(shape.pH == pH) & (shape.ionic_strength == ionic)].sort_values("r")
        m = msd[(msd.pH == pH) & (msd.ionic_strength == ionic)].sort_values("r")
        if s.empty or m.empty:
            continue
        dd0 = m.msd_rot / m.msd_rot.iloc[-3:].mean()
        out[pH] = {
            "r": r,
            "w": w,
            "Swall": np.interp(r, s.r.values, s.Swall.values),
            "dd0": np.clip(np.interp(r, m.r.values, dd0.values), 1e-3, None),
        }
    return out


def excess(prof, A, A_bulk):
    """Adsorption-excess average of A(r) and the partition length K_H (both from bulk)."""
    r, w = prof["r"], prof["w"]
    bw = np.exp(-w)
    K_H = np.trapezoid(bw - 1.0, r)  # surface partition length [A]
    A_bound = A_bulk + np.trapezoid((A - A_bulk) * bw, r) / K_H
    return K_H, A_bound


def specific_area(phi_ob: float, d43_um: float) -> float:
    """a = 6 phi / d43 [1/A] for a dispersion of spheres (surface area per volume)."""
    return 6.0 * phi_ob / (d43_um * 1e4)  # d43 um -> A


def f_bound_henry(a: float, K_H: float) -> float:
    """Bound fraction in the dilute (Henry) limit: PA-concentration-independent."""
    aK = a * K_H
    return aK / (1.0 + aK)


def gvdw_coverage(K_H: float, n_B: float, d2: float, z: float, ionic: float) -> float:
    """Adsorbed-layer coverage theta = Gamma*d^2 from Nordholm's gvdW 2D isotherm.

    Solve the implicit hard-disk (Eq 7.43) / Lennard-Jones (Eq 7.54) isotherm

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


def f_bound_gvdw(a: float, K_H: float, pa_mM: float, d_A: float, z: float, ionic: float) -> float:
    """Bound fraction from the gvdW isotherm: crowding + lateral interactions cap it.

    f_bound = Gamma a / n_B = (theta/d^2) a / n_B; reduces to Henry (a K_H) as theta->0.
    """
    n_B = pa_mM * 1e-3 * 6.02214076e23 / 1e27  # mM -> 1/A^3
    d2 = d_A**2
    theta = gvdw_coverage(K_H, n_B, d2, z, ionic)
    Gamma = theta / d2
    return Gamma * a / n_B


def f_bound_langmuir(a: float, K_H: float, pa_mM: float, site_mM: float, n: float) -> float:
    """Saturation-corrected bound fraction at finite PA loading.

    Henry gives the intrinsic per-PA affinity; a finite pool of n-oleosin sites caps
    the bound amount. Solve the site balance for bound PA b (mM):
        b = S_tot/n * (K c_free) / (1 + K c_free),   c_free = pa - b,
    with the intrinsic K set so the low-coverage slope matches Henry (f0 = a K_H).
    """
    f0 = f_bound_henry(a, K_H)  # low-coverage bound fraction
    smax = site_mM / n  # max bindable PA [mM]
    K = f0 / (1.0 - f0) / max(smax, 1e-9)  # intrinsic constant reproducing Henry slope
    # fixed-point on free PA
    b = 0.0
    for _ in range(200):
        c = pa_mM - b
        b_new = smax * K * c / (1.0 + K * c)
        if abs(b_new - b) < 1e-9:
            break
        b = b_new
    return b / pa_mM


# Measured 31P DeltaR = R1rho - R1 (s^-1) at OB 50%: (free PA_aq, OB-PA). Figure 9G.
NMR_DELTA_R = {4.0: (162.0, 559.0), 6.0: (130.0, 430.0), 10.0: (5.1, 52.0)}

# Bulk phytate valence z(pH) at I=0.15 M (charge.csv, INO selection, large-r mean).
PHYTATE_Z = {4.0: 5.94, 6.0: 7.08, 7.0: 8.67, 10.0: 10.54}


def site_density_for_f(a: float, K_H: float, f_target: float, pa_mM: float, n: float) -> float:
    """Oleosin site conc. (mM) whose saturation caps the bound fraction at f_target.

    Only meaningful when the unsaturated Henry fraction exceeds f_target (sites can
    only *reduce* binding). Bisect on site_mM until Langmuir matches f_target.
    """
    if f_bound_henry(a, K_H) <= f_target:
        return float("nan")  # geometry-limited: no site density reaches f_target
    lo, hi = 1e-4, 1e4
    for _ in range(200):
        mid = (lo * hi) ** 0.5  # geometric bisection (spans decades)
        if f_bound_langmuir(a, K_H, pa_mM, mid, n) < f_target:
            lo = mid
        else:
            hi = mid
    return (lo * hi) ** 0.5


def match_nmr(intensive, pa_mM, n, ionic=0.15, phi=0.5, d43_um=3.0, d_core=10.6):
    """Does the gvdW isotherm reproduce the NMR pH series — in particular pH4 > pH6?

    f_bound now comes from the charge/salt-dependent gvdW 2D isotherm (Eq 7.44/7.54):
    K_H from the PMF, hard-core d from PyGRPY, lateral a2 the coverage-based Forsman
    correction shared with the zeta and ITC models. Compared against the NMR-required
    f_bound from two-site exchange, f_req = (DeltaR_obs - DeltaR_free)/(DeltaR_bound - DeltaR_free).
    """
    a = specific_area(phi, d43_um)
    print(
        f"\n# gvdW isotherm vs NMR pH series (I={ionic}M, {phi * 100:.0f}% OB, d43={d43_um}um, "
        f"d_core={d_core}A, PA={pa_mM}mM)"
    )
    print(
        f"{'pH':>4} {'z':>5} {'K_H':>5} {'alpha_lat':>9} {'f_Henry':>8} {'f_gvdW':>7} | "
        f"{'f_req(DRb=3e4)':>14} {'f_req(1e5)':>10}"
    )
    print("-" * 78)
    fg = {}
    for pH in (4.0, 6.0, 7.0):
        if pH not in intensive or intensive[pH][0] <= 0:
            continue
        K_H = intensive[pH][0]
        z = PHYTATE_Z[pH]
        f_henry = f_bound_henry(a, K_H)
        f_gvdw = f_bound_gvdw(a, K_H, pa_mM, d_core, z, ionic)
        fg[pH] = f_gvdw
        n_B = pa_mM * 1e-3 * 6.02214076e23 / 1e27
        theta = gvdw_coverage(K_H, n_B, d_core**2, z, ionic)  # for display of converged alpha
        alpha = gv.lateral_alpha(ionic, z, theta / d_core**2)
        row = f"{pH:>4} {z:>5.1f} {K_H:>5.0f} {alpha:>9.2f} {f_henry:>8.2%} {f_gvdw:>7.2%} | "
        if pH in NMR_DELTA_R:
            free, obs = NMR_DELTA_R[pH]
            row += f"{(obs - free) / (3e4 - free):>14.2%} {(obs - free) / (1e5 - free):>10.2%}"
        else:
            row += f"{'--':>14} {'--':>10}"
        print(row)

    print("\n# ordering check (the pH4-vs-pH6 problem):")
    for tag, d in (
        ("f_Henry", {p: f_bound_henry(a, intensive[p][0]) for p in fg}),
        ("f_gvdW ", fg),
        (
            "f_req(3e4)",
            {
                p: (NMR_DELTA_R[p][1] - NMR_DELTA_R[p][0]) / (3e4 - NMR_DELTA_R[p][0])
                for p in fg
                if p in NMR_DELTA_R
            },
        ),
    ):
        if 4.0 in d and 6.0 in d:
            order = "pH4 > pH6" if d[4.0] > d[6.0] else "pH6 > pH4"
            print(f"   {tag:11s}: pH4={d[4.0]:.2%}  pH6={d[6.0]:.2%}  ->  {order}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--ionic", type=float, default=0.15, help="ionic strength (M); use the shallow-well 0.15"
    )
    p.add_argument("--pa-mM", type=float, default=7.6, help="total phytate in the NMR sample (mM)")
    p.add_argument("--n-oleosin", type=float, default=2.5, help="ITC stoichiometry, oleosin per PA")
    p.add_argument(
        "--site-mM", type=float, default=20.0, help="accessible oleosin sites (mM), rough"
    )
    args = p.parse_args()

    profs = load_profiles(args.ionic)

    print(f"# I = {args.ionic} M — per-bound-molecule (intensive, PA-independent)\n")
    print(f"{'pH':>3} {'well_kT':>8} {'K_H(A)':>8} {'<S>':>9} {'<D/D0>':>11} {'<tau/tau0>':>11}")
    intensive = {}
    for pH in sorted(profs):
        pr = profs[pH]
        K_H, S_b = excess(pr, pr["Swall"], 0.0)
        _, dd0_b = excess(pr, pr["dd0"], 1.0)
        _, tau_b = excess(pr, 1.0 / pr["dd0"], 1.0)
        intensive[pH] = (K_H, S_b, tau_b)
        well = pr["w"].min()
        undef = "  (no binding)" if K_H < 0 else ""
        print(f"{pH:>3} {well:>8.2f} {K_H:>8.1f} {S_b:>+9.4f} {dd0_b:>11.2f} {tau_b:>11.2f}{undef}")

    # OB-loading scenarios: the 10% and 50% w/w cream NMR samples. phi ~ w/w (oil ~0.92 g/mL).
    print("\n# bound fraction f_bound and predicted NMR observables")
    print("# native cream d43 ~ a few um at prep pH; sweep to show sensitivity\n")
    scenarios = [("10% OB", 0.10), ("50% OB", 0.50)]
    d43s = [1.0, 3.0, 11.0]  # um
    for pH in [x for x in (4.0, 6.0, 7.0) if x in intensive and intensive[x][0] > 0]:
        K_H, S_b, tau_b = intensive[pH]
        print(
            f"pH {pH:g}:  K_H = {K_H:.0f} A,  S_bound = {S_b:+.4f},  tau/tau0|local = {tau_b:.2f}"
        )
        for label, phi in scenarios:
            for d43 in d43s:
                a = specific_area(phi, d43)
                fH = f_bound_henry(a, K_H)
                fL = f_bound_langmuir(a, K_H, args.pa_mM, args.site_mM, args.n_oleosin)
                S_obs = fH * S_b  # A_free(S) = 0
                print(
                    f"   {label:6s} d43={d43:4.1f}um  a={a:.2e}/A  "
                    f"f_bound: Henry={fH:6.2%}  Langmuir({args.pa_mM:g}mM)={fL:6.2%}  "
                    f"<S>_obs={S_obs:+.5f}"
                )
        print()

    match_nmr(intensive, args.pa_mM, args.n_oleosin)


if __name__ == "__main__":
    main()
