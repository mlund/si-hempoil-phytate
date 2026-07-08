#!/usr/bin/env python3
"""Charge-regulated Poisson-Boltzmann (CR-PB) interpretation of the OB data.

A flat Gouy-Chapman surface on which phytate (valence z = -|z|) adsorbs to Langmuir
sites that regulate the surface charge. Each adsorbed anion lowers the local potential
and suppresses the next (regulation feedback), so the self-consistent coverage can settle
*past* neutralisation -- overcharging -- which inverts the diffuse-layer sign and
reproduces the pH-4 zeta charge reversal (+37.7 -> -4.4 mV). The overcharge is
self-limiting, so no correlation term is needed.

The MC feeds two quantities, reflecting the brush *dipole* (cationic N-terminus at the
outer binding edge, anionic anchor at the wall):

  * binding-plane potential psi_beta -- the MC local potential at the brush periphery
    (potential_profiles.csv at r ~ R_BETA). Positive at the binding edge (+52/+5 mV,
    pH 4/6, low salt) even where the net monopole is ~0 or negative, so it *drives
    binding* at pH 6 where a net-charge Grahame surface would predict none; strongly
    negative at pH 10 (-89 mV) -> no binding.
  * net surface charge sigma0 -- from the measured bare (0-PA) zeta via the Grahame
    equation. It sets the diffuse-layer zeta and the overcharging,
    sigma_eff = sigma0 - |z| Gamma.

So the dipole drives *binding*, the monopole + regulation drive *zeta*.

One parameter: K, the intrinsic Langmuir affinity (Gamma_max is a geometric
close-packing prior, 1/(11 A)^2). K is the *non-electrostatic* part of the binding --
the electrostatics are already carried by psi_beta and sigma_eff, so K is not the MC
Henry constant K_H = int(e^{-w}-1)dr (which packs the full PMF and would double-count
them). K reflects the pH-set peptide/phytate ionisation and is nearly salt-independent,
so it is fit per pH to each zeta-vs-[PA] series and held across salt. (eps_ne, an
equivalent contact affinity in the exponent, is degenerate with K and kept off by
default.)

The CR-PB has two uses: (1) computing the zeta-potential / overcharge from the MC, and
(2) providing an isotherm to fit the ITC thermogram (K, area, dH free -- see fit_itc).
NMR is not taken through the CR-PB: a bound fraction at a single salt is read from the MC
PMF directly (K_H = int(e^{-w}-1)dr -> Nordholm gvdW isotherm, nmr_f), parameter-free.
The ITC heat is out of scope: at the physical grafting density the pH-4 coverage is
overcharge-limited to ~0.03 PA/oleosin while the thermogram needs ~10x more, so the heat
reflects inter-droplet cross-linking / aggregation, not monolayer adsorption.
"""

from __future__ import annotations

import os
import sys
from functools import cache
from pathlib import Path

import numpy as np
from scipy import constants
from scipy.optimize import brentq, curve_fit, minimize_scalar

# Reuse the MC-grounded machinery (kappa, Grahame sigma, valence, bare anchors).
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mc_shared as gv  # noqa: E402

KT_MV = gv.KT_MV
LB = gv.LB
N_A = constants.N_A

# Binding plane = brush periphery / PMF-minimum region, where the MC dipolar
# potential is read (the contact position R_CONTACT).
R_BETA = gv.R_CONTACT  # A
# Bare zeta was measured at low salt; sigma0 (a surface property, ~salt-independent
# in e/A^2) is fixed there and carried to every ionic strength via Grahame.
I_ANCHOR = 0.005
# Close-packing saturation of adsorbed phytate: one anion per (11 A)^2 (~1.3x the
# hard-disk area), the geometric ceiling on the regulated charge.
GAMMA_MAX = 1.0 / 11.0**2  # PA / A^2
DILUTION = gv.DILUTION  # emulsion diluted 1000x before the zeta measurement
PHS = gv.PHS  # (4, 6, 10)
I_ITC = 0.015  # ITC cell ionic strength (matches fit_thermogram)

# Salt-robust non-electrostatic contact affinity [kT] (specific phosphate-Lys/Arg
# ion-pairing / H-bonding). It adds to the Langmuir binding free energy and is
# therefore EXACTLY degenerate with K (both enter as K*exp(eps_ne)) -- it cannot
# close the NMR gap in the pure charge-regulated model (see module docstring). Kept
# at 0 (unidentifiable); overridable via $EPS_NE only to demonstrate the degeneracy.
EPS_NE = float(os.environ.get("EPS_NE", "0.0"))  # kT

# Physical grafting area per oleosin [A^2]: the MC brush (50 chains / 130x130 A^2) and
# the measured OB oleosin surface density (~1 umol/g) agree at ~338 A^2. Used to pin
# the ITC capacity to the real grafting density rather than fitting it.
AREA_OLEOSIN = 130.0 * 130.0 / 50.0  # A^2 per oleosin


def _y0(sigma: float, ionic: float) -> float:
    """Grahame reduced surface potential y0 = e psi_d / kT for charge `sigma` [e/A^2]."""
    return 2.0 * np.arcsinh(2.0 * np.pi * LB * sigma / gv.kappa(ionic))


def _psi_d(sigma: float, ionic: float) -> float:
    """Diffuse-layer surface potential [mV] from Grahame."""
    return _y0(sigma, ionic) * KT_MV


@cache
def sigma0(pH: float) -> float:
    """Net brush surface charge [e/A^2] from the measured bare (0-PA) zeta (monopole)."""
    return gv.sigma_from_zeta(gv.ZETA_BARE[pH], I_ANCHOR)


@cache
def _mc_psi_beta(pH: float, ionic: float) -> float:
    """MC local potential [mV] at the binding plane R_BETA (the dipolar signature).

    Read from potential_profiles.csv (total mean potential the phytate feels) and
    log-interpolated in ionic strength between the simulated values, so `ionic`
    need not be a simulated one. Positive at pH 4/6 (cationic binding edge),
    strongly negative at pH 10.
    """
    lo, hi, t = gv._ionic_bracket(ionic)

    def at(ii: float) -> float:
        p = gv._table("potential_profiles.csv")
        d = p[(p.ionic_strength == ii) & (p.pH == pH)].sort_values("r")
        return float(np.interp(R_BETA, d.r.to_numpy(), d.potential_all_mV.to_numpy()))

    return at(lo) if t == 0.0 else (1 - t) * at(lo) + t * at(hi)


class State:
    """Converged charge-regulated state at one (pH, ionic, c_L)."""

    __slots__ = ("gamma", "theta", "sigma_eff", "zeta", "psi_beta")

    def __init__(self, gamma, theta, sigma_eff, zeta, psi_beta):
        self.gamma = gamma
        self.theta = theta
        self.sigma_eff = sigma_eff
        self.zeta = zeta
        self.psi_beta = psi_beta


def crpb_state(
    pH: float,
    ionic: float,
    c_L_M: float,
    K: float,
    eps_ne: float = EPS_NE,
    gamma_max: float = GAMMA_MAX,
):
    """Self-consistent charge-regulated adsorption at free phytate ``c_L_M`` [M].

    Solves Gamma = gamma_max * theta(Gamma), where the Langmuir occupancy
        theta = K c^s / (1 + K c^s),  c^s = c_L exp(z psi_beta / kT + eps_ne)
    uses the regulated binding-plane potential
        psi_beta(Gamma) = psi_beta^MC + [psi_d(sigma_eff) - psi_d(sigma0)],
        sigma_eff = sigma0 - |z| Gamma.
    The MC dipolar baseline drives binding; the adsorbed anion lowers psi_beta and
    self-limits it (overcharge). This mean-field enhancement is too salt-steep to also
    carry the NMR magnitude at 0.15 M, so the NMR bound fraction is read from the MC PMF
    directly (see nmr_f); the CR-PB is used only for the zeta overcharge. Returns a State
    with zeta = psi_d(sigma_eff).
    """
    z = gv.phytate_valence(pH, ionic)  # |z| > 0
    s0 = sigma0(pH)
    psi_mc = _mc_psi_beta(pH, ionic)
    psi_d0 = _psi_d(s0, ionic)

    def theta_of(gamma: float) -> float:
        psi_beta = psi_mc + (_psi_d(s0 - z * gamma, ionic) - psi_d0)
        cs = c_L_M * np.exp(np.clip(z * psi_beta / KT_MV + eps_ne, -300, 300))
        return K * cs / (1.0 + K * cs)

    def residual(gamma: float) -> float:
        return gamma - gamma_max * theta_of(gamma)

    hi = gamma_max * (1 - 1e-9)
    gamma = brentq(residual, 0.0, hi) if residual(hi) > 0 else hi
    theta = gamma / gamma_max
    sigma_eff = s0 - z * gamma
    psi_beta = psi_mc + (_psi_d(sigma_eff, ionic) - psi_d0)
    return State(gamma, theta, sigma_eff, _psi_d(sigma_eff, ionic), psi_beta)


# --- probe 1: zeta ---------------------------------------------------------------
def zeta_pa(
    pH: float, pa_mM: float, K: float, eps_ne: float = EPS_NE, ionic: float = I_ANCHOR
) -> float:
    """Predicted zeta [mV] at total phytate `pa_mM` (1000x diluted for the measurement)."""
    return crpb_state(pH, ionic, pa_mM * 1e-3 / DILUTION, K, eps_ne).zeta


def _zeta_data(pH: float):
    return gv._fig6a(pH)  # (PA_mM, mean_mV, sd_mV), PA>0 low-salt series


def _zeta_chi2(pH: float, K: float, eps_ne: float = EPS_NE) -> float:
    """Error-weighted chi^2 of the zeta-vs-[PA] series at ONE pH."""
    x, y, e = _zeta_data(pH)
    res = np.array([zeta_pa(pH, v, K, eps_ne) - yy for v, yy in zip(x, y, strict=True)])
    sd = np.clip(np.array(e), 0.3, None)
    return float(np.sum((res / sd) ** 2))


def fit_K(pH: float, eps_ne: float = EPS_NE) -> float:
    """Fit the intrinsic K at ONE pH to that pH's zeta-vs-[PA] series.

    K reflects the peptide/PA ionisation, which is set by pH (and nearly
    salt-independent), so it is fit per pH -- not globally -- then held fixed across
    salt. The pH-dependent electrostatics (psi_beta dipole, sigma0 monopole) still
    come from the MC; K adds only the pH-specific contact chemistry.

    When adsorption does not improve the zeta fit (pH 10: the series is flat in [PA],
    so the double layer never inverts and zeta is insensitive to K), K is
    unidentifiable and a free fit runs to the bound; we return 0 (no binding) rather than
    a spurious runaway K.
    """
    sol = minimize_scalar(
        lambda logK: _zeta_chi2(pH, 10**logK, eps_ne),
        bounds=(-6, 6),
        method="bounded",
        options={"xatol": 1e-3},
    )
    K = 10**sol.x
    if _zeta_chi2(pH, K, eps_ne) > _zeta_chi2(pH, 0.0, eps_ne) - 1.0:
        return 0.0  # binding does not improve the fit -> zeta does not constrain K
    return K


def fit_K_by_pH(eps_ne: float = EPS_NE) -> dict[float, float]:
    """Per-pH intrinsic K, each fit to its own zeta-vs-[PA] series."""
    return {pH: fit_K(pH, eps_ne) for pH in PHS}


def eps_degeneracy(pH: float = 4.0, eps_grid=(0.0, 1.0, 2.0, 3.0, 4.0)):
    """Show eps_ne is degenerate with K at one pH: refit K at each eps, read K, the CR-PB
    coverage bound fraction, and chi^2. K*exp(eps_ne), the zeta chi^2 and the CR-PB
    coverage are all invariant -- a salt-robust contact affinity is unidentifiable, so it
    adds nothing; NMR is read from the PMF directly (nmr_f).
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "nmr"))
    import binding_fraction as nb  # noqa: E402

    a = nb.specific_area(0.5, 3.0)
    n_B = 7.6e-3 * N_A / 1e27
    rows = []
    for eps in eps_grid:
        K = fit_K(pH, eps)
        f_cov = crpb_state(pH, 0.15, 7.6e-3, K, eps).gamma * a / n_B  # CR-PB coverage f
        rows.append((eps, K, K * np.exp(eps), f_cov, _zeta_chi2(pH, K, eps)))
    return rows


# --- probe 2: ITC ----------------------------------------------------------------
def fit_itc(pH: float = 4.0, ionic: float = I_ITC, eps_ne: float = EPS_NE):
    """Fit the CR-PB isotherm to the pH-4 thermogram with K, area and dH ALL free.

    A three-parameter fit (like the one-site model) that lets the thermogram choose the
    affinity and capacity. It fits well -- but the returned ``area`` is the diagnostic:
    ~15-20x the physical grafting density (AREA_OLEOSIN), and K-area are strongly
    anti-correlated, so a good fit here is not evidence for monolayer binding (the heat
    is cross-linking; see the write-up). Returns (K, area, dH, rms, corr) with corr the
    K-area fit correlation.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "itc"))
    import extract_thermograms as itc_data  # noqa: E402
    import fit_thermogram as itc_fit  # noqa: E402

    titr = itc_data.load(pH)
    grid = np.logspace(-9, -2, 24)

    def model(_x, logK, area, dH):
        cov = np.array([crpb_state(pH, ionic, c, 10**logK, eps_ne).gamma for c in grid])
        return itc_fit._heat_uJ(titr, itc_fit._bound_from_isotherm(titr, area, grid, cov), dH)

    popt, pcov = curve_fit(
        model,
        titr.mole_ratio,
        titr.heat_uJ,
        p0=[1.8, 5000.0, -35.0],
        bounds=([-3, 50, -1e4], [8, 1e5, 0]),
        maxfev=40000,
    )
    rms = float(np.sqrt(np.mean((titr.heat_uJ - model(titr.mole_ratio, *popt)) ** 2)))
    corr = float(itc_fit.corr_from_cov(pcov)[0, 1])  # logK-area correlation
    return (10 ** popt[0], popt[1], popt[2], rms, corr)


def plot_itc(eps_ne: float, path: Path) -> None:
    """pH-4 ITC thermogram: data (replicate mean +/- SD) with the free-K/area CR-PB fit
    and the one-site model. Both fit the curve well; the point of the write-up is that
    the CR-PB fit's area is unphysical (a good fit is not a good model)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    itc = str(Path(__file__).resolve().parent.parent / "itc")
    sys.path.insert(0, itc)
    import extract_thermograms as itc_data  # noqa: E402
    import fit_thermogram as itc_fit  # noqa: E402

    # Replicate mean +/- SD from the two matched-scheme runs (rep1 has a different scheme).
    reps = [itc_data.load_xls(Path(itc) / f"pH4/pH4_rep{i}.xls", 4.0) for i in (2, 3)]
    n = min(len(r.heat_uJ) for r in reps)
    heats = np.vstack([r.heat_uJ[:n] for r in reps])
    x, y, e = reps[-1].mole_ratio[:n], heats.mean(0), heats.std(0)

    titr = itc_data.load(4.0)
    grid = np.logspace(-9, -2, 24)
    K, area, dH, _, _ = fit_itc(eps_ne=eps_ne)
    cov = np.array([crpb_state(4.0, I_ITC, c, K, eps_ne).gamma for c in grid])
    crpb = itc_fit._heat_uJ(titr, itc_fit._bound_from_isotherm(titr, area, grid, cov), dH)
    Ka, dH1, n1, _ = itc_fit.fit_onesite(titr)
    onesite = itc_fit.onesite_heat(titr, Ka, dH1, n1)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(
        x, y, yerr=e, fmt="o", color="#333333", ms=5, capsize=2, label="experiment (mean ± SD)"
    )
    ax.plot(titr.mole_ratio, crpb, "-", color="#1f4fd8", lw=2, label="charge-regulated (CR-PB)")
    ax.plot(titr.mole_ratio, onesite, "--", color="gray", lw=1.5, label="one-site (independent)")
    ax.axhline(0, color="gray", lw=0.5)
    ax.set_xlabel("mole ratio PA : oleosin")
    ax.set_ylabel(r"heat per injection ($\mu$J)")
    ax.legend(fontsize=9, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=140)


# --- probe 3: NMR ----------------------------------------------------------------
def nmr_f(
    pH: float,
    ionic: float = 0.15,
    pa_mM: float = 7.6,
    phi: float = 0.5,
    d43_um: float = 3.0,
) -> float:
    """NMR bound fraction from the MC PMF directly (not through the CR-PB).

    NMR is a bound fraction at a single salt, so the simulated PMF at that salt is the
    direct, parameter-free input: K_H = int(e^{-w}-1)dr from the MC PMF (with the shared
    non-electrostatic contact well folded in), fed to the Nordholm gvdW isotherm
    (binding_fraction.bound_fraction). No fitted parameter enters.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "nmr"))
    import binding_fraction as nb  # noqa: E402

    a = nb.specific_area(phi, d43_um)  # 1/A
    r_w, w0 = gv._pmf_interp(pH, ionic)  # MC PMF (nonel well included)
    k_henry = float(np.trapezoid(np.expm1(-w0), r_w))  # surface partition length [A]
    z = gv.phytate_valence(pH, ionic)
    return nb.bound_fraction(a, k_henry, pa_mM, 10.6, z, ionic)


def nmr_f_req(pH: float, Rb: float = 3e4) -> float:
    """Bound fraction the two-site fast-exchange model requires from the measured DeltaR.

    Rb is the fully-bound (rigid-limit) rate; f_req = (obs - free)/(Rb - free). The measured
    (free PA_aq, OB-PA) DeltaR pairs are the single copy in binding_fraction.NMR_DELTA_R.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "nmr"))
    import binding_fraction as nb  # noqa: E402

    free, obs = nb.NMR_DELTA_R[pH]
    return (obs - free) / (Rb - free)


# --- sanity limits (minor) -------------------------------------------------------
def _sanity(K_by_pH: dict[float, float], eps_ne: float = EPS_NE) -> list[tuple[str, bool]]:
    """A handful of limit checks so the solver is trustworthy (not a QA suite)."""
    k4, k6 = K_by_pH[4.0], K_by_pH[6.0]
    checks = []
    # c_L -> 0: no coverage, zeta -> bare Grahame at that salt.
    s = crpb_state(4.0, I_ANCHOR, 1e-18, k4, eps_ne)
    bare = _psi_d(sigma0(4.0), I_ANCHOR)
    checks.append(("c->0 -> bare zeta", abs(s.zeta - bare) < 0.1 and s.gamma < 1e-9))
    # Regulation caps coverage below gamma_max: each adsorbed anion lowers psi_beta and
    # suppresses the next, so even c_L -> inf settles at a self-limited ceiling, NOT the
    # geometric gamma_max. (A pure Langmuir would saturate at gamma_max -- the point of
    # a *regulated* model is that it does not.)
    hi = crpb_state(4.0, I_ANCHOR, 1e6, k4, eps_ne)
    checks.append(("regulation caps Gamma < gamma_max", 0 < hi.gamma < GAMMA_MAX))
    # Monotone: more phytate -> more coverage.
    lo_c = crpb_state(4.0, I_ANCHOR, 1e-6, k4, eps_ne)
    hi_c = crpb_state(4.0, I_ANCHOR, 1e-4, k4, eps_ne)
    checks.append(("Gamma monotone in c_L", hi_c.gamma > lo_c.gamma))
    # The pH-4 reversal: enough coverage to invert the sign (sigma_eff, hence zeta, < 0).
    checks.append(("pH4 overcharges (sigma_eff<0)", hi.sigma_eff < 0 and hi.zeta < 0))
    # High salt kills zeta magnitude (Grahame screening).
    lo_salt = abs(_psi_d(sigma0(4.0), 0.005))
    hi_salt = abs(_psi_d(sigma0(4.0), 0.5))
    checks.append(("high salt lowers |zeta|", hi_salt < lo_salt))
    # Self-consistency of the fixed point.
    s6 = crpb_state(6.0, I_ANCHOR, 5e-6, k6, eps_ne)
    checks.append(
        (
            "sigma_eff = sigma0-|z|G",
            abs(s6.sigma_eff - (sigma0(6.0) - gv.phytate_valence(6.0, I_ANCHOR) * s6.gamma))
            < 1e-12,
        )
    )
    return checks


def plot_crpb(K_by_pH: dict[float, float], eps_ne: float, path: Path) -> None:
    """Two-panel CR-PB zeta figure in the measured-vs-model style: experiment as error-bar
    markers, CR-PB model as a smooth absolute-prediction line (per-pH K held fixed across
    salt in panel B). A: zeta vs [PA] (low salt); B: zeta vs [NaCl] at 1.5 mM PA, all pH."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {4.0: "#1f4fd8", 6.0: "#e08a00", 10.0: "#2ca02c"}  # blue / orange / green
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    # A: zeta vs [PA] at low salt, one series per pH; model as a smooth line.
    pa_fine = np.linspace(0.02, 7.7, 200)
    for pH in PHS:
        x, y, e = _zeta_data(pH)
        c, K = colors[pH], K_by_pH[pH]
        ax_a.errorbar(x, y, yerr=e, fmt="x", color=c, ms=6, capsize=2, label=f"pH {pH:g}")
        ax_a.plot(pa_fine, [zeta_pa(pH, v, K, eps_ne) for v in pa_fine], "-", color=c, lw=1.8)
    ax_a.axhline(0, color="gray", lw=0.5)
    ax_a.set_ylim(-75, 40)
    ax_a.set_xlabel("[PA] (mM)")
    ax_a.set_ylabel(r"$\zeta$ (mV)")

    # B: zeta vs [NaCl] at 1.5 mM PA (per-pH K held fixed), all three pH; model as a line.
    # Absolute prediction (not anchored to the measured bare salt curve): zeta comes
    # straight from crpb_state = psi_d(sigma_eff, ionic), so both the bare salt decay
    # and the phytate shift are predicted. The residual (~10 mV at pH 6 high salt) is
    # the flat-Grahame-vs-permeable-brush limitation of the mean-field surface.
    nacl_fine = np.linspace(1.0, 100.0, 200)
    for pH in PHS:
        c, K = colors[pH], K_by_pH[pH]
        xr, yr, er = gv._fig6c(pH, 1.5)
        ax_b.errorbar(xr, yr, yerr=er, fmt="x", color=c, ms=6, capsize=2, label=f"pH {pH:g}")
        model = [
            crpb_state(pH, n / 1000 + 0.003, 1.5e-3 / DILUTION, K, eps_ne).zeta for n in nacl_fine
        ]
        ax_b.plot(nacl_fine, model, "-", color=c, lw=1.8)
    ax_b.axhline(0, color="gray", lw=0.5)
    ax_b.set_xlabel("[NaCl] (mM)")
    ax_b.legend(fontsize=9, frameon=False, ncol=3)

    fig.tight_layout(pad=1)
    fig.savefig(path, dpi=140)


if __name__ == "__main__":
    eps_ne = EPS_NE  # 0: the contact affinity is unidentifiable (degenerate with K)
    # K reflects the pH-set peptide/PA ionisation, so it is fit PER pH (not globally),
    # then held fixed across salt (ionisation is nearly salt-independent).
    K_by_pH = fit_K_by_pH(eps_ne)
    print("CR-PB: intrinsic K fit per pH to each zeta-vs-[PA] series (ionisation is pH-set):")
    for pH in PHS:
        print(
            f"    K(pH {pH:g}) = {K_by_pH[pH]:.3g} /M   (chi^2 = {_zeta_chi2(pH, K_by_pH[pH]):.0f})"
        )
    print(
        f"  (Gamma_max = 1/(11 A)^2 = {GAMMA_MAX:.4f}/A^2, R_beta = {R_BETA:g} A, "
        f"eps_ne = {eps_ne:g} kT; NMR read from the MC PMF directly, not the CR-PB)\n"
    )

    print("non-electrostatic contact affinity is DEGENERATE with K at each pH (pH 4 shown):")
    print(f"  {'eps_ne':>7} {'K(/M)':>9} {'K*exp(eps)':>11} {'f_cover pH4':>11} {'zeta_chi2':>10}")
    for eps, Ke, Kexp, f4, chi2 in eps_degeneracy(4.0):
        print(f"  {eps:>7.1f} {Ke:>9.3g} {Kexp:>11.3g} {f4:>11.3%} {chi2:>10.0f}")

    print("\nsanity limits:")
    for name, ok in _sanity(K_by_pH, eps_ne):
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")

    # sigma0 anchor cross-check: the CR model anchors sigma0 to the *measured* bare zeta
    # (monopole), while feeding the *MC* dipolar psi_beta for binding. The MC's own bare
    # zeta (simulated_zeta.csv, shear-plane rule) reveals why they must be kept separate.
    mcz = gv._table("simulated_zeta.csv")
    mcz = mcz[mcz.ionic_strength == I_ANCHOR]
    print(f"\nsigma0 anchor cross-check (low salt, {I_ANCHOR * 1e3:.0f} mM):")
    print(f"  {'pH':>4} {'measured':>9} {'MC shear-plane':>15} {'note':>12}")
    for pH in PHS:
        meas = gv.ZETA_BARE[pH]
        mc = float(mcz[mcz.pH == pH].zeta_mV.iloc[0])
        note = "sign flip" if (meas < 0) != (mc < 0) else "ok"
        print(f"  {pH:>4g} {meas:>+9.1f} {mc:>+15.1f} {note:>12}")

    print("\nprobe 1 -- zeta vs [PA] (low salt), per-pH K: does CR-PB overcharge at pH 4?")
    print(f"  {'pH':>4} {'[PA]/mM':>8} {'psi_b/mV':>9} {'theta':>7} {'CR zeta':>8} {'data':>10}")
    for pH in PHS:
        x, y, e = _zeta_data(pH)
        for pa_mM, ym, es in zip(x, y, e, strict=True):
            st = crpb_state(pH, I_ANCHOR, pa_mM * 1e-3 / DILUTION, K_by_pH[pH], eps_ne)
            print(
                f"  {pH:>4g} {pa_mM:>8.2f} {st.psi_beta:>+9.1f} {st.theta:>7.3f} "
                f"{st.zeta:>+8.1f} {ym:>+7.1f}+-{es:.1f}"
            )

    # Salt-independence of K: predict zeta vs [NaCl] at 1.5 mM PA with the per-pH K held
    # fixed (fit at low salt). If it tracks the measured NaCl series, ionisation -- hence
    # K -- is indeed nearly salt-independent, as assumed.
    print("\nprobe 1b -- zeta vs [NaCl] at 1.5 mM PA, per-pH K FIXED (salt-independence test):")
    print(f"  {'pH':>4} {'[NaCl]/mM':>10} {'CR zeta':>8} {'data':>10}")
    for pH in PHS:
        xr, yr, er = gv._fig6c(pH, 1.5)
        for nacl, ym, es in zip(xr, yr, er, strict=True):
            ionic = nacl / 1000 + 0.003
            zc = crpb_state(pH, ionic, 1.5e-3 / DILUTION, K_by_pH[pH], eps_ne).zeta
            print(f"  {pH:>4g} {nacl:>10.0f} {zc:>+8.1f} {ym:>+7.1f}+-{es:.1f}")

    print("\nprobe 2 -- NMR bound fraction from the MC PMF DIRECTLY (0.15 M, 7.6 mM):")
    print(f"  {'pH':>4} {'f_sim':>8} {'f_req':>8} {'ratio':>7}")
    for pH in PHS:
        f_sim, f_rq = nmr_f(pH), nmr_f_req(pH)
        ratio = f"{f_rq / f_sim:>5.1f}x" if f_sim > 0 else "  n/a"
        print(f"  {pH:>4g} {f_sim:>8.2%} {f_rq:>8.2%} {ratio:>7}")
    print("  (K_H read from the MC PMF at 0.15 M -> parameter-free; the CR-PB is not used here)")

    print("\nprobe 3 -- ITC (pH 4): the CR-PB isotherm fits well, but the area is unphysical")
    K_itc, area, dH, rms, corr = fit_itc()  # K, area, dH all free
    print(
        f"  free fit (K,area,dH): K={K_itc:.0f}/M  area={area:.0f} A^2  dH={dH:.0f} kJ/mol  "
        f"RMS={rms:.1f} uJ  K-area corr={corr:+.2f}"
    )
    print(
        f"  -> a GOOD fit, but area = {area / AREA_OLEOSIN:.0f}x the physical grafting density "
        f"({AREA_OLEOSIN:.0f} A^2/oleosin, MC = measured 1 umol/g);"
    )
    print(
        "     the heat is dominated by inter-droplet cross-linking (the multiphasic shoulder"
        " the manuscript reports), not monolayer binding -- a good fit is not a good model"
    )

    figure = Path(__file__).parent / "fig_crpb.png"
    plot_crpb(K_by_pH, eps_ne, figure)
    itc_fig = Path(__file__).resolve().parent.parent / "itc" / "fig_itc_crpb.png"
    plot_itc(eps_ne, itc_fig)
    print(f"\nsaved {figure.name}, {itc_fig.name}")
