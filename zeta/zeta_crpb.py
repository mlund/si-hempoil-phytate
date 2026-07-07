#!/usr/bin/env python3
"""Charge-regulated Poisson-Boltzmann (CR-PB) interpretation of the OB data.

An alternative to the gvdW model. The gvdW reconciles ITC, NMR and the pH 6/10
zeta + salt series, but cannot reproduce the pH-4 zeta charge-reversal
(+37.7 -> -4.4 mV): its coverage-based correlation hole self-limits the coverage
below neutralisation, so the mean-field double layer never inverts. Here we ask
whether a *charge-regulated* flat surface -- phytate (z = -|z|) adsorbing on
Langmuir sites that regulate the surface charge -- can capture the reversal while
still fitting ITC and NMR, fed the same MC-derived quantities.

The mechanism the gvdW lacks: an intrinsic site affinity K holds phytate even
against a modestly reversed potential, so the self-consistent coverage settles
*past* neutralisation (sigma_eff < 0) -- overcharging -- rather than stopping at
it. Overcharging is self-limiting because each adsorbed anion lowers the local
potential and suppresses the next (the regulation feedback), so no diverging
correlation term is needed.

MC feeds two distinct quantities (the brush is a *dipole*: cationic N-terminus at
the outer binding edge, anionic anchor at the wall):

  * binding-plane potential psi_beta -- the MC local potential at the brush
    periphery (potential.csv, r ~ R_BETA). It carries the dipole: positive at the
    binding edge (+52/+5 mV, pH 4/6, low salt) even where the net monopole is ~0
    or negative, so it *drives binding* at pH 6 where a net-charge Grahame surface
    would predict none. Strongly negative at pH 10 (-89 mV) -> no binding.
  * net surface charge sigma0 -- from the measured bare (0-PA) zeta via the
    Grahame equation (monopole). It sets the diffuse-layer zeta and the
    overcharging: sigma_eff = sigma0 - |z| Gamma.

So the dipole drives *binding*, the monopole + regulation drive *zeta*. At pH 6
sigma0 < 0, the dipole still binds phytate, and the adsorbed z=-7 anion pushes
sigma_eff (hence zeta) further negative -- the measured trend.

One primary parameter, K (intrinsic Langmuir constant); Gamma_max is pinned by a
geometric close-packing prior 1/(11 A)^2. K reflects the peptide/PA ionisation,
which is set by pH and varies little with salt, so it is fit PER pH to each
zeta-vs-[PA] series -- not globally -- and then held fixed across salt (validated
against the zeta-vs-[NaCl] series) and *transferred* to ITC (fit only area + dH) and
NMR (predict f_bound). A global K would force one chemistry across pH 4/6/10, which
is wrong; the pH-dependent electrostatics (psi_beta dipole, sigma0 monopole) come
from the MC, and K(pH) adds only the pH-specific contact chemistry.

Note what K is NOT: it is not the MC Henry constant K_H = int(e^{-w}-1)dr. That K_H
packs the *full* PMF -- electrostatics + excluded volume + contact -- into one
number; feeding it as the site constant here double-counts the electrostatics (which
PB already carries through psi_beta and sigma_eff) and blows up the low-salt
overcharge (pH-4 zeta -> ~-150 mV). K is instead the *non-electrostatic* part of the
binding (electrostatics handled by PB): essentially K_H with the electrostatic PMF
stripped out, i.e. the same contact content the gvdW packages as EPS_NONEL. That
non-electrostatic K is the SAME physical quantity as eps_ne below -- one written
multiplicatively, the other in the exponent -- which is why they are degenerate. The
CR model therefore has exactly ONE non-electrostatic affinity; the ~2x NMR gap is the
genuine tension between the zeta overcharge and the NMR magnitude that one affinity
(+ PB electrostatics) cannot resolve, not a missing parameter.

The NMR magnitude is ~2x below the value two-site exchange requires, and -- unlike
the gvdW -- this gap CANNOT be closed by a salt-robust non-electrostatic contact
affinity eps_ne. In the pure charge-regulated Langmuir occupancy
    theta = K c^s/(1+K c^s),  c^s = c exp(|z| psi_beta/kT + eps_ne),
K and eps_ne enter only as the product K*exp(eps_ne), so eps_ne is exactly
degenerate with K: refitting K to zeta absorbs any eps_ne and leaves the NMR f
unchanged (see the sweep in __main__). The gvdW's eps DID close its gap because
its coverage saturates at low salt through a *binding-independent* lateral
correlation cap (Forsman a2); that cap makes eps identifiable. But adding the same
cap here destroys the pH-4 reversal -- the cap self-limits the coverage below
neutralisation, exactly the gvdW's failure. So the reversal (needs binding past
neutralisation, no cap) and a closable NMR gap (needs a binding-independent cap)
are mutually exclusive. gvdW and CR-PB are therefore complementary: the gvdW nails
the NMR/ITC/salt magnitude but misses the pH-4 reversal; the CR-PB captures the
reversal but underpredicts the NMR magnitude by ~2x. eps_ne is kept as a knob
(default 0, unidentifiable) only to document this.
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
import zeta_gvdw as gv  # noqa: E402

KT_MV = gv.KT_MV
LB = gv.LB
N_A = constants.N_A

# Binding plane = brush periphery / PMF-minimum region, where the MC dipolar
# potential is read (same position as the gvdW contact well, R_CONTACT).
R_BETA = gv.R_CONTACT  # A
# Bare zeta was measured at low salt; sigma0 (a surface property, ~salt-independent
# in e/A^2) is fixed there and carried to every ionic strength via Grahame.
I_ANCHOR = 0.005
# Close-packing saturation of adsorbed phytate: one anion per (11 A)^2 (~1.3x the
# hard-disk area), the geometric ceiling on the regulated charge.
GAMMA_MAX = 1.0 / 11.0**2  # PA / A^2
DILUTION = gv.DILUTION  # emulsion diluted 1000x before the zeta measurement
PHS = gv.PHS  # (4, 6, 10)
I_ITC = 0.015  # ITC cell ionic strength (matches itc_fit)

# Salt-robust non-electrostatic contact affinity [kT] (specific phosphate-Lys/Arg
# ion-pairing / H-bonding). It adds to the Langmuir binding free energy and is
# therefore EXACTLY degenerate with K (both enter as K*exp(eps_ne)) -- it cannot
# close the NMR gap in the pure charge-regulated model (see module docstring). Kept
# at 0 (unidentifiable); overridable via $EPS_NE only to demonstrate the degeneracy.
EPS_NE = float(os.environ.get("EPS_NE", "0.0"))  # kT


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

    Read from potential.csv (total mean potential the phytate feels) and
    log-interpolated in ionic strength between the simulated values, so `ionic`
    need not be a simulated one. Positive at pH 4/6 (cationic binding edge),
    strongly negative at pH 10.
    """
    lo, hi, t = gv._ionic_bracket(ionic)

    def at(ii: float) -> float:
        p = gv._table("potential.csv")
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
        theta = K c^s / (1 + K c^s),  c^s = c_L exp(+|z| psi_beta / kT + eps_ne)
    uses the regulated binding-plane potential
        psi_beta(Gamma) = psi_beta^MC + [psi_d(sigma_eff) - psi_d(sigma0)],
        sigma_eff = sigma0 - |z| Gamma.
    The dipolar MC baseline drives binding; the adsorbed anion lowers psi_beta and
    self-limits it; the salt-robust contact affinity ``eps_ne`` [kT] deepens binding.
    Returns a State with zeta = psi_d(sigma_eff).
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
    unidentifiable and a free fit runs to the bound; we return 0 (no binding) so the
    transferred NMR/ITC predictions are not driven by a spurious runaway K.
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
    """Show eps_ne is degenerate with K at one pH: refit K at each eps, read K, f, chi^2.

    K*exp(eps_ne), the zeta chi^2 and the NMR f are all invariant, proving a salt-robust
    contact affinity cannot close the NMR gap (see module docstring).
    """
    rows = []
    for eps in eps_grid:
        K = fit_K(pH, eps)
        rows.append((eps, K, K * np.exp(eps), nmr_f(pH, K, eps), _zeta_chi2(pH, K, eps)))
    return rows


# --- probe 2: ITC ----------------------------------------------------------------
def fit_itc(K: float, eps_ne: float = EPS_NE, pH: float = 4.0, ionic: float = I_ITC):
    """Transfer (K, eps_ne) to the pH-4 thermogram; fit only area + dH. Returns (area, dH, rms)."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "itc"))
    import itc_data  # noqa: E402
    import itc_fit  # noqa: E402

    titr = itc_data.load(pH)
    grid = np.logspace(-9, -2, 24)
    cov = np.array([crpb_state(pH, ionic, c, K, eps_ne).gamma for c in grid])

    def model(_x, area, dH):
        return itc_fit._heat_uJ(titr, itc_fit._gvdw_bound(titr, area, grid, cov), dH)

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


# --- probe 3: NMR ----------------------------------------------------------------
# Measured 31P DeltaR (s^-1): (free PA_aq, OB-PA); f_req = (obs-free)/(Rb-free).
NMR_DELTA_R = {4.0: (162.0, 559.0), 6.0: (130.0, 430.0), 10.0: (5.1, 52.0)}


def nmr_f(
    pH: float,
    K: float,
    eps_ne: float = EPS_NE,
    ionic: float = 0.15,
    pa_mM: float = 7.6,
    phi: float = 0.5,
    d43_um: float = 3.0,
) -> float:
    """Transfer (K, eps_ne) to the NMR: bound fraction f = Gamma a / n_B at 0.15 M, 7.6 mM."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "nmr"))
    from nmr_binding_bridge import specific_area  # noqa: E402

    a = specific_area(phi, d43_um)  # 1/A
    gamma = crpb_state(pH, ionic, pa_mM * 1e-3, K, eps_ne).gamma
    n_B = pa_mM * 1e-3 * N_A / 1e27  # 1/A^3
    return gamma * a / n_B


def nmr_f_req(pH: float, Rb: float = 3e4) -> float:
    free, obs = NMR_DELTA_R[pH]
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
    """Overcharge curve (zeta vs [PA], CR vs gvdW vs data) + the three-probe summary."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _, (zeta_6a_gvdw, _) = gv.evaluate_global()
    colors = {4.0: "#1f4fd8", 6.0: "#e08a00", 10.0: "#2ca02c"}
    fig, ax = plt.subplots(figsize=(7, 5))
    pa_fine = np.geomspace(0.05, 8, 60)
    for pH in PHS:
        c, K = colors[pH], K_by_pH[pH]
        x, y, e = _zeta_data(pH)
        ax.errorbar(x, y, yerr=e, fmt="x", color=c, ms=6, capsize=2, label=f"pH {pH:g} (data)")
        ax.plot(pa_fine, [zeta_pa(pH, v, K, eps_ne) for v in pa_fine], "-", color=c, lw=2)
        ax.plot(x, [zeta_6a_gvdw(pH, v) for v in x], "^", mfc="white", mec=c, mew=1.2, ms=7)
    ax.axhline(0, color="gray", lw=0.5)
    ax.set_xlabel("[PA] (mM)")
    ax.set_ylabel(r"$\zeta$ (mV)")
    ax.set_ylim(-75, 45)
    ks = ", ".join(f"K{int(p)}={K_by_pH[p]:.2g}" for p in PHS)
    ax.set_title(
        f"CR-PB (lines) vs gvdW (open triangles) vs data (x).\n"
        f"per-pH {ks} /M (eps_ne={eps_ne:.1f} kT); CR-PB overcharges at pH 4, gvdW cannot"
    )
    ax.legend(fontsize=8, frameon=False, ncol=3)
    fig.tight_layout()
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
        f"eps_ne = {eps_ne:g} kT)\n"
    )

    print("non-electrostatic contact affinity is DEGENERATE with K at each pH (pH 4 shown):")
    print(f"  {'eps_ne':>7} {'K(/M)':>9} {'K*exp(eps)':>11} {'f_NMR pH4':>10} {'zeta_chi2':>10}")
    for eps, Ke, Kexp, f4, chi2 in eps_degeneracy(4.0):
        print(f"  {eps:>7.1f} {Ke:>9.3g} {Kexp:>11.3g} {f4:>10.3%} {chi2:>10.0f}")

    print("\nsanity limits:")
    for name, ok in _sanity(K_by_pH, eps_ne):
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")

    # sigma0 anchor cross-check: the CR model anchors sigma0 to the *measured* bare zeta
    # (monopole), while feeding the *MC* dipolar psi_beta for binding. The MC's own bare
    # zeta (mc/zeta.csv, shear-plane rule) reveals why they must be kept separate.
    mcz = gv._table("zeta.csv")
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

    print("\nprobe 2 -- ITC (pH 4, transfer K(pH 4), fit area + dH):")
    area, dH, rms = fit_itc(K_by_pH[4.0], eps_ne)
    print(f"  CR-PB : area={area:.0f} A^2/oleosin  dH={dH:.1f} kJ/mol  RMS={rms:.1f} uJ")
    print("  gvdW  : area=2945 A^2/oleosin  dH=-33.7 kJ/mol  RMS=7.4 uJ  (for reference)")

    print("\nprobe 3 -- NMR bound fraction (0.15 M, 7.6 mM, transfer per-pH K):")
    print(f"  {'pH':>4} {'f_CR':>8} {'f_req':>8} {'ratio':>7}")
    for pH in PHS:
        f_cr, f_rq = nmr_f(pH, K_by_pH[pH], eps_ne), nmr_f_req(pH)
        ratio = f"{f_rq / f_cr:>5.1f}x" if f_cr > 0 else "  n/a"  # K(pH10)=0: no binding
        print(f"  {pH:>4g} {f_cr:>8.2%} {f_rq:>8.2%} {ratio:>7}")
    print("  (pH 10: zeta gives no binding signal -> K=0; f_req 0.16% is near the NMR floor)")

    figure = Path(__file__).parent / "fig_crpb.png"
    plot_crpb(K_by_pH, eps_ne, figure)
    print(f"\nsaved {figure.name}")
