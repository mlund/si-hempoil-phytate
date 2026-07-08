"""Sanity and unit tests for the charge-regulated Poisson-Boltzmann (CR-PB) model.

These run from the committed CSV layer (no ``mc/runs/`` needed): ``uv run pytest``.
They lock the physical invariants of the model (Grahame consistency, the dilute and
saturation limits, monotonicity, overcharging, the K/eps_ne degeneracy) and regress the
reported paper numbers (per-pH K, the NMR bound fractions, the ITC "good fit / unphysical
area" result).
"""

import math

import charge_regulation as cr
import mc_shared as mcs
import pytest

ANCHOR = cr.I_ANCHOR  # low-salt reference ionic strength (0.005 M)
GAMMA_MAX = cr.GAMMA_MAX
K4 = 64.0  # representative pH-4 intrinsic affinity (cf. fit_K)


# --- Debye / Grahame double-layer utilities --------------------------------------


def test_debye_length_scales_as_sqrt_ionic():
    # kappa ~ sqrt(I): a 10x salt increase multiplies kappa by sqrt(10).
    assert mcs.kappa(0.15) / mcs.kappa(0.015) == pytest.approx(math.sqrt(10), rel=1e-6)


def test_debye_length_magnitude_at_physiological_salt():
    # 1/kappa at 0.15 M 1:1 salt is the standard aqueous Debye length, ~7.9 A.
    debye_A = 1.0 / mcs.kappa(0.15)
    assert 7.0 < debye_A < 8.5


@pytest.mark.parametrize("zeta_mV", [-80.0, -30.0, -5.0, 5.0, 30.0, 80.0])
@pytest.mark.parametrize("ionic", [0.005, 0.05, 0.15])
def test_grahame_roundtrip(zeta_mV, ionic):
    # psi_d(sigma_from_zeta(zeta)) == zeta -- Grahame forward/inverse consistency.
    sigma = mcs.sigma_from_zeta(zeta_mV, ionic)
    assert cr._psi_d(sigma, ionic) == pytest.approx(zeta_mV, rel=1e-6, abs=1e-6)


def test_sigma_from_zeta_sign_and_monotone():
    assert mcs.sigma_from_zeta(20.0, 0.05) > 0
    assert mcs.sigma_from_zeta(-20.0, 0.05) < 0
    assert mcs.sigma_from_zeta(40.0, 0.05) > mcs.sigma_from_zeta(20.0, 0.05)


# --- CR-PB state invariants ------------------------------------------------------


def test_dilute_limit_is_bare_zeta():
    # c_L -> 0: no coverage, zeta -> the bare Grahame value at that salt.
    s = cr.crpb_state(4.0, ANCHOR, 1e-18, K4)
    bare = cr._psi_d(cr.sigma0(4.0), ANCHOR)
    assert s.gamma < 1e-9
    assert s.zeta == pytest.approx(bare, abs=0.1)


@pytest.mark.parametrize("c_L", [1e-8, 1e-6, 1e-4, 1e-2, 1e2])
def test_coverage_within_bounds(c_L):
    s = cr.crpb_state(4.0, ANCHOR, c_L, K4)
    assert 0.0 <= s.theta <= 1.0
    assert 0.0 <= s.gamma <= GAMMA_MAX


def test_coverage_monotone_in_concentration():
    lo = cr.crpb_state(4.0, ANCHOR, 1e-6, K4)
    hi = cr.crpb_state(4.0, ANCHOR, 1e-4, K4)
    assert hi.gamma > lo.gamma


def test_self_consistent_fixed_point():
    # The returned state satisfies both closure relations exactly.
    s = cr.crpb_state(6.0, ANCHOR, 5e-6, 2.9e4)
    z = mcs.phytate_valence(6.0, ANCHOR)
    assert s.sigma_eff == pytest.approx(cr.sigma0(6.0) - z * s.gamma, abs=1e-12)
    assert s.zeta == pytest.approx(cr._psi_d(s.sigma_eff, ANCHOR), rel=1e-9)


def test_regulation_caps_coverage_below_gamma_max():
    # Even at huge c_L the regulated coverage self-limits below the geometric ceiling
    # (a pure Langmuir would saturate at gamma_max; regulation does not).
    s = cr.crpb_state(4.0, ANCHOR, 1e6, K4)
    assert 0.0 < s.gamma < GAMMA_MAX


def test_ph4_overcharge_reverses_zeta():
    # Bare pH-4 zeta is positive; enough phytate inverts sigma_eff and hence zeta.
    bare = cr._psi_d(cr.sigma0(4.0), ANCHOR)
    assert bare > 0
    s = cr.crpb_state(4.0, ANCHOR, 1e6, K4)
    assert s.sigma_eff < 0
    assert s.zeta < 0


@pytest.mark.parametrize("pH", [4.0, 6.0, 10.0])
def test_high_salt_screens_zeta_magnitude(pH):
    lo = abs(cr._psi_d(cr.sigma0(pH), 0.005))
    hi = abs(cr._psi_d(cr.sigma0(pH), 0.5))
    assert hi < lo


@pytest.mark.parametrize("eps", [1.0, 2.5])
def test_eps_ne_degenerate_with_K(eps):
    # K and eps_ne enter only as the product K*exp(eps_ne): (K, eps) == (K*exp(eps), 0).
    a = cr.crpb_state(4.0, ANCHOR, 1e-5, K4, eps_ne=eps)
    b = cr.crpb_state(4.0, ANCHOR, 1e-5, K4 * math.exp(eps), eps_ne=0.0)
    assert a.gamma == pytest.approx(b.gamma, rel=1e-6)
    assert a.zeta == pytest.approx(b.zeta, rel=1e-6)


# --- regression: fitted parameters and probe outputs (paper numbers) -------------


def test_fit_K_per_pH_regression():
    K = cr.fit_K_by_pH()
    assert K[4.0] == pytest.approx(64.0, rel=0.25)
    assert K[6.0] == pytest.approx(2.9e4, rel=0.3)
    assert K[10.0] == 0.0  # pH 10 is repulsive -> no direct binding


@pytest.mark.parametrize("pH,f_expected", [(4.0, 0.0130), (6.0, 0.0105), (10.0, 0.0011)])
def test_nmr_bound_fraction_regression(pH, f_expected):
    assert cr.nmr_f(pH) == pytest.approx(f_expected, rel=0.1)


@pytest.mark.parametrize("pH", [4.0, 6.0])
def test_nmr_matches_two_site_requirement(pH):
    # The simulated bound fraction matches the fraction the two-site exchange requires.
    assert cr.nmr_f(pH) == pytest.approx(cr.nmr_f_req(pH), rel=0.2)


def test_itc_fit_is_good_but_area_unphysical():
    # The charge-regulated isotherm fits the pH-4 thermogram well, yet the fitted area is
    # ~20x the physical grafting density and K/area are near-degenerate -- a good fit that
    # is not evidence for monolayer binding.
    _K, area, _dH, rms, corr = cr.fit_itc()
    assert area / cr.AREA_OLEOSIN > 10.0  # >> physical 338 A^2/oleosin
    assert corr < -0.8  # K and area nearly degenerate
    assert rms < 10.0  # still fits the thermogram


# --- Gouy-Chapman limits ---------------------------------------------------------


@pytest.mark.parametrize("pH", [4.0, 6.0, 10.0])
@pytest.mark.parametrize("ionic", [0.005, 0.05, 0.15])
def test_no_binding_reduces_to_plain_gouy_chapman(pH, ionic):
    # With no adsorption (K = 0) the surface charge is unregulated, so the model collapses
    # to the plain Gouy-Chapman diffuse layer: zeta = psi_d(sigma0), at any [PA].
    s = cr.crpb_state(pH, ionic, 1e-3, K=0.0)
    assert s.gamma == 0.0
    assert s.sigma_eff == pytest.approx(cr.sigma0(pH))
    assert s.zeta == pytest.approx(cr._psi_d(cr.sigma0(pH), ionic), rel=1e-9)


@pytest.mark.parametrize("ionic", [0.005, 0.05, 0.15])
def test_low_charge_recovers_debye_huckel(ionic):
    # Gouy-Chapman linearises to Debye-Hueckel at low charge: psi_d -> 4 pi lB sigma / kappa.
    sigma = 1e-5  # e/A^2, deep in the linear regime
    linear_mV = 4.0 * math.pi * cr.LB * sigma / mcs.kappa(ionic) * cr.KT_MV
    assert cr._psi_d(sigma, ionic) == pytest.approx(linear_mV, rel=1e-3)
