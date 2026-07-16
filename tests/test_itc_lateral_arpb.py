"""Mass-balance and fit checks for the LC-AR-PB ITC bridge."""

import extract_thermograms as itc_data
import fit_lateral_arpb as fit_lc
import fit_thermogram as itc_fit
import numpy as np
import pytest


@pytest.fixture(scope="module")
def isotherm():
    return fit_lc.build_isotherm(4.0, fit_lc.I_ITC, True)


def test_tabulated_coverage_is_bounded_and_monotone(isotherm):
    assert np.all(isotherm.gamma >= 0.0)
    assert np.all(isotherm.gamma <= fit_lc.lc.GAMMA_MAX)
    assert np.all(np.diff(isotherm.gamma) >= 0.0)


def test_isotherm_depends_on_product_K_times_concentration(isotherm):
    first = isotherm.coverage(1e-5, K=20.0)
    second = isotherm.coverage(2e-6, K=100.0)
    assert first == pytest.approx(second)


def test_bound_profile_satisfies_ligand_mass_balance(isotherm):
    titration = itc_data.load(4.0)
    K, area = 100.0, 5000.0
    bound, free = fit_lc.bound_profile(titration, isotherm, K, area)
    oleosin, ligand_total = itc_fit._dilute(titration)
    gamma = isotherm.coverage(free, K)
    assert np.allclose(bound + free, ligand_total, rtol=1e-10, atol=1e-14)
    assert np.allclose(bound, oleosin * area * gamma, rtol=1e-7, atol=1e-12)


def test_bound_amount_is_bounded_by_total_ligand(isotherm):
    titration = itc_data.load(4.0)
    bound, free = fit_lc.bound_profile(titration, isotherm, K=100.0, area_A2=5000.0)
    _oleosin, ligand_total = itc_fit._dilute(titration)
    assert np.all(bound >= 0.0)
    assert np.all(free >= 0.0)
    assert np.all(bound <= ligand_total)


def test_fit_converges_and_reproduces_thermogram():
    titration = itc_data.load(4.0)
    result = fit_lc.fit_lateral_arpb(titration)
    assert result.success
    assert result.K > 0.0
    assert result.area_A2 > 0.0
    assert result.dH_kJ_mol < 0.0
    assert result.rms_uJ < 10.0
    assert np.all(np.isfinite(result.correlation))
