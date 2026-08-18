"""Mass-balance and fit checks for the LC-AR-PB ITC bridge."""

import dataclasses

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


def test_first_injection_is_excluded_from_residuals_but_still_dilutes_the_cell():
    titration = itc_data.load(4.0)
    mask = titration.fit_mask
    assert not mask[0]
    assert np.all(mask[itc_data.SKIP_INJECTIONS :])

    # The skipped injection still delivers phytate, so the cell composition at every
    # later injection must reflect it.  Check the second injection against the
    # two-step overflow balance written out by hand.
    oleosin, ligand_total = itc_fit._dilute(titration)
    volume = titration.vcell_uL * 1e-6
    first, second = titration.vinj_uL[:2] * 1e-6
    after_first = titration.csyr_M * first / (volume + first)
    after_second = (after_first * volume + titration.csyr_M * second) / (volume + second)
    assert ligand_total[1] == pytest.approx(after_second, rel=1e-12)

    # Had the first injection not happened, less phytate would be present.
    alone = titration.csyr_M * second / (volume + second)
    assert ligand_total[1] > alone
    assert np.all(np.diff(oleosin) < 0.0)


def test_fit_ignores_a_corrupted_first_injection_heat():
    titration = itc_data.load(4.0)
    corrupted = dataclasses.replace(titration, heat_uJ=titration.heat_uJ.copy())
    corrupted.heat_uJ[0] += 1000.0
    assert itc_fit.fit_onesite(corrupted) == pytest.approx(itc_fit.fit_onesite(titration))


def test_onesite_reproduces_nanoanalyze_fit():
    """Dropping the first injection matches NanoAnalyze's own independent-model fit."""
    titration = itc_data.load(4.0)
    Ka, dH, n, _rms = itc_fit.fit_onesite(titration)
    reference = titration.na_fit
    assert Ka == pytest.approx(reference["Ka"], rel=0.05)
    assert n == pytest.approx(reference["n"], rel=0.05)
    assert dH == pytest.approx(reference["dH"], rel=0.05)
