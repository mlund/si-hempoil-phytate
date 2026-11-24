"""Physics invariants for the peptide-field plus lateral-correlation zeta model."""

import lateral_arpb as lc
import mc_shared as mcs
import numpy as np
import pytest


def test_binding_field_is_peptide_only_column():
    profiles = mcs._table("potential_profiles.csv")
    rows = profiles[(profiles.pH == 4.0) & (profiles.ionic_strength == 0.005)].sort_values("r")
    expected = np.interp(lc.R_BETA, rows.r, rows.potential_peptide_mV)
    total = np.interp(lc.R_BETA, rows.r, rows.potential_all_mV)
    assert lc.peptide_potential(4.0, 0.005) == pytest.approx(expected)
    assert lc.peptide_potential(4.0, 0.005) != pytest.approx(total)


def test_valence_is_mc_charge_at_binding_coordinate():
    profiles = mcs._table("window_profiles.csv")
    rows = profiles[(profiles.pH == 6.0) & (profiles.ionic_strength == 0.005)].sort_values("r")
    expected = -np.interp(lc.R_BETA, rows.r, rows.q_ino)
    assert lc.phytate_valence(6.0, 0.005) == pytest.approx(expected)


def test_profile_inputs_log_interpolate_in_salt():
    ionic = np.sqrt(0.005 * 0.03)
    expected = 0.5 * (lc.peptide_potential(4.0, 0.005) + lc.peptide_potential(4.0, 0.03))
    assert lc.peptide_potential(4.0, ionic) == pytest.approx(expected)


@pytest.mark.parametrize("pH", lc.PHS)
def test_dilute_limit_is_bare_zeta(pH):
    state = lc.lateral_state(pH, lc.I_ANCHOR, 0.0, K=100.0)
    assert state.gamma == 0.0
    assert state.zeta == pytest.approx(mcs.ZETA_BARE[pH])


def test_pb_increment_has_zero_value_and_slope_at_bare_state():
    assert lc._incremental_pb_free_energy(0.0, 4.0, 0.005, 7.0) == pytest.approx(0.0)
    step = 1e-8
    slope = (
        lc._incremental_pb_free_energy(step, 4.0, 0.005, 7.0)
        - lc._incremental_pb_free_energy(-step, 4.0, 0.005, 7.0)
    ) / (2.0 * step)
    assert slope == pytest.approx(0.0, abs=1e-6)


@pytest.mark.parametrize("ionic", [0.005, 0.03, 0.15])
def test_forsman_residual_does_not_add_full_yukawa_repulsion(ionic):
    gamma = 0.25 * lc.GAMMA_MAX
    z = 7.0
    assert lc._forsman_pair_integral(ionic, z, gamma) <= lc._yukawa_pair_integral(ionic, z)
    assert lc.residual_correlation_free_energy(gamma, ionic, z) <= 0.0


def test_coverage_is_bounded_and_lowers_surface_charge():
    state = lc.lateral_state(4.0, 0.005, 5e-6, K=100.0)
    assert 0.0 <= state.gamma <= lc.GAMMA_MAX
    assert 0.0 <= state.theta <= 1.0
    assert state.sigma_eff == pytest.approx(lc.sigma0(4.0) - state.valence * state.gamma)


def test_effective_binding_potential_adds_only_incremental_pb_response():
    pH, ionic = 4.0, 0.005
    state = lc.lateral_state(pH, ionic, 5e-6, K=100.0)
    expected = (
        lc.peptide_potential(pH, ionic)
        + lc._psi_d(state.sigma_eff, ionic)
        - lc._psi_d(lc.sigma0(pH), ionic)
    )
    assert state.psi_binding_effective == pytest.approx(expected)


def test_finite_activity_has_nonzero_interior_coverage():
    state = lc.lateral_state(6.0, 0.005, 1.5e-7, K=50.0)
    assert 0.0 < state.gamma < lc.GAMMA_MAX


def test_global_solution_does_not_exceed_trial_free_energies():
    pH, ionic, concentration, K = 4.0, 0.005, 5e-6, 100.0
    state = lc.lateral_state(pH, ionic, concentration, K)
    equilibrium = lc.surface_grand_potential(state.gamma, pH, ionic, concentration, K)
    trials = [
        lc.surface_grand_potential(gamma, pH, ionic, concentration, K)
        for gamma in np.linspace(0.0, lc.GAMMA_MAX, 101)
    ]
    assert equilibrium <= min(trials) + 1e-10


def test_correlations_are_a_residual_switch_not_a_second_pb_term():
    pH, ionic, concentration, K = 4.0, 0.005, 5e-6, 100.0
    gamma = 0.2 * lc.GAMMA_MAX
    with_corr = lc.surface_grand_potential(gamma, pH, ionic, concentration, K)
    without_corr = lc.surface_grand_potential(
        gamma,
        pH,
        ionic,
        concentration,
        K,
        include_correlations=False,
    )
    expected = lc.residual_correlation_free_energy(
        gamma,
        ionic,
        lc.phytate_valence(pH, ionic),
    )
    assert with_corr - without_corr == pytest.approx(expected)
