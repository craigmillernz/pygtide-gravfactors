import numpy as np
import pytest

from gravfactors import compute


def test_scalar_input_gives_one_row_of_sixteen():
    result = compute(45.0, 10.0, 100.0)
    assert result.amplitude_factor.shape == (16,)
    assert result.phase_lead_deg.shape == (16,)
    assert result.freq_min_cpd.shape == (16,)
    assert len(result.names) == 16


def test_array_input_gives_points_by_constituents():
    result = compute([0.0, 45.0, -60.0], [0.0, 10.0, 20.0])
    assert result.amplitude_factor.shape == (3, 16)
    assert result.phase_lead_deg.shape == (3, 16)


def test_array_row_matches_scalar_call():
    many = compute([10.0, 45.0, 80.0], [0.0, 0.0, 0.0])
    one = compute(45.0, 0.0)
    np.testing.assert_allclose(many.amplitude_factor[1], one.amplitude_factor)
    np.testing.assert_allclose(many.phase_lead_deg[1], one.phase_lead_deg)


def test_amplitude_factors_are_physically_plausible():
    result = compute(45.0, 0.0)
    factors = dict(zip(result.names, result.amplitude_factor))
    assert 1.15 < factors["M2"] < 1.17
    assert 1.14 < factors["O1"] < 1.17
    assert 1.12 < factors["K1"] < 1.15  # NDFW resonance lowers K1 relative to O1
    assert factors["K1"] < factors["O1"]
    assert 1.05 < factors["M3"] < 1.09
    assert 1.15 < factors["Mf"] < 1.17


def test_elastic_model_has_zero_phase_and_differs_from_inelastic():
    elastic = compute(45.0, 0.0, model="elastic")
    inelastic = compute(45.0, 0.0, model="inelastic")
    np.testing.assert_array_equal(elastic.phase_lead_deg, 0.0)
    assert np.any(inelastic.phase_lead_deg != 0.0)
    assert not np.allclose(elastic.amplitude_factor, inelastic.amplitude_factor)


def test_phase_lag_is_negative_of_phase_lead():
    result = compute(45.0, 0.0)
    np.testing.assert_allclose(result.phase_lag_deg, -result.phase_lead_deg)


def test_latitude_dependence_is_present_but_small():
    equator = compute(0.0, 0.0).amplitude_factor
    pole = compute(90.0, 0.0).amplitude_factor
    difference = np.abs(equator - pole)
    assert difference.max() > 1e-4
    assert difference.max() < 0.01


def test_degree_three_has_no_latitude_dependence():
    equator = compute(0.0, 0.0)
    pole = compute(90.0, 0.0)
    assert equator.amplitude_factor[-1] == pole.amplitude_factor[-1]


def test_longitude_and_height_do_not_change_the_result():
    base = compute(30.0, 0.0)
    moved = compute(30.0, 123.0, 4500.0)
    np.testing.assert_array_equal(base.amplitude_factor, moved.amplitude_factor)
    np.testing.assert_array_equal(base.phase_lead_deg, moved.phase_lead_deg)


def test_longitude_and_height_are_optional():
    np.testing.assert_array_equal(
        compute(30.0).amplitude_factor, compute(30.0, 0.0, 0.0).amplitude_factor
    )


def test_bands_are_reported_per_constituent():
    result = compute(45.0, 0.0)
    assert np.all(result.freq_min_cpd < result.freq_cpd)
    assert np.all(result.freq_cpd < result.freq_max_cpd)
    assert result.species[0] == "long-period" and result.species[-1] == "terdiurnal"


def test_records_contain_one_row_per_point_and_constituent():
    result = compute([10.0, 20.0], [0.0, 5.0])
    records = result.to_records()
    assert len(records) == 32
    assert {"point", "latitude", "longitude", "species", "constituent", "freq_min_cpd",
            "freq_max_cpd", "amplitude_factor", "phase_lead_deg", "phase_lag_deg"} <= set(records[0])


def test_pygtide_wavegroups_for_single_point():
    groups = compute(45.0, 0.0).to_pygtide_wavegroups()
    assert groups.shape == (16, 4)
    assert np.all(groups[:, 0] < groups[:, 1])
    assert np.all(groups[1:, 0] >= groups[:-1, 1])
    assert np.all(groups[:, 2] > 0)


def test_pygtide_wavegroups_for_multiple_points_need_an_index():
    result = compute([10.0, 45.0], [0.0, 0.0])
    with pytest.raises(ValueError, match="point"):
        result.to_pygtide_wavegroups()
    groups = result.to_pygtide_wavegroups(point=1)
    np.testing.assert_allclose(groups, compute(45.0, 0.0).to_pygtide_wavegroups())


@pytest.mark.parametrize("latitude", [91.0, -90.5, np.nan])
def test_invalid_latitude_is_rejected(latitude):
    with pytest.raises(ValueError, match="latitude"):
        compute(latitude, 0.0)


def test_mismatched_input_lengths_are_rejected():
    with pytest.raises(ValueError, match="broadcast"):
        compute([10.0, 20.0], [0.0, 1.0, 2.0])


def test_unknown_model_is_rejected():
    with pytest.raises(ValueError, match="model"):
        compute(45.0, 0.0, model="viscous")
