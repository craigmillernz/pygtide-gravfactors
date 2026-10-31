import pytest

from gravfactors.love import love_numbers

# IERS Conventions 2010, Ch. 7, Table 7.2 (h(0) real, imaginary) at diurnal frequencies.
TABLE_7_2_DIURNAL = {
    0.85461: (0.6039, -0.0027),  # 2Q1
    0.89080: (0.6036, -0.0026),  # Q1
    0.92700: (0.6028, -0.0025),  # O1
    0.99454: (0.5817, -0.0011),  # P1
    1.00000: (0.5236, 0.0030),  # K1
    1.03619: (0.6108, -0.0030),  # J1
}

# Table 7.2 long-period rows: (frequency in cpsd, h(0) real, imaginary).
TABLE_7_2_ZONAL = [
    (0.005461, 0.6182, -0.0054),  # Ssa
    (0.036193, 0.6126, -0.0041),  # Mm
    (0.073002, 0.6109, -0.0037),  # Mf
]

SIDEREAL = 1.002737909


@pytest.mark.parametrize("cpsd", list(TABLE_7_2_DIURNAL))
def test_diurnal_h_matches_iers_table_7_2(cpsd):
    # Table 7.2 also contains ocean-loading and exact-vs-formula corrections (up to
    # ~3e-3 at K1) that the pure body-tide resonance formula does not.
    real, imag = TABLE_7_2_DIURNAL[cpsd]
    h0 = love_numbers(cpsd * SIDEREAL, 1, 2, "inelastic").h0
    assert h0.real == pytest.approx(real, abs=3e-3)
    assert h0.imag == pytest.approx(imag, abs=3e-3)


@pytest.mark.parametrize("cpsd, real, imag", TABLE_7_2_ZONAL)
def test_zonal_h_matches_iers_table_7_2(cpsd, real, imag):
    h0 = love_numbers(cpsd * SIDEREAL, 0, 2, "inelastic").h0
    assert h0.real == pytest.approx(real, abs=3e-4)
    assert h0.imag == pytest.approx(imag, abs=3e-4)


def test_semidiurnal_nominal_values():
    love = love_numbers(1.9322736, 2, 2, "inelastic")
    assert love.h0 == pytest.approx(0.6078 - 0.0022j)
    assert love.k == pytest.approx(0.30102 - 0.00130j)
    assert love.h2 == pytest.approx(-0.0006)


def test_degree_three_is_real_and_constant():
    love = love_numbers(2.8984104, 3, 3, "inelastic")
    assert love.h0 == pytest.approx(0.292)
    assert love.k == pytest.approx(0.093)
    assert love.h0.imag == 0.0 and love.k.imag == 0.0


def test_ndfw_resonance_lowers_k1_love_number():
    o1 = love_numbers(0.9295357, 1, 2, "inelastic")
    k1 = love_numbers(1.0027379, 1, 2, "inelastic")
    assert k1.h0.real < o1.h0.real - 0.05
    assert k1.k.real < o1.k.real - 0.02


@pytest.mark.parametrize("freq, order", [(0.0732020, 0), (0.9295357, 1), (1.0027379, 1), (1.9322736, 2)])
def test_elastic_love_numbers_are_real(freq, order):
    love = love_numbers(freq, order, 2, "elastic")
    assert love.h0.imag == 0.0
    assert love.k.imag == 0.0


def test_elastic_differs_from_inelastic_in_phase_only_slightly_in_amplitude():
    elastic = love_numbers(1.9322736, 2, 2, "elastic")
    inelastic = love_numbers(1.9322736, 2, 2, "inelastic")
    assert abs(elastic.h0.real - inelastic.h0.real) < 0.02
    assert elastic.k.real == pytest.approx(0.29801)


def test_unknown_model_is_rejected():
    with pytest.raises(ValueError, match="model"):
        love_numbers(1.9322736, 2, 2, "viscous")
