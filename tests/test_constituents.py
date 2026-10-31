import pytest

from gravfactors.constituents import BANDS, CONSTITUENTS

KNOWN_FREQ_CPD = {
    "Sa": 0.0027379,
    "Ssa": 0.0054758,
    "Mm": 0.0362916,
    "Mf": 0.0732020,
    "2Q1": 0.8569527,
    "Q1": 0.8932441,
    "O1": 0.9295357,
    "P1": 0.9972621,
    "K1": 1.0027379,
    "J1": 1.0390293,
    "OO1": 1.0759398,
    "N2": 1.8959821,
    "M2": 1.9322736,
    "S2": 2.0,
    "K2": 2.0054758,
    "M3": 2.8984104,
}


def test_there_are_sixteen_unique_constituents():
    names = [c.name for c in CONSTITUENTS]
    assert len(names) == 16
    assert len(set(names)) == 16


@pytest.mark.parametrize("constituent", CONSTITUENTS, ids=lambda c: c.name)
def test_frequencies_match_known_values(constituent):
    assert constituent.freq_cpd == pytest.approx(KNOWN_FREQ_CPD[constituent.name], abs=2e-6)


def test_constituents_are_sorted_by_frequency():
    freqs = [c.freq_cpd for c in CONSTITUENTS]
    assert freqs == sorted(freqs)


def test_each_frequency_lies_in_its_own_band():
    for constituent, (low, high) in zip(CONSTITUENTS, BANDS):
        assert low < constituent.freq_cpd < high


def test_bands_are_increasing_and_do_not_overlap():
    for (_, previous_high), (low, _) in zip(BANDS, BANDS[1:]):
        assert low >= previous_high


def test_bands_do_not_cross_species_boundaries():
    for constituent, (low, high) in zip(CONSTITUENTS, BANDS):
        assert low >= constituent.order - 0.5 - 1e-9
        assert high <= constituent.order + 0.5 + 1e-9


def test_species_labels():
    species = {c.name: c.species for c in CONSTITUENTS}
    assert species["Mf"] == "long-period"
    assert species["K1"] == "diurnal"
    assert species["M2"] == "semidiurnal"
    assert species["M3"] == "terdiurnal"
