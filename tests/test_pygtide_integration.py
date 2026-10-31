"""End-to-end check: PyGTide output scaled by our wave groups reproduces our factors."""

import numpy as np
import pytest

from gravfactors import compute

pygtide = pytest.importorskip("pygtide")

LAT, LON, HEIGHT = 45.0, 10.0, 0.0
WAVES = {"M2": 1.9322736, "S2": 2.0, "O1": 0.9295357, "K1": 1.0027379}
START, HOURS = "2020-01-01", 24 * 365


def _series(wavegroups):
    tide = pygtide.pygtide()
    tide.set_wavegroup(wavegroups)
    tide.predict(LAT, LON, HEIGHT, START, HOURS, 3600, tidalcompo=0)
    frame = tide.results()
    time = (frame["UTC"] - frame["UTC"].iloc[0]).dt.total_seconds().to_numpy() / 86400.0
    return time, frame["Tide [nm/s**2]"].to_numpy()


def _fit(time, signal):
    """Harmonic projection; constituents are separated by whole cycles per year."""
    amplitudes = []
    for freq in WAVES.values():
        phase = 2 * np.pi * freq * time
        cos_part = 2.0 * np.mean(signal * np.cos(phase))
        sin_part = 2.0 * np.mean(signal * np.sin(phase))
        amplitudes.append(cos_part + 1j * sin_part)
    return np.array(amplitudes)


@pytest.mark.parametrize("model", ["inelastic", "elastic"])
def test_wavegroups_scale_pygtide_by_our_factors(model):
    result = compute(LAT, LON, HEIGHT, model=model)
    groups = result.to_pygtide_wavegroups()
    # PyGTide normalises each group by its main wave, so the reference run must use
    # identical group boundaries with unit factor and zero phase.
    unit_groups = groups.copy()
    unit_groups[:, 2:] = [1.0, 0.0]
    time, reference = _series(unit_groups)
    _, scaled = _series(groups)
    ratio = _fit(time, scaled) / _fit(time, reference)
    index = [result.names.index(name) for name in WAVES]
    expected = result.amplitude_factor[index]
    assert np.abs(ratio) == pytest.approx(expected, abs=2e-3)
    # cos(wt + phi) with phase lead phi gives a fitted cos+i*sin ratio of exp(-i*phi)
    lead = np.degrees(-np.angle(ratio))
    assert lead == pytest.approx(result.phase_lead_deg[index], abs=0.05)


def test_positive_wavegroup_phase_is_a_lead_in_pygtide():
    groups = compute(LAT, LON, HEIGHT).to_pygtide_wavegroups()
    groups[:, 2:] = [1.0, 0.0]
    time, reference = _series(groups)
    groups[:, 3] = 30.0
    _, shifted = _series(groups)
    lead = np.degrees(-np.angle(_fit(time, shifted) / _fit(time, reference)))
    assert lead == pytest.approx(30.0, abs=0.2)
