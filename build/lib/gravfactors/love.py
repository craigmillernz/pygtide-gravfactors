"""Complex degree-2 and degree-3 Love numbers for one tidal frequency.

Inelastic model: IERS Conventions (2010) Ch. 7 (h) and Ch. 6 (k), including the
nearly diurnal free wobble resonance in the diurnal band and mantle anelasticity.

Elastic model: derived here, because IERS tabulates no elastic h. The mantle
dispersion term of Ch. 7 Eq. 4a is removed from the real part of h (evaluated
between the tidal period and a 1 s reference period), imaginary parts are set
to zero, and the elastic k of Ch. 6 Table 6.3 is used with the same wobble
resonance shape as the inelastic model.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from . import _data as d

MODELS = ("elastic", "inelastic")
_SECONDS_PER_DAY = 86400.0


@dataclass(frozen=True)
class LoveNumbers:
    """h(0) and h(2) (latitude term) displacement and k potential Love numbers."""

    h0: complex
    h2: complex
    k: complex


def _resonance(freq_cpd: float, coefficients: tuple[complex, tuple[complex, ...]]) -> complex:
    l0, terms = coefficients
    sigma = freq_cpd / d.SIDEREAL_RATIO
    return l0 + sum(term / (sigma - s) for term, s in zip(terms, d.RESONANCE_FREQS))


def _dispersion_ratio(freq_cpd: float) -> float:
    """X = (f_m / f)**alpha, with f_m the frequency of the 200 s reference period."""
    f_m = _SECONDS_PER_DAY / d.REFERENCE_PERIOD_S
    return (f_m / freq_cpd) ** d.ANELASTIC_ALPHA


def _cot_alpha() -> float:
    return 1.0 / math.tan(d.ANELASTIC_ALPHA * math.pi / 2.0)


def _zonal_h(freq_cpd: float) -> complex:
    """IERS Ch. 7 Eq. 4a: anelastic h20 across the zonal band."""
    x = _dispersion_ratio(freq_cpd)
    return d.H20_REFERENCE + d.H20_DISPERSION * complex(_cot_alpha() * (1.0 - x), x)


def _inelastic(freq_cpd: float, order: int) -> LoveNumbers:
    if order == 0:
        return LoveNumbers(_zonal_h(freq_cpd), d.H2_LONG_PERIOD, d.K_ANELASTIC[0])
    if order == 1:
        return LoveNumbers(
            _resonance(freq_cpd, d.H0_DIURNAL),
            _resonance(freq_cpd, d.H2_DIURNAL),
            _resonance(freq_cpd, d.K21_DIURNAL),
        )
    return LoveNumbers(d.H_SEMIDIURNAL, d.H2_SEMIDIURNAL, d.K_ANELASTIC[2])


def _elastic(freq_cpd: float, order: int) -> LoveNumbers:
    anelastic = _inelastic(freq_cpd, order)
    x = _dispersion_ratio(freq_cpd)
    x_elastic = _dispersion_ratio(_SECONDS_PER_DAY / d.ELASTIC_PERIOD_S)
    removed = d.H20_DISPERSION * _cot_alpha() * (x_elastic - x)
    k_offset = d.K_ANELASTIC[order].real - d.K_ELASTIC[order]
    return LoveNumbers(
        complex(anelastic.h0.real - removed, 0.0),
        complex(anelastic.h2.real, 0.0),
        complex(anelastic.k.real - k_offset, 0.0),
    )


def love_numbers(freq_cpd: float, order: int, degree: int, model: str) -> LoveNumbers:
    """Love numbers for a tide of frequency ``freq_cpd`` (cycles per solar day)."""
    if model not in MODELS:
        raise ValueError(f"model must be one of {MODELS}, got {model!r}")
    if degree == 3:
        return LoveNumbers(complex(d.H3), 0j, complex(d.K3))
    if model == "elastic":
        return _elastic(freq_cpd, order)
    return _inelastic(freq_cpd, order)
