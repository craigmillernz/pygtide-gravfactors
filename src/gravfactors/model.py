"""Gravimetric amplitude factors and phase leads for the 16 constituents."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .constituents import BANDS, CONSTITUENTS
from .love import MODELS, love_numbers

logger = logging.getLogger(__name__)

# WGS84 first eccentricity squared, used to convert geodetic to geocentric latitude.
_E2 = 0.00669437999014


@dataclass(frozen=True)
class GravimetricFactors:
    """Result of :func:`compute`.

    ``amplitude_factor``, ``phase_lead_deg`` and ``phase_lag_deg`` have shape
    ``(16,)`` for a single point and ``(n_points, 16)`` otherwise. Frequencies,
    names and species always have shape ``(16,)``.
    """

    model: str
    latitude: NDArray[np.float64]
    longitude: NDArray[np.float64]
    height: NDArray[np.float64]
    names: tuple[str, ...]
    species: tuple[str, ...]
    freq_cpd: NDArray[np.float64]
    freq_min_cpd: NDArray[np.float64]
    freq_max_cpd: NDArray[np.float64]
    amplitude_factor: NDArray[np.float64]
    phase_lead_deg: NDArray[np.float64]

    @property
    def phase_lag_deg(self) -> NDArray[np.float64]:
        return -self.phase_lead_deg

    @property
    def n_points(self) -> int:
        return self.latitude.size

    def to_records(self) -> list[dict[str, Any]]:
        """One dict per (point, constituent); ``pandas.DataFrame(records)`` works directly."""
        amplitude = np.atleast_2d(self.amplitude_factor)
        lead = np.atleast_2d(self.phase_lead_deg)
        return [
            {
                "point": i,
                "latitude": float(self.latitude[i]),
                "longitude": float(self.longitude[i]),
                "species": self.species[j],
                "constituent": self.names[j],
                "freq_min_cpd": float(self.freq_min_cpd[j]),
                "freq_max_cpd": float(self.freq_max_cpd[j]),
                "amplitude_factor": float(amplitude[i, j]),
                "phase_lead_deg": float(lead[i, j]),
                "phase_lag_deg": float(-lead[i, j]),
            }
            for i in range(self.n_points)
            for j in range(len(self.names))
        ]

    def to_pygtide_wavegroups(self, point: int | None = None) -> NDArray[np.float64]:
        """(16, 4) array for ``PyGTide.set_wavegroup``: freq from, freq to, factor, phase lead."""
        if point is None:
            if self.n_points != 1:
                raise ValueError(
                    f"result holds {self.n_points} points; pass point=<index> to select one"
                )
            point = 0
        amplitude = np.atleast_2d(self.amplitude_factor)[point]
        lead = np.atleast_2d(self.phase_lead_deg)[point]
        return np.column_stack([self.freq_min_cpd, self.freq_max_cpd, amplitude, lead])


def _prepare_inputs(
    latitude: ArrayLike, longitude: ArrayLike, height: ArrayLike
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], bool]:
    try:
        lat, lon, hgt = np.broadcast_arrays(
            *(np.asarray(v, dtype=float) for v in (latitude, longitude, height))
        )
    except ValueError as error:
        raise ValueError(f"latitude, longitude and height cannot broadcast together: {error}")
    for name, values in (("latitude", lat), ("longitude", lon), ("height", hgt)):
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{name} must be finite")
    if np.any(np.abs(lat) > 90.0):
        raise ValueError("latitude must be within [-90, 90] degrees")
    return lat.ravel(), lon.ravel(), hgt.ravel(), lat.ndim == 0


def _geocentric_p2(latitude_deg: NDArray[np.float64]) -> NDArray[np.float64]:
    """P2(sin phi') with phi' the geocentric latitude of a point on the WGS84 surface."""
    phi = np.radians(latitude_deg)
    phi_geocentric = np.arctan2((1.0 - _E2) * np.sin(phi), np.cos(phi))
    return (3.0 * np.sin(phi_geocentric) ** 2 - 1.0) / 2.0


def compute(
    latitude: ArrayLike,
    longitude: ArrayLike = 0.0,
    height: ArrayLike = 0.0,
    model: str = "inelastic",
) -> GravimetricFactors:
    """Gravimetric factors and phase leads for 16 constituents at one or many points.

    Parameters
    ----------
    latitude, longitude, height
        Geodetic latitude and longitude in degrees, ellipsoidal height in metres.
        Scalars or broadcastable arrays. The model is axisymmetric and the
        factors do not depend on longitude or height; they are accepted so that
        calls mirror the PyGTide station arguments.
    model
        ``"inelastic"`` (complex Love numbers, nonzero phase) or ``"elastic"``.
    """
    if model not in MODELS:
        raise ValueError(f"model must be one of {MODELS}, got {model!r}")
    lat, lon, hgt, scalar = _prepare_inputs(latitude, longitude, height)
    p2 = _geocentric_p2(lat)
    delta = np.empty((lat.size, len(CONSTITUENTS)), dtype=complex)
    for j, constituent in enumerate(CONSTITUENTS):
        love = love_numbers(constituent.freq_cpd, constituent.order, constituent.degree, model)
        n = constituent.degree
        h = love.h0 + love.h2 * p2
        delta[:, j] = 1.0 + (2.0 / n) * h - ((n + 1.0) / n) * love.k
    amplitude = np.abs(delta)
    phase_lead = np.degrees(np.angle(delta))
    logger.debug("computed %s factors for %d point(s)", model, lat.size)
    if scalar:
        amplitude, phase_lead = amplitude[0], phase_lead[0]
    return GravimetricFactors(
        model=model,
        latitude=lat,
        longitude=lon,
        height=hgt,
        names=tuple(c.name for c in CONSTITUENTS),
        species=tuple(c.species for c in CONSTITUENTS),
        freq_cpd=np.array([c.freq_cpd for c in CONSTITUENTS]),
        freq_min_cpd=np.array([low for low, _ in BANDS]),
        freq_max_cpd=np.array([high for _, high in BANDS]),
        amplitude_factor=amplitude,
        phase_lead_deg=phase_lead,
    )
