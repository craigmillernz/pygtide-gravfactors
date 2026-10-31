"""The 16 tidal constituents and their contiguous wave-group frequency bands."""

from __future__ import annotations

from dataclasses import dataclass

# Fundamental Doodson frequencies in cycles per solar day.
_S = 0.0366009877  # mean longitude of the Moon
_H = 0.0027379093  # mean longitude of the Sun
_P = 0.0003094519  # lunar perigee
_N = 0.0001470950  # negative of the lunar node, N'
_TAU = 1.0027379093 - _S  # lunar time

_SPECIES_EDGES = {0: (0.0, 0.5), 1: (0.5, 1.5), 2: (1.5, 2.5), 3: (2.5, 3.5)}
SPECIES_NAMES = {0: "long-period", 1: "diurnal", 2: "semidiurnal", 3: "terdiurnal"}


@dataclass(frozen=True)
class Constituent:
    """A tidal constituent with its Doodson multipliers (tau, s, h, p, N')."""

    name: str
    order: int  # m: 0 long-period, 1 diurnal, 2 semidiurnal, 3 terdiurnal
    degree: int
    doodson: tuple[int, int, int, int, int]

    @property
    def freq_cpd(self) -> float:
        tau, s, h, p, n = self.doodson
        return tau * _TAU + s * _S + h * _H + p * _P + n * _N

    @property
    def species(self) -> str:
        return SPECIES_NAMES[self.order]


_DEFINITIONS = (
    Constituent("Sa", 0, 2, (0, 0, 1, 0, 0)),
    Constituent("Ssa", 0, 2, (0, 0, 2, 0, 0)),
    Constituent("Mm", 0, 2, (0, 1, 0, -1, 0)),
    Constituent("Mf", 0, 2, (0, 2, 0, 0, 0)),
    Constituent("2Q1", 1, 2, (1, -3, 0, 2, 0)),
    Constituent("Q1", 1, 2, (1, -2, 0, 1, 0)),
    Constituent("O1", 1, 2, (1, -1, 0, 0, 0)),
    Constituent("P1", 1, 2, (1, 1, -2, 0, 0)),
    Constituent("K1", 1, 2, (1, 1, 0, 0, 0)),
    Constituent("J1", 1, 2, (1, 2, 0, -1, 0)),
    Constituent("OO1", 1, 2, (1, 3, 0, 0, 0)),
    Constituent("N2", 2, 2, (2, -1, 0, 1, 0)),
    Constituent("M2", 2, 2, (2, 0, 0, 0, 0)),
    Constituent("S2", 2, 2, (2, 2, -2, 0, 0)),
    Constituent("K2", 2, 2, (2, 2, 0, 0, 0)),
    Constituent("M3", 3, 3, (3, 0, 0, 0, 0)),
)


def _band_edges() -> list[tuple[float, float]]:
    """Split each species range at the midpoints between neighbouring constituents."""
    edges: list[tuple[float, float]] = []
    for order, (low, high) in _SPECIES_EDGES.items():
        freqs = [c.freq_cpd for c in _DEFINITIONS if c.order == order]
        bounds = [low] + [(a + b) / 2 for a, b in zip(freqs, freqs[1:])] + [high]
        edges.extend(zip(bounds[:-1], bounds[1:]))
    return edges


CONSTITUENTS: tuple[Constituent, ...] = _DEFINITIONS
BANDS: tuple[tuple[float, float], ...] = tuple(_band_edges())
