# pygtide-gravfactors

Gravimetric amplitude factors (δ) and phases for 16 tidal constituents, for an
elastic and an inelastic (anelastic-mantle) Earth, packaged as ready-to-use
[PyGTide](https://github.com/hydrogeoscience/pygtide) wave groups.

Inputs are latitude, longitude and ellipsoidal height, as a single point or arrays.

## Install

```bash
pip install pygtide-gravfactors
```

## Usage

```python
from gravfactors import compute

result = compute(latitude=-36.85, longitude=174.76, height=30.0, model="inelastic")

for name, species, lo, hi, delta, lead in zip(
    result.names, result.species, result.freq_min_cpd, result.freq_max_cpd,
    result.amplitude_factor, result.phase_lead_deg,
):
    print(f"{name:4s} {species:12s} {lo:.5f}-{hi:.5f} cpd  delta={delta:.4f}  phase={lead:+.3f} deg")

# Many points: arrays of shape (n_points, 16)
many = compute([-36.85, 48.0, 0.0], [174.76, 11.0, 30.0])

# Straight into PyGTide: columns are freq from, freq to, amplitude factor, phase lead [deg]
import pygtide
tide = pygtide.pygtide()
tide.set_wavegroup(result.to_pygtide_wavegroups())

# Tidy rows, e.g. pandas.DataFrame(many.to_records())
rows = many.to_records()
```

`compute(latitude, longitude=0.0, height=0.0, model="inelastic")` returns a
`GravimetricFactors` object:

| Field | Meaning |
|---|---|
| `names`, `species` | constituent and band (long-period, diurnal, semidiurnal, terdiurnal) |
| `freq_cpd`, `freq_min_cpd`, `freq_max_cpd` | constituent frequency and its wave-group band, cycles per solar day |
| `amplitude_factor` | δ, shape `(16,)` for one point or `(n_points, 16)` |
| `phase_lead_deg` | phase lead in degrees (PyGTide convention); negative means the tide lags |
| `phase_lag_deg` | `-phase_lead_deg` |

Constituents: Sa, Ssa, Mm, Mf (long-period); 2Q1, Q1, O1, P1, K1, J1, OO1 (diurnal);
N2, M2, S2, K2 (semidiurnal); M3 (terdiurnal). Bands are contiguous and
non-overlapping: each species range is split at the midpoints between neighbouring
constituents.

## Method and limitations

Please read before using these numbers in published work.

* **Source.** Love numbers come from the IERS Conventions (2010), Chapter 7
  (Tables 7.1, 7.2, Eq. 4a) and Chapter 6 (Tables 6.3, 6.4, Eq. 6.10). These derive from
  Mathews et al. (1995, 2002) and approximate, but are **not identical to**, the tables of
  Dehant, Defraigne and Wahr (1999), *Tides for a convective Earth*, J. Geophys. Res. 104(B1),
  1035-1058. The package does not reproduce that paper's two models exactly.
* **Gravimetric factor.** δ = 1 + (2/n)·h − ((n+1)/n)·k with n = 2, or n = 3 for M3. h and k are
  complex in the inelastic model; the amplitude is |δ| and the phase lead is arg(δ).
* **Diurnal band.** Complex h and k include the nearly diurnal free wobble resonance
  (K1 is lowered relative to O1). Table 7.2 of IERS also folds in ocean-loading terms that
  the pure body-tide resonance formula used here does not, so h differs from Table 7.2 by up
  to about 3e-3 at K1.
* **Latitude dependence** comes only through h(2) and the geocentric latitude. IERS gives no
  latitude-dependent k, so k is latitude-independent. The effect on δ is small (below 0.01).
  Longitude and height are accepted but have no effect.
* **Elastic model is derived, not tabulated.** IERS lists no elastic h. It is built by removing
  the mantle dispersion term of Eq. 4a, evaluated between the tidal period and a 1 s reference
  period, from the real part of h, setting imaginary parts to zero, and using the elastic k of
  Table 6.3. For the diurnal and semidiurnal bands this borrows the zonal dispersion term, so
  treat it as an approximation. Phases are exactly zero.
* **Terdiurnal (M3)** uses nominal real values (h3 = 0.292, k3 = 0.093), so its phase is zero.
* Within the semidiurnal band IERS gives no frequency dependence, so N2, M2, S2 and K2 share
  the same inelastic δ and phase.

## PyGTide behaviour

PyGTide normalises each wave group by the Dehant (1987) factor of its main wave, so the
factors here are applied as absolute amplitude factors to each group's main wave, while
satellite waves keep PyGTide's internal relative factors. A positive phase is a phase lead,
verified in `tests/test_pygtide_integration.py`.

## Development

```bash
pip install -e ".[test]" pygtide
pytest
```

## License

MIT
