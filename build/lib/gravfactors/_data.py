"""Numerical constants transcribed from the IERS Conventions (2010).

Sources
-------
Chapter 7 (Solid Earth and ocean tides), Tables 7.1 and 7.2 and Eq. 4a, and
Chapter 6 (Geopotential), Tables 6.3 and 6.4 and Eq. 6.10. Those values derive
from Mathews et al. (1995, 2002); they approximate, but are not identical to,
the tables of Dehant, Defraigne and Wahr (1999).
"""

from __future__ import annotations

# Solar days per sidereal day; converts cycles per solar day to cycles per sidereal day.
SIDEREAL_RATIO = 1.002737909

# Resonance frequencies sigma_alpha in cycles per sidereal day (Ch. 6, Eq. 6.10):
# Chandler wobble, retrograde free core nutation, prograde free core nutation.
RESONANCE_FREQS = (
    complex(-0.0026010, -0.0001361),
    complex(1.0023181, 0.000025),
    complex(0.999026, 0.000780),
)

# Diurnal resonance formula L(sigma) = L0 + sum_a L_a / (sigma - sigma_a).
# Each entry is (L0, (L_1, L_2, L_3)).
H0_DIURNAL = (
    complex(0.60671, -0.2420e-2),
    (
        complex(-0.15777e-2, -0.7630e-4),
        complex(0.18053e-3, -0.6292e-5),
        complex(-0.18616e-5, 0.1379e-6),
    ),
)  # IERS Ch. 7, Table 7.1, h(0)
H2_DIURNAL = (
    complex(-0.615e-3, -0.122e-4),
    (
        complex(0.160e-5, 0.116e-6),
        complex(0.201e-6, 0.279e-8),
        complex(-0.329e-7, -0.217e-8),
    ),
)  # IERS Ch. 7, Table 7.1, h(2)
K21_DIURNAL = (
    complex(0.29954, -0.1412e-2),
    (
        complex(-0.77896e-3, -0.3711e-4),
        complex(0.90963e-4, -0.2963e-5),
        complex(-0.11416e-5, 0.5325e-7),
    ),
)  # IERS Ch. 6, Table 6.4, k(0)21

# Latitude-dependence parameter h(2) for the other bands (IERS Ch. 7, Table 7.2).
H2_LONG_PERIOD = -0.0006
H2_SEMIDIURNAL = -0.0006
H_SEMIDIURNAL = complex(0.6078, -0.0022)  # Ch. 7, Table 7.2 / Step 1 nominal values

# Zonal anelastic dispersion (Ch. 7, Eq. 4a):
# h20 = H20_REFERENCE + H20_DISPERSION * {cot(a*pi/2) [1 - X] + i X}, X = (f_m / f)**a
H20_REFERENCE = 0.5998
H20_DISPERSION = -9.96e-4
ANELASTIC_ALPHA = 0.15
REFERENCE_PERIOD_S = 200.0  # period at which the dispersion term vanishes
ELASTIC_PERIOD_S = 1.0  # reference period of the elastic (1 s PREM) model

# Potential Love numbers k (Ch. 6, Table 6.3), indexed by order m.
K_ANELASTIC = {
    0: complex(0.30190, 0.0),
    1: complex(0.29830, -0.00144),
    2: complex(0.30102, -0.00130),
}
K_ELASTIC = {0: 0.29525, 1: 0.29470, 2: 0.29801}

# Degree-3 nominal values (Ch. 7, Step 1; Ch. 6, Table 6.3): real, all orders.
H3 = 0.292
K3 = 0.093
