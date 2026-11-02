from gravfactors import compute
import pygtide
from pygtide import predict_table, predict_series, plot_series
import matplotlib.pyplot as plt
import numpy as np

result = compute(latitude=-36.85, longitude=174.76, height=30.0, model="inelastic")

for name, species, lo, hi, delta, lead in zip(
    result.names, result.species, result.freq_min_cpd, result.freq_max_cpd,
    result.amplitude_factor, result.phase_lead_deg,
):
    print(f"{name:4s} {species:12s} {lo:.5f}-{hi:.5f} cpd  delta={delta:.4f}  phase={lead:+.3f} deg")

# Many points: arrays of shape (n_points, 16)
many = compute([-36.85, 48.0, 0.0], [174.76, 11.0, 30.0])

# Straight into PyGTide: columns are freq from, freq to, amplitude factor, phase lead [deg]

latitude = -38.6875
longitude = 176.0694
height = 350.0
start = "2026-10-16"
duration_hours = 24
interval_s = 60

# predict_table/plot_series create their own pygtide instance, so they ignore wave groups
# set here. Call predict() on this instance instead.
tide = pygtide.pygtide()


def run():
    tide.predict(latitude, longitude, height, start, duration_hours, interval_s, tidalcompo=0)
    return tide.results()


tide.reset_wavegroup()
out_default = run()

tide.set_wavegroup(result.to_pygtide_wavegroups())
out_factors = run()

# Single wave group covering 0-10 cpd, factor 1.16, phase 0
tide.set_wavegroup(np.array([[0, 10, 1.16, 0]]))
out_single = run()

col = "Tide [nm/s**2]"
print("max |factors - default| [nm/s**2]:", np.abs(out_factors[col] - out_default[col]).max())
print("max |single - default| [nm/s**2]:", np.abs(out_single[col] - out_default[col]).max())

# %%
fig, ax = plt.subplots()
ax.plot(out_default["UTC"], out_default[col], label="default wave group [0, 10, 1, 0]")
ax.plot(out_factors["UTC"], out_factors[col], label="gravfactors wave groups")
ax.plot(out_single["UTC"], out_single[col], label="single group [0, 10, 1.16, 0]")
ax.set_ylabel(col)
ax.legend()

fig2, ax2 = plt.subplots()
ax2.plot(out_single["UTC"], out_single[col] - out_factors[col])
ax2.set_ylabel(f"single group - gravfactors groups, {col}")
plt.show()
# %%
