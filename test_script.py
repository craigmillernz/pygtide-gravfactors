from gravfactors import compute
import pygtide
from pygtide import predict_table, predict_series, plot_series
import matplotlib.pyplot as plt
result = compute(latitude=-36.85, longitude=174.76, height=30.0, model="inelastic")

for name, species, lo, hi, delta, lead in zip(
    result.names, result.species, result.freq_min_cpd, result.freq_max_cpd,
    result.amplitude_factor, result.phase_lead_deg,
):
    print(f"{name:4s} {species:12s} {lo:.5f}-{hi:.5f} cpd  delta={delta:.4f}  phase={lead:+.3f} deg")

# Many points: arrays of shape (n_points, 16)
many = compute([-36.85, 48.0, 0.0], [174.76, 11.0, 30.0])

# Straight into PyGTide: columns are freq from, freq to, amplitude factor, phase lead [deg]

tide = pygtide.pygtide()
tide.set_wavegroup(result.to_pygtide_wavegroups())
latitude=-36.85
longitude=174.76
height=30.0
start='2026-10-16'
duration= 24*30
samprate=600

out = predict_table(latitude, longitude, height,start, duration, samprate )
fig, ax = plt.subplots()
plot_series(latitude, longitude, height,start, duration, samprate, ax=ax)
# Tidy rows, e.g. pandas.DataFrame(many.to_records())
tide.reset_wavegroup()
plot_series(latitude, longitude, height,start, duration, samprate, ax=ax)