import runpy
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from gravfactors import compute


SCRIPT = Path(__file__).resolve().parents[1] / "plot_wavegroups.py"


@pytest.mark.parametrize("model", ["inelastic", "elastic"])
def test_plots_each_wave_group_with_correct_values(model):
    plot_wavegroups = runpy.run_path(str(SCRIPT))["plot_wavegroups"]
    result = compute(np.linspace(-90.0, 90.0, 361), model=model)
    for values in (result.amplitude_factor, result.phase_lead_deg):
        fig = plot_wavegroups(result, values, "Test")
        try:
            assert len(fig.axes) == len(result.names)
            for index, ax in enumerate(fig.axes):
                assert result.names[index] in ax.get_title()
                np.testing.assert_array_equal(ax.lines[0].get_xdata(), result.latitude)
                np.testing.assert_array_equal(ax.lines[0].get_ydata(), values[:, index])
        finally:
            plt.close(fig)


def test_script_saves_both_figures(tmp_path):
    prefix = tmp_path / "wavegroups"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--model", "elastic", "--output-prefix", str(prefix)],
        check=True,
    )
    for suffix in ("amplitude", "phase"):
        assert (tmp_path / f"wavegroups_{suffix}.png").read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
