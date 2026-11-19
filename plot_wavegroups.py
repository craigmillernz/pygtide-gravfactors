"""Plot gravimetric wave-group factors and phases against latitude."""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

from gravfactors import GravimetricFactors, compute


def plot_wavegroups(
    result: GravimetricFactors, values: NDArray[np.float64], ylabel: str
) -> plt.Figure:
    fig, axes = plt.subplots(4, 4, figsize=(14, 10), sharex=True, layout="constrained")
    for index, ax in enumerate(axes.flat):
        ax.plot(result.latitude, values[:, index])
        ax.set_title(f"{result.names[index]} ({result.species[index]})")
        ax.grid(True, alpha=0.3)
    fig.supxlabel("Geodetic latitude (degrees)")
    fig.supylabel(ylabel)
    fig.suptitle(f"{ylabel} by latitude ({result.model} model)")
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", choices=("inelastic", "elastic"), default="inelastic")
    parser.add_argument(
        "--output-prefix",
        help="Save PNG figures with this path prefix instead of showing them",
    )
    args = parser.parse_args()

    result = compute(np.linspace(-90.0, 90.0, 361), model=args.model)
    amplitude = plot_wavegroups(result, result.amplitude_factor, "Amplitude factor")
    phase = plot_wavegroups(result, result.phase_lead_deg, "Phase lead (degrees)")

    if args.output_prefix:
        amplitude.savefig(f"{args.output_prefix}_amplitude.png", dpi=150)
        phase.savefig(f"{args.output_prefix}_phase.png", dpi=150)
        plt.close("all")
    else:
        plt.show()


if __name__ == "__main__":
    main()
