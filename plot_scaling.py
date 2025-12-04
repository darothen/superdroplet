#!/usr/bin/env -S uv run --script

# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "matplotlib",
#   "pandas",
#   "seaborn",
# ]
# ///
"""Plot the scaling performance of Superdroplet model implementations."""

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import seaborn as sns


def create_timing_data() -> pd.DataFrame:
    """Create a DataFrame from the timing data in the README table."""
    # Data extracted from README.md table
    data = {
        "Language": [],
        "N_droplets": [],
        "Time_sec": [],
    }

    # Number of superdroplets (powers of 2)
    n_droplets = [2**15, 2**17, 2**19, 2**21]

    # Timing data for each language (in seconds)
    timing_data = {
        "Python w/ Numba": [3.9, 15.8, 175.2, 874.5],
        "Julia": [0.9, 4.4, 43.8, 373.5],
        "C++": [0.7, 2.9, 19.6, 155.3],
        "Fortran": [1.9, 11.5, 68.2, 347.9],
        "Rust": [1.0, 4.2, 38.2, 245.9],
    }

    # Build the DataFrame
    for language, times in timing_data.items():
        for n, time in zip(n_droplets, times):
            data["Language"].append(language)
            data["N_droplets"].append(n)
            data["Time_sec"].append(time)

    return pd.DataFrame(data)


def plot_scaling_performance(output_file: str = "scaling.png") -> None:
    """Create and save a scaling performance visualization."""
    # Create the data
    df = create_timing_data()

    # Set up the plot style
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.3)

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot each language with different colors and markers
    markers = ["o", "s", "^", "D", "v"]
    for (language, group), marker in zip(df.groupby("Language"), markers):
        ax.plot(
            group["N_droplets"],
            group["Time_sec"],
            marker=marker,
            markersize=8,
            linewidth=2,
            label=language,
            alpha=0.8,
        )

    # Set log scale for both axes
    ax.set_xscale("log")
    ax.set_yscale("log")

    # Set axis limits
    ax.set_xlim(20_000, 5_000_000)
    ax.set_ylim(0.1, 1_500)

    # Format x-axis ticks
    ax.xaxis.set_major_locator(ticker.LogLocator(base=10, numticks=10))
    ax.xaxis.set_major_formatter(
        ticker.FuncFormatter(
            lambda x, p: f"{int(x):,}" if x >= 1 else f"1/{int(1 / x)}"
        )
    )
    ax.xaxis.set_minor_locator(ticker.LogLocator(base=10, subs="auto", numticks=10))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())

    # Format y-axis ticks
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=10))
    ax.yaxis.set_major_formatter(
        ticker.FuncFormatter(
            lambda x, p: f"{int(x):,}" if x >= 1 else f"1/{int(1 / x)}"
        )
    )
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10, subs="auto", numticks=10))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())

    # Configure tick marks to be visible and point outward
    ax.tick_params(axis="both", which="major", direction="out", length=6, width=1.5)
    ax.tick_params(axis="both", which="minor", direction="out", length=3, width=1)

    # Add grid for major ticks only
    ax.grid(True, which="major", alpha=0.3, linestyle="--", linewidth=0.8)

    # Labels and title
    ax.set_xlabel("Number of Superdroplets", fontsize=12, fontweight="bold")
    ax.set_ylabel("Time (seconds)", fontsize=12, fontweight="bold")
    ax.set_title(
        "Superdroplet Model Scaling Performance", fontsize=14, fontweight="bold", pad=20
    )

    # Legend
    ax.legend(
        frameon=False,
        fancybox=False,
        shadow=False,
        loc="upper left",
        fontsize=11,
    )

    # Remove top and right spines while keeping ticks on left and bottom
    sns.despine(fig=fig, offset=5, trim=False)

    plt.show()

    # Save figure
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"Scaling performance plot saved to: {output_file}")


if __name__ == "__main__":
    plot_scaling_performance()
