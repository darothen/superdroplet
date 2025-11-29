#!/usr/bin/env uv run
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "matplotlib",
#   "numpy",
#   "seaborn",
# ]
# ///
"""Plot the output of a Superdroplet model run."""

import argparse
from pathlib import Path
from typing import Optional

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from cycler import cycler


def colortext_legend(
    text_color_map: dict[str, str],
    ax: Optional[plt.Axes] = None,
    text_labels: Optional[dict[str, str]] = None,
    **kwargs: dict,
) -> Optional["Legend"]:
    """Add a custom-built legend to a plot with colored text.

    All items in the legend have colored text corresponding to colored elements
    in the figure.

    Args:
        text_color_map: A mapping from labels to colors which will be used to
            construct the patch elements to include in the legend.
        ax: The axes instance on which to add the legend. If not passed,
            then the legend will simply be created and returned.
        text_labels: A mapping from labels to longer descriptions to be used
            in their place in the legend.
        **kwargs: Additional arguments to pass to the legend function.

    Returns:
        The legend object, for further customization.

    """

    legend_elems = []
    labels = []
    for text, color in text_color_map.items():
        legend_elems.append(
            (mpatches.Rectangle((0, 0), 1, 1, facecolor="none", edgecolor="none"), text)
        )
        labels.append(text)

    # Set up a legend with colored-text elements
    elems, texts = zip(*legend_elems)
    if ax is not None:
        ax.legend(elems, texts, **kwargs)
        leg = ax.get_legend()
    else:
        leg = plt.legend(elems, texts, **kwargs)

    # Change the label color
    for label in leg.get_texts():
        old_label = label.get_text()
        if text_labels is None:
            new_label = old_label
        else:
            try:
                new_label = text_labels[old_label]
            except KeyError:
                new_label = old_label

        plt.setp(label, color=text_color_map[old_label])
        label.set_text(new_label)

    return leg


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    parser.add_argument("output_file", type=Path, help="The output file to plot.")
    parser.add_argument("-m", "--model-name", help="The name of the model to plot.")
    return parser.parse_args()


def plot_output(output_file: Path, model_name: Optional[str] = None) -> None:
    data = np.loadtxt(output_file, delimiter=",", dtype=float)
    times_sec = data[1:, 0]
    bin_radii_um = data[0, 1:]
    bin_counts = data[1:, 1:]

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)

    labels = []
    colors = []
    styles = cycler(color=[f"C{i}" for i in range(len(times_sec))])
    for i, (time_sec, bin_counts, style) in enumerate(
        zip(times_sec, bin_counts, styles)
    ):
        if i == 0:
            label = "Initial"
            lw = 1.5
            color = "k"
        else:
            label = f"{time_sec:g} sec"
            lw = 1.0
            color = style["color"]
        ax.plot(bin_radii_um, bin_counts, label=label, lw=lw, color=color)
        labels.append(label)
        colors.append(color)

    _ = colortext_legend(
        text_color_map={label: color for label, color in zip(labels, colors)},
        ax=ax,
        frameon=False,
        bbox_to_anchor=(0.8, 0.5),
        loc="center left",
    )
    ax.set_xlabel("Bin Radius (µm)")
    ax.set_ylabel("Droplet Count")
    ax.set_title("Droplet Size Distribution", loc="left")
    if model_name:
        ax.set_title(model_name, loc="right")
    ax.set_xscale("log")

    ax.set_ylim(0)
    ax.set_xlim(0)
    sns.despine(fig=fig, offset=5)

    plt.show()


if __name__ == "__main__":
    args = parse_args()

    plot_output(args.output_file, args.model_name)
