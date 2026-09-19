#!/usr/bin/env python3
"""Plot the audited GaN B3-to-B1 VCNEB enthalpy path.

The line connects the computed images without interpolation or smoothing.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--output-stem", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    analysis = json.loads(args.analysis.read_text(encoding="utf-8"))
    energies = [float(value) for value in analysis["relative_enthalpies_eV_per_GaN"]]
    structures = analysis["structures"]
    images = list(range(len(energies)))
    saddle = int(analysis["highest_image_index"])
    literature = float(analysis["literature_barrier_eV_per_GaN"])

    args.output_stem.parent.mkdir(parents=True, exist_ok=True)
    source_path = args.output_stem.with_name(args.output_stem.name + "_source_data.csv")
    with source_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "image_index",
                "relative_enthalpy_eV_per_GaN",
                "space_group_number_symprec_1e-4_A",
                "space_group_symbol_symprec_1e-4_A",
                "volume_A3",
            ],
        )
        writer.writeheader()
        for image, energy, structure in zip(images, energies, structures, strict=True):
            symmetry = structure["symmetry"]["0.0001"]
            writer.writerow(
                {
                    "image_index": image,
                    "relative_enthalpy_eV_per_GaN": f"{energy:.12f}",
                    "space_group_number_symprec_1e-4_A": symmetry["number"],
                    "space_group_symbol_symprec_1e-4_A": symmetry["international"],
                    "volume_A3": f"{float(structure['volume_A3']):.12f}",
                }
            )

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 7,
            "axes.labelsize": 8,
            "axes.titlesize": 8,
            "axes.linewidth": 0.8,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )

    fig, ax = plt.subplots(figsize=(3.50, 2.55), constrained_layout=True)
    path_color = "#244A73"
    accent_color = "#B64949"
    endpoint_color = "#2F7D6D"

    ax.plot(
        images,
        energies,
        color=path_color,
        linewidth=1.35,
        marker="o",
        markersize=3.1,
        markerfacecolor="white",
        markeredgecolor=path_color,
        markeredgewidth=0.75,
        zorder=2,
    )
    ax.scatter(
        [images[0], images[-1]],
        [energies[0], energies[-1]],
        s=28,
        color=endpoint_color,
        edgecolor="white",
        linewidth=0.7,
        zorder=4,
    )
    ax.scatter(
        [saddle],
        [energies[saddle]],
        s=34,
        color=accent_color,
        edgecolor="white",
        linewidth=0.7,
        zorder=5,
    )

    ax.axhline(literature, color=accent_color, linewidth=0.9, linestyle=(0, (4, 3)), alpha=0.85)
    ax.text(
        27.7,
        literature + 0.025,
        f"Qian et al. ≈ {literature:.2f} eV/GaN",
        color=accent_color,
        ha="right",
        va="bottom",
        fontsize=6.6,
    )
    ax.annotate(
        f"image {saddle}\n{energies[saddle]:.3f} eV/GaN",
        xy=(saddle, energies[saddle]),
        xytext=(16.0, 0.88),
        arrowprops={"arrowstyle": "-", "color": accent_color, "linewidth": 0.8},
        color=accent_color,
        ha="left",
        va="center",
        fontsize=6.8,
    )
    ax.text(0.0, energies[0] + 0.045, "B3", color=endpoint_color, ha="left", va="bottom")
    ax.text(28.0, energies[-1] + 0.045, "B1", color=endpoint_color, ha="right", va="bottom")
    ax.text(
        0.02,
        0.96,
        r"45.0 GPa  |  29 images  |  final fmax = 0.099 eV Å$^{-1}$",
        transform=ax.transAxes,
        ha="left",
        va="top",
        color="#4B5563",
        fontsize=6.4,
    )

    ax.set_xlim(-0.6, 28.6)
    ax.set_ylim(-0.055, 1.055)
    ax.set_xticks([0, 5, 10, 15, 20, 25, 28])
    ax.set_xlabel("Image index")
    ax.set_ylabel("Relative enthalpy (eV/GaN)")
    ax.set_title("GaN B3→B1 converged VCNEB path", loc="left", pad=8, fontweight="bold")
    ax.grid(axis="y", color="#D7DCE2", linewidth=0.55, alpha=0.8)
    ax.set_axisbelow(True)

    for suffix, kwargs in {
        ".png": {"dpi": 300},
        ".tiff": {"dpi": 600},
        ".svg": {},
        ".pdf": {},
    }.items():
        fig.savefig(args.output_stem.with_suffix(suffix), bbox_inches="tight", **kwargs)
    plt.close(fig)


if __name__ == "__main__":
    main()
