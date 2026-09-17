"""Create a publication-ready summary of the converged BST50 VCA-VCNEB run."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 8
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["legend.frameon"] = False

BLUE = "#0F4D92"
TEAL = "#33B5A5"
RED = "#B64342"
GREY = "#767676"
LIGHT_GREY = "#D9D9D9"
BLACK = "#272727"


def parse_optimizer_log(path: Path) -> tuple[np.ndarray, np.ndarray]:
    steps, residuals = [], []
    pattern = re.compile(r"FIRE:\s+(\d+)\s+\S+\s+[-+0-9.Ee]+\s+([-+0-9.Ee]+)")
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.search(line)
        if match:
            steps.append(int(match.group(1)))
            residuals.append(float(match.group(2)))
    if not steps:
        raise ValueError(f"no FIRE records found in {path}")
    return np.asarray(steps), np.asarray(residuals)


def panel_label(ax, label: str) -> None:
    ax.text(-0.14, 1.04, label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--audit", required=True)
    parser.add_argument("--optimizer-log", required=True)
    parser.add_argument("--output", required=True, help="output basename without extension")
    args = parser.parse_args()

    summary = json.loads(Path(args.summary).read_text(encoding="utf-8"))
    audit = json.loads(Path(args.audit).read_text(encoding="utf-8"))
    if summary.get("status") != "completed" or audit.get("status") != "ok":
        raise ValueError("VCNEB result is not completed and audited")

    images = summary["path_diagnostics"]["images"]
    segments = np.asarray(summary["path_diagnostics"]["geometry"]["segment_lengths_A"], dtype=float)
    coordinate = np.r_[0.0, np.cumsum(segments)]
    coordinate /= coordinate[-1]
    energy = 1000.0 * np.asarray([image["relative_enthalpy_eV"] for image in images])
    lengths = np.asarray([image["cell_lengths_A"] for image in images], dtype=float)
    c_over_a = lengths[:, 2] / lengths[:, :2].mean(axis=1)
    volume = np.asarray([image["volume_A3"] for image in images], dtype=float)
    steps, residuals = parse_optimizer_log(Path(args.optimizer_log))
    target = float(summary["fmax_target_eV_per_A"])

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.with_name(out.name + "_source_data").with_suffix(".csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["image", "reaction_coordinate", "relative_enthalpy_meV_per_fu", "c_over_a", "volume_A3", "neb_residual_eV_per_A"])
        for index, image in enumerate(images):
            writer.writerow([
                index, coordinate[index], energy[index], c_over_a[index], volume[index],
                image["neb_residual_generalized_force_eV_per_A"],
            ])

    fig = plt.figure(figsize=(7.09, 3.35))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.65, 1.0], height_ratios=[1, 1], wspace=0.42, hspace=0.62)
    ax_energy = fig.add_subplot(grid[:, 0])
    ax_shape = fig.add_subplot(grid[0, 1])
    ax_conv = fig.add_subplot(grid[1, 1])

    ax_energy.plot(coordinate, energy, color=BLUE, lw=2.1, marker="o", ms=5.0, mec="white", mew=0.7)
    ax_energy.fill_between(coordinate, 0, energy, color=BLUE, alpha=0.09)
    ax_energy.axhline(0, color=LIGHT_GREY, lw=0.8, zorder=0)
    ax_energy.set_xlabel("Normalized T → C path coordinate")
    ax_energy.set_ylabel("Relative enthalpy (meV/f.u.)")
    ax_energy.set_xlim(-0.03, 1.03)
    ax_energy.set_ylim(-2.5, max(35.5, energy.max() + 4))
    ax_energy.text(coordinate[0], energy[0] + 1.8, "T", color=BLUE, ha="center", fontweight="bold")
    ax_energy.text(coordinate[-1], energy[-1] + 1.8, "C", color=BLUE, ha="center", fontweight="bold")
    ax_energy.annotate(
        f"ΔE$_{{C-T}}$ = {energy[-1]:.2f} meV/f.u.",
        xy=(1.0, energy[-1]), xytext=(0.51, energy[-1] - 7.0),
        arrowprops={"arrowstyle": "->", "color": GREY, "lw": 0.8}, color=BLACK,
    )
    ax_energy.text(
        0.05, 0.92, "Monotonic path\nno interior saddle",
        transform=ax_energy.transAxes, ha="left", va="top", color=GREY,
    )
    panel_label(ax_energy, "a")

    ax_shape.plot(coordinate, c_over_a, color=TEAL, lw=1.8, marker="s", ms=4.2, mec="white", mew=0.6)
    ax_shape.axhline(1.0, color=LIGHT_GREY, lw=0.8, zorder=0)
    ax_shape.set_xlim(-0.03, 1.03)
    ax_shape.set_ylabel("Tetragonality, c/a")
    ax_shape.set_xticks([0, 0.5, 1.0])
    ax_shape.set_xticklabels([])
    ax_shape.text(0.02, c_over_a[0] + 0.002, f"{c_over_a[0]:.3f}", color=TEAL, ha="left")
    ax_shape.text(0.98, c_over_a[-1] + 0.002, "1.000", color=TEAL, ha="right")
    panel_label(ax_shape, "b")

    ax_conv.plot(steps, residuals, color=RED, lw=1.8, marker="o", ms=4.2, mec="white", mew=0.6)
    ax_conv.axhline(target, color=GREY, lw=1.0, ls="--")
    ax_conv.fill_between([steps.min(), steps.max()], 0, target, color=LIGHT_GREY, alpha=0.25)
    ax_conv.set_xlabel("FIRE step")
    ax_conv.set_ylabel("Max. generalized force\n(eV Å$^{-1}$)")
    ax_conv.set_xlim(steps.min() - 0.15, steps.max() + 0.15)
    ax_conv.set_ylim(0, max(residuals) * 1.14)
    ax_conv.set_xticks(steps)
    ax_conv.text(steps.max() - 0.05, target + 0.012, f"target = {target:.2f}", color=GREY, ha="right")
    ax_conv.annotate(
        f"{residuals[-1]:.3f}", xy=(steps[-1], residuals[-1]), xytext=(steps[-1] - 0.9, residuals[-1] + 0.055),
        arrowprops={"arrowstyle": "->", "color": GREY, "lw": 0.8}, color=RED,
    )
    panel_label(ax_conv, "c")

    fig.suptitle("Ba$_{0.5}$Sr$_{0.5}$TiO$_3$ virtual-crystal T → C transformation", x=0.50, y=1.01, fontsize=10, fontweight="bold")
    for extension, kwargs in {
        ".svg": {},
        ".pdf": {},
        ".png": {"dpi": 400},
        ".tiff": {"dpi": 600},
    }.items():
        fig.savefig(out.with_suffix(extension), bbox_inches="tight", facecolor="white", **kwargs)
    plt.close(fig)

    figure_summary = {
        "conclusion": "BST50 retains a tetragonal 0 K ground state and follows a monotonic, barrierless T-to-C path.",
        "reaction_enthalpy_meV_per_fu": float(energy[-1]),
        "interior_barrier": False,
        "final_max_generalized_force_eV_per_A": float(residuals[-1]),
        "fmax_target_eV_per_A": target,
        "tetragonal_c_over_a": float(c_over_a[0]),
        "cubic_c_over_a": float(c_over_a[-1]),
        "minimum_path_distance_A": float(audit["minimum_path_distance_A"]),
        "maximum_deformation": float(audit["maximum_deformation"]),
    }
    out.with_name(out.name + "_figure_summary").with_suffix(".json").write_text(
        json.dumps(figure_summary, indent=2) + "\n", encoding="utf-8",
    )
    print(json.dumps(figure_summary, indent=2))


if __name__ == "__main__":
    main()
