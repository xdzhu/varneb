#!/usr/bin/env python3
"""Plot the unsmoothed enthalpy profile of one evaluated VCNEB snapshot."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from ase.units import GPa


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--output-stem", type=Path, required=True)
    parser.add_argument("--pressure-gpa", type=float, required=True)
    parser.add_argument("--formula-units", type=int, required=True)
    parser.add_argument("--formula-label", default="f.u.")
    parser.add_argument("--step", type=int, required=True)
    parser.add_argument("--fmax", type=float, required=True)
    parser.add_argument("--fmax-target", type=float, default=0.10)
    parser.add_argument("--reference-barrier", type=float)
    parser.add_argument("--reference-label", default="Published barrier")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.formula_units < 1:
        raise ValueError("formula-units must be positive")
    images = read(args.snapshot, index=":")
    if len(images) < 3 or any(image.calc is None for image in images):
        raise ValueError("snapshot must contain evaluated endpoints and interior images")
    if len({len(image) for image in images}) != 1:
        raise ValueError("all images must have the same atom count")

    energy = np.asarray([image.get_potential_energy() for image in images], dtype=float)
    volume = np.asarray([image.get_volume() for image in images], dtype=float)
    if not np.isfinite(energy).all() or not np.isfinite(volume).all():
        raise ValueError("snapshot contains non-finite energy or volume")
    enthalpy = energy + args.pressure_gpa * GPa * volume
    relative = (enthalpy - enthalpy[0]) / args.formula_units
    indices = np.arange(len(images))
    peak = int(np.argmax(relative[1:-1])) + 1
    barrier = float(relative[peak])
    status = "converged" if args.fmax < args.fmax_target else "unconverged"

    args.output_stem.parent.mkdir(parents=True, exist_ok=True)
    source = args.output_stem.with_name(args.output_stem.name + "_source_data.csv")
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["image_index", "energy_eV_per_cell", "volume_A3", "enthalpy_eV_per_cell", "relative_enthalpy_eV_per_formula_unit"])
        for row in zip(indices, energy, volume, enthalpy, relative, strict=True):
            writer.writerow([int(row[0]), *(f"{float(value):.12f}" for value in row[1:])])

    provenance = {
        "snapshot": Path(os.path.relpath(args.snapshot.resolve(),
                                         start=args.output_stem.parent.resolve())).as_posix(),
        "snapshot_sha256": hashlib.sha256(args.snapshot.read_bytes()).hexdigest(),
        "step": args.step,
        "status": status,
        "fmax_eV_per_A": args.fmax,
        "fmax_target_eV_per_A": args.fmax_target,
        "n_images_total": len(images),
        "n_images_interior": len(images) - 2,
        "n_atoms_per_image": len(images[0]),
        "formula_units_per_cell": args.formula_units,
        "pressure_gpa": args.pressure_gpa,
        "peak_image": peak,
        "barrier_eV_per_formula_unit": barrier,
        "reaction_enthalpy_eV_per_formula_unit": float(relative[-1]),
        "reference_barrier_eV_per_formula_unit": args.reference_barrier,
        "reference_label": args.reference_label if args.reference_barrier is not None else None,
        "source_data": source.name,
    }
    args.output_stem.with_name(args.output_stem.name + "_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )

    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 10,
        "axes.labelsize": 10.5,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "axes.linewidth": 0.8,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    })
    fig, ax = plt.subplots(figsize=(4.7, 3.2), constrained_layout=True)
    ax.plot(indices, relative, color="#27547A", lw=1.8, marker="o", ms=3.6,
            mfc="white", mec="#27547A", mew=0.9, zorder=3)
    ax.scatter([peak], [barrier], s=42, color="#B24742", edgecolor="white", lw=0.7, zorder=5)
    ax.axhline(0, color="#AAB4BF", lw=0.75, zorder=1)
    if args.reference_barrier is not None:
        ax.axhline(args.reference_barrier, color="#866C4D", lw=1.1,
                   linestyle=(0, (4, 3)), zorder=2)
        ax.text(len(images) - 1, args.reference_barrier + 0.008,
                f"{args.reference_label}: {args.reference_barrier:.2f}",
                ha="right", va="bottom", color="#725B3E", fontsize=9)
    ax.annotate(f"{barrier:.3f} eV/{args.formula_label}", xy=(peak, barrier),
                xytext=(peak + 1.2, barrier + 0.010),
                ha="left", va="bottom", color="#A43E3A", fontsize=9.5,
                arrowprops={"arrowstyle": "-", "color": "#A43E3A", "lw": 0.85})
    ax.text(0.02, 0.97, f"Step {args.step} ({status})",
            transform=ax.transAxes, ha="left", va="top", fontsize=9, color="#4A5562")
    ax.text(0, relative[0] + 0.015, "B4", ha="left", va="bottom", color="#27547A")
    ax.text(len(images) - 1, relative[-1] + 0.015, "B1", ha="right", va="bottom", color="#27547A")
    ax.set(xlim=(-0.6, len(images) - 0.4),
           ylim=(min(-0.085, float(relative.min()) - 0.025),
                 max(0.44, float(relative.max()) + 0.10)),
           xlabel="Image index", ylabel=f"Relative enthalpy (eV/{args.formula_label})")
    ax.tick_params(direction="in", top=True, right=True, length=4, width=0.8)
    ax.grid(axis="y", color="#E2E6EA", lw=0.6)
    ax.set_axisbelow(True)
    for suffix, kwargs in {
        ".png": {"dpi": 300},
        ".tiff": {"dpi": 600},
        ".svg": {},
        ".pdf": {},
    }.items():
        fig.savefig(args.output_stem.with_suffix(suffix), bbox_inches="tight", **kwargs)
    plt.close(fig)
    print(json.dumps({"barrier_eV_per_formula_unit": barrier, "peak_image": peak,
                      "reaction_eV_per_formula_unit": float(relative[-1]),
                      "source_data": str(source)}, indent=2))


if __name__ == "__main__":
    main()
