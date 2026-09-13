"""Plot archived VC-NEB per-image metrics without re-running a calculator."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable


def read_metrics(path: Path) -> list[dict[str, float | bool]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows: list[dict[str, float | bool]] = []
        for raw in csv.DictReader(handle):
            rows.append(
                {
                    "reaction_coordinate": float(raw["reaction_coordinate"]),
                    "relative_enthalpy_eV": float(raw["relative_enthalpy_eV"]),
                    "volume_A3": float(raw["volume_A3"]),
                    "a_A": float(raw["a_A"]),
                    "b_A": float(raw["b_A"]),
                    "c_A": float(raw["c_A"]),
                    "max_true_generalized_force_eV_per_A": float(
                        raw["max_true_generalized_force_eV_per_A"]
                    ),
                    "neb_residual_generalized_force_eV_per_A": float(
                        raw["neb_residual_generalized_force_eV_per_A"]
                    )
                    if raw["neb_residual_generalized_force_eV_per_A"]
                    else float("nan"),
                    "is_climbing_image": raw["is_climbing_image"].lower() == "true",
                }
            )
    if len(rows) < 2:
        raise ValueError(f"{path} must contain at least two image rows")
    return rows


def plot_metrics(inputs: Iterable[tuple[str, Path]], output: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - optional plotting extra
        raise RuntimeError("plotting requires matplotlib; install the plot extra") from exc

    series = [(label, read_metrics(path)) for label, path in inputs]
    plt.style.use("seaborn-v0_8-whitegrid")
    figure, axes = plt.subplots(2, 2, figsize=(10.0, 7.0), constrained_layout=True)
    energy_axis, cell_axis, volume_axis, force_axis = axes.flat
    palette = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for series_index, (label, rows) in enumerate(series):
        color = palette[series_index % len(palette)]
        coordinate = [float(row["reaction_coordinate"]) for row in rows]
        energy = [float(row["relative_enthalpy_eV"]) for row in rows]
        volume = [float(row["volume_A3"]) for row in rows]
        cell_lengths = {
            name: [float(row[name]) for row in rows] for name in ("a_A", "b_A", "c_A")
        }
        true_force = [
            float(row["max_true_generalized_force_eV_per_A"]) for row in rows
        ]
        residual = [
            float(row["neb_residual_generalized_force_eV_per_A"]) for row in rows
        ]
        energy_axis.plot(coordinate, energy, marker="o", color=color, label=label)
        cell_axis.plot(coordinate, cell_lengths["a_A"], marker="o", color=color, linestyle="-", label=f"{label}: a")
        cell_axis.plot(coordinate, cell_lengths["b_A"], marker="s", color=color, linestyle="--", label=f"{label}: b")
        cell_axis.plot(coordinate, cell_lengths["c_A"], marker="^", color=color, linestyle=":", label=f"{label}: c")
        volume_axis.plot(coordinate, volume, marker="o", color=color, label=label)
        force_axis.plot(coordinate, true_force, marker="o", color=color, label=f"{label}: true")
        if any(value == value for value in residual):
            force_axis.plot(
                coordinate,
                residual,
                marker="x",
                color=color,
                linestyle="--",
                label=f"{label}: NEB",
            )
        climbing = [
            (x, y)
            for x, y, row in zip(coordinate, energy, rows)
            if bool(row["is_climbing_image"])
        ]
        if climbing:
            energy_axis.scatter(
                [point[0] for point in climbing],
                [point[1] for point in climbing],
                facecolors="none",
                edgecolors="black",
                s=80,
                linewidths=1.2,
                zorder=5,
            )

    energy_axis.set(xlabel="Normalized reaction coordinate", ylabel="Relative enthalpy (eV)")
    cell_axis.set(xlabel="Normalized reaction coordinate", ylabel="Cell length (A)")
    volume_axis.set(xlabel="Normalized reaction coordinate", ylabel="Volume (A^3)")
    force_axis.set(
        xlabel="Normalized reaction coordinate",
        ylabel="Maximum generalized force (eV/A)",
        yscale="log",
    )
    for axis in axes.flat:
        axis.legend(fontsize=8)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=220)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="metrics CSV files; labels default to each file stem",
    )
    parser.add_argument("--label", action="append", help="label for each input CSV")
    parser.add_argument("--output", type=Path, required=True, help="PNG or SVG output")
    args = parser.parse_args()
    if args.label is not None and len(args.label) != len(args.inputs):
        parser.error("--label must be supplied once per input CSV")
    labels = args.label or [path.stem for path in args.inputs]
    try:
        plot_metrics(zip(labels, args.inputs), args.output)
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        parser.error(str(exc))
    print(f"wrote={args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
