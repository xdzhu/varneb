#!/usr/bin/env python3
"""Extract and compare latest and lowest-force CdSe VCNEB path snapshots."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from ase.io import Trajectory, write


LOG_PATTERN = re.compile(
    r"^CheckedFIRE:\s+(?P<step>\d+)\s+\S+\s+"
    r"(?P<energy>[-+0-9.eE]+)\s+(?P<fmax>[-+0-9.eE]+)\s*$"
)


@dataclass(frozen=True)
class LogRecord:
    step: int
    energy_eV: float
    fmax_eV_per_A: float


@dataclass(frozen=True)
class RouteInput:
    key: str
    title: str
    trajectory: Path
    logfile: Path


def parse_log(path: Path) -> list[LogRecord]:
    records: list[LogRecord] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        match = LOG_PATTERN.match(line)
        if match:
            records.append(
                LogRecord(
                    step=int(match.group("step")),
                    energy_eV=float(match.group("energy")),
                    fmax_eV_per_A=float(match.group("fmax")),
                )
            )
    if not records:
        raise ValueError(f"No CheckedFIRE records found in {path}")
    if [record.step for record in records] != list(range(records[-1].step + 1)):
        raise ValueError(f"Optimizer steps are not contiguous in {path}")
    return records


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_chain(trajectory: Trajectory, step: int, n_images: int) -> list:
    start = step * n_images
    stop = start + n_images
    if stop > len(trajectory):
        raise ValueError(
            f"Trajectory has {len(trajectory)} frames, but step {step} requires {stop}"
        )
    return [trajectory[index] for index in range(start, stop)]


def extended_coordinate(images: list) -> np.ndarray:
    reference = np.asarray(images[0].cell.array, dtype=float)
    cell_scale = abs(float(np.linalg.det(reference))) ** (1.0 / 3.0)
    vectors = []
    for atoms in images:
        fractional = atoms.get_scaled_positions(wrap=False)
        deformation = np.linalg.solve(reference, atoms.cell.array).T
        vectors.append(
            np.concatenate(
                [
                    (fractional @ reference).reshape(-1),
                    (cell_scale * (deformation - np.eye(3))).reshape(-1),
                ]
            )
        )
    coordinate = [0.0]
    for left, right in zip(vectors[:-1], vectors[1:]):
        coordinate.append(coordinate[-1] + float(np.linalg.norm(right - left)))
    coordinate = np.asarray(coordinate)
    if coordinate[-1] <= 0.0:
        raise ValueError("Path has zero extended-coordinate length")
    return coordinate / coordinate[-1]


def energy_profile(images: list) -> tuple[np.ndarray, np.ndarray]:
    energies = np.asarray([atoms.get_potential_energy() for atoms in images])
    n_atoms = len(images[0])
    relative_mev_per_atom = 1000.0 * (energies - energies[0]) / n_atoms
    return energies, relative_mev_per_atom


def snapshot_record(
    route: RouteInput,
    kind: str,
    record: LogRecord,
    images: list,
) -> tuple[dict, list[dict]]:
    coordinate = extended_coordinate(images)
    energies, relative = energy_profile(images)
    rows = []
    for image_index, (atoms, x, energy, delta) in enumerate(
        zip(images, coordinate, energies, relative)
    ):
        rows.append(
            {
                "route": route.key,
                "snapshot": kind,
                "optimizer_step": record.step,
                "fmax_eV_per_A": record.fmax_eV_per_A,
                "image_index": image_index,
                "path_coordinate_normalized": float(x),
                "energy_eV_total": float(energy),
                "relative_energy_meV_per_atom": float(delta),
                "volume_A3": float(atoms.get_volume()),
            }
        )
    summary = {
        "optimizer_step": record.step,
        "fmax_eV_per_A": record.fmax_eV_per_A,
        "barrier_meV_per_atom": float(relative.max()),
        "barrier_image_index": int(np.argmax(relative)),
        "reaction_energy_meV_per_atom": float(relative[-1]),
    }
    return summary, rows


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 7,
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "axes.spines.right": False,
            "axes.spines.top": False,
        }
    )


def plot_figure(route_data: list[dict], output_stem: Path) -> None:
    configure_matplotlib()
    latest_color = "#D55E00"
    best_color = "#0072B2"
    threshold_color = "#6E6E6E"

    fig = plt.figure(figsize=(7.2, 5.0), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[3.2, 1.25], hspace=0.08, wspace=0.18)
    energy_axes = [fig.add_subplot(grid[0, column]) for column in range(2)]
    force_axes = [fig.add_subplot(grid[1, column]) for column in range(2)]

    all_relative = [
        row["relative_energy_meV_per_atom"]
        for route in route_data
        for snapshot in route["snapshots"].values()
        for row in snapshot["rows"]
    ]
    y_min = min(all_relative)
    y_max = max(all_relative)
    padding = 0.06 * max(1.0, y_max - y_min)

    for column, route in enumerate(route_data):
        ax = energy_axes[column]
        best = route["snapshots"]["lowest_force"]
        latest = route["snapshots"]["latest"]
        for snapshot, color, linestyle, marker, zorder, label in (
            (best, best_color, "--", "o", 2, "Lowest-force snapshot"),
            (latest, latest_color, "-", "s", 3, "Latest snapshot"),
        ):
            rows = snapshot["rows"]
            ax.plot(
                [row["path_coordinate_normalized"] for row in rows],
                [row["relative_energy_meV_per_atom"] for row in rows],
                color=color,
                linestyle=linestyle,
                linewidth=1.5,
                marker=marker,
                markersize=3.2,
                markerfacecolor="white" if marker == "o" else color,
                markeredgewidth=0.8,
                label=(
                    f"{label}: step {snapshot['summary']['optimizer_step']}, "
                    f"$f_{{max}}$={snapshot['summary']['fmax_eV_per_A']:.3f} eV Å$^{{-1}}$"
                ),
                zorder=zorder,
            )
        ax.axhline(0.0, color="#B8B8B8", linewidth=0.7, zorder=0)
        ax.set_title(route["title"], loc="left", fontweight="bold")
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(y_min - padding, y_max + padding)
        ax.set_xlabel("Normalized extended path coordinate")
        if column == 0:
            ax.set_ylabel("Relative energy (meV atom$^{-1}$)")
        else:
            ax.tick_params(labelleft=False)
        ax.legend(loc="lower left", handlelength=2.5)
        if best["summary"]["optimizer_step"] == latest["summary"]["optimizer_step"]:
            ax.text(
                0.98,
                0.96,
                "Latest = lowest-force snapshot",
                ha="right",
                va="top",
                transform=ax.transAxes,
                color="#404040",
                fontsize=7,
            )
        ax.text(
            -0.13,
            1.05,
            chr(ord("a") + column),
            transform=ax.transAxes,
            fontweight="bold",
            fontsize=10,
        )

        fax = force_axes[column]
        records = route["log_records"]
        steps = [record.step for record in records]
        forces = [record.fmax_eV_per_A for record in records]
        fax.plot(steps, forces, color="#333333", linewidth=1.2, marker="o", markersize=2.3)
        fax.axhline(0.10, color=threshold_color, linewidth=0.9, linestyle=":", label="Target 0.10")
        for snapshot, color, marker in ((best, best_color, "o"), (latest, latest_color, "s")):
            summary = snapshot["summary"]
            fax.scatter(
                summary["optimizer_step"],
                summary["fmax_eV_per_A"],
                s=28,
                color=color,
                marker=marker,
                edgecolor="white",
                linewidth=0.6,
                zorder=4,
            )
        fax.set_xlim(-1, max(steps) + 1)
        fax.set_ylim(0.0, max(forces) * 1.08)
        fax.set_xlabel("VCNEB optimization step")
        if column == 0:
            fax.set_ylabel("$f_{max}$ (eV Å$^{-1}$)")
        else:
            fax.tick_params(labelleft=False)
        fax.legend(loc="upper right")
        fax.text(
            -0.13,
            1.05,
            chr(ord("c") + column),
            transform=fax.transAxes,
            fontweight="bold",
            fontsize=10,
        )

    fig.suptitle(
        "CdSe rock-salt → wurtzite VCNEB: latest versus lowest-force paths",
        fontsize=10,
        fontweight="bold",
    )
    fig.savefig(output_stem.with_suffix(".png"), dpi=400, bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell-trajectory", type=Path, required=True)
    parser.add_argument("--cell-log", type=Path, required=True)
    parser.add_argument("--atomic-trajectory", type=Path, required=True)
    parser.add_argument("--atomic-log", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-images", type=int, default=17)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    routes = [
        RouteInput("cell_mapping", "Cell-mapped initialization", args.cell_trajectory, args.cell_log),
        RouteInput(
            "atomic_mapping",
            "Atomic-mapped initialization",
            args.atomic_trajectory,
            args.atomic_log,
        ),
    ]
    route_data = []
    source_rows: list[dict] = []
    force_rows: list[dict] = []
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "n_images_total": args.n_images,
        "energy_normalization": "meV per atom relative to image 0 of each snapshot",
        "path_coordinate": "cumulative VARNEB extended-coordinate distance normalized to [0, 1]",
        "routes": {},
    }

    for route in routes:
        records = parse_log(route.logfile)
        latest_record = records[-1]
        best_record = min(records, key=lambda record: (record.fmax_eV_per_A, record.step))
        trajectory = Trajectory(str(route.trajectory), mode="r")
        if len(trajectory) % args.n_images:
            raise ValueError(f"Incomplete image group in {route.trajectory}")
        available_steps = len(trajectory) // args.n_images
        if latest_record.step >= available_steps:
            raise ValueError(
                f"Latest logged step {latest_record.step} is unavailable in {route.trajectory}"
            )

        snapshots = {}
        route_metadata = {
            "trajectory_file": route.trajectory.name,
            "trajectory_sha256": file_sha256(route.trajectory),
            "log_file": route.logfile.name,
            "log_sha256": file_sha256(route.logfile),
            "available_complete_chains": available_steps,
            "snapshots": {},
        }
        for kind, record in (("lowest_force", best_record), ("latest", latest_record)):
            images = read_chain(trajectory, record.step, args.n_images)
            summary, rows = snapshot_record(route, kind, record, images)
            snapshots[kind] = {"summary": summary, "rows": rows}
            source_rows.extend(rows)
            route_metadata["snapshots"][kind] = summary
            output_traj = args.output_dir / (
                f"{route.key}_{kind}_step{record.step:04d}.traj"
            )
            write(str(output_traj), images)
            route_metadata["snapshots"][kind]["extracted_trajectory"] = output_traj.name
            route_metadata["snapshots"][kind]["extracted_sha256"] = file_sha256(output_traj)

        best_step = best_record.step
        latest_step = latest_record.step
        for record in records:
            force_rows.append(
                {
                    "route": route.key,
                    "optimizer_step": record.step,
                    "optimizer_energy_eV": record.energy_eV,
                    "fmax_eV_per_A": record.fmax_eV_per_A,
                    "is_lowest_force_snapshot": record.step == best_step,
                    "is_latest_snapshot": record.step == latest_step,
                }
            )
        route_data.append(
            {
                "key": route.key,
                "title": route.title,
                "snapshots": snapshots,
                "log_records": records,
            }
        )
        metadata["routes"][route.key] = route_metadata
        trajectory.close()

    source_csv = args.output_dir / "cdse_path_snapshot_source_data.csv"
    with source_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(source_rows[0]))
        writer.writeheader()
        writer.writerows(source_rows)
    force_csv = args.output_dir / "cdse_force_history_source_data.csv"
    with force_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(force_rows[0]))
        writer.writeheader()
        writer.writerows(force_rows)

    output_stem = args.output_dir / "cdse_latest_vs_lowest_force_paths"
    plot_figure(route_data, output_stem)
    metadata["figure_files"] = [
        output_stem.with_suffix(suffix).name for suffix in (".png", ".svg", ".pdf", ".tiff")
    ]
    metadata["source_data_files"] = [source_csv.name, force_csv.name]
    (args.output_dir / "cdse_path_snapshot_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
