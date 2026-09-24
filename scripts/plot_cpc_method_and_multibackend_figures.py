"""Build the CPC architecture and GaN multi-backend validation figures.

Figure contract
---------------
Core conclusion: VARNEB separates the path controller from static calculators,
and five independently converged first-principles backends recover the same
dominant GaN B4-to-B1 barrier topology at 45.7 GPa.
Archetypes: schematic-led composite (architecture) and quantitative grid
(backend validation).
Backend: Python/matplotlib only.
Output: 183-mm editable SVG/PDF plus 600-dpi PNG and a source-data CSV.
Reviewer risk: backend curves use distinct basis/pseudopotential contracts;
the figure therefore compares topology and barrier scale, not total energies.
Only force-converged production paths are plotted.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams.update(
    {
        "font.size": 8.5,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.linewidth": 0.9,
        "axes.spines.top": True,
        "axes.spines.right": True,
        "legend.frameon": True,
        "legend.framealpha": 0.84,
        "legend.facecolor": "white",
        "legend.edgecolor": "#7A7A7A",
    }
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EVIDENCE = ROOT / "paper" / "VARNEB_CPC" / "evidence" / "gan_45p7_multibackend_vcneb_20260924.json"
DEFAULT_VASP = ROOT / "outputs" / "neb_literature_benchmarks" / "gan_b4_to_b1_tetragonal_literature_comparison_source_data.csv"
DEFAULT_OUTPUT = ROOT / "paper" / "VARNEB_CPC" / "figures"

COLORS = {
    "ABACUS": "#0F4D92",
    "VASP": "#B64342",
    "QE": "#42949E",
    "ABINIT": "#9A4D8E",
    "CP2K": "#A77722",
    "literature": "#272727",
    "core": "#DCE9F6",
    "worker": "#DDF3DE",
    "evidence": "#F0E0D0",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument("--vasp-source", type=Path, default=DEFAULT_VASP)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def read_vasp_and_literature(path: Path) -> tuple[list[float], list[float], list[float]]:
    vasp: list[tuple[int, float]] = []
    literature: list[tuple[float, float]] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row["dataset"].startswith("VARNEB"):
                vasp.append((int(row["image_index"]), float(row["relative_energy"])))
            elif row["dataset"].startswith("Qian"):
                literature.append((float(row["normalized_coordinate"]), float(row["relative_energy"])))
    vasp.sort()
    literature.sort()
    if len(vasp) != 29:
        raise ValueError(f"expected 29 VASP images, found {len(vasp)}")
    return [index / 28.0 for index, _ in vasp], [value for _, value in vasp], [value for _, value in literature]


def build_series(payload: dict, vasp_source: Path) -> tuple[dict[str, tuple[np.ndarray, np.ndarray]], tuple[np.ndarray, np.ndarray]]:
    vasp_x, vasp_y, literature_y = read_vasp_and_literature(vasp_source)
    series: dict[str, tuple[np.ndarray, np.ndarray]] = {
        "VASP": (np.asarray(vasp_x), np.asarray(vasp_y)),
    }
    for key, label in (("abacus", "ABACUS"), ("qe", "QE"), ("abinit", "ABINIT"), ("cp2k", "CP2K")):
        values = np.asarray(payload["backends"][key]["relative_enthalpy_eV_per_cell"], dtype=float) / 2.0
        series[label] = (np.linspace(0.0, 1.0, len(values)), values)
    literature = (np.linspace(0.0, 1.0, len(literature_y)), np.asarray(literature_y))
    return series, literature


def write_source_data(path: Path, payload: dict, series: dict[str, tuple[np.ndarray, np.ndarray]], literature: tuple[np.ndarray, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "backend",
                "status",
                "image_index",
                "normalized_image_index",
                "relative_enthalpy_eV_per_GaN",
                "barrier_eV_per_GaN",
                "final_fmax_eV_per_A",
                "provenance",
            ),
        )
        writer.writeheader()
        for label, (coordinate, values) in series.items():
            record = payload["backends"][label.lower()]
            for index, (x_value, energy) in enumerate(zip(coordinate, values)):
                writer.writerow(
                    {
                        "backend": label,
                        "status": record["status"],
                        "image_index": index,
                        "normalized_image_index": f"{x_value:.12g}",
                        "relative_enthalpy_eV_per_GaN": f"{energy:.12g}",
                        "barrier_eV_per_GaN": f"{record['barrier_eV_per_GaN']:.12g}",
                        "final_fmax_eV_per_A": f"{record['final_max_generalized_force_eV_per_A']:.12g}",
                        "provenance": f"Slurm job {record['job_id']}",
                    }
                )
        for index, (x_value, energy) in enumerate(zip(*literature), start=1):
            writer.writerow(
                {
                    "backend": "Qian et al.",
                    "status": "digitized reference",
                    "image_index": index,
                    "normalized_image_index": f"{x_value:.12g}",
                    "relative_enthalpy_eV_per_GaN": f"{energy:.12g}",
                    "barrier_eV_per_GaN": payload["literature"]["barrier_eV_per_GaN"],
                    "final_fmax_eV_per_A": "",
                    "provenance": "Qian et al. (2013), Fig. 4; approximate digitization",
                }
            )


def _save(figure: plt.Figure, base: Path) -> tuple[Path, ...]:
    base.parent.mkdir(parents=True, exist_ok=True)
    outputs = []
    for suffix, kwargs in ((".svg", {}), (".pdf", {}), (".png", {"dpi": 600})):
        output = base.with_suffix(suffix)
        figure.savefig(output, bbox_inches="tight", **kwargs)
        if suffix == ".svg":
            output.write_text(
                "\n".join(line.rstrip() for line in output.read_text(encoding="utf-8").splitlines()) + "\n",
                encoding="utf-8",
            )
        outputs.append(output)
    plt.close(figure)
    return tuple(outputs)


def plot_architecture(output: Path) -> tuple[Path, ...]:
    figure, axis = plt.subplots(figsize=(7.2, 2.55))
    axis.set_xlim(0, 10)
    axis.set_ylim(0, 4.1)
    axis.axis("off")

    def box(x: float, y: float, w: float, h: float, text: str, color: str, *, size: float = 8) -> None:
        axis.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.08", facecolor=color, edgecolor="#4D4D4D", linewidth=0.9))
        axis.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=size)

    def arrow(x0: float, y0: float, x1: float, y1: float) -> None:
        axis.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=10, linewidth=0.95, color="#4D4D4D"))

    box(0.08, 2.45, 1.72, 1.05, "Endpoint identity\n+ path preflight", COLORS["evidence"])
    box(2.12, 2.45, 2.20, 1.05, "VCNEB manager\npath force + convergence", COLORS["core"])
    box(4.65, 2.45, 1.57, 1.05, "Active interior\nimage queue", COLORS["worker"])
    box(6.55, 2.45, 3.35, 1.05, "Isolated workers\nstatic $E, F, \\sigma$ per active image", COLORS["worker"])
    arrow(1.80, 2.98, 2.12, 2.98)
    arrow(4.32, 2.98, 4.65, 2.98)
    arrow(6.22, 2.98, 6.55, 2.98)

    box(0.08, 0.88, 1.72, 0.90, "Endpoint + input\nprovenance", COLORS["evidence"], size=7.8)
    box(2.12, 0.88, 2.20, 0.90, "Independent update policy\nFIRE / block / staged", "#E8E0F3", size=7.6)
    box(4.65, 0.88, 1.57, 0.90, "Cache + complete\nsnapshots", COLORS["evidence"], size=7.6)
    box(6.55, 0.88, 3.35, 0.90, "ASE calculator contract\nABACUS | VASP | QE | CP2K | ABINIT | LAMMPS", "#F2F2F2", size=7.3)
    arrow(0.94, 2.45, 0.94, 1.78)
    arrow(3.22, 1.78, 3.22, 2.45)
    arrow(5.43, 2.45, 5.43, 1.78)
    arrow(8.23, 2.45, 8.23, 1.78)

    box(0.08, 0.14, 9.82, 0.50, "Saved chain + provenance  $\\longrightarrow$  independent audit + $\\Gamma$-mode interpretation", COLORS["evidence"], size=8)
    arrow(5.43, 0.88, 5.43, 0.64)
    figure.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.06)
    return _save(figure, output)


def plot_multibackend(output: Path, payload: dict, series: dict[str, tuple[np.ndarray, np.ndarray]], literature: tuple[np.ndarray, np.ndarray]) -> tuple[Path, ...]:
    figure = plt.figure(figsize=(7.2, 3.95))
    grid = figure.add_gridspec(2, 2, width_ratios=(1.58, 1), wspace=0.43, hspace=0.42)
    path_axis = figure.add_subplot(grid[:, 0])
    barrier_axis = figure.add_subplot(grid[0, 1])
    force_axis = figure.add_subplot(grid[1, 1])
    figure.subplots_adjust(left=0.11, right=0.975, top=0.92, bottom=0.15)

    labels = ("ABACUS", "VASP", "QE", "ABINIT", "CP2K")
    for label in labels:
        coordinate, values = series[label]
        path_axis.plot(coordinate, values, "o-", ms=3.1, lw=1.55, color=COLORS[label], label=label)
    path_axis.plot(literature[0], literature[1], linestyle="--", lw=1.25, color=COLORS["literature"], label="Qian et al.")
    path_axis.set(xlabel="Normalized image index", ylabel="Relative enthalpy (eV/GaN)", xlim=(-0.025, 1.025), ylim=(-0.078, 0.39))
    path_axis.legend(fontsize=8, ncol=1, loc="upper left", handlelength=1.8, borderpad=0.4)

    y_positions = np.arange(len(labels))[::-1]
    barrier_axis.axvline(payload["literature"]["barrier_eV_per_GaN"], color=COLORS["literature"], linestyle="--", lw=1.1, zorder=0)
    force_axis.axvline(payload["contract"]["fmax_target_eV_per_A"], color=COLORS["literature"], linestyle="--", lw=1.1, zorder=0)

    for label, y_value in zip(labels, y_positions, strict=True):
        record = payload["backends"][label.lower()]
        barrier = record["barrier_eV_per_GaN"]
        force = record["final_max_generalized_force_eV_per_A"]
        for axis, value, annotation_x in ((barrier_axis, barrier, 0.365), (force_axis, force, 0.108)):
            axis.plot([0, value], [y_value, y_value], color=COLORS[label], lw=1.25, alpha=0.65, zorder=1)
            axis.scatter([value], [y_value], color=COLORS[label], edgecolor="white", linewidth=0.7, s=45, zorder=3)
            axis.text(annotation_x, y_value, f"{value:.3f}", ha="left", va="center", fontsize=8, color="#272727")

    for axis, xlabel, xlim in (
        (barrier_axis, "Barrier (eV/GaN)", (0, 0.48)),
        (force_axis, r"Final $f_{\max}$ (eV/$\AA$)", (0, 0.15)),
    ):
        axis.set(xlim=xlim, ylim=(-0.6, 4.6), xlabel=xlabel)
        axis.set_yticks(y_positions, labels)
        axis.tick_params(direction="out", length=3.2, width=0.8, top=True, right=True)

    for label, axis in zip(("(a)", "(b)", "(c)"), (path_axis, barrier_axis, force_axis)):
        axis.text(-0.15, 1.04, label, transform=axis.transAxes, fontweight="normal", fontsize=10.5, va="bottom", ha="left")
        axis.tick_params(direction="out", length=3.2, width=0.8, top=True, right=True)
    return _save(figure, output)


def main() -> int:
    args = parse_args()
    payload = json.loads(args.evidence.read_text(encoding="utf-8"))
    series, literature = build_series(payload, args.vasp_source)
    source_data = args.output_dir / "gan_multibackend_validation_source_data.csv"
    write_source_data(source_data, payload, series, literature)
    outputs = [source_data]
    outputs.extend(plot_architecture(args.output_dir / "varneb_architecture"))
    outputs.extend(plot_multibackend(args.output_dir / "gan_multibackend_validation", payload, series, literature))
    for output in outputs:
        print(f"wrote={output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
