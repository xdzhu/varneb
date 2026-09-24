"""Create a traceable BTO/HfO2 VC-NEB validation figure from completed summaries.

The script is deliberately calculator-free: it reads only archived, completed
``vcneb_summary.json`` files.  It writes source-data CSV files alongside
editable SVG, PDF, and preview PNG figures.  Literature values are explicit
arguments so a different like-for-like reference can be substituted without
changing the plotting implementation.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


HFO2_FORMULA_UNITS_PER_CELL = 4
KCAL_PER_MOL_PER_EV = 23.060547830619


@dataclass(frozen=True)
class PathSeries:
    """One completed path, normalized by its total geometric arc length."""

    label: str
    material: str
    image_count: int
    coordinate: tuple[float, ...]
    relative_enthalpy_eV: tuple[float, ...]
    volume_A3: tuple[float, ...]
    climbing_image_index: int | None

    @property
    def barrier_eV(self) -> float:
        return max(self.relative_enthalpy_eV)

    @property
    def formula_units(self) -> int:
        if self.material == "BaTiO3":
            return 1
        if self.material == "HfO2":
            return HFO2_FORMULA_UNITS_PER_CELL
        raise ValueError(f"unknown formula-unit count for {self.material}")

    @property
    def relative_enthalpy_per_formula_unit_eV(self) -> tuple[float, ...]:
        return tuple(value / self.formula_units for value in self.relative_enthalpy_eV)


def _coordinates(summary: dict, n_images: int) -> tuple[float, ...]:
    geometry = (summary.get("path_diagnostics") or {}).get("geometry") or {}
    lengths = [float(value) for value in geometry.get("segment_lengths_A") or []]
    if len(lengths) != n_images - 1:
        raise ValueError("segment_lengths_A does not match the number of images")
    coordinate = [0.0]
    for length in lengths:
        coordinate.append(coordinate[-1] + length)
    if coordinate[-1] <= 0.0:
        raise ValueError("path has zero total geometric length")
    return tuple(value / coordinate[-1] for value in coordinate)


def read_completed_summary(path: Path, *, label: str, material: str) -> PathSeries:
    """Read one completed VC-NEB summary and retain only plot-ready evidence."""

    try:
        summary = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read {path}: {exc}") from exc
    if summary.get("status") != "completed":
        raise ValueError(f"{path} is not a completed VC-NEB summary")
    diagnostics = summary.get("path_diagnostics") or {}
    images = diagnostics.get("images") or []
    n_images = int(summary.get("n_images", diagnostics.get("n_images", len(images))))
    if n_images < 2 or len(images) != n_images:
        raise ValueError(f"{path} has inconsistent image diagnostics")
    ordered = sorted(images, key=lambda item: int(item["image_index"]))
    if [int(item["image_index"]) for item in ordered] != list(range(n_images)):
        raise ValueError(f"{path} image indices are not contiguous")
    relative_enthalpy = tuple(float(item["relative_enthalpy_eV"]) for item in ordered)
    volume = tuple(float(item["volume_A3"]) for item in ordered)
    climbing = [
        int(item["image_index"])
        for item in ordered
        if bool(item.get("is_climbing_image", False))
    ]
    if len(climbing) > 1:
        raise ValueError(f"{path} has more than one climbing image")
    return PathSeries(
        label=label,
        material=material,
        image_count=n_images,
        coordinate=_coordinates(summary, n_images),
        relative_enthalpy_eV=relative_enthalpy,
        volume_A3=volume,
        climbing_image_index=climbing[0] if climbing else None,
    )


def write_source_data(
    output_dir: Path,
    *,
    bto_paths: Iterable[PathSeries],
    hfo2_paths: Iterable[PathSeries],
    hfo2_literature_barrier_eV_per_fu: float,
    bto_literature_barrier_kcal_per_mol: float,
) -> tuple[Path, Path]:
    """Write long-form path and barrier source-data tables."""

    output_dir.mkdir(parents=True, exist_ok=True)
    path_data = output_dir / "vcneb_path_comparison_source_data.csv"
    with path_data.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "material",
                "path_label",
                "total_images",
                "image_index",
                "normalized_reaction_coordinate",
                "relative_enthalpy_eV_per_cell",
                "relative_enthalpy_eV_per_formula_unit",
                "volume_A3",
                "is_climbing_image",
            ),
        )
        writer.writeheader()
        for series in (*tuple(bto_paths), *tuple(hfo2_paths)):
            for index, (coordinate, enthalpy, enthalpy_per_fu, volume) in enumerate(
                zip(series.coordinate, series.relative_enthalpy_eV, series.relative_enthalpy_per_formula_unit_eV, series.volume_A3)
            ):
                writer.writerow(
                    {
                        "material": series.material,
                        "path_label": series.label,
                        "total_images": series.image_count,
                        "image_index": index,
                        "normalized_reaction_coordinate": f"{coordinate:.12g}",
                        "relative_enthalpy_eV_per_cell": f"{enthalpy:.12g}",
                        "relative_enthalpy_eV_per_formula_unit": f"{enthalpy_per_fu:.12g}",
                        "volume_A3": f"{volume:.12g}",
                        "is_climbing_image": index == series.climbing_image_index,
                    }
                )

    barrier_data = output_dir / "vcneb_barrier_literature_comparison.csv"
    with barrier_data.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "material",
                "method",
                "barrier_eV_per_formula_unit",
                "barrier_kcal_per_mol_per_formula_unit",
                "reference_note",
            ),
        )
        writer.writeheader()
        for series in bto_paths:
            writer.writerow(
                {
                    "material": "BaTiO3",
                    "method": series.label,
                    "barrier_eV_per_formula_unit": f"{series.barrier_eV:.12g}",
                    "barrier_kcal_per_mol_per_formula_unit": f"{series.barrier_eV * KCAL_PER_MOL_PER_EV:.12g}",
                    "reference_note": "ABACUS PBE 100 Ry, 10 au DZP",
                }
            )
        writer.writerow(
            {
                "material": "BaTiO3",
                "method": "Literature restrained NEB",
                "barrier_eV_per_formula_unit": f"{bto_literature_barrier_kcal_per_mol / KCAL_PER_MOL_PER_EV:.12g}",
                "barrier_kcal_per_mol_per_formula_unit": f"{bto_literature_barrier_kcal_per_mol:.12g}",
                "reference_note": "PBEsol restrained NEB; not a like-for-like VC-NEB benchmark",
            }
        )
        for series in hfo2_paths:
            writer.writerow(
                {
                    "material": "HfO2",
                    "method": series.label,
                    "barrier_eV_per_formula_unit": f"{series.barrier_eV / HFO2_FORMULA_UNITS_PER_CELL:.12g}",
                    "barrier_kcal_per_mol_per_formula_unit": f"{series.barrier_eV * KCAL_PER_MOL_PER_EV / HFO2_FORMULA_UNITS_PER_CELL:.12g}",
                    "reference_note": "ABACUS PBE 100 Ry, 10 au DZP",
                }
            )
        writer.writerow(
            {
                "material": "HfO2",
                "method": "Literature VC-NEB",
                "barrier_eV_per_formula_unit": f"{hfo2_literature_barrier_eV_per_fu:.12g}",
                "barrier_kcal_per_mol_per_formula_unit": f"{hfo2_literature_barrier_eV_per_fu * KCAL_PER_MOL_PER_EV:.12g}",
                "reference_note": "LDA/QE/USPEX, 40 images, RMS force threshold 0.025 eV/A",
            }
        )
    return path_data, barrier_data


def plot_validation_figure(
    output_base: Path,
    *,
    bto_paths: tuple[PathSeries, ...],
    hfo2_paths: tuple[PathSeries, ...],
    hfo2_literature_barrier_eV_per_fu: float,
    bto_literature_barrier_kcal_per_mol: float,
) -> tuple[Path, ...]:
    """Render the four-panel evidence figure with editable vector text."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
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
    colors = {
        "n5": "#9A4D8E",
        "n7": "#0F4D92",
        "n9": "#42949E",
        "ci": "#B64342",
        "linear": "#767676",
        "literature": "#272727",
    }

    figure, axes = plt.subplots(2, 2, figsize=(7.205, 4.85))
    bto_ax, hfo2_ax, barrier_ax, volume_ax = axes.flat
    figure.subplots_adjust(left=0.115, right=0.975, top=0.91, bottom=0.205, wspace=0.34, hspace=0.39)
    hfo2_handles = []

    def draw_paths(axis, paths, *, collect_handles: bool = False) -> None:
        for series in paths:
            key = "ci" if "CI" in series.label else "linear" if "linear" in series.label.lower() else f"n{series.image_count}"
            values = series.relative_enthalpy_per_formula_unit_eV
            line, = axis.plot(series.coordinate, values, marker="o", markersize=3.5, linewidth=1.55, color=colors[key], label=series.label)
            if collect_handles:
                hfo2_handles.append(line)
            if series.climbing_image_index is not None:
                index = series.climbing_image_index
                axis.plot(series.coordinate[index], values[index], marker="o", markersize=7.0, markerfacecolor="white", markeredgecolor="#272727", markeredgewidth=1.0, linestyle="None", zorder=6)
        axis.set(xlim=(-0.03, 1.03), xlabel="Normalized path coordinate", ylabel="Relative enthalpy (eV/f.u.)")
        axis.set_xticks([0, 0.25, 0.5, 0.75, 1.0])

    draw_paths(bto_ax, bto_paths)
    bto_ax.legend(loc="upper left", fontsize=8, handlelength=1.7, borderpad=0.35)
    draw_paths(hfo2_ax, hfo2_paths, collect_handles=True)
    literature_line = hfo2_ax.axhline(hfo2_literature_barrier_eV_per_fu, linestyle=":", linewidth=1.15, color=colors["literature"], zorder=0)
    hfo2_handles.append(literature_line)

    hfo2_barriers = [series.barrier_eV / HFO2_FORMULA_UNITS_PER_CELL for series in hfo2_paths]
    def short_label(series: PathSeries) -> str:
        if "CI" in series.label:
            return "CI n7"
        if "linear" in series.label.lower():
            return "linear n7"
        return f"n{series.image_count}"

    hfo2_labels = [short_label(series) for series in hfo2_paths] + ["lit."]
    values_eV_per_fu = hfo2_barriers + [hfo2_literature_barrier_eV_per_fu]
    values_meV_per_fu = [value * 1000.0 for value in values_eV_per_fu]
    bar_colors = [colors["n7"], colors["n9"], colors["ci"], colors["linear"], colors["literature"]]
    positions = list(range(len(values_meV_per_fu)))
    bars = barrier_ax.bar(positions, values_meV_per_fu, color=bar_colors[: len(values_meV_per_fu)], edgecolor="#272727", linewidth=0.45)
    for bar, value in zip(bars, values_meV_per_fu):
        barrier_ax.text(bar.get_x() + bar.get_width() / 2, value + 0.9, f"{value:.1f}", ha="center", va="bottom", fontsize=8)
    barrier_ax.set(xticks=positions, xticklabels=hfo2_labels, ylabel="Barrier (meV/f.u.)")
    barrier_ax.set_ylim(0, max(values_meV_per_fu) * 1.2)

    for series in hfo2_paths:
        key = "ci" if "CI" in series.label else "linear" if "linear" in series.label.lower() else f"n{series.image_count}"
        normalized_volume = [(value / series.volume_A3[0] - 1.0) * 100.0 for value in series.volume_A3]
        volume_ax.plot(series.coordinate, normalized_volume, marker="o", markersize=3.3, linewidth=1.45, color=colors[key], label=series.label)
    volume_ax.axhline(0.0, color="#767676", linewidth=0.65)
    volume_ax.set(xlim=(-0.03, 1.03), xlabel="Normalized path coordinate", ylabel="$\\Delta V/V_0$ (%)")
    volume_ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])

    figure.legend(
        hfo2_handles,
        [series.label for series in hfo2_paths] + ["literature"],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.018),
        ncol=3,
        fontsize=8,
        frameon=True,
        framealpha=0.84,
        handlelength=1.7,
        columnspacing=1.1,
    )

    for label, axis in zip(("(a)", "(b)", "(c)", "(d)"), axes.flat):
        axis.text(-0.15, 1.035, label, transform=axis.transAxes, fontweight="normal", fontsize=10.5, va="bottom", ha="left")
        axis.tick_params(direction="out", length=3.2, width=0.8, top=True, right=True)

    output_base.parent.mkdir(parents=True, exist_ok=True)
    paths = []
    for suffix, kwargs in (
        (".svg", {}),
        (".pdf", {}),
        (".png", {"dpi": 600}),
    ):
        output = output_base.with_suffix(suffix)
        figure.savefig(output, bbox_inches="tight", **kwargs)
        if suffix == ".svg":
            # Matplotlib intentionally wraps long SVG path commands with a
            # trailing space.  Normalize generated line endings so the source
            # artifact passes the repository whitespace gate unchanged.
            output.write_text(
                "\n".join(line.rstrip() for line in output.read_text(encoding="utf-8").splitlines())
                + "\n",
                encoding="utf-8",
            )
        paths.append(output)
    plt.close(figure)
    return tuple(paths)


def _path_argument(value: str) -> tuple[str, Path]:
    try:
        label, raw_path = value.split("=", 1)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("path values must use LABEL=SUMMARY_JSON") from exc
    if not label:
        raise argparse.ArgumentTypeError("path label must not be empty")
    return label, Path(raw_path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bto", action="append", required=True, type=_path_argument, metavar="LABEL=SUMMARY_JSON")
    parser.add_argument("--hfo2", action="append", required=True, type=_path_argument, metavar="LABEL=SUMMARY_JSON")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--hfo2-literature-barrier-eV-per-fu", type=float, default=0.032)
    parser.add_argument("--bto-literature-barrier-kcal-per-mol", type=float, default=2.1)
    args = parser.parse_args()
    try:
        bto_paths = tuple(read_completed_summary(path, label=label, material="BaTiO3") for label, path in args.bto)
        hfo2_paths = tuple(read_completed_summary(path, label=label, material="HfO2") for label, path in args.hfo2)
        if not any("CI" in series.label for series in hfo2_paths):
            raise ValueError("HfO2 inputs must include one explicitly labeled CI path")
        source_data = write_source_data(
            args.output_dir,
            bto_paths=bto_paths,
            hfo2_paths=hfo2_paths,
            hfo2_literature_barrier_eV_per_fu=args.hfo2_literature_barrier_eV_per_fu,
            bto_literature_barrier_kcal_per_mol=args.bto_literature_barrier_kcal_per_mol,
        )
        figure_paths = plot_validation_figure(
            args.output_dir / "vcneb_material_validation",
            bto_paths=bto_paths,
            hfo2_paths=hfo2_paths,
            hfo2_literature_barrier_eV_per_fu=args.hfo2_literature_barrier_eV_per_fu,
            bto_literature_barrier_kcal_per_mol=args.bto_literature_barrier_kcal_per_mol,
        )
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    for path in (*source_data, *figure_paths):
        print(f"wrote={path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
