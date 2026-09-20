"""Plot one literature comparison figure for every converged VARNEB path.

The script reads archived ``vcneb_summary.json`` files, converts energies to
the normalization used by the corresponding paper, and writes one PNG, SVG,
PDF, and long-form source-data CSV per physical path.  Literature curves are
either copied from author-provided source data (HfO2), or explicitly marked as
figure-digitized approximations (GaN and CdSe).  BaTiO3 has only a published
barrier/endpoint-rise value in the evidence currently archived by the project.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


KCAL_PER_MOL_PER_EV = 23.060547830619


@dataclass(frozen=True)
class Curve:
    label: str
    coordinate: tuple[float, ...]
    energy: tuple[float, ...]
    unit: str
    provenance: str
    source_type: str
    x_definition: str
    image_index: tuple[int, ...]
    notes: str = ""

    @property
    def barrier(self) -> float:
        return max(self.energy)


@dataclass(frozen=True)
class FigureSpec:
    stem: str
    title: str
    ylabel: str
    current: Curve
    literature: Curve | None
    literature_barrier: float
    literature_label: str
    comparison_note: str


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _normalized_arc_coordinate(summary: dict, n_images: int) -> tuple[float, ...]:
    geometry = (summary.get("path_diagnostics") or {}).get("geometry") or {}
    lengths = [float(value) for value in geometry.get("segment_lengths_A") or []]
    if len(lengths) != n_images - 1 or sum(lengths) <= 0.0:
        return tuple(index / (n_images - 1) for index in range(n_images))
    cumulative = [0.0]
    for length in lengths:
        cumulative.append(cumulative[-1] + length)
    return tuple(value / cumulative[-1] for value in cumulative)


def read_completed_summary(
    path: Path,
    *,
    label: str,
    denominator: float,
    multiplier: float,
    unit: str,
    provenance: str,
    notes: str,
) -> Curve:
    summary = json.loads(path.read_text(encoding="utf-8"))
    if summary.get("status") != "completed":
        raise ValueError(f"{path} is not a completed summary")
    diagnostics = summary.get("path_diagnostics") or {}
    images = sorted(diagnostics.get("images") or [], key=lambda item: int(item["image_index"]))
    n_images = int(summary.get("n_images", len(images)))
    if len(images) != n_images or n_images < 2:
        raise ValueError(f"{path} has inconsistent path diagnostics")
    indices = tuple(int(item["image_index"]) for item in images)
    if indices != tuple(range(n_images)):
        raise ValueError(f"{path} image indices are not contiguous")
    energy = tuple(
        float(item["relative_enthalpy_eV"]) * multiplier / denominator for item in images
    )
    return Curve(
        label=label,
        coordinate=_normalized_arc_coordinate(summary, n_images),
        energy=energy,
        unit=unit,
        provenance=provenance,
        source_type="completed vcneb_summary.json",
        x_definition="normalized VARNEB geometric arc length",
        image_index=indices,
        notes=notes,
    )


def literature_curve(
    *,
    label: str,
    values: Sequence[float],
    unit: str,
    provenance: str,
    source_type: str,
    notes: str,
) -> Curve:
    n_values = len(values)
    return Curve(
        label=label,
        coordinate=tuple(index / (n_values - 1) for index in range(n_values)),
        energy=tuple(float(value) for value in values),
        unit=unit,
        provenance=provenance,
        source_type=source_type,
        x_definition="normalized published image index / digitized path distance",
        image_index=tuple(range(1, n_values + 1)),
        notes=notes,
    )


def _relative_per_formula_unit(energies_eV_per_cell: Sequence[float], n_formula_units: int) -> list[float]:
    reference = float(energies_eV_per_cell[0])
    return [(float(value) - reference) / n_formula_units for value in energies_eV_per_cell]


def write_source_data(path: Path, spec: FigureSpec) -> None:
    fields = (
        "dataset",
        "source_type",
        "image_index",
        "normalized_coordinate",
        "relative_energy",
        "unit",
        "x_definition",
        "provenance",
        "notes",
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        curves = [spec.current] + ([spec.literature] if spec.literature is not None else [])
        for curve in curves:
            for image_index, coordinate, energy in zip(
                curve.image_index, curve.coordinate, curve.energy
            ):
                writer.writerow(
                    {
                        "dataset": curve.label,
                        "source_type": curve.source_type,
                        "image_index": image_index,
                        "normalized_coordinate": f"{coordinate:.12g}",
                        "relative_energy": f"{energy:.12g}",
                        "unit": curve.unit,
                        "x_definition": curve.x_definition,
                        "provenance": curve.provenance,
                        "notes": curve.notes,
                    }
                )


def draw_figure(output_dir: Path, spec: FigureSpec) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
            "font.size": 8,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )
    figure, axis = plt.subplots(figsize=(3.54, 2.95), constrained_layout=False)
    figure.subplots_adjust(left=0.17, right=0.97, top=0.87, bottom=0.25)
    axis.plot(
        spec.current.coordinate,
        spec.current.energy,
        color="#0F4D92",
        marker="o",
        markersize=3.4,
        linewidth=1.45,
        label=spec.current.label,
        zorder=3,
    )
    if spec.literature is not None:
        axis.plot(
            spec.literature.coordinate,
            spec.literature.energy,
            color="#C45A23",
            marker="s",
            markersize=2.7,
            linewidth=1.2,
            linestyle="--",
            label=spec.literature.label,
            zorder=2,
        )
    else:
        axis.axhline(
            spec.literature_barrier,
            color="#C45A23",
            linewidth=1.1,
            linestyle="--",
            label=spec.literature_label,
            zorder=1,
        )
    axis.axhline(0.0, color="#777777", linewidth=0.55, zorder=0)
    axis.set(
        xlim=(-0.025, 1.025),
        xlabel="Normalized path coordinate",
        ylabel=spec.ylabel,
        title=spec.title,
    )
    axis.xaxis.set_major_locator(MaxNLocator(5))
    axis.yaxis.set_major_locator(MaxNLocator(6))
    axis.tick_params(direction="out", length=3.0, width=0.7)
    axis.legend(loc="best", fontsize=6.5, handlelength=2.5)
    barrier_text = (
        f"Barrier: VARNEB {spec.current.barrier:.3g}; "
        f"literature {spec.literature_barrier:.3g} {spec.current.unit}"
    )
    figure.text(0.17, 0.105, barrier_text, ha="left", va="bottom", fontsize=6.5)
    figure.text(0.17, 0.055, spec.comparison_note, ha="left", va="bottom", fontsize=5.8, color="#555555")
    outputs: list[Path] = []
    for extension, dpi in (("png", 600), ("svg", None), ("pdf", None)):
        target = output_dir / f"{spec.stem}.{extension}"
        figure.savefig(target, dpi=dpi, bbox_inches="tight")
        outputs.append(target)
    plt.close(figure)
    return outputs


def build_specs(args: argparse.Namespace) -> tuple[FigureSpec, ...]:
    bto = read_completed_summary(
        args.bto,
        label="VARNEB (ABACUS/PBE)",
        denominator=1,
        multiplier=1.0,
        unit="eV/f.u.",
        provenance=str(args.bto),
        notes="7 total images; accepted under the project default 0.10 eV/A criterion",
    )
    hfo2_tpo = read_completed_summary(
        args.hfo2_tpo,
        label="VARNEB (VASP/PBE)",
        denominator=4,
        multiplier=1.0,
        unit="eV/f.u.",
        provenance=(
            f"{args.hfo2_tpo}; origin hf:/public/home/iai806/abacus/agent-runs/"
            "20260920-varneb-feasible-stage/hfo2_tpo_symprec3e3_candidate_staged_r1/"
            "vcneb_summary.json"
        ),
        notes="20 total images; 4 HfO2 formula units per cell",
    )
    hfo2_pom = read_completed_summary(
        args.hfo2_pom,
        label="VARNEB (VASP/PBE)",
        denominator=4,
        multiplier=1.0,
        unit="eV/f.u.",
        provenance=(
            f"{args.hfo2_pom}; origin hf:/public/home/iai806/abacus/agent-runs/"
            "20260920-varneb-feasible-stage/hfo2_pom_symprec3e3_candidate_staged_r2/"
            "vcneb_summary.json"
        ),
        notes="20 total images; 4 HfO2 formula units per cell",
    )
    gan_tet = read_completed_summary(
        args.gan_tetragonal,
        label="VARNEB (VASP/PBE)",
        denominator=2,
        multiplier=1.0,
        unit="eV/GaN",
        provenance=(
            f"{args.gan_tetragonal}; origin hf:/public/home/iai806/abacus/agent-runs/"
            "20260918-varneb-vasp-hf/gan_b4_b1_vasp_pbe_paw_qian/"
            "vcneb_hf_n29_w9_mpi32/vcneb_summary.json"
        ),
        notes="29 total images; 2 GaN formula units per cell; tetragonal mapping",
    )
    gan_hex = read_completed_summary(
        args.gan_hexagonal,
        label="VARNEB (VASP/PBE)",
        denominator=2,
        multiplier=1.0,
        unit="eV/GaN",
        provenance=(
            f"{args.gan_hexagonal}; origin hf:/public/home/iai806/abacus/agent-runs/"
            "20260920-varneb-accelerated-resumes/gan_hex_block01_r1/vcneb_summary.json"
        ),
        notes="29 total images; 2 GaN formula units per cell; hexagonal mapping",
    )
    gan_b3 = read_completed_summary(
        args.gan_b3,
        label="VARNEB (VASP/PBE)",
        denominator=4,
        multiplier=1.0,
        unit="eV/GaN",
        provenance=(
            f"{args.gan_b3}; origin hf:/public/home/iai806/abacus/agent-runs/"
            "20260918-varneb-gan-suite/cases/b3_b1_diagonal/vcneb_n29_w9_mpi32/"
            "vcneb_summary.json"
        ),
        notes="29 total images; 4 GaN formula units per cell; single diagonal mapping",
    )
    cdse = read_completed_summary(
        args.cdse_cell,
        label="VARNEB (VASP/PBE)",
        denominator=8,
        multiplier=1000.0,
        unit="meV/atom",
        provenance=(
            f"{args.cdse_cell}; origin hf:/public/home/iai806/abacus/agent-runs/"
            "20260920-varneb-accelerated-resumes/cdse_cell_block01/vcneb_summary.json"
        ),
        notes="17 total images; 8 atoms per cell; cell-mapping branch",
    )

    hfo2_tpo_lit = literature_curve(
        label="Liu et al. PRL 2023 source data",
        values=_relative_per_formula_unit(
            [-121.4266, -121.4245, -121.4202, -121.4166, -121.4168,
             -121.4218, -121.4324, -121.4491, -121.4713, -121.4984,
             -121.5293, -121.5629, -121.5975, -121.6319, -121.6645,
             -121.6939, -121.7185, -121.7377, -121.7499, -121.7542],
            4,
        ),
        unit="eV/f.u.",
        provenance="PRL 130, 226801 (2023), Fig. 2(a) author workbook, sheet T-pca21",
        source_type="author-provided workbook values",
        notes="20 images; raw cell energies divided by 4 HfO2 formula units",
    )
    hfo2_pom_lit = literature_curve(
        label="Liu et al. PRL 2023 source data",
        values=_relative_per_formula_unit(
            [-121.7542, -121.7416, -121.7309, -121.7359, -121.6946,
             -121.6347, -121.5306, -121.4012, -121.3339, -121.5419,
             -121.7559, -121.8845, -121.9478, -121.9995, -122.0464,
             -122.0419, -122.0478, -122.0554, -122.0743, -122.0916],
            4,
        ),
        unit="eV/f.u.",
        provenance="PRL 130, 226801 (2023), Fig. 2(a) author workbook, sheet pca21-M",
        source_type="author-provided workbook values",
        notes="20 images; raw cell energies divided by 4 HfO2 formula units",
    )
    gan_tet_lit = literature_curve(
        label="Qian et al. Fig. 4 (digitized)",
        values=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.005, 0.025,
                0.080, 0.165, 0.250, 0.315, 0.342, 0.330, 0.280, 0.225,
                0.170, 0.135, 0.100, 0.060, 0.025, 0.005, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0],
        unit="eV/GaN",
        provenance="Qian et al., Comput. Phys. Commun. 184, 2111-2118 (2013), Fig. 4",
        source_type="approximate manual digitization of published figure",
        notes="Figure resolution limits point precision; paper states 0.34 eV/formula",
    )
    gan_hex_lit = literature_curve(
        label="Qian et al. Fig. 6 (digitized)",
        values=[0.0, 0.015, 0.055, 0.105, 0.150, 0.190, 0.225, 0.255,
                0.280, 0.297, 0.307, 0.311, 0.315, 0.325, 0.332, 0.342,
                0.353, 0.367, 0.380, 0.388, 0.390, 0.380, 0.355, 0.335,
                0.292, 0.240, 0.150, 0.058, 0.0],
        unit="eV/GaN",
        provenance="Qian et al., Comput. Phys. Commun. 184, 2111-2118 (2013), Fig. 6",
        source_type="approximate manual digitization of published figure",
        notes="Figure resolution limits point precision; paper states 0.39 eV/formula",
    )
    gan_b3_lit = literature_curve(
        label="Qian et al. Fig. 8 (digitized)",
        values=[0.0, 0.10, 0.23, 0.38, 0.50, 0.575, 0.415, 0.30, 0.16,
                0.045, 0.0, 0.13, 0.275, 0.425, 0.57, 0.405, 0.295,
                0.165, 0.075, 0.015, 0.0, 0.16, 0.28, 0.415, 0.575,
                0.38, 0.25, 0.125, 0.0],
        unit="eV/GaN",
        provenance="Qian et al., Comput. Phys. Commun. 184, 2111-2118 (2013), Fig. 8(a)",
        source_type="approximate manual digitization of published figure",
        notes="Paper reports three peaks near 0.57 eV/formula at images 6, 15, and 25",
    )
    cdse_lit = literature_curve(
        label="Sheppard et al. Fig. 11 (digitized)",
        values=[0.0, 2.4, -20.0, -68.0, -110.0, -132.0, -145.0,
                -132.0, -138.0, -78.0, -100.0, -137.0, -148.0],
        unit="meV/atom",
        provenance="Sheppard et al., J. Chem. Phys. 136, 074103 (2012), Fig. 11",
        source_type="approximate manual digitization of published figure",
        notes="Composite atom-dominated a-c then cell-dominated c-e route; paper uses PW91",
    )

    return (
        FigureSpec(
            stem="batio3_t_to_c_literature_comparison",
            title="BaTiO$_3$ T$\\rightarrow$C",
            ylabel="Relative enthalpy (eV/f.u.)",
            current=bto,
            literature=None,
            literature_barrier=2.1 / KCAL_PER_MOL_PER_EV,
            literature_label="PBEsol restrained NEB (barrier only)",
            comparison_note="Literature archive provides a 2.1 kcal mol$^{-1}$ endpoint rise, not path points.",
        ),
        FigureSpec(
            stem="hfo2_t_to_po_literature_comparison",
            title="HfO$_2$ T$\\rightarrow$PO",
            ylabel="Relative enthalpy (eV/f.u.)",
            current=hfo2_tpo,
            literature=hfo2_tpo_lit,
            literature_barrier=hfo2_tpo_lit.barrier,
            literature_label=hfo2_tpo_lit.label,
            comparison_note="VARNEB: VASP/PBE; literature: author workbook underlying PRL Fig. 2(a).",
        ),
        FigureSpec(
            stem="hfo2_po_to_m_literature_comparison",
            title="HfO$_2$ PO$\\rightarrow$M",
            ylabel="Relative enthalpy (eV/f.u.)",
            current=hfo2_pom,
            literature=hfo2_pom_lit,
            literature_barrier=hfo2_pom_lit.barrier,
            literature_label=hfo2_pom_lit.label,
            comparison_note="VARNEB: VASP/PBE; literature: author workbook underlying PRL Fig. 2(a).",
        ),
        FigureSpec(
            stem="gan_b4_to_b1_tetragonal_literature_comparison",
            title="GaN B4$\\rightarrow$B1: tetragonal path",
            ylabel="Relative enthalpy (eV/GaN)",
            current=gan_tet,
            literature=gan_tet_lit,
            literature_barrier=0.34,
            literature_label=gan_tet_lit.label,
            comparison_note="Both at 45.7 GPa; literature curve is approximate Fig. 4 digitization.",
        ),
        FigureSpec(
            stem="gan_b4_to_b1_hexagonal_literature_comparison",
            title="GaN B4$\\rightarrow$B1: hexagonal path",
            ylabel="Relative enthalpy (eV/GaN)",
            current=gan_hex,
            literature=gan_hex_lit,
            literature_barrier=0.39,
            literature_label=gan_hex_lit.label,
            comparison_note="Both at 45.7 GPa; literature curve is approximate Fig. 6 digitization.",
        ),
        FigureSpec(
            stem="gan_b3_to_b1_literature_comparison",
            title="GaN B3$\\rightarrow$B1",
            ylabel="Relative enthalpy (eV/GaN)",
            current=gan_b3,
            literature=gan_b3_lit,
            literature_barrier=0.57,
            literature_label=gan_b3_lit.label,
            comparison_note="45.0 GPa; mechanism mismatch is visible rather than hidden by barrier-only comparison.",
        ),
        FigureSpec(
            stem="cdse_rs_to_wz_cell_mapping_literature_comparison",
            title="CdSe rock-salt$\\rightarrow$wurtzite",
            ylabel="Relative energy (meV/atom)",
            current=cdse,
            literature=cdse_lit,
            literature_barrier=2.4,
            literature_label=cdse_lit.label,
            comparison_note="VARNEB PBE vs literature PW91; each path coordinate is normalized independently.",
        ),
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/neb_literature_benchmarks"))
    parser.add_argument(
        "--bto",
        type=Path,
        default=Path("outputs/batio3_t_to_c_pbe100_dzp10au/bto_tetragonal_to_cubic_n7_static_audit_summary.json"),
    )
    parser.add_argument("--hfo2-tpo", type=Path, default=Path("tmp/plot_data/hfo2/tpo_vasp_summary.json"))
    parser.add_argument("--hfo2-pom", type=Path, default=Path("tmp/plot_data/hfo2/pom_vasp_summary.json"))
    parser.add_argument("--gan-tetragonal", type=Path, default=Path("tmp/plot_data/gan/b4_legacy_vasp_summary.json"))
    parser.add_argument("--gan-hexagonal", type=Path, default=Path("tmp/plot_data/gan/hex_vasp_summary.json"))
    parser.add_argument("--gan-b3", type=Path, default=Path("tmp/plot_data/gan/b3_vasp_summary.json"))
    parser.add_argument("--cdse-cell", type=Path, default=Path("tmp/plot_data/cdse/cell_vasp_summary.json"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    specs = build_specs(args)
    outputs: list[Path] = []
    for spec in specs:
        source_path = args.output_dir / f"{spec.stem}_source_data.csv"
        write_source_data(source_path, spec)
        outputs.append(source_path)
        outputs.extend(draw_figure(args.output_dir, spec))
    input_paths = {
        "bto": args.bto,
        "hfo2_tpo": args.hfo2_tpo,
        "hfo2_pom": args.hfo2_pom,
        "gan_tetragonal": args.gan_tetragonal,
        "gan_hexagonal": args.gan_hexagonal,
        "gan_b3": args.gan_b3,
        "cdse_cell": args.cdse_cell,
    }
    manifest = {
        "description": "Per-path VARNEB/literature benchmark figures",
        "inputs": {
            name: {"path": str(path), "sha256": _sha256(path)}
            for name, path in input_paths.items()
        },
        "figures": [
            {
                "stem": spec.stem,
                "current_barrier": spec.current.barrier,
                "literature_barrier": spec.literature_barrier,
                "unit": spec.current.unit,
                "comparison_note": spec.comparison_note,
            }
            for spec in specs
        ],
        "outputs": [str(path) for path in outputs],
    }
    (args.output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
