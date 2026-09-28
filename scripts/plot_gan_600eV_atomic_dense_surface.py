"""Draw the validated central GaN 600-eV path-adapted enthalpy cut.

Figure contract: the GaN B4→B1 barrier has a measured transverse enthalpy
response around the variable-cell path.  Panel (a) is a two-coordinate,
atomic-only frozen central cut with all 90 measured grid locations; panel (b)
places its domain on the full 29-image barrier profile.  The only plotted
interpolant is linear in path arc fraction and quadratic in atomic q, and is
refused unless both documented holdout gates pass. The along-s gate was in
the submitted manifest; the q gate was fixed before inspecting new outputs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from ase.units import GPa
from matplotlib.lines import Line2D

from scripts.prepare_gan_600eV_atomic_tube_dense import CENTRAL_IMAGES, TRANSVERSE_Q_A


ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = ROOT / "benchmarks/numerical_integrity"
DEFAULT_REPORT = BENCHMARKS / "gan_600eV_atomic_tube_dense_20260928.json"
DEFAULT_TRAJECTORY = (
    ROOT / "paper/VARNEB_CPC/evidence/gan_45p7_final_chains_20260927"
    / "gan_vasp_45p7_final_chain.traj"
)
DEFAULT_PREFIX = ROOT / "paper/VARNEB_CPC/figures/gan_600eV_atomic_dense_surface"
INK = "#263440"
BLUE = "#285D80"
ORANGE = "#C47D45"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "font.size": 10.0, "axes.labelsize": 10.6,
    "xtick.labelsize": 9.4, "ytick.labelsize": 9.4,
    "axes.linewidth": 0.9,
    "axes.spines.top": True, "axes.spines.right": True,
})


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(report_path: Path, trajectory: Path) -> tuple[dict, list]:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if (report.get("status") != "GaN_600eV_central_dense_atomic_tube_raw_audited"
            or report.get("central_contour_gate_pass") is not True
            or report.get("s_holdout_gate_pass") is not True
            or report.get("q_halfstep_gate_pass") is not True
            or report.get("pressure_GPa") != 45.7
            or report.get("formula_units_per_cell") != 2
            or report.get("central_image_indices") != list(CENTRAL_IMAGES)
            or report.get("q_atom_A") != list(TRANSVERSE_Q_A)
            or report.get("n_new_raw_audited_statics") != 44
            or report.get("n_reused_offpath_statics") != 28
            or report["source_sha256"].get("trajectory") != sha256(trajectory)):
        raise ValueError("GaN dense-grid raw audit or prospective contour gate failed")
    frames = read(trajectory, index=":")
    if len(frames) != 29:
        raise ValueError("GaN production trajectory is not 29 total images")
    path_h = np.array([
        frame.get_potential_energy() + 45.7 * GPa * frame.get_volume()
        for frame in frames
    ])
    if not np.allclose(path_h[list(CENTRAL_IMAGES)],
                       report["path_enthalpy_eV_per_cell"], atol=1e-8, rtol=0):
        raise ValueError("dense-grid q=0 enthalpy differs from the production chain")
    matrix = np.asarray(report["excess_enthalpy_meV_per_GaN"], dtype=float)
    if matrix.shape != (18, 5) or not np.isfinite(matrix).all() or not np.allclose(matrix[:, 2], 0):
        raise ValueError("incomplete or nonfinite central grid")
    return report, frames


def model_surface(report: dict, full_h: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    arc = np.asarray(report["arc_fraction_s"], dtype=float)
    matrix = np.asarray(report["excess_enthalpy_meV_per_GaN"], dtype=float)
    width = 0.05
    slope = (matrix[:, 4] - matrix[:, 0]) / (2 * width)
    curvature_half = (matrix[:, 4] + matrix[:, 0]) / (2 * width**2)
    s_dense = np.linspace(arc[0], arc[-1], 401)
    q_dense = np.linspace(-width, width, 201)
    s_mesh, q_mesh = np.meshgrid(s_dense, q_dense)
    h0 = np.interp(s_dense, arc, np.asarray(report["path_enthalpy_eV_per_cell"]))
    slope_dense = np.interp(s_dense, arc, slope)
    curvature_dense = np.interp(s_dense, arc, curvature_half)
    surface = ((h0[None, :] - full_h[0]) * 500
               + slope_dense[None, :] * q_mesh
               + curvature_dense[None, :] * q_mesh**2)
    return s_mesh, q_mesh, surface


def source_rows(report: dict, frames: list, trajectory_sha: str) -> list[dict]:
    all_cases = {}
    first = json.loads(
        (BENCHMARKS / "gan_600eV_atomic_tube_central_20260928.json").read_text()
    )
    refinement = json.loads(
        (BENCHMARKS / "gan_600eV_atomic_tube_refinement_20260928.json").read_text()
    )
    if (sha256(BENCHMARKS / "gan_600eV_atomic_tube_central_20260928.json")
            != report["source_sha256"]["first_audit"]
            or sha256(BENCHMARKS / "gan_600eV_atomic_tube_refinement_20260928.json")
            != report["source_sha256"]["refinement_audit"]):
        raise ValueError("prior audited raw cases changed")
    for origin, cases in (("prior_static", first["cases"] + refinement["new_cases"]),
                          ("new_static", report["new_cases"])):
        for case in cases:
            key = (int(case["image_index"]), round(float(case["q_atom_A"]), 6))
            if key in all_cases:
                raise ValueError("duplicate raw surface point")
            all_cases[key] = (origin, case["outcar_sha256"])
    matrix = np.asarray(report["excess_enthalpy_meV_per_GaN"])
    arc = report["arc_fraction_s"]
    full_h = np.array([
        frame.get_potential_energy() + 45.7 * GPa * frame.get_volume()
        for frame in frames
    ])
    rows = []
    for row_index, index in enumerate(CENTRAL_IMAGES):
        for col_index, q in enumerate(TRANSVERSE_Q_A):
            source, digest = ("path_trajectory", trajectory_sha) if q == 0.0 else all_cases[(index, q)]
            rows.append({
                "image_index": index,
                "arc_fraction_s": arc[row_index],
                "q_atom_A": q,
                "measured_H_minus_B4_meV_per_GaN": float((full_h[index] - full_h[0]) * 500
                                                       + matrix[row_index, col_index]),
                "offpath_excess_meV_per_GaN": float(matrix[row_index, col_index]),
                "source_kind": source,
                "source_sha256": digest,
                "prospective_s_holdout": index in (10, 16),
            })
    if len(rows) != 90:
        raise ValueError("not all 90 grid source rows were assembled")
    return rows


def draw(report: dict, frames: list) -> plt.Figure:
    path_h = np.asarray([
        frame.get_potential_energy() + 45.7 * GPa * frame.get_volume()
        for frame in frames
    ])
    path_mev = (path_h - path_h[0]) * 500
    full_arc = np.asarray(json.loads(
        (BENCHMARKS / "gan_600eV_atomic_tube_refinement_20260928.json").read_text()
    )["arc_fraction_s"], dtype=float)
    if not np.allclose(full_arc[list(CENTRAL_IMAGES)], report["arc_fraction_s"], atol=1e-12):
        raise ValueError("full-path arc fraction changed")
    s_mesh, q_mesh, surface = model_surface(report, path_h)
    fig = plt.figure(figsize=(7.2, 4.4), facecolor="white")
    ax_a = fig.add_axes((0.105, 0.315, 0.525, 0.555))
    ax_b = fig.add_axes((0.710, 0.315, 0.250, 0.555))
    cax = fig.add_axes((0.155, 0.125, 0.435, 0.025))
    fig.text(0.105, 0.935, "(a)", fontsize=13, fontweight="normal", color=INK)
    fig.text(0.710, 0.935, "(b)", fontsize=13, fontweight="normal", color=INK)
    levels = np.linspace(np.floor(surface.min() / 25) * 25,
                         np.ceil(surface.max() / 25) * 25, 15)
    contour = ax_a.contourf(s_mesh, q_mesh, surface, levels=levels,
                            cmap="cividis", extend="neither")
    ax_a.contour(s_mesh, q_mesh, surface, levels=levels[::2],
                 colors=INK, alpha=0.28, linewidths=0.55)
    s = np.asarray(report["arc_fraction_s"])
    q = np.asarray(TRANSVERSE_Q_A)
    grid_s, grid_q = np.meshgrid(s, q, indexing="ij")
    ax_a.scatter(grid_s.ravel(), grid_q.ravel(), s=14, marker="o",
                 facecolors="white", edgecolors=INK, linewidths=0.55, zorder=4)
    ax_a.plot(s, np.zeros_like(s), color="white", lw=2.8, zorder=5)
    ax_a.plot(s, np.zeros_like(s), color=INK, lw=1.45, zorder=6)
    peak_local = 15 - CENTRAL_IMAGES[0]
    ax_a.scatter(s[peak_local], 0, s=135, marker="*", color=ORANGE,
                 edgecolor="white", linewidth=0.7, zorder=7)
    ax_a.set_xlim(s[0] - 0.004, s[-1] + 0.004)
    ax_a.set_ylim(-0.054, 0.054)
    ax_a.set_yticks([-0.05, -0.025, 0, 0.025, 0.05])
    ax_a.set_xlabel(r"VCNEB arc fraction $s$")
    ax_a.set_ylabel(r"Atomic transverse coordinate $q_{\perp}$ (Å)")
    ax_a.tick_params(direction="in", top=True, right=True, length=4, pad=4)
    cbar = fig.colorbar(contour, cax=cax, orientation="horizontal")
    cbar.ax.xaxis.set_label_position("top")
    cbar.set_label(r"$H-H_{\mathrm{B4}}$ (meV/GaN)", fontsize=9.6, labelpad=4)
    cbar.ax.tick_params(labelsize=8.8, direction="in", length=3)
    fig.legend(handles=[
        Line2D([], [], marker="o", linestyle="none", markerfacecolor="white",
               markeredgecolor=INK, markersize=5, label="DFT grid (90)"),
        Line2D([], [], color=INK, lw=1.5, label="VCNEB centerline"),
    ], loc="center", bbox_to_anchor=(0.370, 0.928), ncol=2,
        fontsize=8.6, frameon=True, facecolor="white",
        edgecolor="#8B969C", framealpha=0.84, borderpad=0.35)

    ax_b.axvspan(s[0], s[-1], color="#DCE8EE", alpha=0.75, zorder=0)
    ax_b.plot(full_arc, path_mev, color=BLUE, lw=1.65, zorder=2)
    ax_b.scatter(full_arc, path_mev, s=19, color=BLUE,
                 edgecolors="white", linewidths=0.4, zorder=3)
    ax_b.scatter(full_arc[15], path_mev[15], s=155, marker="*",
                 color=ORANGE, edgecolors="white", linewidths=0.7, zorder=4)
    ax_b.set_xlim(0, 1)
    ax_b.set_ylim(min(-18, path_mev.min() - 8), path_mev.max() + 23)
    ax_b.set_xlabel(r"VCNEB arc fraction $s$")
    ax_b.set_ylabel(r"$H-H_{\mathrm{B4}}$ (meV/GaN)")
    ax_b.tick_params(direction="in", top=True, right=True, length=4, pad=4)
    fig.legend(handles=[
        Line2D([], [], color=BLUE, lw=1.7, marker="o", markersize=4,
               label="29-image path"),
        Line2D([], [], color="#DCE8EE", lw=7, label="2D domain"),
    ], loc="center", bbox_to_anchor=(0.835, 0.160), fontsize=8.4,
        frameon=True, facecolor="white",
        edgecolor="#8B969C", framealpha=0.84, borderpad=0.35)
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--trajectory", type=Path, default=DEFAULT_TRAJECTORY)
    parser.add_argument("--output-prefix", type=Path, default=DEFAULT_PREFIX)
    args = parser.parse_args()
    prefix = args.output_prefix
    paths = {
        suffix: prefix.with_name(prefix.name + suffix)
        for suffix in ("_source_data.csv", "_qa.json", ".pdf", ".svg", ".tiff", ".png")
    }
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("dense GaN figure output exists; use a fresh prefix")
    report, frames = load(args.report, args.trajectory)
    rows = source_rows(report, frames, sha256(args.trajectory))
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with paths["_source_data.csv"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    fig = draw(report, frames)
    fig.savefig(paths[".pdf"])
    fig.savefig(paths[".svg"])
    svg = paths[".svg"]
    svg.write_text(
        "\n".join(line.rstrip() for line in svg.read_text(encoding="utf-8").splitlines())
        + "\n", encoding="utf-8"
    )
    fig.savefig(paths[".tiff"], dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(paths[".png"], dpi=200)
    plt.close(fig)
    qa = {
        "status": "GaN_600eV_central_dense_atomic_enthalpy_contour_holdout_validated",
        "core_conclusion": "The 600-eV GaN VCNEB barrier has a measured atomic-transverse enthalpy response across its central path segment.",
        "archetype": "quantitative grid; central 2D contour hero plus full-path context",
        "backend": "Python/matplotlib",
        "final_size_mm": [182.88, 111.76],
        "panel_map": {
            "a": "Measured central 18-by-5 frozen atomic-only grid with a prospectively validated linear-s/quadratic-q enthalpy interpolant",
            "b": "Complete 29-image variable-cell enthalpy barrier and exact central chart domain",
        },
        "electronic_contract": "Original VASP 600 eV/PBE/Ga_d+N/Gamma8x8x6, unchanged at every point",
        "n_path_centers": 18,
        "n_offpath_statics": 72,
        "maximum_prospective_s_holdout_error_meV_per_GaN": report[
            "maximum_prospective_s_holdout_error_meV_per_GaN"],
        "maximum_q_halfstep_error_meV_per_GaN": report["maximum_q_halfstep_error_meV_per_GaN"],
        "interpolation": "Linear in s for path energy and outer-pair odd/even coefficients; quadratic in q between ±0.05 Å only",
        "claim_limit": report["claim_limit"],
        "TS_certified": False,
        "whole_path_2D_surface_certified": False,
        "source_sha256": {
            "raw_audit": sha256(args.report),
            "trajectory": sha256(args.trajectory),
            "source_csv": sha256(paths["_source_data.csv"]),
            "plotter": sha256(Path(__file__)),
        },
        "reviewer_risk": "The atomic transverse direction has endpoint seams; same-600-eV strain energy/stress inconsistency prevents strict TS certification.",
        "image_integrity": "All 90 measured locations are marked; interpolation is confined to the central sampled rectangle.",
        "raster_export_policy": "PDF/SVG/PNG tracked; 600-dpi TIFF generated locally and repository-ignored.",
    }
    paths["_qa.json"].write_text(json.dumps(qa, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": qa["status"], "source_rows": len(rows)}))


if __name__ == "__main__":
    main()
