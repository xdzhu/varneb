"""Plot the 27-point raw-audited BTO Qy=0 restricted conditional sheet.

The contour interpolates only measured same-contract conditional relaxations.
It does not minimize over the omitted unstable Qy branch and is not a global
T-to-C free-energy surface or barrier.
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
from matplotlib.lines import Line2D
from scipy.interpolate import PchipInterpolator

from scripts.plot_bto_qy0_restricted_sheet import evidence


ROOT = Path(__file__).resolve().parents[1]
INK, RED, GOLD = "#263746", "#CA5940", "#E5A547"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> tuple[dict, list[dict], list[dict], float]:
    report = json.loads(path.read_text(encoding="utf-8"))
    old_qa, old_points, images, t_qz, _ = evidence()
    qz = np.asarray(report.get("grid_qz_sqrt_amu_A"), dtype=float)
    qx = np.asarray(report.get("grid_qx_sqrt_amu_A"), dtype=float)
    energy = np.asarray(report.get("energy_minus_C_meV_per_BTO"), dtype=float)
    if (report.get("status") != "BTO_Qy0_restricted_27_node_sheet_raw_audited"
            or report.get("n_reused_measured_nodes") != 12
            or report.get("n_new_raw_audited_nodes") != 15
            or report.get("n_total_measured_nodes") != 27
            or report.get("source_sha256", {}).get("measured_csv")
               != sha256(ROOT / "paper/VARNEB_CPC/figures/bto_qy0_restricted_even_mode_sheet_points_source_data.csv")
            or not report.get("calculator_id", "").startswith("abacus-bto-pbe100-dzp10au-4x4x4")
            or qz.shape != (9,) or qx.shape != (3,) or energy.shape != (9, 3)
            or not np.isfinite(energy).all()
            or not np.allclose(np.diff(qz), 0.15, atol=1e-12, rtol=0)
            or not np.allclose(np.diff(qx), 0.15, atol=1e-12, rtol=0)
            or not np.isclose(energy[0, 0], old_points[0]["observed_meV_per_BTO"],
                              atol=1e-8, rtol=0)
            or old_qa["prospective_holdout_nodes"] != 2):
        raise ValueError("BTO 9x3 measured sheet or prior endpoint/path contract invalid")
    return report, old_points, images, t_qz


def surface(report: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    qz = np.asarray(report["grid_qz_sqrt_amu_A"], dtype=float)
    qx = np.asarray(report["grid_qx_sqrt_amu_A"], dtype=float)
    values = np.asarray(report["energy_minus_C_meV_per_BTO"], dtype=float)
    zz = np.linspace(qz[0], qz[-1], 361)
    xx = np.linspace(qx[0], qx[-1], 151)
    along_x = PchipInterpolator(qx, values, axis=1)(xx)
    interpolated = PchipInterpolator(qz, along_x, axis=0)(zz)
    if (not np.isfinite(interpolated).all()
            or not np.allclose(PchipInterpolator(qz, values, axis=0)(qz), values,
                               atol=1e-10, rtol=0)):
        raise ValueError("restricted-sheet display interpolation failed measured-node check")
    return zz, xx, interpolated


def draw(report: dict, old_points: list[dict], images: list[dict], t_qz: float) -> plt.Figure:
    zz, xx, interpolated = surface(report)
    qz = np.asarray(report["grid_qz_sqrt_amu_A"], dtype=float)
    qx = np.asarray(report["grid_qx_sqrt_amu_A"], dtype=float)
    mesh_z, mesh_x = np.meshgrid(zz, xx, indexing="ij")
    grid_z, grid_x = np.meshgrid(qz, qx, indexing="ij")
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "font.size": 10.2, "axes.labelsize": 11.0,
        "xtick.labelsize": 9.8, "ytick.labelsize": 9.8,
        "axes.linewidth": 0.9, "axes.spines.top": True, "axes.spines.right": True,
    })
    fig = plt.figure(figsize=(8.5, 4.7), facecolor="white")
    ax_a = fig.add_axes((0.09, 0.28, 0.55, 0.58))
    ax_b = fig.add_axes((0.735, 0.28, 0.23, 0.58))
    cax = fig.add_axes((0.16, 0.105, 0.43, 0.028))
    fig.text(0.09, 0.925, "(a)", fontsize=13, fontweight="normal", color=INK)
    fig.text(0.735, 0.925, "(b)", fontsize=13, fontweight="normal", color=INK)
    levels = np.arange(-105, 5, 5)
    contour = ax_a.contourf(mesh_z, mesh_x, interpolated, levels=levels,
                            cmap="cividis_r", extend="both")
    ax_a.contour(mesh_z, mesh_x, interpolated, levels=np.arange(-90, 1, 15),
                 colors="white", linewidths=0.6, alpha=0.55)
    ax_a.scatter(grid_z.ravel(), grid_x.ravel(), s=33, marker="o",
                 facecolors="white", edgecolors=INK, linewidths=0.75, zorder=6)
    path_z = [row["q_parallel_sqrt_amu_A"] for row in images]
    path_x = [row["q_transverse_sqrt_amu_A"] for row in images]
    ax_a.plot(path_z, path_x, color="white", lw=3.0, zorder=7)
    ax_a.plot(path_z, path_x, color=RED, lw=1.55, zorder=8)
    ax_a.scatter([0], [0], s=85, marker="s", facecolor="white",
                 edgecolor=RED, linewidths=1.2, zorder=9)
    ax_a.scatter([t_qz], [0], s=165, marker="*", facecolor=GOLD,
                 edgecolor=RED, linewidths=0.9, zorder=9)
    ax_a.set_xlim(-0.015, max(1.225, t_qz + 0.012))
    ax_a.set_ylim(-0.012, 0.312)
    ax_a.set_xlabel(r"$Q_z$ ($\sqrt{\mathrm{amu}}$ Å)")
    ax_a.set_ylabel(r"$Q_x$ ($\sqrt{\mathrm{amu}}$ Å)")
    ax_a.tick_params(direction="in", top=True, right=True, length=4)
    colorbar = fig.colorbar(contour, cax=cax, orientation="horizontal")
    colorbar.ax.xaxis.set_label_position("top")
    colorbar.set_label(r"Restricted $E-E_C$ (meV/BTO)", fontsize=10, labelpad=4)
    colorbar.ax.tick_params(labelsize=9.1, direction="in", length=3)
    fig.legend(handles=[
        Line2D([], [], marker="o", linestyle="none", markerfacecolor="white",
               markeredgecolor=INK, markersize=6, label="27 DFT nodes"),
        Line2D([], [], color=RED, lw=1.5, label="VCNEB projection"),
    ], loc="center", bbox_to_anchor=(0.395, 0.929), ncol=2, fontsize=9.2,
        frameon=True, facecolor="white", edgecolor="#8B969C", framealpha=0.86)
    ax_b.axhspan(-2, 2, color="#DDEDEA", alpha=0.72, zorder=0)
    ax_b.axhline(0, color=INK, lw=0.9, zorder=1)
    colors = {0.0: "#315D85", 0.15: "#D38A40", 0.3: "#885A86"}
    for x, color in colors.items():
        rows = [row for row in report["new"] if np.isclose(row["q_sqrt_amu_A"][1], x)]
        ax_b.scatter([row["q_sqrt_amu_A"][0] for row in rows],
                     [row["DFT_minus_model_meV_per_BTO"] for row in rows],
                     s=44, marker="o", facecolor=color, edgecolor="white",
                     linewidths=0.6, zorder=3, label=rf"$Q_x={x:g}$")
    prior_blind = [row for row in old_points if row["class"] == "prospective"]
    ax_b.scatter([row["qz"] for row in prior_blind],
                 [row["observed_minus_predicted_meV_per_BTO"] for row in prior_blind],
                 s=68, marker="D", facecolor=GOLD, edgecolor=INK,
                 linewidths=0.75, zorder=4, label="prior blind")
    residual = [abs(row["DFT_minus_model_meV_per_BTO"]) for row in report["new"]]
    ax_b.set_ylim(-max(2.6, max(residual) + 0.5), max(2.6, max(residual) + 0.5))
    ax_b.set_xlim(-0.04, 1.24)
    ax_b.set_xlabel(r"$Q_z$ ($\sqrt{\mathrm{amu}}$ Å)")
    ax_b.set_ylabel("DFT − nine-node model (meV/BTO)")
    ax_b.tick_params(direction="in", top=True, right=True, length=4)
    ax_b.legend(loc="best", fontsize=8.0, frameon=True, facecolor="white",
                edgecolor="#8B969C", framealpha=0.86, borderpad=0.3)
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()
    outputs = {suffix: args.output_prefix.with_name(args.output_prefix.name + suffix)
               for suffix in (".pdf", ".svg", ".png", "_qa.json", "_source_data.csv")}
    if any(path.exists() for path in outputs.values()):
        raise FileExistsError("use a fresh restricted-sheet figure prefix")
    report, old_points, images, t_qz = load(args.report)
    fig = draw(report, old_points, images, t_qz)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".pdf", ".svg", ".png"):
        fig.savefig(outputs[suffix], dpi=350 if suffix == ".png" else None,
                    facecolor="white")
    plt.close(fig)
    with outputs["_source_data.csv"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("Qz_sqrt_amu_A", "Qx_sqrt_amu_A", "E_minus_C_meV_per_BTO",
                         "source"))
        grid = np.asarray(report["energy_minus_C_meV_per_BTO"])
        for i, z in enumerate(report["grid_qz_sqrt_amu_A"]):
            for j, x in enumerate(report["grid_qx_sqrt_amu_A"]):
                writer.writerow((z, x, grid[i, j], "audited_DFT"))
    qa = {
        "status": "BTO_Qy0_restricted_27_node_measured_sheet_figure",
        "claim_limit": report["claim_limit"],
        "n_measured_DFT_nodes": 27,
        "interpolation": "separable shape-preserving cubic Hermite within measured grid only",
        "new_node_model_max_abs_error_meV_per_BTO": report[
            "new_node_model_max_abs_error_meV_per_BTO"],
        "new_node_model_gate_pass": report["new_node_model_gate_pass"],
        "source_sha256": {
            "raw_audit": sha256(args.report),
            "source_data": sha256(outputs["_source_data.csv"]),
            "plotter": sha256(Path(__file__)),
        },
    }
    outputs["_qa.json"].write_text(json.dumps(qa, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": qa["status"], "n_measured_DFT_nodes": 27}))


if __name__ == "__main__":
    main()
