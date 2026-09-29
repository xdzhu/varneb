"""Draw the 162-measured-point GaN central frozen transverse enthalpy cut.

The display contour is shape-preserving interpolation of audited static DFT
samples, not a relaxed PES or a stationary-TS certificate. The right panel is
the independently audited whole-path 45.7-GPa enthalpy trace.
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
from scipy.interpolate import PchipInterpolator


ROOT = Path(__file__).resolve().parents[1]
INK, BLUE, ORANGE = "#283844", "#315E82", "#D47B42"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(report_path: Path, refinement_path: Path, trajectory: Path) -> tuple[dict, np.ndarray, np.ndarray]:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    refinement = json.loads(refinement_path.read_text(encoding="utf-8"))
    frames = read(trajectory, index=":")
    indices = report.get("central_image_indices")
    s = np.asarray(report.get("arc_fraction_s"), dtype=float)
    q = np.asarray(report.get("q_atom_A"), dtype=float)
    values = np.asarray(report.get("excess_enthalpy_meV_per_GaN"), dtype=float)
    if (report.get("status") != "GaN_600eV_central_18x9_raw_audited"
            or report.get("pressure_GPa") != 45.7
            or report.get("n_reused_measured_points") != 90
            or report.get("n_new_raw_audited_statics") != 72
            or report.get("prospective_gate_pass") is not True
            or report.get("source_sha256", {}).get("trajectory") != sha256(trajectory)
            or report.get("source_sha256", {}).get("dense_report")
               != sha256(ROOT / "benchmarks/numerical_integrity/gan_600eV_atomic_tube_dense_20260928.json")
            or indices != list(range(5, 23)) or len(frames) != 29
            or s.shape != (18,) or q.shape != (9,) or values.shape != (18, 9)
            or not np.all(np.diff(s) > 0) or not np.all(np.diff(q) > 0)
            or not np.isfinite(values).all()
            or not np.allclose(values[:, 4], 0, atol=1e-8, rtol=0)
            or not np.allclose(s, np.asarray(refinement["arc_fraction_s"])[indices],
                               atol=1e-12, rtol=0)):
        raise ValueError("18x9 GaN grid or its original 600-eV contract is not auditable")
    pressure = 45.7 * GPa
    enthalpy = np.asarray([frame.get_potential_energy() + pressure * frame.get_volume()
                           for frame in frames], dtype=float)
    if not np.allclose(enthalpy, refinement["path_enthalpy_eV_per_cell"],
                       atol=1e-8, rtol=0):
        raise ValueError("whole-path reference enthalpy changed")
    return report, np.asarray(refinement["arc_fraction_s"], dtype=float), enthalpy


def interpolated_surface(report: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    s = np.asarray(report["arc_fraction_s"], dtype=float)
    q = np.asarray(report["q_atom_A"], dtype=float)
    measured = np.asarray(report["excess_enthalpy_meV_per_GaN"], dtype=float)
    s_dense = np.linspace(s[0], s[-1], 401)
    q_dense = np.linspace(q[0], q[-1], 241)
    along_q = PchipInterpolator(q, measured, axis=1)(q_dense)
    surface = PchipInterpolator(s, along_q, axis=0)(s_dense)
    if (not np.isfinite(surface).all()
            or not np.allclose(PchipInterpolator(s, measured, axis=0)(s), measured,
                               atol=1e-10, rtol=0)
            or np.max(np.abs(PchipInterpolator(q_dense, surface.T, axis=0)([0.0]))) > 1e-8):
        raise ValueError("shape-preserving display interpolation failed node/baseline check")
    return s_dense, q_dense, surface


def draw(report: dict, full_arc: np.ndarray, full_h: np.ndarray) -> plt.Figure:
    s_dense, q_dense, z = interpolated_surface(report)
    s = np.asarray(report["arc_fraction_s"], dtype=float)
    q = np.asarray(report["q_atom_A"], dtype=float)
    path = (full_h - full_h[0]) * 500.0
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "font.size": 10.2, "axes.labelsize": 11.0,
        "xtick.labelsize": 9.6, "ytick.labelsize": 9.6,
        "axes.linewidth": 0.9, "axes.spines.top": True, "axes.spines.right": True,
    })
    fig = plt.figure(figsize=(8.4, 4.9), facecolor="white")
    ax_a = fig.add_axes((0.095, 0.31, 0.540, 0.57))
    ax_b = fig.add_axes((0.735, 0.31, 0.225, 0.57))
    cax = fig.add_axes((0.16, 0.13, 0.410, 0.028))
    fig.text(0.095, 0.93, "(a)", fontsize=13, fontweight="normal", color=INK)
    fig.text(0.735, 0.93, "(b)", fontsize=13, fontweight="normal", color=INK)
    ss, qq = np.meshgrid(s_dense, q_dense, indexing="ij")
    levels = np.r_[[-1.0, -0.5, 0.0], np.arange(2.0, 24.0, 2.0)]
    contour = ax_a.contourf(ss, qq, z, levels=levels, cmap="cividis", extend="both")
    ax_a.contour(ss, qq, z, levels=levels[4::2], colors=INK,
                 linewidths=0.45, alpha=0.3)
    grid_s, grid_q = np.meshgrid(s, q, indexing="ij")
    ax_a.scatter(grid_s.ravel(), grid_q.ravel(), s=9, marker="o",
                 facecolors="white", edgecolors=INK, linewidths=0.5, zorder=6)
    measured = np.asarray(report["excess_enthalpy_meV_per_GaN"], dtype=float)
    lower = np.argwhere(measured < -0.05)
    ax_a.scatter([s[i] for i, _ in lower], [q[j] for _, j in lower],
                 s=50, marker="D", facecolors=ORANGE, edgecolors="white",
                 linewidths=0.8, zorder=9)
    ax_a.plot(s, np.zeros_like(s), color="white", lw=2.8, zorder=7)
    ax_a.plot(s, np.zeros_like(s), color=INK, lw=1.5, zorder=8)
    ax_a.scatter(s[10], 0, s=135, marker="*", color=ORANGE,
                 edgecolors="white", linewidths=0.7, zorder=9)
    ax_a.set_xlim(s[0] - 0.004, s[-1] + 0.004)
    ax_a.set_ylim(-0.054, 0.054)
    ax_a.set_yticks([-0.05, -0.025, 0, 0.025, 0.05])
    ax_a.set_xlabel(r"VCNEB arc fraction $s$")
    ax_a.set_ylabel(r"Atomic transverse coordinate $q_{\perp}$ (Å)")
    ax_a.tick_params(direction="in", top=True, right=True, length=4)
    colorbar = fig.colorbar(contour, cax=cax, orientation="horizontal")
    colorbar.ax.xaxis.set_label_position("top")
    colorbar.set_label(r"$H(s,q_{\perp})-H(s,0)$ (meV/GaN)", fontsize=10, labelpad=4)
    colorbar.ax.tick_params(labelsize=9.2, direction="in", length=3)
    fig.legend(handles=[
        Line2D([], [], marker="o", linestyle="none", markerfacecolor="white",
               markeredgecolor=INK, markersize=5, label="162 DFT points"),
        Line2D([], [], color=INK, lw=1.5, label="VCNEB centerline"),
    ], loc="center", bbox_to_anchor=(0.4, 0.935), ncol=2, fontsize=9.2,
        frameon=True, facecolor="white", edgecolor="#89959D", framealpha=0.86)
    ax_b.axvspan(s[0], s[-1], color="#DCE8EE", alpha=0.73, zorder=0)
    ax_b.plot(full_arc, path, color=BLUE, lw=1.8, zorder=2)
    ax_b.scatter(full_arc, path, s=21, color=BLUE, edgecolors="white",
                 linewidths=0.5, zorder=3)
    ax_b.scatter(full_arc[15], path[15], s=155, marker="*", color=ORANGE,
                 edgecolors="white", linewidths=0.7, zorder=4)
    ax_b.set_xlim(0, 1)
    ax_b.set_ylim(min(-18, path.min() - 8), path.max() + 23)
    ax_b.set_xlabel(r"VCNEB arc fraction $s$")
    ax_b.set_ylabel(r"$H-H_{\mathrm{B4}}$ (meV/GaN)")
    ax_b.tick_params(direction="in", top=True, right=True, length=4)
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--refinement", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()
    outputs = {suffix: args.output_prefix.with_name(args.output_prefix.name + suffix)
               for suffix in (".pdf", ".svg", ".png", "_qa.json", "_source_data.csv")}
    if any(path.exists() for path in outputs.values()):
        raise FileExistsError("use a new GaN 18x9 figure prefix")
    report, full_arc, full_h = load(args.report, args.refinement, args.trajectory)
    fig = draw(report, full_arc, full_h)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".pdf", ".svg", ".png"):
        fig.savefig(outputs[suffix], dpi=350 if suffix == ".png" else None,
                    facecolor="white")
    plt.close(fig)
    with outputs["_source_data.csv"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("image_index", "arc_fraction_s", "q_atom_A",
                         "excess_enthalpy_meV_per_GaN", "kind"))
        for row, index in enumerate(report["central_image_indices"]):
            for column, q in enumerate(report["q_atom_A"]):
                writer.writerow((index, report["arc_fraction_s"][row], q,
                                 report["excess_enthalpy_meV_per_GaN"][row][column],
                                 "audited_DFT"))
    qa = {
        "status": "GaN_600eV_central_18x9_frozen_transverse_figure",
        "claim_limit": report["claim_limit"],
        "n_audited_DFT_points": 162,
        "n_measured_lower_than_centerline_by_0p05_meV": int(np.count_nonzero(
            np.asarray(report["excess_enthalpy_meV_per_GaN"]) < -0.05)),
        "minimum_measured_excess_meV_per_GaN": float(np.min(
            report["excess_enthalpy_meV_per_GaN"])),
        "interpolation": "separable shape-preserving cubic Hermite; only 162 marked nodes are DFT",
        "maximum_prospective_error_meV_per_GaN": report[
            "prospective_max_abs_error_meV_per_GaN"],
        "source_sha256": {
            "audit": sha256(args.report), "refinement": sha256(args.refinement),
            "trajectory": sha256(args.trajectory),
            "source_data": sha256(outputs["_source_data.csv"]),
            "plotter": sha256(Path(__file__)),
        },
    }
    outputs["_qa.json"].write_text(json.dumps(qa, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": qa["status"], "n_audited_DFT_points": 162}))


if __name__ == "__main__":
    main()
