"""Plot the 289-point GaN local frozen joint-coordinate enthalpy cut.

The two-axis contour is only a measured-and-interpolated local slice at the
original 600-eV VASP/45.7-GPa contract. Panel (b) prospectively tests the
previous 9x9 interpolant at 208 independently calculated midpoint nodes.
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
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np
from scipy.interpolate import CloughTocher2DInterpolator

from scripts.plot_gan_600eV_local_joint_cut import source_rows
from scripts.prepare_gan_600eV_ts_2d_refinement import coordinate_key


INK = "#263440"
BLUE = "#285D80"
ORANGE = "#D47B42"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_rows(pilot_path: Path, refined5_path: Path, dense9_path: Path,
              dense17_path: Path) -> tuple[list[dict], dict]:
    pilot = json.loads(pilot_path.read_text(encoding="utf-8"))
    refined5 = json.loads(refined5_path.read_text(encoding="utf-8"))
    dense9 = json.loads(dense9_path.read_text(encoding="utf-8"))
    dense17 = json.loads(dense17_path.read_text(encoding="utf-8"))
    if (dense17.get("status") != "GaN_600eV_local_joint_17x17_refinement_raw_audited"
            or dense17.get("pressure_GPa") != 45.7
            or dense17.get("n_total_DFT_points") != 289
            or dense17.get("n_new_DFT_points") != 208
            or len(dense17.get("cases", [])) != 208
            or dense17.get("source_sha256", {}).get("pilot_audit") != sha256(pilot_path)
            or dense17.get("source_sha256", {}).get("refinement5_audit") != sha256(refined5_path)
            or dense17.get("source_sha256", {}).get("dense9_audit") != sha256(dense9_path)):
        raise ValueError("17x17 raw audit is not tied to the preceding 81-point cut")
    rows = source_rows(pilot, sha256(pilot_path), refined5, sha256(refined5_path),
                       dense9, sha256(dense9_path))
    fields = list(rows[0])
    for case in dense17["cases"]:
        rows.append({
            "case": case["case"], "kind": "dense17_midpoint",
            "q_u_A": float(case["q_u_A"]), "q_v_A": float(case["q_v_A"]),
            "measured_delta_H_meV_per_GaN": float(case["delta_enthalpy_meV_per_GaN"]),
            "quadratic_model_delta_H_meV_per_GaN": "",
            "residual_meV_per_GaN": float(case["prospective_error_meV_per_GaN"]),
            "raw_OUTCAR_sha256": case["outcar_sha256"],
            "source_report_sha256": sha256(dense17_path),
        })
    u_values = sorted({r["q_u_A"] for r in rows})
    v_values = sorted({r["q_v_A"] for r in rows})
    errors = np.asarray([r["prospective_error_meV_per_GaN"]
                         for r in dense17["cases"]], dtype=float)
    if (len(rows) != 289 or len({coordinate_key(r["q_u_A"], r["q_v_A"])
                                     for r in rows}) != 289
            or len(u_values) != 17 or len(v_values) != 17
            or not np.allclose(u_values, np.linspace(-0.02, 0.02, 17), atol=1e-10, rtol=0)
            or not np.allclose(v_values, np.linspace(-0.0125, 0.0125, 17), atol=1e-10, rtol=0)
            or not np.isfinite([r["measured_delta_H_meV_per_GaN"] for r in rows]).all()
            or not np.isfinite(errors).all()
            or not np.isclose(np.max(np.abs(errors)),
                              dense17["prior_9x9_prospective_max_abs_error_meV_per_GaN"],
                              atol=1e-9, rtol=0)
            or not np.isclose(np.sqrt(np.mean(errors**2)),
                              dense17["prior_9x9_prospective_rms_error_meV_per_GaN"],
                              atol=1e-9, rtol=0)
            or any(set(row) != set(fields) for row in rows)):
        raise ValueError("17x17 source rows do not form a complete measured grid")
    rows.sort(key=lambda r: (r["q_v_A"], r["q_u_A"]))
    return rows, dense17


def surface(rows: list[dict]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xy = np.asarray([(r["q_u_A"], r["q_v_A"]) for r in rows], dtype=float)
    energy = np.asarray([r["measured_delta_H_meV_per_GaN"] for r in rows], dtype=float)
    u = np.linspace(-0.02, 0.02, 401)
    v = np.linspace(-0.0125, 0.0125, 251)
    uu, vv = np.meshgrid(u, v)
    interpolation = CloughTocher2DInterpolator(xy, energy)
    z = np.asarray(interpolation(uu, vv), dtype=float)
    if (not np.isfinite(z).all()
            or not np.allclose(interpolation(xy), energy, atol=1e-9, rtol=0)):
        raise ValueError("17x17 display interpolation fails measured-node test")
    return uu, vv, z


def draw(rows: list[dict], dense17: dict) -> plt.Figure:
    uu, vv, z = surface(rows)
    levels = [-0.50, -0.40, -0.30, -0.20, -0.10, -0.05, 0.0,
              0.025, 0.050, 0.075, 0.10, 0.15]
    if float(np.min(z)) < levels[0] or float(np.max(z)) > levels[-1]:
        raise ValueError("measured cut exceeds the fixed cross-figure color range")
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "font.size": 10.2, "axes.labelsize": 10.5,
        "xtick.labelsize": 9.2, "ytick.labelsize": 9.2,
        "axes.linewidth": 0.9, "axes.spines.top": True, "axes.spines.right": True,
    })
    fig = plt.figure(figsize=(7.2, 4.1), facecolor="white")
    ax_a = fig.add_axes((0.105, 0.32, 0.435, 0.55))
    ax_b = fig.add_axes((0.650, 0.32, 0.295, 0.55))
    cax = fig.add_axes((0.16, 0.095, 0.32, 0.027))
    fig.text(0.105, 0.928, "(a)", fontsize=13, fontweight="normal", color=INK)
    fig.text(0.650, 0.928, "(b)", fontsize=13, fontweight="normal", color=INK)
    contour = ax_a.contourf(
        uu, vv, z, levels=levels, cmap="PuOr_r",
        norm=TwoSlopeNorm(vmin=-0.50, vcenter=0, vmax=0.15), extend="neither",
    )
    ax_a.contour(uu, vv, z, levels=[0], colors=[INK], linewidths=0.9,
                 linestyles="--")
    ax_a.scatter([r["q_u_A"] for r in rows], [r["q_v_A"] for r in rows],
                 s=12, marker="o", facecolors="white", edgecolors=INK,
                 linewidths=0.48, alpha=0.94, zorder=5)
    ax_a.scatter([0], [0], s=150, marker="*", color=BLUE,
                 edgecolors="white", linewidths=0.75, zorder=7)
    ax_a.set_xlim(-0.0225, 0.0225)
    ax_a.set_ylim(-0.0144, 0.0144)
    ax_a.set_xticks([-0.02, -0.01, 0, 0.01, 0.02])
    ax_a.set_yticks([-0.0125, 0, 0.0125])
    ax_a.set_xlabel(r"Unstable joint coordinate $q_u$ (Å)")
    ax_a.set_ylabel(r"Stable joint coordinate $q_v$ (Å)")
    ax_a.tick_params(direction="in", top=True, right=True, length=4, pad=4)
    colorbar = fig.colorbar(contour, cax=cax, orientation="horizontal",
                            ticks=[-0.4, -0.2, 0, 0.1])
    colorbar.ax.xaxis.set_label_position("top")
    colorbar.set_label(r"$\Delta H$ (meV/GaN)", fontsize=9.4, labelpad=4)
    colorbar.ax.tick_params(labelsize=8.7, direction="in", length=3)

    new = dense17["cases"]
    predicted = np.asarray([r["prior_9x9_interpolation_meV_per_GaN"] for r in new])
    measured = np.asarray([r["delta_enthalpy_meV_per_GaN"] for r in new])
    limits = (-0.57, 0.20)
    diagonal = np.array(limits)
    band = float(dense17["prior_9x9_prospective_gate"]["max_abs_limit_meV_per_GaN"])
    ax_b.fill_between(diagonal, diagonal - band, diagonal + band,
                      color="#EAEFF2", alpha=0.85, zorder=0)
    ax_b.plot(diagonal, diagonal, color="#64727C", lw=1.2, zorder=1)
    ax_b.scatter(predicted, measured, s=20, marker="o", color=BLUE,
                 edgecolor="white", linewidth=0.35, alpha=0.85, zorder=3)
    ax_b.set_xlim(*limits)
    ax_b.set_ylim(*limits)
    ax_b.set_xticks([-0.5, -0.3, -0.1, 0.1])
    ax_b.set_yticks([-0.5, -0.3, -0.1, 0.1])
    ax_b.set_xlabel(r"Prior 9×9 prediction $\Delta H$ (meV/GaN)")
    ax_b.set_ylabel(r"New DFT $\Delta H$ (meV/GaN)")
    ax_b.tick_params(direction="in", top=True, right=True, length=4, pad=4)
    ax_b.legend(handles=[
        Line2D([], [], marker="o", color="none", markerfacecolor=BLUE,
               markeredgecolor="white", markersize=6, label="208 new DFT statics"),
        Line2D([], [], color="#C8D2D9", lw=8, alpha=0.7,
               label=r"Prospective ±0.02 meV/GaN"),
    ], loc="upper left", fontsize=8.2, frameon=True, facecolor="white",
        edgecolor="#8B969C", framealpha=0.86, borderpad=0.42,
        handlelength=1.2, handletextpad=0.45)
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("pilot-audit", "refinement5-audit", "dense9-audit",
                 "dense17-audit", "output-prefix"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    prefix = args.output_prefix
    suffixes = ("_source_data.csv", ".pdf", ".svg", ".tiff", ".png", "_qa.json")
    paths = {suffix: prefix.with_name(prefix.name + suffix) for suffix in suffixes}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("17x17 figure output exists; use a fresh prefix")
    rows, dense17 = load_rows(args.pilot_audit, args.refinement5_audit,
                              args.dense9_audit, args.dense17_audit)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with paths["_source_data.csv"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    fig = draw(rows, dense17)
    fig.savefig(paths[".pdf"])
    fig.savefig(paths[".svg"])
    fig.savefig(paths[".tiff"], dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(paths[".png"], dpi=250)
    plt.close(fig)
    qa = {
        "status": "GaN_600eV_local_joint_17x17_DFT_refinement_visualized",
        "core_conclusion": "The local frozen joint-coordinate cut is resolved by 289 measured DFT nodes; 208 new nodes independently test the previous 9x9 interpolation.",
        "archetype": "quantitative 289-point contour plus prospective midpoint validation",
        "backend": "Python/matplotlib", "final_size_mm": [182.88, 104.14],
        "pressure_GPa": 45.7,
        "electronic_contract": "Original VASP 600 eV/PBE/Ga_d+N/Γ8×8×6",
        "n_measured_DFT_points": len(rows),
        "n_new_midpoint_DFT_points": 208,
        "prospective_max_abs_error_meV_per_GaN": dense17["prior_9x9_prospective_max_abs_error_meV_per_GaN"],
        "prospective_rms_error_meV_per_GaN": dense17["prior_9x9_prospective_rms_error_meV_per_GaN"],
        "prospective_gate": dense17["prior_9x9_prospective_gate"],
        "interpolation": "Clough-Tocher within the measured 17x17 rectangle; pixels are not DFT nodes",
        "claim_limit": dense17["claim_limit"],
        "source_sha256": {
            "pilot_audit": sha256(args.pilot_audit),
            "refinement5_audit": sha256(args.refinement5_audit),
            "dense9_audit": sha256(args.dense9_audit),
            "dense17_audit": sha256(args.dense17_audit),
            "source_data": sha256(paths["_source_data.csv"]),
            "plotter": sha256(Path(__file__)),
        },
    }
    paths["_qa.json"].write_text(json.dumps(qa, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: qa[key] for key in (
        "status", "n_measured_DFT_points", "prospective_max_abs_error_meV_per_GaN",
        "prospective_gate",
    )}))


if __name__ == "__main__":
    main()
