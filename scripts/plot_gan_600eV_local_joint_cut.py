"""Plot the raw-audited local GaN 600-eV joint-mode enthalpy cut.

Figure contract: within the predeclared small coordinate rectangle, one
coupled atom–strain direction lowers enthalpy and its orthogonal direction
raises it. A local quadratic model passes four axial half-step holdouts, but
the center is not a certified full-variable-cell TS and this is not a
whole-path or finite-temperature free-energy surface. With the optional full
5x5 audit, panel (a) interpolates the 25 measured DFT enthalpies inside their
rectangle; otherwise it shows the original quadratic model. Panel (b)
validates the original quadratic prediction against the measured statics.
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
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from scipy.interpolate import CloughTocher2DInterpolator

from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = ROOT / "benchmarks/numerical_integrity/gan_600eV_ts_2d_pilot_20260928.json"
DEFAULT_HESSIAN_ROOT = ROOT / "benchmarks/numerical_integrity/gan_600eV_ts_hessian_0p02_20260928"
DEFAULT_PREFIX = ROOT / "paper/VARNEB_CPC/figures/gan_600eV_local_joint_cut"
INK = "#263440"
BLUE = "#285D80"
ORANGE = "#C47D45"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "font.size": 10.0, "axes.labelsize": 10.6,
    "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
    "axes.linewidth": 0.9,
    "axes.spines.top": True, "axes.spines.right": True,
})


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(report_path: Path, hessian_root: Path) -> dict:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    hessian_audit = hessian_root / "audit.json"
    hessian_npz = hessian_root / "joint_hessian.npz"
    gates = report.get("predeclared_gate", {})
    if (report.get("status") != "GaN_600eV_local_2D_frozen_enthalpy_pilot_raw_audited"
            or report.get("pressure_GPa") != 45.7
            or report.get("formula_units_per_cell") != 2
            or report["source_sha256"].get("hessian_audit") != sha256(hessian_audit)
            or report["source_sha256"].get("hessian_npz") != sha256(hessian_npz)
            or gates.get("negative_u_positive_v_at_both_scales") is not True
            or gates.get("holdout_model_error_below_0p2_meV_per_GaN") is not True
            or len(report.get("cases", [])) != 12
            or report["quadratic_model_max_axial_holdout_error_meV_per_GaN"] > 0.2
            or any(case.get("input_sha256", {}).get(name) != digest
                   for case in report["cases"]
                   for name, digest in PRODUCTION_INPUT_SHA256.items())):
        raise ValueError("local GaN cut lacks same-600-eV raw audit or holdout gate")
    coordinates = {(float(case["q_u_A"]), float(case["q_v_A"]))
                   for case in report["cases"]}
    if len(coordinates) != 12 or not np.all(np.isfinite(
            [*report["center_gradient_uv_eV_per_A"],
             *report["hessian_eigenvalues_uv_eV_per_A2"]])):
        raise ValueError("local joint-coordinate grid changed")
    return report


def model_meV_per_GaN(report: dict, q_u: np.ndarray, q_v: np.ndarray) -> np.ndarray:
    gradient = np.asarray(report["center_gradient_uv_eV_per_A"], dtype=float)
    eig = np.asarray(report["hessian_eigenvalues_uv_eV_per_A2"], dtype=float)
    return 500 * (gradient[0] * q_u + gradient[1] * q_v
                  + 0.5 * (eig[0] * q_u**2 + eig[1] * q_v**2))


def source_rows(report: dict, report_sha: str, refinement: dict | None = None,
                refinement_sha: str | None = None, dense9: dict | None = None,
                dense9_sha: str | None = None) -> list[dict]:
    rows = [{
        "case": "center", "kind": "center", "q_u_A": 0.0, "q_v_A": 0.0,
        "measured_delta_H_meV_per_GaN": 0.0,
        "quadratic_model_delta_H_meV_per_GaN": 0.0,
        "residual_meV_per_GaN": 0.0,
        "raw_OUTCAR_sha256": report["source_sha256"]["center_OUTCAR"],
        "source_report_sha256": report_sha,
    }]
    for case in report["cases"]:
        q_u, q_v = float(case["q_u_A"]), float(case["q_v_A"])
        predicted = float(model_meV_per_GaN(report, q_u, q_v))
        measured = float(case["delta_enthalpy_meV_per_GaN"])
        residual = measured - predicted
        if not np.isclose(residual, case["quadratic_model_residual_meV_per_GaN"],
                          atol=1e-8, rtol=0):
            raise ValueError(f"model and raw audit disagree at {case['case']}")
        rows.append({
            "case": case["case"],
            "kind": "axial_holdout" if case["case"].startswith("hold_") else "grid",
            "q_u_A": q_u, "q_v_A": q_v,
            "measured_delta_H_meV_per_GaN": measured,
            "quadratic_model_delta_H_meV_per_GaN": predicted,
            "residual_meV_per_GaN": residual,
            "raw_OUTCAR_sha256": case["outcar_sha256"],
            "source_report_sha256": report_sha,
        })
    if ([row["kind"] for row in rows].count("grid") != 8
            or [row["kind"] for row in rows].count("axial_holdout") != 4):
        raise ValueError("local grid or independent axial holdouts missing")
    if refinement is not None:
        if (refinement.get("status") != "GaN_600eV_local_joint_5x5_refinement_raw_audited"
                or refinement.get("n_total_measured_coordinates") != 25
                or refinement.get("source_sha256", {}).get("original_pilot_audit") != report_sha
                or len(refinement.get("cases", [])) != 12
                or not refinement_sha):
            raise ValueError("local 5x5 refinement is not tied to the pilot")
        for case in refinement["cases"]:
            q_u, q_v = float(case["q_u_A"]), float(case["q_v_A"])
            predicted = float(model_meV_per_GaN(report, q_u, q_v))
            measured = float(case["delta_enthalpy_meV_per_GaN"])
            residual = measured - predicted
            if not np.isclose(residual, case["quadratic_model_residual_meV_per_GaN"],
                              atol=1e-8, rtol=0):
                raise ValueError(f"refinement/model discrepancy at {case['case']}")
            rows.append({
                "case": case["case"], "kind": "offaxis_refinement",
                "q_u_A": q_u, "q_v_A": q_v,
                "measured_delta_H_meV_per_GaN": measured,
                "quadratic_model_delta_H_meV_per_GaN": predicted,
                "residual_meV_per_GaN": residual,
                "raw_OUTCAR_sha256": case["outcar_sha256"],
                "source_report_sha256": refinement_sha,
            })
        if (len({(row["q_u_A"], row["q_v_A"]) for row in rows}) != 25
                or [row["kind"] for row in rows].count("offaxis_refinement") != 12):
            raise ValueError("the local 5x5 DFT coordinates are incomplete")
    if dense9 is not None:
        if (refinement is None or not dense9_sha
                or dense9.get("status") != "GaN_600eV_local_joint_9x9_refinement_raw_audited"
                or dense9.get("n_total_DFT_points") != 81
                or len(dense9.get("cases", [])) != 56
                or dense9.get("source_sha256", {}).get("pilot_audit") != report_sha
                or dense9.get("source_sha256", {}).get("refinement5_audit") != refinement_sha):
            raise ValueError("9x9 data are not tied to both preceding audits")
        for case in dense9["cases"]:
            q_u, q_v = float(case["q_u_A"]), float(case["q_v_A"])
            predicted = float(model_meV_per_GaN(report, q_u, q_v))
            rows.append({
                "case": case["case"], "kind": "dense9_refinement",
                "q_u_A": q_u, "q_v_A": q_v,
                "measured_delta_H_meV_per_GaN": float(case["delta_enthalpy_meV_per_GaN"]),
                "quadratic_model_delta_H_meV_per_GaN": predicted,
                "residual_meV_per_GaN": float(case["delta_enthalpy_meV_per_GaN"] - predicted),
                "raw_OUTCAR_sha256": case["outcar_sha256"],
                "source_report_sha256": dense9_sha,
            })
        if len({(row["q_u_A"], row["q_v_A"]) for row in rows}) != 81:
            raise ValueError("the 9x9 DFT coordinates are incomplete or repeated")
    return rows


def measured_interpolant(rows: list[dict], u_mesh: np.ndarray,
                         v_mesh: np.ndarray) -> np.ndarray:
    if len(rows) not in (25, 81) or len({(r["q_u_A"], r["q_v_A"]) for r in rows}) != len(rows):
        raise ValueError("interpolated DFT cut requires 25 or 81 unique measured coordinates")
    xy = np.asarray([(r["q_u_A"], r["q_v_A"]) for r in rows], dtype=float)
    if len({float(x) for x in xy[:, 0]}) != int(np.sqrt(len(rows))) or len({float(y) for y in xy[:, 1]}) != int(np.sqrt(len(rows))):
        raise ValueError("interpolated DFT cut requires a complete regular grid")
    energy = np.asarray([r["measured_delta_H_meV_per_GaN"] for r in rows], dtype=float)
    interpolator = CloughTocher2DInterpolator(xy, energy)
    z = np.asarray(interpolator(u_mesh, v_mesh), dtype=float)
    if not np.isfinite(z).all():
        raise ValueError("interpolated DFT cut has an uncovered coordinate")
    if not np.allclose(interpolator(xy), energy, atol=1e-9, rtol=0):
        raise ValueError("interpolated DFT cut does not reproduce samples")
    return z


def draw(report: dict, rows: list[dict]) -> plt.Figure:
    q_u = np.linspace(-0.02, 0.02, 241)
    q_v = np.linspace(-0.0125, 0.0125, 151)
    u_mesh, v_mesh = np.meshgrid(q_u, q_v)
    new_points = [row for row in rows if row["kind"] == "offaxis_refinement"]
    dense_points = [row for row in rows if row["kind"] == "dense9_refinement"]
    z = (measured_interpolant(rows, u_mesh, v_mesh) if new_points
         else model_meV_per_GaN(report, u_mesh, v_mesh))
    levels = [-0.50, -0.40, -0.30, -0.20, -0.10, -0.05, 0.0,
              0.025, 0.050, 0.075, 0.10, 0.15]
    if float(np.min(z)) < levels[0] or float(np.max(z)) > levels[-1]:
        raise ValueError("local enthalpy cut would be clipped by figure levels")
    fig = plt.figure(figsize=(7.2, 4.1), facecolor="white")
    ax_a = fig.add_axes((0.105, 0.320, 0.420, 0.550))
    ax_b = fig.add_axes((0.640, 0.320, 0.300, 0.550))
    cax = fig.add_axes((0.158, 0.095, 0.315, 0.027))
    fig.text(0.105, 0.928, "(a)", fontsize=13, fontweight="normal", color=INK)
    fig.text(0.640, 0.928, "(b)", fontsize=13, fontweight="normal", color=INK)

    norm = TwoSlopeNorm(vmin=-0.50, vcenter=0, vmax=0.15)
    contour = ax_a.contourf(u_mesh, v_mesh, z, levels=levels,
                            cmap="PuOr_r", norm=norm, extend="neither")
    ax_a.contour(u_mesh, v_mesh, z, levels=[0], colors=[INK], linewidths=0.9,
                 linestyles="--")
    grid = [row for row in rows if row["kind"] == "grid"]
    holdouts = [row for row in rows if row["kind"] == "axial_holdout"]
    if dense_points:
        measured = [row for row in rows if row["kind"] != "center"]
        ax_a.scatter([row["q_u_A"] for row in measured],
                     [row["q_v_A"] for row in measured], s=33, marker="o",
                     facecolors="white", edgecolors=INK, linewidths=0.7, zorder=6)
    else:
        ax_a.scatter([row["q_u_A"] for row in grid], [row["q_v_A"] for row in grid],
                     s=62, marker="o", facecolors="white", edgecolors=INK,
                     linewidths=1.25, zorder=5)
        ax_a.scatter([row["q_u_A"] for row in holdouts],
                     [row["q_v_A"] for row in holdouts], s=64, marker="s",
                     facecolors=ORANGE, edgecolors="white", linewidths=0.8, zorder=6)
        if new_points:
            ax_a.scatter([row["q_u_A"] for row in new_points],
                         [row["q_v_A"] for row in new_points], s=58, marker="^",
                         facecolors="#2B8F8A", edgecolors="white", linewidths=0.7, zorder=6)
    ax_a.scatter([0], [0], s=155, marker="*", color=BLUE,
                 edgecolors="white", linewidths=0.7, zorder=7)
    ax_a.set_xlim(-0.023, 0.023)
    ax_a.set_ylim(-0.0146, 0.0146)
    ax_a.set_xticks([-0.02, -0.01, 0, 0.01, 0.02])
    ax_a.set_yticks([-0.0125, 0, 0.0125])
    ax_a.set_xlabel(r"Unstable joint coordinate $q_u$ (Å)")
    ax_a.set_ylabel(r"Stable joint coordinate $q_v$ (Å)")
    ax_a.tick_params(axis="both", direction="in", top=True, right=True, length=4, pad=4)
    colorbar = fig.colorbar(contour, cax=cax, orientation="horizontal",
                            ticks=[-0.4, -0.2, 0, 0.1])
    colorbar.ax.xaxis.set_label_position("top")
    colorbar.set_label(r"$\Delta H$ (meV/GaN)", fontsize=9.4, labelpad=4)
    colorbar.ax.tick_params(labelsize=8.7, direction="in", length=3)

    limits = (-0.57, 0.20)
    diagonal = np.array(limits)
    ax_b.fill_between(diagonal, diagonal - 0.2, diagonal + 0.2,
                      color="#EAEFF2", alpha=0.75, zorder=0)
    ax_b.plot(diagonal, diagonal, color="#64727C", lw=1.25, zorder=1)
    if dense_points:
        measured = [row for row in rows if row["kind"] != "center"]
        ax_b.scatter([row["quadratic_model_delta_H_meV_per_GaN"] for row in measured],
                     [row["measured_delta_H_meV_per_GaN"] for row in measured],
                     s=33, marker="o", color=BLUE, edgecolor="white",
                     linewidth=0.55, zorder=3)
    else:
        for kind, marker, color, size in (("grid", "o", BLUE, 58),
                                         ("axial_holdout", "s", ORANGE, 66),
                                         ("offaxis_refinement", "^", "#2B8F8A", 58)):
            subset = [row for row in rows if row["kind"] == kind]
            ax_b.scatter([row["quadratic_model_delta_H_meV_per_GaN"] for row in subset],
                         [row["measured_delta_H_meV_per_GaN"] for row in subset],
                         s=size, marker=marker, color=color, edgecolor="white",
                         linewidth=0.8, zorder=3)
    ax_b.scatter([0], [0], s=110, marker="*", color=INK,
                 edgecolor="white", linewidth=0.7, zorder=4)
    ax_b.set_xlim(*limits)
    ax_b.set_ylim(*limits)
    ax_b.set_xticks([-0.5, -0.3, -0.1, 0.1])
    ax_b.set_yticks([-0.5, -0.3, -0.1, 0.1])
    ax_b.set_xlabel(r"Quadratic model $\Delta H$ (meV/GaN)")
    ax_b.set_ylabel(r"DFT $\Delta H$ (meV/GaN)")
    ax_b.tick_params(axis="both", direction="in", top=True, right=True, length=4, pad=4)
    handles = ([Line2D([], [], marker="o", color="none", markerfacecolor=BLUE,
                       markeredgecolor="white", markersize=7, label="DFT statics (80)")]
               if dense_points else [
        Line2D([], [], marker="o", color="none", markerfacecolor="white",
               markeredgecolor=INK, markeredgewidth=1.1, markersize=7,
               label="Grid (8)"),
        Line2D([], [], marker="s", color="none", markerfacecolor=ORANGE,
               markeredgecolor="white", markersize=7, label="Holdout (4)"),
    ])
    if new_points and not dense_points:
        handles.append(Line2D([], [], marker="^", color="none", markerfacecolor="#2B8F8A",
                              markeredgecolor="white", markersize=7,
                              label="New DFT (12)"))
    handles.append(Line2D([], [], marker="*", color="none", markerfacecolor=BLUE,
                          markeredgecolor="white", markersize=10, label="Center"))
    ax_b.legend(handles=handles, loc="upper left", fontsize=8.9,
                frameon=True, facecolor="white", edgecolor="#8B969C",
                framealpha=0.86, borderpad=0.45, handlelength=1.0,
                handletextpad=0.45)
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--hessian-root", type=Path, default=DEFAULT_HESSIAN_ROOT)
    parser.add_argument("--refinement", type=Path,
                        help="optional raw-audited twelve-point 5x5 completion")
    parser.add_argument("--dense9", type=Path,
                        help="optional raw-audited 56-point 9x9 completion; requires --refinement")
    parser.add_argument("--output-prefix", type=Path, default=DEFAULT_PREFIX)
    args = parser.parse_args()
    prefix = args.output_prefix
    paths = {suffix: prefix.with_name(prefix.name + suffix) for suffix in
             ("_source_data.csv", ".pdf", ".svg", ".tiff", ".png", "_qa.json")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("local joint-cut figure output exists; use a fresh prefix")
    report = load(args.report, args.hessian_root)
    refinement = json.loads(args.refinement.read_text(encoding="utf-8")) if args.refinement else None
    dense9 = json.loads(args.dense9.read_text(encoding="utf-8")) if args.dense9 else None
    rows = source_rows(report, sha256(args.report), refinement,
                       sha256(args.refinement) if args.refinement else None,
                       dense9, sha256(args.dense9) if args.dense9 else None)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with paths["_source_data.csv"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    fig = draw(report, rows)
    fig.savefig(paths[".pdf"])
    fig.savefig(paths[".svg"])
    svg_path = paths[".svg"]
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_path.read_text(encoding="utf-8").splitlines())
        + "\n", encoding="utf-8"
    )
    fig.savefig(paths[".tiff"], dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(paths[".png"], dpi=200)
    plt.close(fig)
    hessian = json.loads((args.hessian_root / "audit.json").read_text(encoding="utf-8"))
    qa = {
        "status": ("GaN_600eV_local_joint_9x9_DFT_refinement_visualized" if dense9 else
                   "GaN_600eV_local_joint_5x5_DFT_refinement_visualized" if refinement else
                   "GaN_600eV_local_joint_quadratic_cut_axial_holdout_validated"),
        "core_conclusion": (
            "One local joint direction lowers 600-eV GaN enthalpy and an orthogonal "
            "direction raises it; a quadratic cut predicts axial holdouts locally."
        ),
        "archetype": ("quantitative 81-point DFT-interpolated contour with model parity validation"
                      if dense9 else "quantitative 25-point DFT-interpolated contour with model parity validation"
                      if refinement else "quantitative grid with model-contour hero and DFT parity validation"),
        "backend": "Python/matplotlib", "final_size_mm": [182.88, 104.14],
        "panel_map": {
            "a": ("Clough-Tocher display interpolation of 81 actual DFT enthalpies within measured bounds"
                  if dense9 else "Clough-Tocher display interpolation of 25 actual DFT enthalpies within measured bounds"
                  if refinement else "Local quadratic enthalpy cut only inside measured joint-coordinate bounds, with eight grid statics, four axial holdouts, and center"),
            "b": ("Eighty actual VASP statics versus the original quadratic prediction; shaded ±0.2-meV/GaN original axial gate band"
                  if dense9 else "Twenty-four actual VASP statics versus the original quadratic prediction; shaded ±0.2-meV/GaN original axial gate band"
                  if refinement else "Twelve actual VASP statics versus Hessian-based quadratic prediction; shaded ±0.2-meV/GaN model-error band"),
        },
        "pressure_GPa": 45.7,
        "electronic_contract": "Original VASP 600 eV/PBE/Ga_d+N/Γ8×8×6",
        "model_gradient_eV_per_A": report["center_gradient_uv_eV_per_A"],
        "model_eigenvalues_eV_per_A2": report["hessian_eigenvalues_uv_eV_per_A2"],
        "sampled_domain_A": {"q_u": [-0.02, 0.02], "q_v": [-0.0125, 0.0125]},
        "n_grid_statics": 8, "n_axial_holdout_statics": 4,
        "n_new_offaxis_refinement_statics": 12 if refinement else 0,
        "n_new_dense9_statics": 56 if dense9 else 0,
        "maximum_grid_model_error_meV_per_GaN": report[
            "quadratic_model_max_grid_error_meV_per_GaN"],
        "maximum_axial_holdout_error_meV_per_GaN": report[
            "quadratic_model_max_axial_holdout_error_meV_per_GaN"],
        **({"maximum_new_offaxis_model_error_meV_per_GaN": refinement[
            "new_point_quadratic_model_max_abs_error_meV_per_GaN"]} if refinement else {}),
        **({"prospective_5x5_max_abs_error_meV_per_GaN": dense9[
            "prior_5x5_prospective_max_abs_error_meV_per_GaN"]} if dense9 else {}),
        "predeclared_axial_holdout_gate_meV_per_GaN": 0.2,
        "energy_force_gradient_mismatch_eV_per_A": hessian[
            "energy_gradient_max_abs_difference_eV_per_A"],
        "predeclared_local_gate_pass": report["predeclared_gate"],
        "TS_certified": False, "whole_path_2D_surface_certified": False,
        "model_scope": "Frozen local joint atom-strain enthalpy slice; no orthogonal relaxation, no FES",
        "image_integrity": ("The contour is an explicitly labeled interpolation of raw-audited DFT enthalpies only within measured q bounds; measured nodes are overlaid."
                            if refinement else "Analytic quadratic model is drawn only within measured q bounds; actual statics are overlaid."),
        "reviewer_risk": ("The 81 measured points validate a local frozen enthalpy cut, not orthogonal minimization, a full-variable-cell TS, or a whole-path PES."
                          if dense9 else "The extra off-axis DFT points test this local quadratic cut but neither certify a full-variable-cell TS nor a whole-path PES."
                          if refinement else "Axial half-step holdouts do not independently validate every off-axis interior point or certify a full-variable-cell TS."),
        "raster_export_policy": "PDF/SVG/PNG are tracked; 600-dpi TIFF is generated locally and repository-ignored.",
        "caption_draft": (
            "Local frozen joint-mode enthalpy slice at the GaN B4→B1 highest-image candidate "
            "(45.7 GPa, VASP/PBE, 600 eV). (a) Clough–Tocher display interpolation "
            "of 81 distinct raw-audited DFT coordinates on a complete 9×9 grid; circles "
            "mark the 80 off-center statics and a star marks the center. The source data "
            "identify the 13-point pilot, 12-point 5×5 refinement, and 56-point 9×9 "
            "refinement. (b) Direct DFT enthalpy changes against the original local "
            "quadratic model; the shaded ±0.2-meV/GaN band was predeclared only for "
            "the four axial holdouts. The 56 new points prospectively test the prior "
            "5×5 interpolation. This local frozen slice neither certifies a stationary "
            "transition state nor describes the whole path."
            if dense9 else
            "Local joint-mode enthalpy cut at the GaN B4→B1 highest-image candidate "
            "(45.7 GPa, VASP/PBE, 600 eV). (a) Clough–Tocher display interpolation "
            "within the measured rectangle, with 25 actual DFT coordinates: "
            "the center, eight original grid statics, four original axial holdouts, "
            "and twelve independent off-axis refinement statics. (b) Direct DFT "
            "enthalpy changes against the original model; the shaded ±0.2-meV/GaN "
            "band was predeclared for the original axial holdouts, not retroactively "
            "for the new off-axis points. This frozen local cut neither certifies a "
            "stationary transition state nor describes the whole path."
            if refinement else
            "Local joint-mode enthalpy cut at the GaN B4→B1 highest-image candidate "
            "(45.7 GPa, VASP/PBE, 600 eV). (a) Quadratic model from the audited "
            "atom–strain Hessian and center gradient, restricted to the measured "
            "coordinate rectangle; symbols show eight grid statics, four axial "
            "half-step holdouts, and the center. (b) Direct DFT enthalpy changes "
            "against model predictions; the band is ±0.2 meV/GaN. The largest "
            "axial holdout error is 0.050 meV/GaN. This frozen local cut neither "
            "certifies a true transition state nor describes the whole path."
        ),
        "source_sha256": {
            "pilot_raw_audit": sha256(args.report),
            **({"refinement_raw_audit": sha256(args.refinement)} if refinement else {}),
            **({"dense9_raw_audit": sha256(args.dense9)} if dense9 else {}),
            "hessian_raw_audit": sha256(args.hessian_root / "audit.json"),
            "hessian_npz": sha256(args.hessian_root / "joint_hessian.npz"),
            "plot_script": sha256(Path(__file__)),
            "source_csv": sha256(paths["_source_data.csv"]),
        },
        "exports_sha256": {suffix.lstrip("."): sha256(paths[suffix])
                           for suffix in (".pdf", ".svg", ".tiff", ".png")},
    }
    paths["_qa.json"].write_text(json.dumps(qa, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": qa["status"],
                      "max_axial_holdout_error_meV_per_GaN": qa[
                          "maximum_axial_holdout_error_meV_per_GaN"],
                      "TS_certified": qa["TS_certified"]}))


if __name__ == "__main__":
    main()
