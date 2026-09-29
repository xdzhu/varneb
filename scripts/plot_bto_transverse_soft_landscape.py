"""Exploratory BTO soft/soft frozen-C-cell contour with the VCNEB shadow.

Figure contract: the real T-to-C VCNEB path projects onto the longitudinal
soft axis, while its cell strain and omitted stable-mode motion keep it off
the two-dimensional frozen cubic-cell surface. Panel (a) shows 25, 41, 59,
81, or 289 audited DFT samples and an explicitly interpolated display contour; (b) the seven actual
variable-cell image energies; (c) their off-plane atomic residual. This is a
quantitative-grid figure, not a conditional PES or an activation-barrier map.
Python/matplotlib is the existing project plotting workflow. Editable SVG and
PDF plus a PNG preview and source-data CSV are exported together.
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
import numpy as np
from scipy.interpolate import CloughTocher2DInterpolator


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "font.size": 9, "axes.labelsize": 10, "xtick.labelsize": 9,
    "ytick.labelsize": 9, "legend.fontsize": 8.5,
    "axes.linewidth": 0.8,
})


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _source(assembled: dict, path_audit: dict, path_report: dict,
            path_report_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if (assembled.get("kind") != "BTO_transverse_soft_frozen_C_cell_25_raw_DFT_points_not_conditional_PES_or_MEP"
            or assembled.get("status") != "assembled_from_25_independently_audited_real_DFT_points"
            or assembled.get("n_reused_DFT_points") != 5
            or assembled.get("n_new_DFT_points") != 20
            or len(assembled.get("samples", [])) != 25
            or path_audit.get("kind") != "offline_transverse_soft_plane_audit_not_PES"
            or path_audit.get("input_sha256", {}).get("report") != _sha256(path_report_path)
            or path_report.get("kind") != "posthoc_q1q2_projection_not_PES"
            or path_audit.get("axis_labels") != assembled.get("axis_labels")
            or path_audit.get("axis_mode_weights") != assembled.get("axis_mode_weights")
            or len(path_audit.get("images", [])) != 7
            or len(path_report.get("images", [])) != 7):
        raise ValueError("figure inputs lack one complete, aligned soft-plane provenance chain")
    q1 = sorted({float(sample["q_parallel_sqrt_amu_A"]) for sample in assembled["samples"]})
    q2 = sorted({float(sample["q_transverse_sqrt_amu_A"]) for sample in assembled["samples"]})
    if len(q1) != 5 or len(q2) != 5:
        raise ValueError("expected a full five-by-five measured grid")
    energy = np.full((5, 5), np.nan)
    for sample in assembled["samples"]:
        i, j = sample["grid_index_q1_q2"]
        if (not (0 <= i < 5 and 0 <= j < 5) or np.isfinite(energy[j, i])
                or abs(q1[i] - sample["q_parallel_sqrt_amu_A"]) > 1e-12
                or abs(q2[j] - sample["q_transverse_sqrt_amu_A"]) > 1e-12):
            raise ValueError("measured grid has a duplicate or misplaced sample")
        energy[j, i] = float(sample["energy_minus_C_eV_per_BTO"])
    if not np.all(np.isfinite(energy)):
        raise ValueError("measured DFT grid is incomplete")
    path_q = np.asarray([[row["q_parallel_sqrt_amu_A"], row["q_transverse_sqrt_amu_A"]]
                         for row in path_audit["images"]], dtype=float)
    path_residual = np.asarray([row["off_plane_atomic_norm_sqrt_amu_A"]
                                for row in path_audit["images"]], dtype=float)
    path_e = np.asarray([row["relative_energy_eV_per_formula_unit"]
                         for row in path_report["images"]], dtype=float)
    if (not all(np.all(np.isfinite(item)) for item in (path_q, path_residual, path_e))
            or np.max(np.abs(path_q[:, 1])) > 1e-8
            or abs(path_e[0]) > 1e-9
            or len({row["index"] for row in path_audit["images"]}) != 7):
        raise ValueError("archived VCNEB path projection or energy is inconsistent")
    return np.asarray(q1), np.asarray(q2), energy, np.c_[path_q, path_e, path_residual]


def _interpolation_audit(points: np.ndarray, values: np.ndarray,
                         q1: np.ndarray, q2: np.ndarray) -> tuple[np.ndarray, dict]:
    dense_q1 = np.linspace(q1[0], q1[-1], 401)
    dense_q2 = np.linspace(q2[0], q2[-1], 201)
    gx, gy = np.meshgrid(dense_q1, dense_q2)
    surface = np.asarray(CloughTocher2DInterpolator(points, values)(gx, gy), dtype=float)
    if not np.all(np.isfinite(surface)):
        raise ValueError("display interpolator is undefined within the measured rectangle")
    errors = []
    interior = np.flatnonzero((points[:, 0] > q1[0] + 1e-12)
                             & (points[:, 0] < q1[-1] - 1e-12)
                             & (points[:, 1] > q2[0] + 1e-12)
                             & (points[:, 1] < q2[-1] - 1e-12))
    for index in interior:
        retain = np.arange(len(points)) != index
        predicted = float(CloughTocher2DInterpolator(points[retain], values[retain])(points[index]))
        if not np.isfinite(predicted):
            raise ValueError("interior leave-one-out interpolation is undefined")
        errors.append(1000.0 * (predicted - values[index]))
    audit = {
        "n_measured_DFT_samples": len(points),
        "n_interior_leave_one_out": len(interior),
        "interior_leave_one_out_max_abs_error_meV_per_BTO": float(np.max(np.abs(errors))),
        "interior_leave_one_out_rms_error_meV_per_BTO": float(np.sqrt(np.mean(np.square(errors)))),
        "interior_leave_one_out_signed_errors_meV_per_BTO": errors,
        "display_pixel_grid": [401, 201],
        "sampled_energy_range_eV_per_BTO": [float(np.min(values)), float(np.max(values))],
        "interpolated_energy_range_eV_per_BTO": [float(np.min(surface)), float(np.max(surface))],
        "interpolated_downward_overshoot_meV_per_BTO": float(max(0.0, 1000.0 * (np.min(values) - np.min(surface)))),
        "claim_limit": "interpolation is visual only; interior leave-one-out checks do not certify the boundary or any extremum",
    }
    return surface, audit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembled", type=Path, required=True)
    parser.add_argument("--path-audit", type=Path, required=True)
    parser.add_argument("--path-report", type=Path, required=True)
    parser.add_argument("--analysis41", type=Path,
                        help="optional independently scored 16-center result; plots 41 real points")
    parser.add_argument("--analysis59", type=Path,
                        help="optional independently scored 18-edge result; requires --analysis41")
    parser.add_argument("--analysis81", type=Path,
                        help="optional raw-audited 22-node completion; requires --analysis59")
    parser.add_argument("--analysis289", type=Path,
                        help="optional raw-audited nested 17x17 completion; requires --analysis81")
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--allow-exploratory-contours", action="store_true")
    args = parser.parse_args()
    if not args.allow_exploratory_contours:
        raise ValueError("exploratory interpolated contours require --allow-exploratory-contours")
    if args.analysis59 is not None and args.analysis41 is None:
        raise ValueError("59-point plot requires its 41-point provenance source")
    if args.analysis81 is not None and args.analysis59 is None:
        raise ValueError("81-point plot requires its 59-point provenance source")
    if args.analysis289 is not None and args.analysis81 is None:
        raise ValueError("289-point plot requires its 81-point provenance source")
    outputs = [args.output_prefix.with_suffix(suffix) for suffix in (".png", ".pdf", ".svg")]
    outputs += [Path(str(args.output_prefix) + suffix) for suffix in ("_source_data.csv", "_qa.json")]
    if any(path.exists() for path in outputs):
        raise FileExistsError("refusing to overwrite a previous figure or source-data file")
    assembled, path_audit, path_report = map(_load, (args.assembled, args.path_audit, args.path_report))
    q1, q2, energy, path = _source(assembled, path_audit, path_report, args.path_report)
    xx, yy = np.meshgrid(q1, q2)
    measured_points = np.c_[xx.ravel(), yy.ravel()]
    measured_energies = energy.ravel()
    analysis41 = None
    if args.analysis41 is not None:
        analysis41 = _load(args.analysis41)
        if (analysis41.get("kind") != "BTO_transverse_soft_frozen_C_cell_41_real_DFT_points_not_conditional_PES_or_MEP"
                or analysis41.get("status") != "sixteen_independent_center_holdouts_scored"
                or analysis41.get("n_real_DFT_points") != 41
                or analysis41.get("source_sha256", {}).get("assembled25") != _sha256(args.assembled)
                or analysis41.get("axis_labels") != assembled.get("axis_labels")
                or len(analysis41.get("samples", [])) != 41):
            raise ValueError("41-point analysis is not tied to this audited 25-point base")
        measured_points = np.asarray([[sample["q_parallel_sqrt_amu_A"],
                                       sample["q_transverse_sqrt_amu_A"]]
                                      for sample in analysis41["samples"]], dtype=float)
        measured_energies = np.asarray([sample["energy_minus_C_eV_per_BTO"]
                                        for sample in analysis41["samples"]], dtype=float)
        if (len({tuple(point) for point in measured_points}) != 41
                or not np.all(np.isfinite(measured_energies))):
            raise ValueError("41-point measured source has duplicate or nonfinite values")
    analysis59 = None
    if args.analysis59 is not None:
        analysis59 = _load(args.analysis59)
        if (analysis59.get("kind") != "BTO_transverse_soft_frozen_C_cell_59_real_DFT_points_not_conditional_PES_or_MEP"
                or analysis59.get("status") != "eighteen_independent_adaptive_edge_holdouts_scored"
                or analysis59.get("n_real_DFT_points") != 59
                or analysis59.get("source_sha256", {}).get("analysis41") != _sha256(args.analysis41)
                or analysis59.get("source_sha256", {}).get("assembled25") != _sha256(args.assembled)
                or analysis59.get("axis_labels") != assembled.get("axis_labels")
                or len(analysis59.get("samples", [])) != 59):
            raise ValueError("59-point analysis is not tied to the audited 25/41-point sources")
        measured_points = np.asarray([[sample["q_parallel_sqrt_amu_A"],
                                       sample["q_transverse_sqrt_amu_A"]]
                                      for sample in analysis59["samples"]], dtype=float)
        measured_energies = np.asarray([sample["energy_minus_C_eV_per_BTO"]
                                        for sample in analysis59["samples"]], dtype=float)
        if (len({tuple(point) for point in measured_points}) != 59
                or not np.all(np.isfinite(measured_energies))):
            raise ValueError("59-point measured source has duplicate or nonfinite values")
    analysis81 = None
    if args.analysis81 is not None:
        analysis81 = _load(args.analysis81)
        if (analysis81.get("kind") != "BTO_transverse_soft_frozen_C_cell_81_real_DFT_points_not_conditional_PES_or_MEP"
                or analysis81.get("status") != "22_new_nodes_raw_energy_force_stress_audited"
                or analysis81.get("n_total_DFT_points") != 81
                or analysis81.get("source_sha256", {}).get("prior_59_audit") != _sha256(args.analysis59)
                or analysis81.get("axis_labels") != assembled.get("axis_labels")
                or len(analysis81.get("samples", [])) != 81):
            raise ValueError("81-point result is not tied to its audited 59-point source")
        measured_points = np.asarray([[sample["q_parallel_sqrt_amu_A"],
                                       sample["q_transverse_sqrt_amu_A"]]
                                      for sample in analysis81["samples"]], dtype=float)
        measured_energies = np.asarray([sample["energy_minus_C_eV_per_BTO"]
                                        for sample in analysis81["samples"]], dtype=float)
        if (len({tuple(point) for point in measured_points}) != 81
                or not np.all(np.isfinite(measured_energies))):
            raise ValueError("81-point measured source has duplicate or nonfinite values")
    analysis289 = None
    if args.analysis289 is not None:
        analysis289 = _load(args.analysis289)
        if (analysis289.get("status") != "BTO_frozen_soft_mode_17x17_raw_audited"
                or analysis289.get("n_reused_measured_points") != 81
                or analysis289.get("n_new_raw_audited_statics") != 208
                or analysis289.get("source_sha256", {}).get("old_audit")
                   != _sha256(args.analysis81)):
            raise ValueError("289-point result is not tied to the audited 81-point source")
        fine_q1 = np.asarray(analysis289["q1"], dtype=float)
        fine_q2 = np.asarray(analysis289["q2"], dtype=float)
        measured_energies = np.asarray(
            analysis289["energy_minus_C_eV_per_BTO"], dtype=float).T.ravel()
        fine_x, fine_y = np.meshgrid(fine_q1, fine_q2)
        measured_points = np.c_[fine_x.ravel(), fine_y.ravel()]
        if (len(measured_points) != 289 or not np.isfinite(measured_energies).all()
                or not np.allclose(fine_q1[::4], q1, atol=1e-12, rtol=0)
                or not np.allclose(fine_q2[::4], q2, atol=1e-12, rtol=0)):
            raise ValueError("289-point nested measured grid is incomplete")
    surface, interpolation = _interpolation_audit(measured_points, measured_energies, q1, q2)
    if analysis41 is not None:
        interpolation["independent_16_center_prediction_errors"] = analysis41["holdout_errors"]
    if analysis59 is not None:
        interpolation["independent_18_adaptive_edge_prediction_errors"] = analysis59["holdout_errors"]
    if analysis81 is not None:
        interpolation["independent_22_missing_node_max_abs_error_meV_per_BTO"] = analysis81[
            "prospective_22_point_max_abs_error_meV_per_BTO"
        ]
        interpolation["independent_22_missing_node_rms_error_meV_per_BTO"] = analysis81[
            "prospective_22_point_rms_error_meV_per_BTO"
        ]
    if analysis289 is not None:
        interpolation["prospective_208_nested_node_max_abs_error_meV_per_BTO"] = analysis289[
            "prospective_max_abs_error_meV_per_BTO"
        ]
        interpolation["prospective_208_nested_node_rms_error_meV_per_BTO"] = analysis289[
            "prospective_rms_error_meV_per_BTO"
        ]
    gx, gy = np.meshgrid(np.linspace(q1[0], q1[-1], surface.shape[1]),
                         np.linspace(q2[0], q2[-1], surface.shape[0]))
    span = max(abs(np.min(surface)), abs(np.max(surface)))
    norm = TwoSlopeNorm(vmin=-span, vcenter=0.0, vmax=span)
    # Fixed axes avoid constrained-layout/colorbar interactions that shrink
    # the two diagnostic panels and let their y labels collide with panel (a).
    fig = plt.figure(figsize=(10.4, 5.0))
    ax = fig.add_axes((0.09, 0.20, 0.50, 0.70))
    ax_e = fig.add_axes((0.70, 0.57, 0.26, 0.31))
    ax_r = fig.add_axes((0.70, 0.19, 0.26, 0.31), sharex=ax_e)
    cax = fig.add_axes((0.13, 0.075, 0.42, 0.035))
    filled = ax.contourf(gx, gy, surface, levels=np.linspace(-span, span, 25),
                         cmap="RdBu_r", norm=norm, extend="both")
    ax.contour(gx, gy, surface, levels=np.linspace(-span, span, 13),
               colors="#3e4c59", linewidths=0.45, alpha=0.7)
    ax.scatter(measured_points[:, 0], measured_points[:, 1],
               s=8 if analysis289 is not None else 15,
               facecolors="white", edgecolors="#30485c",
               linewidths=0.65, zorder=4,
               label=f"{len(measured_points)} computed DFT points")
    ax.plot(path[:, 0], path[:, 1], color="#171717", lw=1.35, zorder=5,
            label="VCNEB projection only")
    ax.scatter(path[:, 0], path[:, 1], s=29, facecolors="#ffcf70",
               edgecolors="#171717", linewidths=0.8, zorder=6)
    ax.annotate("C", (path[-1, 0], path[-1, 1]), xytext=(-12, 11),
                textcoords="offset points", fontsize=9,
                bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.8})
    ax.annotate("T projection", (path[0, 0], path[0, 1]), xytext=(-65, 11),
                textcoords="offset points", fontsize=8,
                bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.8})
    ax.set_xlabel(r"$Q_{\parallel}$ ($\sqrt{\mathrm{amu}}$ Å)")
    ax.set_ylabel(r"$Q_{\perp}$ ($\sqrt{\mathrm{amu}}$ Å)")
    ax.set_xlim(q1[0] - 0.02, q1[-1] + 0.04)
    ax.set_ylim(q2[0], q2[-1])
    ax.legend(loc="upper center", ncol=2, frameon=True, fancybox=True,
              facecolor="white", edgecolor="#79818a", framealpha=0.85)
    colorbar = fig.colorbar(filled, cax=cax, orientation="horizontal")
    colorbar.set_ticks(np.linspace(-0.06, 0.06, 5))
    colorbar.set_label(r"Frozen $E-E_{\mathrm{C}}$ (eV/BTO)")
    indices = np.arange(len(path))
    ax_e.plot(indices, path[:, 2], "o-", lw=1.8, ms=4.6, color="#1f5a99")
    ax_e.axhline(0.0, color="#666666", lw=0.8, ls="--")
    ax_e.set_ylabel(r"$E-E_{\mathrm{T}}$ (eV/BTO)")
    ax_e.tick_params(labelbottom=False)
    ax_r.plot(indices, path[:, 3], "o-", lw=1.8, ms=4.6, color="#a5513e")
    ax_r.set_ylabel(r"Off-plane norm ($\sqrt{\mathrm{amu}}$ Å)")
    ax_r.set_xlabel("VCNEB image (T → C)")
    ax_r.set_xticks(indices)
    for panel in (ax, ax_e, ax_r):
        panel.tick_params(direction="in", top=True, right=True, width=0.8, length=4)
        for spine in panel.spines.values():
            spine.set_visible(True)
    for position, label in (((0.09, 0.92), "(a)"), ((0.70, 0.90), "(b)"),
                            ((0.70, 0.52), "(c)")):
        fig.text(*position, label, fontsize=12, fontweight="normal")
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for output in outputs[:3]:
        fig.savefig(output, dpi=400, bbox_inches="tight", facecolor="white")
    svg_path = args.output_prefix.with_suffix(".svg")
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_path.read_text(encoding="utf-8").splitlines())
        + "\n", encoding="utf-8"
    )
    plt.close(fig)
    with outputs[3].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["record", "q_parallel_sqrt_amu_A", "q_transverse_sqrt_amu_A",
                         "E_minus_C_eV_per_BTO", "VCNEB_E_minus_T_eV_per_BTO",
                         "VCNEB_offplane_atomic_residual_sqrt_amu_A"])
        for point, value in zip(measured_points, measured_energies):
            writer.writerow(["frozen_DFT", point[0], point[1], value, "", ""])
        for row in path:
            writer.writerow(["VCNEB_image", row[0], row[1], "", row[2], row[3]])
    qa = {
        "kind": "BTO_exploratory_transverse_soft_frozen_C_cell_contour_not_conditional_PES_or_MEP",
        "source_sha256": {"assembled": _sha256(args.assembled),
                          "path_audit": _sha256(args.path_audit),
                          "path_report": _sha256(args.path_report),
                          **({"analysis41": _sha256(args.analysis41)}
                             if args.analysis41 is not None else {}),
                          **({"analysis59": _sha256(args.analysis59)}
                             if args.analysis59 is not None else {}),
                          **({"analysis81": _sha256(args.analysis81)}
                             if args.analysis81 is not None else {}),
                          **({"analysis289": _sha256(args.analysis289)}
                             if args.analysis289 is not None else {})},
        "interpolation": interpolation,
        "figure_contract": {
            "core_conclusion": "The variable-cell T-to-C path projects onto the longitudinal soft mode but leaves the frozen cubic-cell two-mode plane.",
            "archetype": "quantitative grid with a hero contour and two path diagnostics",
            "backend": "Python/matplotlib", "size_inches": [10.4, 5.0],
            "panels": {"a": "measured frozen DFT points plus display interpolation and path shadow",
                       "b": "seven actual variable-cell image energies",
                       "c": "seven atomic off-plane residuals"},
            "review_risk": "interpolated pixels, projected path and relaxed barrier are distinct objects",
        },
        "source_data_csv_sha256": _sha256(outputs[3]),
    }
    outputs[4].write_text(json.dumps(qa, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"outputs": [str(path) for path in outputs], "interpolation": interpolation},
                     indent=2))


if __name__ == "__main__":
    main()
