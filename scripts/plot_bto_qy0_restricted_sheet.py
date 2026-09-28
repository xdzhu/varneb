"""Plot the audited BTO Qy=0 local sheet and prospective model holdouts.

The shaded field is an even-mode polynomial trained on nine DFT nodes. It is
never extrapolated outside the convex hull of measured nodes plus the audited
T endpoint. Every DFT location used to discuss the model is drawn explicitly.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from scipy.spatial import Delaunay

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.bto_q1q2_reference import (
    bto_transverse_soft_plane, load_bto_q1q2_reference, validate_bto_gamma_source,
)
from examples.preflight_bto_transverse_soft_conditional import remaining_soft_y_direction
from scripts.plan_bto_qy0_even_mode_validation import even_basis, plan as reconstruct_plan

DATA = ROOT / "outputs" / "batio3_t_to_c_pbe100_dzp10au"
BENCH = ROOT / "benchmarks" / "numerical_integrity"
FIGURES = ROOT / "paper" / "VARNEB_CPC" / "figures"
BASENAME = "bto_qy0_restricted_even_mode_sheet"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def same_structure_up_to_translation(atoms, record: dict) -> bool:
    if len(atoms) != record["n_atoms"] or atoms.get_chemical_symbols() != record["species_order"]:
        return False
    if np.max(np.abs(atoms.cell.array - np.asarray(record["cell_A"]))) > 1e-8:
        return False
    delta = np.asarray(record["fractional_positions_wrapped"]) - atoms.get_scaled_positions(wrap=True)
    delta = (delta + 0.5) % 1.0 - 0.5
    relative = (delta - delta[0] + 0.5) % 1.0 - 0.5
    return bool(np.max(np.abs(relative)) <= 1e-8)


def evidence() -> tuple[dict, list[dict], list[dict], float, float]:
    measured_path = BENCH / "bto_qy0_ten_measured_nodes_20260928.json"
    old_path = BENCH / "bto_qy0_three_holdout_gate_20260928.json"
    plan_path = BENCH / "bto_qy0_even_mode_two_holdout_plan_20260928.json"
    gate_path = BENCH / "bto_qy0_even_mode_two_holdout_gate_20260928.json"
    endpoint_path = DATA / "bto_endpoint_completed_static_audit_job27778739.json"
    endpoint_summary_path = DATA / "bto_endpoint_completed_static_summary_job27778739.json"
    path_projection_path = DATA / "bto_transverse_soft_plane_offline_audit_2026-09-26.json"
    trajectory_path = DATA / "bto_tetragonal_to_cubic_n7_static_audit_step3.traj"
    manifest_path = DATA / "result_manifest_hf_27777454.json"
    sources = [measured_path, old_path, plan_path, gate_path, endpoint_path,
               endpoint_summary_path, path_projection_path, trajectory_path, manifest_path,
               DATA / "bto_t_to_c_q1q2_strain_gradient_preflight_2026-09-25.json",
               DATA / "cubic_CONTCAR", DATA / "bto_cubic_gamma_force_constants.npz",
               DATA / "bto_cubic_phonopy_gamma_eigenpairs.npz",
               DATA / "bto_cubic_gamma_phonon_provenance.json",
               DATA / "bto_cubic_gamma_FORCE_SETS",
               DATA / "bto_cubic_phonopy_gamma_eigenpairs_provenance.json"]
    hashes = {str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path) for path in sources}
    measured, old, frozen, gate = (load_json(path) for path in
                                   (measured_path, old_path, plan_path, gate_path))
    reconstructed = reconstruct_plan(measured, old)
    if (gate["status"] != "two_prospective_local_gates_passed_not_global_PES_certificate"
            or len(gate["points"]) != 2 or not all(point["all_gates_pass"] for point in gate["points"])
            or gate["frozen_plan_sha256"] != sha256(plan_path)
            or gate["gate_script_sha256"] != sha256(ROOT / "scripts" / "audit_bto_qy0_even_mode_holdouts.py")
            or any(frozen.get(key) != value for key, value in reconstructed.items())):
        raise ValueError("frozen model and prospectively audited gate disagree")
    endpoint, endpoint_summary = load_json(endpoint_path), load_json(endpoint_summary_path)
    manifest, path_projection = load_json(manifest_path), load_json(path_projection_path)
    if (endpoint["status"] != "passed"
            or endpoint["static_summary_sha256"] != sha256(endpoint_summary_path)
            or endpoint_summary["calculator_parameters"] != manifest["calculator_parameters"]
            or path_projection["input_sha256"]["trajectory"] != sha256(trajectory_path)):
        raise ValueError("endpoint or path source violates the BTO calculator/trajectory contract")
    cubic = [point for point in manifest["points"] if abs(point["q1"]) < 1e-12
             and abs(point["q2"]) < 1e-12]
    if len(cubic) != 1 or abs(cubic[0]["energy_eV"] - endpoint["endpoints"][1]["energy_eV"]) > 1e-6:
        raise ValueError("endpoint and sheet cubic energy references disagree")
    images = read(str(trajectory_path), index=":")
    if len(images) != 7 or len(path_projection["images"]) != 7:
        raise ValueError("seven-image T-to-C path is missing")
    for atoms, record in zip((images[0], images[-1]),
                             (endpoint_summary["endpoint_structures"]["initial"],
                              endpoint_summary["endpoint_structures"]["final"])):
        if not same_structure_up_to_translation(atoms, record):
            raise ValueError("path endpoint and independently audited static endpoint differ")
    validate_bto_gamma_source(DATA / "bto_cubic_gamma_phonon_provenance.json",
                              DATA / "bto_cubic_gamma_FORCE_SETS",
                              DATA / "bto_cubic_phonopy_gamma_eigenpairs_provenance.json")
    reference = load_bto_q1q2_reference(
        DATA / "bto_t_to_c_q1q2_strain_gradient_preflight_2026-09-25.json",
        DATA / "cubic_CONTCAR", DATA / "bto_cubic_gamma_force_constants.npz",
        DATA / "bto_cubic_phonopy_gamma_eigenpairs.npz")
    plane = bto_transverse_soft_plane(reference, strain_metric_weights_amu_A2=np.ones(6))
    third = remaining_soft_y_direction(reference.modes, plane)
    centered, _ = reference.chart.remove_mass_weighted_translation(
        reference.chart.from_atoms(images[0]), reference.modes.masses_amu)
    t_q = plane.project(centered)
    t_qy = float(third @ (plane.metric_weights * centered))
    archived_t = path_projection["images"][0]
    if (abs(t_q[0] - archived_t["q_parallel_sqrt_amu_A"]) > 1e-8
            or abs(t_q[1] - archived_t["q_transverse_sqrt_amu_A"]) > 1e-8
            or abs(t_qy) > 1e-8
            or abs(path_projection["images"][-1]["q_parallel_sqrt_amu_A"]) > 1e-8):
        raise ValueError("T/C endpoint is not in the Qy=0 coordinate slice")
    coefficients = np.asarray(frozen["coefficients_eV_per_BTO"])
    predict = lambda q: float(even_basis(np.asarray(q))[0] @ coefficients)
    points = []
    for row in measured["points"]:
        points.append({"qz": row["q"][0], "qx": row["q"][1],
                       "class": "fit_node" if row["role"] == "grid_node" else "retrospective",
                       "observed_meV_per_BTO": 1000 * row["energy_minus_c_eV_per_BTO"],
                       "predicted_meV_per_BTO": 1000 * predict(row["q"])})
    for row in old["points"]:
        points.append({"qz": row["q_sqrt_amu_A"][0], "qx": row["q_sqrt_amu_A"][1],
                       "class": "retrospective", "observed_meV_per_BTO":
                       1000 * row["audited_energy_minus_C_eV_per_BTO"],
                       "predicted_meV_per_BTO": 1000 * predict(row["q_sqrt_amu_A"])})
    for row in gate["points"]:
        points.append({"qz": row["q_sqrt_amu_A"][0], "qx": row["q_sqrt_amu_A"][1],
                       "class": "prospective", "observed_meV_per_BTO":
                       1000 * row["audited_energy_minus_C_eV_per_BTO"],
                       "predicted_meV_per_BTO":
                       1000 * row["frozen_model_prediction_energy_minus_C_eV_per_BTO"]})
    t_energy = endpoint["endpoints"][0]["energy_eV"] - endpoint["endpoints"][1]["energy_eV"]
    points.append({"qz": float(t_q[0]), "qx": float(t_q[1]), "class": "T_endpoint_retrospective",
                   "observed_meV_per_BTO": 1000 * t_energy,
                   "predicted_meV_per_BTO": 1000 * predict(t_q)})
    for point in points:
        point["observed_minus_predicted_meV_per_BTO"] = (
            point["observed_meV_per_BTO"] - point["predicted_meV_per_BTO"])
    qa = {"status": "restricted_Qy0_local_model_two_prospective_gates_pass_not_global_PES",
          "evidence_sources_sha256": hashes, "model_training_nodes": 9,
          "prospective_holdout_nodes": 2,
          "maximum_absolute_prospective_error_meV_per_BTO":
          max(abs(point["observed_minus_predicted_meV_per_BTO"])
              for point in points if point["class"] == "prospective"),
          "gate_meV_per_BTO": frozen["energy_absolute_error_gate_meV_per_BTO"],
          "T_qz_qx_qy_sqrt_amu_A": [float(t_q[0]), float(t_q[1]), t_qy],
          "T_atomic_off_two_mode_plane_sqrt_amu_A":
          float(path_projection["images"][0]["off_plane_atomic_norm_sqrt_amu_A"]),
          "T_strain_voigt": path_projection["images"][0]["symmetric_strain_voigt"],
          "T_endpoint_observed_minus_model_meV_per_BTO":
          points[-1]["observed_minus_predicted_meV_per_BTO"],
          "figure_scope": "Qy=0 symmetry-restricted variable-cell local sheet; fitted model only inside measured convex hull; not unrestricted PES, finite-temperature FES or activation barrier",
          "limitations": "Two blind points do not certify all branches or restricted Hessian positivity throughout the domain. T is an independently optimized endpoint, not a model-training node."}
    return qa, points, path_projection["images"], float(t_q[0]), float(t_qy)


def render() -> None:
    qa, points, images, t_qz, _ = evidence()
    frozen = load_json(BENCH / "bto_qy0_even_mode_two_holdout_plan_20260928.json")
    coefficients = np.asarray(frozen["coefficients_eV_per_BTO"])
    fit_nodes = np.asarray([[point["qz"], point["qx"]] for point in points
                            if point["class"] == "fit_node"])
    hull = Delaunay(np.vstack((fit_nodes, [t_qz, 0.0])))
    qz_axis = np.linspace(0, t_qz, 362)
    qx_axis = np.linspace(0, 0.3, 91)
    qz, qx = np.meshgrid(qz_axis, qx_axis)
    coords = np.column_stack((qz.ravel(), qx.ravel()))
    inside = hull.find_simplex(coords, tol=1e-10) >= 0
    model = (1000 * even_basis(coords) @ coefficients).reshape(qz.shape)
    model = np.ma.masked_where(~inside.reshape(qz.shape), model)

    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 9.3,
                         "axes.linewidth": 0.85, "axes.spines.top": True,
                         "axes.spines.right": True, "legend.frameon": True})
    fig = plt.figure(figsize=(7.2, 4.35), layout="constrained")
    grid = fig.add_gridspec(2, 2, width_ratios=[1, 0.028], height_ratios=[1.55, 1.0],
                            hspace=0.055, wspace=0.04)
    ax = fig.add_subplot(grid[0, 0])
    residual_ax = fig.add_subplot(grid[1, 0], sharex=ax)
    color_ax = fig.add_subplot(grid[0, 1])
    levels = np.arange(-100, 6, 5)
    contour = ax.contourf(qz, qx, model, levels=levels, cmap="cividis_r", extend="both")
    ax.contour(qz, qx, model, levels=np.arange(-90, 1, 15), colors="white",
               alpha=0.65, linewidths=0.55)
    fig.colorbar(contour, cax=color_ax, label=r"$E-E_C$ (meV/BTO)")
    ax.plot([row["q_parallel_sqrt_amu_A"] for row in images],
            [row["q_transverse_sqrt_amu_A"] for row in images],
            color="#ce533d", marker="o", markersize=3.8, linewidth=1.5,
            label="7-image VCNEB projection", zorder=5)
    ax.scatter(fit_nodes[:, 0], fit_nodes[:, 1], s=26, color="white",
               edgecolor="#25344a", linewidth=0.9, zorder=6, label="9 fit nodes")
    retrospective = [row for row in points if row["class"] == "retrospective"]
    ax.scatter([row["qz"] for row in retrospective],
               [row["qx"] for row in retrospective], s=36, facecolor="white",
               edgecolor="#687382", linewidth=1.0, zorder=7,
               label="4 retrospective checks")
    blind = [row for row in points if row["class"] == "prospective"]
    ax.scatter([row["qz"] for row in blind], [row["qx"] for row in blind], s=58,
               marker="D", facecolor="#e5a547", edgecolor="#583519", linewidth=0.9,
               zorder=7, label="2 blind holdouts")
    ax.scatter([0.0], [0.0], marker="s", s=56, facecolor="white",
               edgecolor="#a21d28", linewidth=1.2, zorder=8, label="C endpoint")
    ax.scatter([t_qz], [0.0], marker="*", s=130, facecolor="#e5a547",
               edgecolor="#a21d28", linewidth=1.0, zorder=8, label="T endpoint")
    ax.set_ylabel(r"$Q_x$ ($\sqrt{\mathrm{amu}}\,\mathrm{\AA}$)", labelpad=5)
    ax.set_ylim(-0.023, 0.325)
    ax.text(-0.06, 1.03, "(a)", transform=ax.transAxes, va="bottom", fontsize=11,
            fontweight="normal", color="#182335")

    residual_ax.axhspan(-2, 2, color="#dbebea", alpha=0.53, zorder=0)
    residual_ax.axhline(0, color="#43556b", linewidth=0.8)
    for limit in (-2, 2):
        residual_ax.axhline(limit, color="#638282", linestyle=(0, (3, 3)), linewidth=0.8)
    residual_ax.scatter([row["qz"] for row in retrospective],
                        [row["observed_minus_predicted_meV_per_BTO"] for row in retrospective],
                        marker="o", s=36, facecolor="white", edgecolor="#687382",
                        linewidth=1.0, zorder=4, label="retrospective")
    residual_ax.scatter([row["qz"] for row in blind],
                        [row["observed_minus_predicted_meV_per_BTO"] for row in blind],
                        marker="D", s=55, facecolor="#e5a547", edgecolor="#583519",
                        linewidth=0.9, zorder=5, label="prospective")
    endpoint = points[-1]
    residual_ax.scatter([endpoint["qz"]], [endpoint["observed_minus_predicted_meV_per_BTO"]],
                        marker="*", s=120, facecolor="#e5a547", edgecolor="#a21d28",
                        linewidth=1.0, zorder=6, label="T endpoint (retrospective)")
    residual_ax.set_ylabel("DFT − model\n(meV/BTO)", labelpad=5)
    residual_ax.set_xlabel(r"$Q_z$ ($\sqrt{\mathrm{amu}}\,\mathrm{\AA}$)")
    residual_ax.set_ylim(-2.7, 2.7)
    residual_ax.text(-0.06, 1.03, "(b)", transform=residual_ax.transAxes,
                     va="bottom", fontsize=11, fontweight="normal", color="#182335")
    ax.set_xlim(-0.025, 1.225)
    ax.tick_params(top=True, right=True, direction="in", labelbottom=False, labelsize=9)
    residual_ax.tick_params(top=True, right=True, direction="in", labelsize=9)
    color_ax.tick_params(labelsize=8.5)
    handles, labels = ax.get_legend_handles_labels()
    legend = fig.legend(handles, labels, loc="outside lower center", ncol=3,
                        fontsize=8.5, columnspacing=1.1, handletextpad=0.4)
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_alpha(0.88)
    legend.get_frame().set_edgecolor("#9aa7b3")

    FIGURES.mkdir(parents=True, exist_ok=True)
    for suffix, options in (("pdf", {}), ("svg", {}), ("png", {"dpi": 300}),
                            ("tiff", {"dpi": 600})):
        fig.savefig(FIGURES / f"{BASENAME}.{suffix}", facecolor="white", **options)
    plt.close(fig)
    # Matplotlib's multiline SVG paths end lines with a cosmetic space.
    # Normalize only that generated whitespace so the vector artifact passes
    # the repository's diff check without changing any path coordinates.
    svg_path = FIGURES / f"{BASENAME}.svg"
    svg_path.write_text("\n".join(line.rstrip() for line in
                                  svg_path.read_text(encoding="utf-8").splitlines()) + "\n",
                        encoding="utf-8")
    with (FIGURES / f"{BASENAME}_points_source_data.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(points[0]))
        writer.writeheader()
        writer.writerows(points)
    with (FIGURES / f"{BASENAME}_contour_source_data.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("Qz_sqrt_amu_A", "Qx_sqrt_amu_A", "model_E_minus_C_meV_per_BTO"))
        for (z, x), valid, energy in zip(coords, inside, model.filled(np.nan).ravel()):
            if valid:
                writer.writerow((float(z), float(x), float(energy)))
    qa["script_sha256"] = sha256(Path(__file__))
    qa["source_data_sha256"] = {suffix: sha256(FIGURES / f"{BASENAME}_{suffix}_source_data.csv")
                                 for suffix in ("points", "contour")}
    qa["exports_sha256"] = {suffix: sha256(FIGURES / f"{BASENAME}.{suffix}")
                            for suffix in ("pdf", "svg", "png", "tiff")}
    qa["contour_grid_inside_measured_hull"] = int(inside.sum())
    qa["display_interpolation"] = "even sixth-order mode polynomial; no new DFT between marked points"
    (FIGURES / f"{BASENAME}_qa.json").write_text(
        json.dumps(qa, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": qa["status"], "max_blind_error_meV_per_BTO":
                      qa["maximum_absolute_prospective_error_meV_per_BTO"],
                      "T_QzQxQy": qa["T_qz_qx_qy_sqrt_amu_A"],
                      "figure": str(FIGURES / f"{BASENAME}.pdf")}, ensure_ascii=False))


if __name__ == "__main__":
    render()
