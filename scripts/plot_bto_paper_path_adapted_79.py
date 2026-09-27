"""Plot the audited 79-point frozen BTO path-adapted Gamma plane.

Q1 combines cubic unstable modes 0--2; Q2 combines stable modes 6--8.
The plotted variable-cell path is a projection, not a MEP on this frozen cut.
The separate BTO two-soft-mode conditional pilot uses a different plane.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
from scipy.interpolate import CloughTocher2DInterpolator, PchipInterpolator


ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "paper" / "VARNEB_CPC" / "figures"
STEM = "bto_frozen_path_adapted_79"
GRID = FIGURES / f"{STEM}_grid_source_data.csv"
PATH = FIGURES / f"{STEM}_path_source_data.csv"
SOURCE_QA = FIGURES / f"{STEM}_source_qa.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_data() -> tuple[list[dict], list[dict], dict]:
    with GRID.open(encoding="utf-8-sig", newline="") as handle:
        nodes = list(csv.DictReader(handle))
    with PATH.open(encoding="utf-8-sig", newline="") as handle:
        images = list(csv.DictReader(handle))
    qa = json.loads(SOURCE_QA.read_text(encoding="utf-8"))
    if (len(nodes) != 79 or len(images) != 7 or qa["n_grid_points"] != 79
            or qa["n_path_images_total"] != 7
            or qa["claim_scope"] != "frozen_cubic_atomic_mode_slice_plus_variable_cell_path_shadow_not_barrier"
            or qa["contour_claim_status"] != "exploratory_not_globally_validated"):
        raise ValueError("79-point frozen-plane or seven-image path contract changed")
    if {row["sample_group"] for row in nodes} != {
        "original_25_DFT", "adaptive_38_DFT", "nonuniform_refinement_16_DFT",
    }:
        raise ValueError("unexpected real-DFT sample groups")
    for index, row in enumerate(images):
        if int(row["image"]) != index:
            raise ValueError("path image order differs from the audited seven-image chain")
        if abs(float(row["eta_xx"]) - float(row["eta_yy"])) > 1e-8:
            raise ValueError("the archived path no longer has eta_xx = eta_yy")
    if abs(float(images[-1]["H_minus_C_eV_per_BTO"])) > 1e-8:
        raise ValueError("the cubic final endpoint is not the declared energy zero")
    return nodes, images, qa


def values(rows: list[dict], key: str) -> np.ndarray:
    array = np.asarray([float(row[key]) for row in rows], dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"non-finite source data: {key}")
    return array


def interpolation_diagnostics(nodes: list[dict]) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    q = np.column_stack([values(nodes, "Q1_sqrt_amu_A"), values(nodes, "Q2_sqrt_amu_A")])
    energy = values(nodes, "E_minus_C_eV_per_BTO")
    if len({tuple(point) for point in q}) != 79:
        raise ValueError("sampled mode coordinates are not unique")
    refinement = np.asarray([row["sample_group"] == "nonuniform_refinement_16_DFT"
                             for row in nodes])
    model63 = CloughTocher2DInterpolator(q[~refinement], energy[~refinement])
    preregistered_errors = 1000 * (np.asarray(model63(q[refinement])) - energy[refinement])
    loo_errors = []
    for index in np.flatnonzero(refinement):
        keep = np.ones(len(nodes), dtype=bool)
        keep[index] = False
        estimate = float(CloughTocher2DInterpolator(q[keep], energy[keep])(q[index]))
        loo_errors.append(1000 * (estimate - energy[index]))
    x = np.linspace(-1.2, 1.2, 401)
    y = np.linspace(-0.3, 0.3, 201)
    xx, yy = np.meshgrid(x, y)
    interpolated = np.asarray(CloughTocher2DInterpolator(q, energy)(xx, yy))
    if not np.all(np.isfinite(interpolated)):
        raise ValueError("the final display interpolation has non-finite pixels")
    diagnostics = {
        "preregistered_16_point_cubic_max_abs_error_meV_per_BTO": float(np.max(np.abs(preregistered_errors))),
        "new_16_point_leave_one_out_cubic_max_abs_error_meV_per_BTO": float(np.max(np.abs(loo_errors))),
        "global_cubic_downward_overshoot_meV_per_BTO": float(1000 * (np.min(energy) - np.min(interpolated))),
        "sampled_minimum_eV_per_BTO": float(np.min(energy)),
        "display_interpolation_minimum_eV_per_BTO": float(np.min(interpolated)),
        "display_pixel_grid": [401, 201],
        "n_real_DFT_points": 79,
        "n_display_pixels_not_DFT": int(interpolated.size),
    }
    return xx, yy, interpolated, diagnostics


def plot(nodes: list[dict], images: list[dict], xx: np.ndarray, yy: np.ndarray,
         interpolated: np.ndarray) -> list[dict]:
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8.8, "axes.labelsize": 9.5, "xtick.labelsize": 8.6,
        "ytick.labelsize": 8.6, "axes.linewidth": 0.85,
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "axes.spines.top": True, "axes.spines.right": True,
    })
    fig = plt.figure(figsize=(7.2, 4.25))  # CPC double-column width: about 183 mm.
    grid = fig.add_gridspec(2, 2, left=0.105, right=0.985, bottom=0.235, top=0.885,
                            width_ratios=[1.72, 1.0], wspace=0.46, hspace=0.38)
    ax_map = fig.add_subplot(grid[:, 0])
    ax_energy = fig.add_subplot(grid[0, 1])
    ax_strain = fig.add_subplot(grid[1, 1], sharex=ax_energy)

    norm = TwoSlopeNorm(vmin=-0.05, vcenter=0.0, vmax=0.08)
    fill = ax_map.contourf(xx, yy, interpolated, levels=np.linspace(-0.05, 0.08, 104),
                           cmap="RdBu_r", norm=norm, extend="both")
    ax_map.contour(xx, yy, interpolated,
                   levels=np.arange(-0.04, 0.061, 0.02), colors="#444B54",
                   linewidths=0.5, alpha=0.8)
    ax_map.scatter(values(nodes, "Q1_sqrt_amu_A"), values(nodes, "Q2_sqrt_amu_A"),
                   s=16, facecolors="white", edgecolors="#405260", linewidths=0.75,
                   zorder=4, label="79 computed DFT points")

    image_index = values(images, "image")
    dense_index = np.linspace(0.0, 6.0, 301)
    guide_fields = {
        key: PchipInterpolator(image_index, values(images, key))(dense_index)
        for key in ("Q1_sqrt_amu_A", "Q2_sqrt_amu_A", "H_minus_C_eV_per_BTO",
                    "eta_xx", "eta_zz")
    }
    ax_map.plot(guide_fields["Q1_sqrt_amu_A"], guide_fields["Q2_sqrt_amu_A"],
                color="#22272B", lw=1.55, zorder=5, label="VCNEB projection only")
    ax_map.scatter(values(images, "Q1_sqrt_amu_A"), values(images, "Q2_sqrt_amu_A"),
                   s=33, facecolors="#F6C65E", edgecolors="#242B30", linewidths=0.9,
                   zorder=6)
    ax_map.text(0.02, -0.045, "C", fontsize=8.7, ha="left", va="center",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 1.5})
    ax_map.text(1.14, 0.253, "T", fontsize=8.7, ha="right", va="center",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 1.5})
    # The T endpoint projects to Q1=1.2043, just beyond the sampled 1.2 grid.
    # Leave a narrow white margin instead of clipping it or extrapolating DFT.
    ax_map.set(xlim=(-1.23, 1.23), ylim=(-0.31, 0.31),
               xlabel=r"$Q_1$ ($\sqrt{\mathrm{amu}}$ Å)",
               ylabel=r"$Q_2$ ($\sqrt{\mathrm{amu}}$ Å)")
    ax_map.legend(loc="upper center", bbox_to_anchor=(0.52, 1.16), ncol=2,
                  fontsize=7.3, frameon=True, facecolor="white", framealpha=0.86,
                  edgecolor="#7D868C", handlelength=1.5, borderpad=0.3)

    ax_energy.plot(dense_index, guide_fields["H_minus_C_eV_per_BTO"],
                   color="#0F4D92", lw=1.8, label="PCHIP guide only")
    ax_energy.scatter(image_index, values(images, "H_minus_C_eV_per_BTO"),
                      s=25, color="#0F4D92", edgecolors="white", linewidths=0.5,
                      zorder=4, label="7 computed images")
    ax_energy.axhline(0, color="#8A8F94", ls="--", lw=0.75, zorder=0)
    ax_energy.set_ylabel(r"$H-H_C$ (eV/BTO)")
    ax_energy.legend(loc="upper left", fontsize=7.1, frameon=True, facecolor="white",
                     framealpha=0.86, edgecolor="#7D868C", borderpad=0.3)
    ax_energy.tick_params(labelbottom=False)

    ax_strain.plot(dense_index, guide_fields["eta_xx"], color="#3278BD", lw=1.6,
                   label=r"$\eta_{xx}=\eta_{yy}$")
    ax_strain.plot(dense_index, guide_fields["eta_zz"], color="#B84445", lw=1.6,
                   label=r"$\eta_{zz}$")
    ax_strain.scatter(image_index, values(images, "eta_xx"), s=20, marker="s",
                      color="#3278BD", zorder=4)
    ax_strain.scatter(image_index, values(images, "eta_zz"), s=22,
                      color="#B84445", zorder=4)
    ax_strain.axhline(0, color="#8A8F94", ls="--", lw=0.75, zorder=0)
    ax_strain.set(xlabel="VCNEB image (T $\\rightarrow$ C)",
                  ylabel="Reference-cell strain")
    ax_strain.legend(loc="upper right", fontsize=7.7, frameon=True, facecolor="white",
                     framealpha=0.86, edgecolor="#7D868C", borderpad=0.3)

    for ax, label in ((ax_map, "a"), (ax_energy, "b"), (ax_strain, "c")):
        ax.tick_params(direction="in", top=True, right=True, width=0.75, length=3.0)
        ax.text(-0.15, 1.045, f"({label})", transform=ax.transAxes, fontsize=10.8,
                fontweight="normal", ha="left", va="bottom")
    for ax in (ax_energy, ax_strain):
        ax.set_xlim(-0.2, 6.2)
        ax.set_xticks(np.arange(7))

    colorbar_ax = fig.add_axes([0.105, 0.095, 0.515, 0.025])
    colorbar = fig.colorbar(fill, cax=colorbar_ax, orientation="horizontal",
                            ticks=[-0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08])
    colorbar.set_label(r"Frozen $E-E_C$ (eV/BTO)", fontsize=8.8)
    colorbar.ax.tick_params(labelsize=8.1, width=0.7, length=2.5)

    paths = []
    for suffix, dpi in ((".svg", 600), (".pdf", 600), (".png", 600)):
        destination = FIGURES / f"{STEM}{suffix}"
        fig.savefig(destination, dpi=dpi, facecolor="white")
        if suffix == ".svg":
            # Matplotlib SVG path wraps have cosmetic trailing spaces.
            destination.write_text(
                "\n".join(line.rstrip() for line in destination.read_text(encoding="utf-8").splitlines()) + "\n",
                encoding="utf-8",
            )
        paths.append(destination)
    plt.close(fig)

    guide_path = FIGURES / f"{STEM}_visual_guide_not_DFT.csv"
    with guide_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["image_parameter_not_DFT", *guide_fields])
        for index, parameter in enumerate(dense_index):
            writer.writerow([parameter, *(guide_fields[key][index] for key in guide_fields)])
    paths.append(guide_path)
    return [{"file": path.name, "sha256": sha256(path)} for path in paths]


def main() -> None:
    nodes, images, source_qa = read_data()
    xx, yy, interpolated, diagnostics = interpolation_diagnostics(nodes)
    generated = plot(nodes, images, xx, yy, interpolated)
    report = {
        "status": "exploratory_frozen_soft_stable_plane_not_conditional_PES_or_barrier",
        "claim": "one unstable plus one stable cubic Gamma mode axis captures the path shadow in a frozen C-cell slice",
        "axis_mode_groups_zero_based": [[0, 1, 2], [6, 7, 8]],
        "calculator": "ABACUS PBE 100 Ry Ba/Ti/O 10 au DZP; electronic k mesh 4x4x4",
        "source_sha256": {path.name: sha256(path) for path in (GRID, PATH, SOURCE_QA)},
        "source_qa_contour_claim_status": source_qa["contour_claim_status"],
        "source_qa_holdout_values_meV_per_BTO": source_qa["smooth_contour_validation"],
        "recomputed_interpolation": diagnostics,
        "maximum_path_atomic_projection_residual_sqrt_amu_A": source_qa[
            "maximum_atomic_projection_residual_sqrt_amu_A"],
        "sampled_Q1_bounds_sqrt_amu_A": [
            float(np.min(values(nodes, "Q1_sqrt_amu_A"))),
            float(np.max(values(nodes, "Q1_sqrt_amu_A"))),
        ],
        "T_endpoint_Q1_outside_sampled_range_sqrt_amu_A": max(
            0.0, float(images[0]["Q1_sqrt_amu_A"]) -
            float(np.max(values(nodes, "Q1_sqrt_amu_A"))),
        ),
        "n_actual_path_images": 7,
        "n_PCHIP_guide_samples_not_DFT": 301,
        "generated": generated,
        "limitations": [
            "The contour is fixed cubic-cell static E, while the path uses variable-cell H.",
            "The path line and smooth guide samples are projections/interpolations, not additional DFT.",
            "This soft-plus-stable plane is different from the two-soft-mode conditional pilot.",
            "The T endpoint Q1 projection exceeds the sampled frozen grid by 0.00433 sqrt(amu) A; no contour extrapolation is drawn there.",
            "No sampled saddle, mode-constrained barrier, global contour error bound, or finite-T free energy is claimed.",
        ],
    }
    qa_path = FIGURES / f"{STEM}_qa.json"
    qa_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"[DONE] {qa_path}")
    print("[DONE] 79 real frozen DFT nodes, 7 real VCNEB images; no new calculation")


if __name__ == "__main__":
    main()
