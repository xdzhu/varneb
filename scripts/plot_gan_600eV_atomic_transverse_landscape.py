"""Plot the measured GaN central transverse response without hiding its scale.

The left-hand contour is H(s, q_atom) - H(s, 0), not the absolute barrier and
not a globally relaxed two-phonon potential-energy surface.  The right-hand
panel supplies the absolute 45.7-GPa VCNEB enthalpy profile for context.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase.units import GPa
from matplotlib.lines import Line2D

from scripts.plot_gan_600eV_atomic_dense_surface import (
    BENCHMARKS,
    CENTRAL_IMAGES,
    DEFAULT_REPORT,
    DEFAULT_TRAJECTORY,
    INK,
    BLUE,
    ORANGE,
    TRANSVERSE_Q_A,
    load,
    model_surface,
    sha256,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PREFIX = (
    ROOT / "paper/VARNEB_CPC/figures/gan_600eV_atomic_transverse_landscape"
)


def transverse_surface(report: dict, full_h: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Remove the path baseline from the independently gated absolute model."""
    s_mesh, q_mesh, absolute = model_surface(report, full_h)
    h0 = np.interp(
        s_mesh[0], report["arc_fraction_s"], report["path_enthalpy_eV_per_cell"]
    )
    transverse = absolute - (h0[None, :] - full_h[0]) * 500
    if not np.isfinite(transverse).all() or np.max(np.abs(transverse[100])) > 1e-9:
        raise ValueError("central transverse baseline is not zero")
    return s_mesh, q_mesh, transverse


def draw(report: dict, frames: list) -> plt.Figure:
    full_h = np.asarray([
        frame.get_potential_energy() + 45.7 * GPa * frame.get_volume()
        for frame in frames
    ])
    path_mev = (full_h - full_h[0]) * 500
    full_arc = np.asarray(json.loads(
        (BENCHMARKS / "gan_600eV_atomic_tube_refinement_20260928.json").read_text()
    )["arc_fraction_s"], dtype=float)
    if not np.allclose(full_arc[list(CENTRAL_IMAGES)], report["arc_fraction_s"], atol=1e-12):
        raise ValueError("full-path arc fraction changed")
    s_mesh, q_mesh, transverse = transverse_surface(report, full_h)
    s = np.asarray(report["arc_fraction_s"], dtype=float)
    q = np.asarray(TRANSVERSE_Q_A, dtype=float)

    fig = plt.figure(figsize=(7.2, 4.5), facecolor="white")
    ax_a = fig.add_axes((0.095, 0.305, 0.535, 0.575))
    ax_b = fig.add_axes((0.715, 0.305, 0.245, 0.575))
    cax = fig.add_axes((0.145, 0.125, 0.445, 0.027))
    fig.text(0.095, 0.935, "(a)", fontsize=13, fontweight="normal", color=INK)
    fig.text(0.715, 0.935, "(b)", fontsize=13, fontweight="normal", color=INK)

    # The fitted local minimum is slightly below the q=0 path at some s.
    # Include that small negative range rather than leaving white gaps that
    # could be mistaken for missing DFT data.
    levels = np.concatenate(([-0.5, 0.0], np.arange(2.0, 22.0 + 0.01, 2.0)))
    contour = ax_a.contourf(s_mesh, q_mesh, transverse, levels=levels,
                            cmap="cividis", extend="max")
    ax_a.contour(s_mesh, q_mesh, transverse, levels=levels[3::2],
                 colors=INK, alpha=0.28, linewidths=0.55)
    grid_s, grid_q = np.meshgrid(s, q, indexing="ij")
    ax_a.scatter(grid_s.ravel(), grid_q.ravel(), s=14, marker="o",
                 facecolors="white", edgecolors=INK, linewidths=0.55, zorder=7)
    ax_a.plot(s, np.zeros_like(s), color="white", lw=2.8, zorder=5)
    ax_a.plot(s, np.zeros_like(s), color=INK, lw=1.45, zorder=6)
    ax_a.scatter(s[15 - CENTRAL_IMAGES[0]], 0, s=135, marker="*",
                 color=ORANGE, edgecolor="white", linewidth=0.7, zorder=8)
    ax_a.set_xlim(s[0] - 0.004, s[-1] + 0.004)
    ax_a.set_ylim(-0.054, 0.054)
    ax_a.set_yticks([-0.05, -0.025, 0, 0.025, 0.05])
    ax_a.set_xlabel(r"VCNEB arc fraction $s$")
    ax_a.set_ylabel(r"Atomic transverse coordinate $q_{\perp}$ (Å)")
    ax_a.tick_params(direction="in", top=True, right=True, length=4, pad=4)
    colorbar = fig.colorbar(contour, cax=cax, orientation="horizontal")
    colorbar.ax.xaxis.set_label_position("top")
    colorbar.set_label(r"$H(s,q_{\perp})-H(s,0)$ (meV/GaN)", fontsize=9.6, labelpad=4)
    colorbar.ax.tick_params(labelsize=8.8, direction="in", length=3)
    fig.legend(handles=[
        Line2D([], [], marker="o", linestyle="none", markerfacecolor="white",
               markeredgecolor=INK, markersize=5, label="DFT grid (90)"),
        Line2D([], [], color=INK, lw=1.5, label="VCNEB centerline"),
    ], loc="center", bbox_to_anchor=(0.370, 0.930), ncol=2,
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
    ], loc="center", bbox_to_anchor=(0.835, 0.155), fontsize=8.4,
        frameon=True, facecolor="white", edgecolor="#8B969C",
        framealpha=0.84, borderpad=0.35)
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--trajectory", type=Path, default=DEFAULT_TRAJECTORY)
    parser.add_argument("--output-prefix", type=Path, default=DEFAULT_PREFIX)
    args = parser.parse_args()
    prefix = args.output_prefix
    paths = {suffix: prefix.with_name(prefix.name + suffix)
             for suffix in ("_qa.json", ".pdf", ".svg", ".tiff", ".png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("GaN transverse figure output exists; use a fresh prefix")
    report, frames = load(args.report, args.trajectory)
    full_h = np.asarray([
        frame.get_potential_energy() + 45.7 * GPa * frame.get_volume()
        for frame in frames
    ])
    _, _, transverse = transverse_surface(report, full_h)
    fig = draw(report, frames)
    prefix.parent.mkdir(parents=True, exist_ok=True)
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
        "status": "GaN_600eV_central_atomic_transverse_landscape",
        "core_conclusion": "The audited GaN path has a small but resolved central atomic-transverse enthalpy response, distinct from its absolute barrier.",
        "backend": "Python/matplotlib",
        "archetype": "quantitative grid with contour hero and path context",
        "electronic_contract": "VASP PBE/Ga_d+N/600 eV/Gamma 8x8x6 at 45.7 GPa",
        "n_measured_coordinates": 90,
        "transverse_energy_min_max_meV_per_GaN": [float(transverse.min()), float(transverse.max())],
        "s_domain": [float(report["arc_fraction_s"][0]), float(report["arc_fraction_s"][-1])],
        "q_atom_domain_A": [-0.05, 0.05],
        "maximum_s_holdout_error_meV_per_GaN": report["maximum_prospective_s_holdout_error_meV_per_GaN"],
        "maximum_q_halfstep_error_meV_per_GaN": report["maximum_q_halfstep_error_meV_per_GaN"],
        "interpolation": "Same audited linear-s/quadratic-q model as the absolute-enthalpy contour; only the path baseline has been subtracted.",
        "claim_limit": "Frozen central (s,q_atom) cut only; not a global two-phonon PES or certified TS.",
        "source_sha256": {
            "raw_audit": sha256(args.report),
            "trajectory": sha256(args.trajectory),
            "measured_csv": sha256(ROOT / "paper/VARNEB_CPC/figures/gan_600eV_atomic_dense_surface_source_data.csv"),
            "plotter": sha256(Path(__file__)),
        },
    }
    paths["_qa.json"].write_text(json.dumps(qa, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": qa["status"], "n_measured_coordinates": 90}))


if __name__ == "__main__":
    main()
