"""Plot audited GaN split VCNEB profile beside its local 289-point mode cut.

Figure contract: the shared candidate remains the highest point of the
converged stitched chain and the path approaches its local negative-curvature
direction. The local cut contains only the shared center as an actual VCNEB
image. Dashed rays in that panel are clipped *linear coordinate connections*
to the nearest outside images, not extra DFT results or a relaxed local MEP.

Archetype: two-panel quantitative grid; Python/matplotlib only. Output includes
editable PDF/SVG, PNG preview, path source CSV, and a claim-limit QA record.
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
from ase.io.trajectory import Trajectory
from ase.units import GPa

from scripts.plot_gan_600eV_local_joint_cut17 import surface
from vcneb import VCNEB
from vcneb.joint_curvature import symmetric_strain_basis


CELL_SCALE_A = 3.3982714330050063
BLUE = "#245F89"
ORANGE = "#C26A3D"
INK = "#24323C"
GRID = "#D6DEE2"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def final_frames(path: Path, n: int) -> list:
    with Trajectory(str(path), "r") as trajectory:
        if len(trajectory) < n or len(trajectory) % n:
            raise ValueError(f"{path.name}: missing complete final chain")
        return [trajectory[i] for i in range(len(trajectory) - n, len(trajectory))]


def project(image, center, modes: np.ndarray) -> tuple[float, float, float]:
    """Project exact finite image displacement on center joint Hessian modes."""
    cell = center.cell.array
    delta = image.get_scaled_positions(wrap=False) - center.get_scaled_positions(wrap=False)
    delta -= np.rint(delta)
    atoms = (delta @ cell).ravel()
    deform = np.linalg.solve(cell, image.cell.array)
    strain = 0.5 * (deform + deform.T) - np.eye(3)
    strain_coords = CELL_SCALE_A * np.einsum("ij,aij->a", strain, symmetric_strain_basis())
    displacement = np.concatenate([atoms, strain_coords])
    q = modes[:, :2].T @ displacement
    residual = np.linalg.norm(displacement - modes[:, :2] @ q)
    return float(q[0]), float(q[1]), float(residual)


def clipped_ray(neighbor: tuple[float, float], limits: tuple[float, float]) -> tuple[float, float]:
    """Clip a center-to-neighbor straight segment to a rectangular plot cut."""
    u, v = neighbor
    u_limit, v_limit = limits
    factors = [1.0]
    if abs(u) > 0:
        factors.append(u_limit / abs(u))
    if abs(v) > 0:
        factors.append(v_limit / abs(v))
    t = min(factors)
    return t * u, t * v


def evidence_rows(audit_path: Path, evidence_dir: Path, hessian_path: Path,
                  surface_csv: Path) -> tuple[list[dict], list[dict], dict]:
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if (audit.get("status") != "raw_per_image_audit_passed"
            or audit.get("pressure_GPa") != 45.7
            or audit.get("n_total_stitched_images") != 29
            or audit.get("n_total_raw_interior_OUTCARs_audited") != 26):
        raise ValueError("two-segment raw audit is missing or incomplete")
    segments = audit["segments"]
    left_path, right_path = evidence_dir / "left.traj", evidence_dir / "right.traj"
    if (sha256(left_path) != segments["left"]["trajectory_sha256"]
            or sha256(right_path) != segments["right"]["trajectory_sha256"]):
        raise ValueError("segment trajectory differs from raw audit")
    left = final_frames(left_path, 16)
    right = final_frames(right_path, 14)
    if (np.max(np.abs(left[-1].cell.array - right[0].cell.array)) > 1e-10
            or np.max(np.abs(left[-1].positions - right[0].positions)) > 1e-10):
        raise ValueError("two segments do not share an exact center image")
    stitched = left + right[1:]
    with np.load(hessian_path, allow_pickle=False) as archive:
        modes = np.asarray(archive["eigenvectors"], dtype=float)
        eigenvalues = np.asarray(archive["eigenvalues"], dtype=float)
    if (modes.shape != (18, 15) or not eigenvalues[0] < 0 < eigenvalues[1]
            or not np.allclose(modes.T @ modes, np.eye(15), atol=1e-8)):
        raise ValueError("local index-one joint-mode basis failed")
    center = stitched[15]
    coords = VCNEB(stitched, pressure=45.7 * GPa, mic=True).reaction_coordinate()
    coords = coords / coords[-1]
    h0 = stitched[0].get_potential_energy() + 45.7 * GPa * stitched[0].get_volume()
    path_rows = []
    for index, (image, arc) in enumerate(zip(stitched, coords)):
        e = float(image.get_potential_energy())
        h = e + 45.7 * GPa * image.get_volume()
        q_u, q_v, residual = project(image, center, modes)
        expected = (segments["left"]["records"][index] if index <= 15 else
                    segments["right"]["records"][index - 15])
        if abs(h - expected["enthalpy_eV_cell"]) > 1e-6:
            raise ValueError(f"image {index} enthalpy differs from raw audit")
        path_rows.append({
            "stitched_index": index,
            "segment": "left" if index <= 15 else "right",
            "arc_fraction": float(arc),
            "enthalpy_meV_per_GaN_relative_B4": 1000 * (h - h0) / 2,
            "q_u_A": q_u, "q_v_A": q_v, "two_mode_residual_A": residual,
            "inside_measured_local_rectangle": abs(q_u) <= 0.02 and abs(q_v) <= 0.0125,
        })
    if [r["stitched_index"] for r in path_rows if r["inside_measured_local_rectangle"]] != [15]:
        raise ValueError("the local-cut occupancy changed; revisit overlay semantics")
    with surface_csv.open(newline="", encoding="utf-8") as handle:
        surface_rows = list(csv.DictReader(handle))
    if (len(surface_rows) != 289 or
            len({(float(r["q_u_A"]), float(r["q_v_A"])) for r in surface_rows}) != 289):
        raise ValueError("289-node measured cut is incomplete")
    for row in surface_rows:
        for key in ("q_u_A", "q_v_A", "measured_delta_H_meV_per_GaN"):
            row[key] = float(row[key])
    return path_rows, surface_rows, audit


def draw(path_rows: list[dict], surface_rows: list[dict]) -> plt.Figure:
    uu, vv, z = surface(surface_rows)
    levels = [-0.50, -0.40, -0.30, -0.20, -0.10, -0.05, 0,
              0.025, 0.050, 0.075, 0.10, 0.15]
    if np.nanmin(z) < levels[0] or np.nanmax(z) > levels[-1]:
        raise ValueError("surface color range exceeded")
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "font.size": 10, "axes.labelsize": 10.5, "xtick.labelsize": 9.2,
        "ytick.labelsize": 9.2, "axes.linewidth": 0.9,
        "axes.spines.top": True, "axes.spines.right": True,
    })
    fig = plt.figure(figsize=(7.25, 3.75), facecolor="white")
    ax_map = fig.add_axes((0.105, 0.305, 0.40, 0.53))
    ax_profile = fig.add_axes((0.625, 0.305, 0.32, 0.53))
    cax = fig.add_axes((0.148, 0.095, 0.31, 0.025))
    fig.text(0.105, 0.925, "(a)", fontsize=13, fontweight="normal", color=INK)
    fig.text(0.625, 0.925, "(b)", fontsize=13, fontweight="normal", color=INK)
    contour = ax_map.contourf(
        uu, vv, z, levels=levels, cmap="PuOr_r",
        norm=TwoSlopeNorm(vmin=-0.50, vcenter=0, vmax=0.15),
    )
    ax_map.contour(uu, vv, z, levels=[0], colors=[INK], linewidths=0.8,
                   linestyles="--")
    ax_map.scatter([r["q_u_A"] for r in surface_rows],
                   [r["q_v_A"] for r in surface_rows], s=11.0,
                   facecolors="white", edgecolors=INK, linewidths=0.42,
                   alpha=0.9, zorder=4)
    for neighbor, color in ((path_rows[14], BLUE), (path_rows[16], ORANGE)):
        end = clipped_ray((neighbor["q_u_A"], neighbor["q_v_A"]),
                          (0.02, 0.0125))
        ax_map.plot([0, end[0]], [0, end[1]], color=color, lw=2.0,
                    linestyle=(0, (4.2, 2.8)), zorder=5)
    ax_map.scatter(0, 0, s=170, marker="*", color=INK,
                   edgecolors="white", linewidths=0.75, zorder=6)
    ax_map.set(xlim=(-0.0208, 0.0208), ylim=(-0.0131, 0.0131),
               xticks=[-0.02, -0.01, 0, 0.01, 0.02],
               yticks=[-0.0125, 0, 0.0125],
               xlabel=r"Unstable joint coordinate $q_u$ (Å)",
               ylabel=r"Stable joint coordinate $q_v$ (Å)")
    ax_map.tick_params(direction="in", top=True, right=True, length=4, pad=4)
    cb = fig.colorbar(contour, cax=cax, orientation="horizontal",
                      ticks=[-0.4, -0.2, 0, 0.1])
    cb.ax.xaxis.set_label_position("top")
    cb.set_label(r"Local frozen-cut $\Delta H$ (meV/GaN)", fontsize=9.2, labelpad=4)
    cb.ax.tick_params(labelsize=8.6, direction="in", length=3)

    arc = [r["arc_fraction"] for r in path_rows]
    h = [r["enthalpy_meV_per_GaN_relative_B4"] for r in path_rows]
    ax_profile.plot(arc[:16], h[:16], color=BLUE, lw=2.0, marker="o", ms=4.1,
                    markeredgecolor="white", markeredgewidth=0.4, zorder=3)
    ax_profile.plot(arc[15:], h[15:], color=ORANGE, lw=2.0, marker="o", ms=4.1,
                    markeredgecolor="white", markeredgewidth=0.4, zorder=3)
    ax_profile.scatter(arc[15], h[15], s=150, marker="*", color=INK,
                       edgecolors="white", linewidths=0.7, zorder=5)
    ax_profile.axhline(0, color=GRID, lw=0.9, zorder=0)
    ax_profile.set(xlim=(-0.025, 1.025), ylim=(-24, 376),
                   xticks=[0, 0.25, 0.50, 0.75, 1.0],
                   yticks=[0, 100, 200, 300],
                   xlabel="Full-chain generalized arc fraction",
                   ylabel=r"$H-H_{\mathrm{B4}}$ (meV/GaN)")
    ax_profile.tick_params(direction="in", top=True, right=True, length=4, pad=4)
    fig.legend(handles=[
        Line2D([0], [0], color=BLUE, lw=2, marker="o", ms=5, label="B4 → candidate"),
        Line2D([0], [0], color=ORANGE, lw=2, marker="o", ms=5, label="Candidate → B1"),
    ], loc="upper left", bbox_to_anchor=(0.625, 0.875), ncol=2, fontsize=7.8,
        frameon=True, facecolor="white", edgecolor="#8A969D", framealpha=0.88,
        handlelength=1.6, handletextpad=0.35, columnspacing=0.55, borderpad=0.28)
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--hessian", type=Path, required=True)
    parser.add_argument("--surface-csv", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()
    prefix = args.output_prefix
    suffixes = (".pdf", ".svg", ".png", "_source_data.csv", "_qa.json")
    paths = {suffix: prefix.with_name(prefix.name + suffix) for suffix in suffixes}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("output exists; choose a new prefix")
    evidence_dir = args.evidence_dir
    rows, grid, audit = evidence_rows(evidence_dir / "raw_audit.json",
                                      evidence_dir, args.hessian, args.surface_csv)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with paths["_source_data.csv"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    fig = draw(rows, grid)
    fig.savefig(paths[".pdf"], bbox_inches="tight")
    fig.savefig(paths[".svg"], bbox_inches="tight")
    svg_lines = paths[".svg"].read_text(encoding="utf-8").splitlines()
    paths[".svg"].write_text(
        "\n".join(line.rstrip() for line in svg_lines) + "\n", encoding="utf-8",
    )
    fig.savefig(paths[".png"], dpi=250, bbox_inches="tight")
    plt.close(fig)
    qa = {
        "status": "GaN_two_segment_29_image_profile_and_local_projection_visualized",
        "core_conclusion": "The same-contract stitched path has a fixed highest center and approaches its local unstable joint mode.",
        "archetype": "two-panel quantitative grid",
        "backend": "Python/matplotlib",
        "pressure_GPa": 45.7,
        "n_measured_local_DFT_nodes": len(grid),
        "n_stitched_VCNEB_images": len(rows),
        "VCNEB_images_inside_local_cut": [r["stitched_index"] for r in rows if r["inside_measured_local_rectangle"]],
        "local_path_line_semantics": "Dashed rays are clipped linear projections from shared center to adjacent converged images; only star is an actual VCNEB image inside the measured cut.",
        "claim_limit": audit["claim_limit"] + "; frozen local cut is not a relaxed full-path PES",
        "forward_barrier_meV_per_GaN": audit["forward_barrier_meV_per_GaN"],
        "reverse_barrier_meV_per_GaN": audit["reverse_barrier_meV_per_GaN"],
        "source_sha256": {
            "raw_audit": sha256(evidence_dir / "raw_audit.json"),
            "left_trajectory": sha256(evidence_dir / "left.traj"),
            "right_trajectory": sha256(evidence_dir / "right.traj"),
            "hessian": sha256(args.hessian),
            "measured_surface_csv": sha256(args.surface_csv),
            "path_source_csv": sha256(paths["_source_data.csv"]),
            "plotter": sha256(Path(__file__)),
        },
    }
    paths["_qa.json"].write_text(json.dumps(qa, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: qa[key] for key in (
        "status", "n_stitched_VCNEB_images", "VCNEB_images_inside_local_cut",
        "forward_barrier_meV_per_GaN", "reverse_barrier_meV_per_GaN")}))


if __name__ == "__main__":
    main()
