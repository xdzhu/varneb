"""Draw a *sampled*, non-interpolated GaN transverse enthalpy diagnostic.

Figure contract: the original 600-eV B4→B1 path has a barrier, while 28
same-input off-path statics resolve a central frozen-atomic transverse response
but do not validate a smooth two-dimensional contour. Panel (a) gives the
full path; panel (b) shows actual DFT coordinates only. This quantitative-grid
figure is exploratory and is not a free-energy surface or TS certificate.
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


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "benchmarks/numerical_integrity"
DEFAULT_FIRST = BASE / "gan_600eV_atomic_tube_central_20260928.json"
DEFAULT_REFINED = BASE / "gan_600eV_atomic_tube_refinement_20260928.json"
DEFAULT_COMPONENTS = BASE / "gan_600eV_atomic_tube_components_20260928.json"
DEFAULT_PREFIX = ROOT / "paper/VARNEB_CPC/figures/gan_600eV_atomic_tube_samples"
INK = "#263440"
BLUE = "#285D80"
ORANGE = "#C47D45"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "font.size": 10.2, "axes.labelsize": 10.8,
    "xtick.labelsize": 9.8, "ytick.labelsize": 9.8,
    "axes.linewidth": 0.9,
    "axes.spines.top": True, "axes.spines.right": True,
})


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(first_path: Path, refined_path: Path, components_path: Path) -> tuple[dict, dict, dict]:
    first = json.loads(first_path.read_text(encoding="utf-8"))
    refined = json.loads(refined_path.read_text(encoding="utf-8"))
    components = json.loads(components_path.read_text(encoding="utf-8"))
    if (first.get("status") != "GaN_600eV_central_atomic_transverse_tube_raw_audited"
            or refined.get("status")
            != "GaN_600eV_central_atomic_tube_targeted_refinement_raw_audited"
            or components.get("status") != "GaN_600eV_atomic_tube_even_odd_force_diagnostic"
            or refined["source_sha256"].get("first_audit") != sha256(first_path)
            or components["source_sha256"].get("first_raw_audit") != sha256(first_path)
            or components["source_sha256"].get("refinement_raw_audit") != sha256(refined_path)
            or refined.get("central_contour_gate_pass") is not False
            or components.get("contour_certified") is not False
            or refined.get("predeclared_LOO_gate_meV_per_GaN") != 1.0
            or len(first.get("cases", [])) != 16
            or len(refined.get("new_cases", [])) != 12):
        raise ValueError("GaN 600-eV sampled-figure provenance or claim gate changed")
    if (first["arc_fraction_s"] != refined["arc_fraction_s"]
            or first["path_enthalpy_eV_per_cell"] != refined["path_enthalpy_eV_per_cell"]):
        raise ValueError("audited 600-eV path changed")
    return first, refined, components


def source_rows(first: dict, refined: dict, first_sha: str, refined_sha: str) -> list[dict]:
    arc = first["arc_fraction_s"]
    h0 = first["path_enthalpy_eV_per_cell"]
    rows = []
    for index, (s, h) in enumerate(zip(arc, h0)):
        rows.append({
            "panel": "a", "kind": "archived_path", "image_index": index,
            "arc_fraction_s": s, "q_atom_A": "",
            "enthalpy_rel_B4_eV_per_GaN": (h - h0[0]) / 2,
            "offpath_excess_meV_per_GaN": "",
            "case": "", "raw_outcar_sha256": "",
            "source_report_sha256": first_sha,
        })
        if 5 <= index <= 22:
            rows.append({
                "panel": "b", "kind": "archived_path_center", "image_index": index,
                "arc_fraction_s": s, "q_atom_A": 0.0,
                "enthalpy_rel_B4_eV_per_GaN": (h - h0[0]) / 2,
                "offpath_excess_meV_per_GaN": 0.0,
                "case": "", "raw_outcar_sha256": "",
                "source_report_sha256": first_sha,
            })
    for cases, report_sha in ((first["cases"], first_sha),
                              (refined["new_cases"], refined_sha)):
        for case in cases:
            index = int(case["image_index"])
            q = float(case["q_atom_A"])
            h = float(case["enthalpy_eV_per_cell"])
            excess = (h - h0[index]) * 500
            if (index < 5 or index > 22 or not np.isclose(
                    excess, case["delta_enthalpy_meV_per_GaN"], atol=1e-8, rtol=0)):
                raise ValueError("off-path case inconsistent with 600-eV q=0 path")
            rows.append({
                "panel": "b", "kind": "offpath_static", "image_index": index,
                "arc_fraction_s": arc[index], "q_atom_A": q,
                "enthalpy_rel_B4_eV_per_GaN": (h - h0[0]) / 2,
                "offpath_excess_meV_per_GaN": excess,
                "case": case["case"], "raw_outcar_sha256": case["outcar_sha256"],
                "source_report_sha256": report_sha,
            })
    if (sum(row["kind"] == "offpath_static" for row in rows) != 28
            or sum(row["kind"] == "archived_path_center" for row in rows) != 18):
        raise ValueError("sampled central chart is incomplete")
    return rows


def draw(rows: list[dict]) -> plt.Figure:
    fig = plt.figure(figsize=(7.2, 4.6), facecolor="white")
    # Equal axis widths keep both panels' path coordinates aligned. The
    # colorbar sits beyond the right spine without shrinking panel (b).
    ax_a = fig.add_axes((0.112, 0.650, 0.726, 0.265))
    ax_b = fig.add_axes((0.112, 0.155, 0.726, 0.360), sharex=ax_a)
    cax = fig.add_axes((0.872, 0.174, 0.021, 0.324))
    fig.text(0.112, 0.944, "(a)", fontsize=13, fontweight="normal", color=INK)
    fig.text(0.112, 0.537, "(b)", fontsize=13, fontweight="normal", color=INK)

    path = sorted((row for row in rows if row["kind"] == "archived_path"),
                  key=lambda row: row["image_index"])
    s = np.array([row["arc_fraction_s"] for row in path], dtype=float)
    h = np.array([row["enthalpy_rel_B4_eV_per_GaN"] for row in path], dtype=float)
    peak = int(np.argmax(h))
    ax_a.axvspan(s[5], s[22], color="#E8EDF0", alpha=0.62, zorder=0)
    ax_a.plot(s, h, color=BLUE, lw=2.25, marker="o", markersize=3.3,
              markerfacecolor="white", markeredgewidth=0.8, zorder=2)
    ax_a.scatter([s[peak]], [h[peak]], s=68, color=ORANGE, edgecolor="white",
                 linewidth=0.8, zorder=4)
    ax_a.annotate(f"{h[peak]:.3f} eV/GaN", (s[peak], h[peak]),
                  xytext=(8, 13), textcoords="offset points", fontsize=9.5,
                  color=INK, ha="left")
    ax_a.set_xlim(-0.018, 1.018)
    ax_a.set_ylim(-0.05, 0.40)
    ax_a.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
    ax_a.set_ylabel(r"$H(s)-H_{\mathrm{B4}}$ (eV/GaN)")
    ax_a.tick_params(axis="both", which="both", direction="in", top=True,
                     right=True, length=4.0, pad=5)
    ax_a.tick_params(axis="x", labelbottom=False)

    center = sorted((row for row in rows if row["kind"] == "archived_path_center"),
                    key=lambda row: row["image_index"])
    ax_b.plot([row["arc_fraction_s"] for row in center],
              [0.0] * len(center), color="#697983", lw=1.45, zorder=1)
    ax_b.scatter([row["arc_fraction_s"] for row in center],
                 [0.0] * len(center), s=18, color="#697983", zorder=2)
    offpath = [row for row in rows if row["kind"] == "offpath_static"]
    cmap = plt.get_cmap("cividis")
    norm = plt.Normalize(vmin=0.0, vmax=21.0)
    scatter = None
    for amplitude, marker, size in ((0.05, "o", 81), (0.025, "s", 74)):
        subset = [row for row in offpath if np.isclose(abs(row["q_atom_A"]), amplitude)]
        scatter = ax_b.scatter(
            [row["arc_fraction_s"] for row in subset],
            [row["q_atom_A"] for row in subset],
            c=[row["offpath_excess_meV_per_GaN"] for row in subset],
            cmap=cmap, norm=norm, marker=marker, s=size,
            edgecolors="white", linewidths=0.95, zorder=4,
        )
    assert scatter is not None
    cbar = fig.colorbar(scatter, cax=cax, ticks=[0, 5, 10, 15, 20])
    cbar.set_label(r"$H(s,q)-H(s,0)$ (meV/GaN)", fontsize=9.8, labelpad=8)
    cbar.ax.tick_params(labelsize=9.2, direction="in", length=3)
    ax_b.set_ylim(-0.067, 0.070)
    ax_b.set_yticks([-0.05, -0.025, 0.0, 0.025, 0.05],
                    ["−0.050", "−0.025", "0", "+0.025", "+0.050"])
    ax_b.set_xticks(np.linspace(0, 1, 6))
    ax_b.set_xlabel("Path arc fraction $s$")
    ax_b.set_ylabel(r"Frozen atomic $q$ (Å)")
    ax_b.tick_params(axis="both", which="both", direction="in", top=True,
                     right=True, length=4.0, pad=5)
    handles = [
        Line2D([], [], color="#697983", marker="o", markersize=5, lw=1.3,
               label="$q=0$"),
        Line2D([], [], color=INK, marker="o", markersize=7, lw=0,
               markerfacecolor="white", label="$\pm0.05$ Å"),
        Line2D([], [], color=INK, marker="s", markersize=7, lw=0,
               markerfacecolor="white", label="$\pm0.025$ Å"),
    ]
    ax_b.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.005, 0.985),
                fontsize=8.9, frameon=True, facecolor="white",
                edgecolor="#8B969C", framealpha=0.86, borderpad=0.45,
                handlelength=0.9, handletextpad=0.45)
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-audit", type=Path, default=DEFAULT_FIRST)
    parser.add_argument("--refinement-audit", type=Path, default=DEFAULT_REFINED)
    parser.add_argument("--components", type=Path, default=DEFAULT_COMPONENTS)
    parser.add_argument("--output-prefix", type=Path, default=DEFAULT_PREFIX)
    args = parser.parse_args()
    first, refined, components = load(args.first_audit, args.refinement_audit,
                                      args.components)
    prefix = args.output_prefix
    paths = {suffix: prefix.with_name(prefix.name + suffix)
             for suffix in ("_source_data.csv", ".pdf", ".svg", ".tiff", ".png", "_qa.json")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("sampled figure output exists; use a new prefix")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    rows = source_rows(first, refined, sha256(args.first_audit),
                       sha256(args.refinement_audit))
    with paths["_source_data.csv"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    fig = draw(rows)
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
    path_h = np.array(first["path_enthalpy_eV_per_cell"], dtype=float)
    qa = {
        "status": "GaN_600eV_atomic_tube_discrete_samples_exploratory_figure",
        "core_conclusion": "Same-input GaN path and central off-path statics show a transverse response, not a certified smooth 2D contour.",
        "archetype": "quantitative grid; shared x axis, lower sampled-point hero",
        "backend": "Python/matplotlib",
        "final_size_mm": [182.88, 116.84],
        "panel_map": {
            "a": "29-image archived 600-eV variable-cell enthalpy path at 45.7 GPa",
            "b": "18 central q=0 path points and 28 raw off-path DFT statics; no surface interpolation",
        },
        "n_full_amplitude_offpath_statics": 24,
        "n_half_amplitude_offpath_statics": 4,
        "forward_barrier_eV_per_GaN": float((max(path_h) - path_h[0]) / 2),
        "reverse_barrier_eV_per_GaN": float((max(path_h) - path_h[-1]) / 2),
        "path_peak_image": int(np.argmax(path_h)),
        "max_linear_LOO_error_meV_per_GaN": refined[
            "maximum_leave_one_out_error_meV_per_GaN"],
        "predeclared_LOO_gate_meV_per_GaN": refined[
            "predeclared_LOO_gate_meV_per_GaN"],
        "contour_certified": False,
        "image_integrity": "Raw sample locations only; no smoothing, contour, masking, or rescaling of enthalpy values.",
        "raster_export_policy": (
            "PDF/SVG/PNG are tracked; the reproducible 600-dpi TIFF is generated "
            "locally and ignored by the repository's figure policy."
        ),
        "reviewer_risk": "Do not interpret discrete point colors as an interpolated PES or a TS certificate.",
        "caption_draft": (
            "Exploratory GaN central transverse sampling at 45.7 GPa. "
            "(a) Original 600-eV VCNEB enthalpy per GaN formula unit; shading marks "
            "the central atomic-coordinate domain. (b) Archived path centerline "
            "and 28 same-input frozen-atomic off-path static points, colored by "
            "excess enthalpy relative to the corresponding path image. Symbols "
            "mark actual calculations only; no two-dimensional contour is drawn. "
            "The maximum linear leave-one-anchor-out error is 2.329 meV/GaN, "
            "above the predeclared 1.0-meV/GaN gate."
        ),
        "source_sha256": {
            "first_raw_audit": sha256(args.first_audit),
            "refinement_raw_audit": sha256(args.refinement_audit),
            "components": sha256(args.components),
            "plot_script": sha256(Path(__file__)),
            "source_csv": sha256(paths["_source_data.csv"]),
        },
        "exports_sha256": {suffix.lstrip("."): sha256(paths[suffix])
                           for suffix in (".pdf", ".svg", ".tiff", ".png")},
    }
    paths["_qa.json"].write_text(json.dumps(qa, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": qa["status"], "peak_image": qa["path_peak_image"],
                      "forward_barrier_eV_per_GaN": qa["forward_barrier_eV_per_GaN"],
                      "contour_certified": qa["contour_certified"]}))


if __name__ == "__main__":
    main()
