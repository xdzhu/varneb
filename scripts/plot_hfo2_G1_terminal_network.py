"""E044 editorial update: actual 6/1(refined)/69/45 bands, zero new DFT."""
from __future__ import annotations

import argparse
import csv
import importlib
import json
from pathlib import Path

import numpy as np
from ase.io import read

from scripts.analyze_hfo2_chain_observations import analyze_observation
from scripts.analyze_hfo2_parent_patterns import rotated_t_triplet
from scripts.audit_hfo2_static_replica import sha256
from scripts.plot_hfo2_network_update import verify_material_replay

CASE = Path("benchmarks/hfo2_channels/20261008/reversing_peak_sampling_20261009")


def build_data(repository):
    root = repository / CASE
    specification = root / "network_specification.json"
    material = root / "terminal_network_local.json"
    adapter = importlib.import_module("benchmarks.hfo2_channels.20261008.reversing_peak_sampling_20261009.check_terminal_network")
    actual = adapter.audit(specification)
    differences = verify_material_replay(json.loads(material.read_text()), actual)
    if (len(actual["channels"]) != 4 or actual["existing_image_records"] != 40
            or [c["snapshot_step"] for c in actual["channels"]] != [6, 1, 69, 45]
            or not all(c["ordinary_residual_passed"] for c in actual["channels"])
            or actual["full_G1_passed"] or actual["TS_certified"]):
        raise ValueError("the four original ordinary terminal bands are required")
    spec = json.loads(specification.read_text())
    t = read(root / spec["T_reference"], format="vasp")
    parent, patterns, _ = rotated_t_triplet(t)
    with np.load(root / spec["T_Gamma_reference"]) as gamma:
        preserving = analyze_observation(root / spec["channels"][2]["observation"], t, parent, patterns, gamma)
    if preserving["source_observation_sha256"] != actual["channels"][2]["source_observation_sha256"]:
        raise ValueError("tangent and material audit source mismatch")
    tangent = {r["image_index"]: r.get("true_tangential_euclidean_eV_A") for r in preserving["images"]}
    rows = []
    for c in actual["channels"]:
        for r in c["rows"]:
            rows.append({"channel": c["name"], "source_job_id": c["source_job_id"], "snapshot_step": c["snapshot_step"],
                         **{k: v for k, v in r.items() if k != "structure_audit"},
                         "preserving_true_tangent_eV_A": tangent[r["source_image_index"]] if c["name"] == actual["channels"][2]["name"] else None,
                         **{f"space_group_symprec_{s['symprec_A']}_A": s["symbol"] for s in r["structure_audit"]["symmetry_sweep"]}})
    return {"channels": actual["channels"], "preserving": preserving, "rows": rows,
            "material_report_sha256": sha256(material), "specification_sha256": sha256(specification),
            "physical_analysis_source_sha256": actual["analysis_source_sha256"],
            "runtime_metadata_differences_retained": differences,
            "existing_image_records_including_cached_duplicates": 40, "new_DFT_calls": 0,
            "full_G1_passed": False, "TS_certified": False, "smooth_MEP_interpolation_used": False}


def make_figure(data):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 10, "axes.labelsize": 11, "legend.fontsize": 9.3,
        "axes.linewidth": .9, "svg.fonttype": "none", "pdf.fonttype": 42})
    fig, grid = plt.subplots(2, 2, figsize=(183/25.4, 177/25.4), sharex=True)
    fig.subplots_adjust(left=.13, right=.985, bottom=.10, top=.815, wspace=.42, hspace=.43)
    axes = grid.ravel()
    colors = ("#256B91", "#757D85", "#23877E", "#C97738")
    markers = ("o", "s", "D", "^")
    labels = (r"PO$^+\to$T (10 total)", r"PO$^+\to$M (12 total)",
              "Preserving flip (9 total)", "Reversing flip (9 total)")
    handles = []
    for j, c in enumerate(data["channels"]):
        rows = c["rows"]
        options = {"color": colors[j], "marker": markers[j], "markersize": 4., "linewidth": 1.5}
        x = [r["s_normalized"] for r in rows]
        h, = axes[0 if j < 3 else 1].plot(x, [r["energy_relative_PO_meV_fu"] for r in rows], **options)
        handles.append(h)
        axes[2].plot(x[1:-1], [r["NEB_max_vector_eV_A"] for r in rows[1:-1]], **options)
    center = data["channels"][2]["rows"][4]
    axes[0].annotate("Pbcn image 4", xy=(center["s_normalized"], center["energy_relative_PO_meV_fu"]),
        xytext=(.34, -53), ha="center", fontsize=9.4,
        arrowprops={"arrowstyle": "-", "color": colors[2], "linewidth": .8},
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": .85, "pad": 1.})
    axes[0].set_ylim(-85, 130)
    axes[1].set_ylim(-15, 425)
    for a in axes[:2]:
        a.set_ylabel(r"$E-E_{\rm PO^+}$ (meV/f.u.)")
        a.axhline(0, color="#D4D9DE", linewidth=.7, zorder=0)
    axes[2].set_ylabel(r"NEB max vector (eV/$\mathrm{\AA}$)")
    axes[2].set_ylim(0, .115)
    rows = data["preserving"]["images"][1:-1]
    arc = data["preserving"]["images"][-1]["extended_reaction_coordinate_A"]
    x = [r["extended_reaction_coordinate_A"]/arc for r in rows]
    for key, label, color, marker, style in (
        ("true_tangential_euclidean_eV_A", r"$|F_{\rm tangent}^{\rm physical}|$", "#333B43", "o", "-"),
        ("NEB_atomic_block_max_vector_eV_A", "NEB atomic", colors[2], "D", "--"),
        ("NEB_cell_block_max_vector_eV_A", "NEB cell", "#8BAA6D", "s", ":")):
        axes[3].plot(x, [abs(r[key]) for r in rows], label=label, color=color, marker=marker,
                     markersize=4., linewidth=1.4, linestyle=style)
    axes[3].set_ylabel(r"Preserving: force (eV/$\mathrm{\AA}$)")
    axes[3].set_ylim(0, .36)
    axes[3].legend(loc="upper center", bbox_to_anchor=(.5, 1.32), ncols=1, fontsize=8.8,
        frameon=True, facecolor="white", edgecolor="#9CA3AA", framealpha=.85, labelspacing=.25, borderpad=.3)
    for a in axes[2:]:
        a.axhline(.10, color="#9AA2AA", linewidth=1., linestyle="--", zorder=0)
        a.set_xlabel("Normalized generalized arc")
    for i, a in enumerate(axes):
        a.set_xlim(-.025, 1.025)
        a.set_xticks([0, .25, .5, .75, 1.])
        a.tick_params(direction="in", top=True, right=True, length=4)
        for s in a.spines.values():
            s.set_visible(True)
        a.text(-.19, 1.035, f"({chr(97+i)})", transform=a.transAxes, fontsize=13, fontweight="normal", va="bottom")
    fig.legend(handles, labels, ncols=2, loc="upper center", bbox_to_anchor=(.53, .985),
        frameon=True, facecolor="white", edgecolor="#9CA3AA", framealpha=.85, handlelength=1.9,
        columnspacing=1.1, labelspacing=.55)
    fig.align_ylabels(axes[::2])
    fig.align_ylabels(axes[1::2])
    fig.canvas.draw()
    return fig


def export(data, output):
    if output.exists():
        raise FileExistsError("refusing existing terminal figure export")
    fig = make_figure(data)
    output.mkdir(parents=True)
    stem = output/"hfo2_G1_terminal"
    for ext in ("svg", "pdf", "png", "tiff"):
        extra = {"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tiff" else {}
        fig.savefig(stem.with_suffix("."+ext), dpi=600 if ext == "tiff" else 300, facecolor="white", **extra)
    source = output/"source_data.csv"
    with source.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(data["rows"][0]))
        writer.writeheader()
        writer.writerows(data["rows"])
    positions = [a.get_position().bounds for a in fig.axes]
    qa = {k: v for k, v in data.items() if k not in ("rows", "channels", "preserving")}
    qa.update({"plot_script_sha256": sha256(Path(__file__)), "source_data_sha256": sha256(source),
        "panel_positions_normalized": positions, "panel_labels": ["(a)", "(b)", "(c)", "(d)"],
        "panel_label_weight": "normal", "panel_titles_empty": all(not a.get_title() for a in fig.axes),
        "all_spines_visible": all(s.get_visible() for a in fig.axes for s in a.spines.values()),
        "row_alignment_verified": all(abs(positions[i][1]-positions[i+1][1]) < 1e-12 for i in (0, 2)),
        "column_alignment_verified": all(abs(positions[i][0]-positions[i+2][0]) < 1e-12 for i in (0, 1)),
        "figure_dimensions_mm": (fig.get_size_inches()*25.4).tolist(),
        "shared_legend_frame_alpha": fig.legends[0].get_frame().get_alpha(),
        "statistics": "40 records including cached/reused endpoints; no independent replicates or error bars",
        "visual_review": "pending; actual pixel inspection required",
        "dated_bundle": "four_terminal_bands_6_1_69_45",
        "file_sha256": {p.name: sha256(p) for p in sorted(output.glob(stem.name+".*"))}})
    (output/"qa.json").write_text(json.dumps(qa, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return qa


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    qa = export(build_data(a.repository), a.output)
    print(json.dumps({"records": qa["existing_image_records_including_cached_duplicates"], "new_DFT_calls": 0}))


if __name__ == "__main__":
    main()
