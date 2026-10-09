"""Plot the frozen 6/39/69/32 G1 stage, without promoting a sampled peak to TS.

Historical figure adapters are deliberately left unchanged. This dated bundle
joins actual samples only and replays source-bound descriptors before drawing.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from ase.io import read

from scripts.analyze_hfo2_chain_observations import analyze_observation
from scripts.analyze_hfo2_network_progress import analyze_progress
from scripts.analyze_hfo2_parent_patterns import rotated_t_triplet
from scripts.audit_hfo2_static_replica import sha256
from scripts.plot_hfo2_converged_band import assert_replay_equal
from scripts.plot_hfo2_network_update import verify_material_replay


def build_data(specification: Path, material_report: Path, chain_report: Path) -> dict:
    recorded = json.loads(material_report.read_text())
    actual = analyze_progress(specification)
    differences = verify_material_replay(recorded, actual)
    names = ["PO_to_T", "PO_to_M", "PO_flip_T_pattern_preserving",
             "PO_flip_T_pattern_reversing"]
    channels = actual["channels"]
    if ([c["name"] for c in channels] != names
            or [c["snapshot_step"] for c in channels] != [6, 39, 69, 32]
            or [c["ordinary_residual_passed"] for c in channels] != [True, True, True, False]):
        raise ValueError("the frozen 6/39/69/32 three-pass stage is required")
    if actual["full_G1_passed"] or actual["TS_certified"] or actual["new_DFT_calls"] != 0:
        raise ValueError("this observation must not certify unmeasured gates")

    # Network rows deliberately do not expose physical tangential forces. Replay
    # the independent chain descriptor, rather than infer them from a polyline.
    specification_data = json.loads(specification.read_text())
    root = specification.parent
    t_path = root / specification_data["T_reference"]
    gamma_path = root / specification_data["T_Gamma_reference"]
    archived = json.loads(chain_report.read_text())
    if (archived["analysis_script_sha256"] != sha256(Path(analyze_observation.__code__.co_filename))
            or archived["T_reference_sha256"] != sha256(t_path)
            or archived["T_Gamma_source_sha256"] != sha256(gamma_path)):
        raise ValueError("stale chain analysis or reference hash")
    t = read(t_path, format="vasp")
    parent, patterns, _ = rotated_t_triplet(t)
    preserving = archived["observations"][0]
    folder = root / specification_data["channels"][2]["observation"]
    with np.load(gamma_path) as gamma:
        replayed = analyze_observation(folder, t, parent, patterns, gamma)
    assert_replay_equal(preserving, replayed)
    if preserving["source_observation_sha256"] != channels[2]["source_observation_sha256"]:
        raise ValueError("the chain and network snapshots differ")
    tangent = {r["image_index"]: r.get("true_tangential_euclidean_eV_A")
               for r in preserving["images"]}
    rows = []
    for channel in channels:
        for r in channel["rows"]:
            rows.append({"channel": channel["name"], "source_job_id": channel["source_job_id"],
                         "snapshot_step": channel["snapshot_step"],
                         **{k: v for k, v in r.items() if k != "structure_audit"},
                         "preserving_true_tangent_eV_A": tangent[r["source_image_index"]]
                         if channel["name"] == names[2] else None,
                         **{f"space_group_symprec_{sym['symprec_A']}_A": sym["symbol"]
                            for sym in r["structure_audit"]["symmetry_sweep"]}})
    return {"channels": channels, "preserving": preserving, "rows": rows,
            "material_report_sha256": sha256(material_report),
            "chain_report_sha256": sha256(chain_report),
            "specification_sha256": actual["specification_sha256"],
            "analysis_source_sha256": actual["analysis_source_sha256"],
            "physical_analysis_source_sha256": actual["physical_analysis_source_sha256"],
            "runtime_metadata_differences_retained": differences,
            "existing_image_records_including_cached_duplicates": len(rows),
            "new_DFT_calls": 0, "full_G1_passed": False, "TS_certified": False,
            "smooth_MEP_interpolation_used": False}


def make_figure(data):
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 10, "axes.labelsize": 11, "legend.fontsize": 9.3,
                         "axes.linewidth": .9, "svg.fonttype": "none", "pdf.fonttype": 42})
    fig, grid = plt.subplots(2, 2, figsize=(183 / 25.4, 177 / 25.4), sharex=True)
    fig.subplots_adjust(left=.13, right=.985, bottom=.10, top=.815, wspace=.42, hspace=.43)
    axes = grid.ravel()
    colors = ("#256B91", "#757D85", "#23877E", "#C97738")
    markers = ("o", "s", "D", "^")
    labels = (r"PO$^+\to$T (step 6, pass)", r"PO$^+\to$M (step 39, pass)",
              "Preserving flip (step 69, pass)", "Reversing flip (step 32, not pass)")
    handles = []
    for j, channel in enumerate(data["channels"]):
        rows = channel["rows"]
        options = {"color": colors[j], "marker": markers[j], "markersize": 4.,
                   "linewidth": 1.5, "linestyle": "--" if j == 3 else "-"}
        x = [r["s_normalized"] for r in rows]
        handle, = axes[0 if j < 3 else 1].plot(
            x, [r["energy_relative_PO_meV_fu"] for r in rows], **options)
        handles.append(handle)
        axes[2].plot(x[1:-1], [r["NEB_max_vector_eV_A"] for r in rows[1:-1]], **options)
    axes[0].set_ylim(-85, 130)
    center = data["channels"][2]["rows"][4]
    axes[0].annotate("Pbcn image 4", xy=(center["s_normalized"], center["energy_relative_PO_meV_fu"]),
                     xytext=(.30, -55), ha="center", fontsize=9.4,
                     arrowprops={"arrowstyle": "-", "color": colors[2], "linewidth": .8},
                     bbox={"facecolor": "white", "edgecolor": "none", "alpha": .85, "pad": 1.})
    axes[1].set_ylim(-15, 430)
    for axis in axes[:2]:
        axis.set_ylabel(r"$E-E_{\rm PO^+}$ (meV/f.u.)")
        axis.axhline(0, color="#D4D9DE", linewidth=.7, zorder=0)
    axes[2].set_ylabel(r"NEB max vector (eV/$\mathrm{\AA}$)")
    axes[2].set_ylim(0, .23)
    rows = data["preserving"]["images"][1:-1]
    total_arc = data["preserving"]["images"][-1]["extended_reaction_coordinate_A"]
    x = [r["extended_reaction_coordinate_A"] / total_arc for r in rows]
    for key, label, color, marker, style in (
        ("true_tangential_euclidean_eV_A", r"$|F_{\rm tangent}^{\rm physical}|$", "#333B43", "o", "-"),
        ("NEB_atomic_block_max_vector_eV_A", "NEB atomic", colors[2], "D", "--"),
        ("NEB_cell_block_max_vector_eV_A", "NEB cell", "#8BAA6D", "s", ":"),
    ):
        axes[3].plot(x, [abs(r[key]) for r in rows], label=label, color=color,
                     marker=marker, markersize=4., linewidth=1.4, linestyle=style)
    axes[3].set_ylabel(r"Preserving: force (eV/$\mathrm{\AA}$)")
    axes[3].set_ylim(0, .36)
    axes[3].legend(loc="upper center", bbox_to_anchor=(.5, 1.32), ncols=1,
                   fontsize=8.8, frameon=True, facecolor="white", edgecolor="#9CA3AA",
                   framealpha=.85, labelspacing=.25, borderpad=.3)
    for axis in axes[2:]:
        axis.axhline(.10, color="#9AA2AA", linewidth=1., linestyle="--", zorder=0)
        axis.set_xlabel("Normalized source generalized arc")
    for i, axis in enumerate(axes):
        axis.set_xlim(-.025, 1.025)
        axis.set_xticks([0, .25, .5, .75, 1.])
        axis.tick_params(direction="in", top=True, right=True, length=4)
        for spine in axis.spines.values():
            spine.set_visible(True)
        axis.text(-.19, 1.035, f"({chr(97+i)})", transform=axis.transAxes,
                  fontsize=13, fontweight="normal", va="bottom")
    fig.legend(handles, labels, ncols=2, loc="upper center", bbox_to_anchor=(.53, .985),
               frameon=True, facecolor="white", edgecolor="#9CA3AA", framealpha=.85,
               handlelength=1.9, columnspacing=1.1, labelspacing=.55)
    fig.align_ylabels(axes[::2])
    fig.align_ylabels(axes[1::2])
    fig.canvas.draw()
    return fig


def export(data, output: Path) -> dict:
    if output.exists():
        raise FileExistsError("refusing an existing figure bundle")
    fig = make_figure(data)
    output.mkdir(parents=True)
    stem = output / "hfo2_G1_preserving_pass"
    for ext in ("svg", "pdf", "png", "tiff"):
        kwargs = {"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tiff" else {}
        fig.savefig(stem.with_suffix("." + ext), dpi=600 if ext == "tiff" else 300,
                    facecolor="white", **kwargs)
    csv_path = output / "source_data.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(data["rows"][0]))
        writer.writeheader()
        writer.writerows(data["rows"])
    positions = [a.get_position().bounds for a in fig.axes]
    qa = {k: v for k, v in data.items() if k not in ("rows", "channels", "preserving")}
    qa.update({"plot_script_sha256": sha256(Path(__file__)), "source_data_sha256": sha256(csv_path),
               "panel_positions_normalized": positions, "panel_labels": ["(a)", "(b)", "(c)", "(d)"],
               "panel_label_weight": "normal", "panel_titles_empty": all(not a.get_title() for a in fig.axes),
               "all_spines_visible": all(s.get_visible() for a in fig.axes for s in a.spines.values()),
               "row_alignment_verified": all(abs(positions[i][1]-positions[i+1][1]) < 1e-12 for i in (0, 2)),
               "column_alignment_verified": all(abs(positions[i][0]-positions[i+2][0]) < 1e-12 for i in (0, 1)),
               "figure_dimensions_mm": (fig.get_size_inches()*25.4).tolist(),
               "shared_legend_frame_alpha": fig.legends[0].get_frame().get_alpha(),
               "statistics": "37 records including cached/reused endpoints; no independent replicates or error bars",
               "visual_review": "pending; actual pixel inspection required",
               "file_sha256": {p.name: sha256(p) for p in sorted(output.glob(stem.name + ".*"))}})
    (output / "qa.json").write_text(json.dumps(qa, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return qa


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("specification", "material-report", "chain-report", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("refusing an existing figure bundle")
    data = build_data(args.specification, args.material_report, args.chain_report)
    qa = export(data, args.output)
    print(json.dumps({"output": str(args.output), "records": qa["existing_image_records_including_cached_duplicates"],
                      "new_DFT_calls": 0, "TS_certified": False}))


if __name__ == "__main__":
    main()
