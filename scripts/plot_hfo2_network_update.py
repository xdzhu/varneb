"""Plot dated G1 observations, not final barriers or certified transition states."""

from __future__ import annotations

import argparse
import copy
import csv
import json
from pathlib import Path

import numpy as np

from scripts.analyze_hfo2_network_update import analyze_case
from scripts.audit_hfo2_static_replica import sha256
from scripts.plot_hfo2_converged_band import assert_replay_equal


def verify_material_replay(recorded, replayed):
    """Recheck all physical/source fields; retain runtime-warning/cache metadata.

    Pymatgen's cached spglib calls can emit different counts of the same
    deprecation warning across identical replays. They are preserved on both
    sides and reported separately, never silently discarded or treated as
    a physical descriptor. Package versions are likewise provenance metadata.
    """
    left, right = copy.deepcopy(recorded), copy.deepcopy(replayed)
    differences = []
    def metadata(a, b, key, path):
        x, y = a.pop(key), b.pop(key)
        if x != y:
            differences.append({"field": path + "/" + key, "recorded": x, "replayed": y})
    metadata(left, right, "packages", "")
    metadata(left["reference_T_structure_audit"], right["reference_T_structure_audit"], "warnings", "/reference_T_structure_audit")
    for i, (a, b) in enumerate(zip(left["channels"], right["channels"])):
        for j, (x, y) in enumerate(zip(a["rows"], b["rows"])):
            metadata(x["structure_audit"], y["structure_audit"], "warnings", f"/channels/{i}/rows/{j}/structure_audit")
    assert_replay_equal(left, right)
    return differences


def build_data(specification, material_report):
    recorded = json.loads(Path(material_report).read_text())
    actual = analyze_case(Path(specification))
    metadata_differences = verify_material_replay(recorded, actual)
    # This figure is intentionally tied to the disclosed incomplete four-path set.
    names = [c["name"] for c in actual["channels"]]
    expected = ["PO_to_T", "PO_to_M", "PO_flip_T_pattern_preserving", "PO_flip_T_pattern_reversing"]
    if names != expected or [c["snapshot_step"] for c in actual["channels"]] != [6, 20, 14, 10]:
        raise ValueError("the registered dated G1 snapshots and channel order required")
    if [c["ordinary_residual_passed"] for c in actual["channels"]] != [True, False, False, False]:
        raise ValueError("do not silently promote or substitute the figure's observation stage")
    rows = []
    for c in actual["channels"]:
        for r in c["rows"]:
            rows.append({"channel": c["name"], "source_job_id": c["source_job_id"],
                         "snapshot_step": c["snapshot_step"], "source_direction": c["source_direction"],
                         **{k: v for k, v in r.items() if k != "structure_audit"},
                         "space_group_symprec_0p001_A": r["structure_audit"]["symmetry_sweep"][0]["symbol"],
                         "space_group_symprec_0p01_A": r["structure_audit"]["symmetry_sweep"][1]["symbol"],
                         "space_group_symprec_0p05_A": r["structure_audit"]["symmetry_sweep"][2]["symbol"]})
    return {"channels": actual["channels"], "rows": rows,
            "material_report_sha256": sha256(Path(material_report)),
            "specification_sha256": actual["specification_sha256"],
            "source_T_reference_sha256": actual["T_reference_sha256"],
            "source_gamma_sha256": actual["T_Gamma_sha256"],
            "material_analysis_sha256": actual["analysis_source_sha256"],
            "existing_image_records": len(rows), "new_DFT_calls": 0,
            "runtime_metadata_differences_retained": metadata_differences,
            "TS_or_final_channel_ranking_certified": False,
            "reversing_energy_profile_in_panel_a": False,
            "all_four_energy_profiles_in_source_data": True}


def make_figure(data):
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 10, "axes.labelsize": 11, "legend.fontsize": 9.5,
                         "axes.linewidth": .9, "svg.fonttype": "none", "pdf.fonttype": 42})
    fig, grid = plt.subplots(2, 2, figsize=(183 / 25.4, 168 / 25.4), sharex=True)
    fig.subplots_adjust(left=.13, right=.985, bottom=.105, top=.82, wspace=.42, hspace=.34)
    axes = grid.ravel()
    colors = ("#256B91", "#757D85", "#23877E", "#C97738")
    markers = ("o", "s", "D", "^")
    labels = (r"PO$^+\to$T (ordinary converged)", r"PO$^+\to$M (step 20)",
              "Preserving flip (step 14)", "Reversing flip (step 10)")
    handles = []
    for j, channel in enumerate(data["channels"]):
        rows = channel["rows"]
        x = np.array([r["s_normalized"] for r in rows])
        options = {"color": colors[j], "marker": markers[j], "markersize": 3.9,
                   "linewidth": 1.5, "label": labels[j]}
        if j < 3:
            axes[0].plot(x, [r["energy_relative_PO_meV_fu"] for r in rows], **options)
        internal = rows[1:-1]
        axes[1].plot([r["s_normalized"] for r in internal],
                     [r["NEB_max_vector_eV_A"] for r in internal], **options)
        handle, = axes[2].plot(x, [r["Q_pattern_x_A"] for r in rows], **options)
        handles.append(handle)
        axes[3].plot(x, [100 * r["T_Green_strain_Frobenius_norm"] for r in rows], **options)
    t_energy = data["channels"][0]["rows"][-1]["energy_relative_PO_meV_fu"]
    axes[0].axhline(t_energy, color="#ADB4BB", linewidth=1., linestyle="--", zorder=0)
    axes[0].text(.035, t_energy + 20., r"relaxed T (dashed)", fontsize=9,
                 bbox={"facecolor": "white", "edgecolor": "none", "alpha": .8, "pad": 1.})
    axes[0].set_ylabel(r"$E-E_{\rm PO^+}$ (meV/f.u.)")
    axes[0].set_ylim(-85, 140)
    preserving = data["channels"][2]
    peak = preserving["rows"][preserving["peak_view_image_index"]]
    axes[0].annotate("Pbcn snapshot", xy=(peak["s_normalized"], peak["energy_relative_PO_meV_fu"]),
                     xytext=(.98, 68), textcoords="data", fontsize=9.5, ha="right", va="center",
                     arrowprops={"arrowstyle": "-", "color": colors[2], "linewidth": .8},
                     bbox={"facecolor": "white", "edgecolor": "none", "alpha": .8, "pad": 1.})
    axes[1].set_ylabel(r"NEB max vector (eV/$\mathrm{\AA}$)")
    axes[1].set_ylim(0, .84)
    axes[1].axhline(.10, color="#9AA2AA", linewidth=1., linestyle="--", zorder=0)
    axes[1].text(.04, .34, "ordinary target: 0.10", fontsize=9,
                 bbox={"facecolor": "white", "edgecolor": "none", "alpha": .8, "pad": 1.})
    axes[2].set_ylabel(r"Geometric pattern $Q_x$ ($\mathrm{\AA}$)")
    axes[2].set_ylim(-1.40, 1.10)
    axes[2].axhline(0, color="#D4D9DE", linewidth=.7, zorder=0)
    axes[3].set_ylabel(r"T-referenced $\|\eta\|_{\rm F}$ (%)")
    axes[3].set_ylim(0, 16)
    for i, axis in enumerate(axes):
        axis.set_xlim(-.025, 1.025)
        axis.set_xticks([0, .25, .5, .75, 1.])
        axis.tick_params(direction="in", top=True, right=True, length=4)
        for spine in axis.spines.values():
            spine.set_visible(True)
        axis.text(-.19, 1.035, f"({chr(97+i)})", transform=axis.transAxes,
                  fontsize=13, fontweight="normal", va="bottom")
    for axis in axes[-2:]:
        axis.set_xlabel("Normalized source generalized arc")
    fig.legend(handles, labels, ncols=2, loc="upper center", bbox_to_anchor=(.53, .985),
               frameon=True, facecolor="white", edgecolor="#9CA3AA", framealpha=.85,
               handlelength=1.9, columnspacing=1.3, labelspacing=.55)
    fig.align_ylabels(axes[::2])
    fig.align_ylabels(axes[1::2])
    fig.canvas.draw()
    return fig


def export(data, output):
    output = Path(output)
    if output.exists():
        raise FileExistsError("refusing an existing figure bundle")
    fig = make_figure(data)
    output.mkdir(parents=True)
    stem = output / "hfo2_G1_provisional_mechanisms"
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
    qa = {k: v for k, v in data.items() if k not in ("rows", "channels")}
    qa.update({"plot_script_sha256": sha256(Path(__file__)), "source_data_sha256": sha256(csv_path),
               "panel_positions_normalized": positions, "panel_labels": ["(a)", "(b)", "(c)", "(d)"],
               "panel_label_weight": "normal", "panel_titles_empty": all(not a.get_title() for a in fig.axes),
               "all_spines_visible": all(s.get_visible() for a in fig.axes for s in a.spines.values()),
               "row_alignment_verified": all(abs(positions[i][1]-positions[i+1][1]) < 1e-12 for i in (0, 2)),
               "column_alignment_verified": all(abs(positions[i][0]-positions[i+2][0]) < 1e-12 for i in (0, 1)),
               "figure_dimensions_mm": (fig.get_size_inches()*25.4).tolist(),
               "shared_legend_frame_alpha": fig.legends[0].get_frame().get_alpha(),
               "visual_review": "pending; numeric/layout checks do not replace actual visual inspection",
               "file_sha256": {p.name: sha256(p) for p in sorted(output.glob(stem.name + ".*"))}})
    (output / "qa.json").write_text(json.dumps(qa, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return qa


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--specification", type=Path, required=True)
    parser.add_argument("--material-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("refusing an existing figure bundle")
    data = build_data(args.specification, args.material_report)
    qa = export(data, args.output)
    print(json.dumps({"output": str(args.output), "records": qa["existing_image_records"], "new_DFT_calls": 0}))


if __name__ == "__main__":
    main()
