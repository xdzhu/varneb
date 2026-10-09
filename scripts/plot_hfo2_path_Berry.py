"""Plot the actual E044 R3 measurements; never choose an ambiguous branch."""
from __future__ import annotations
import argparse
import csv
import importlib
import json
from pathlib import Path

from scripts import hfo2_switching_path_polarization as flow
from scripts.analyze_hfo2_channel_network import read_evaluated_observation
from scripts.audit_hfo2_static_replica import sha256
from vcneb import VCNEB


CASE = flow.CASE/"switching_path_polarization_20261009"


def build_data(repository):
    root = repository/CASE
    checker = importlib.import_module("benchmarks.hfo2_channels.20261008.switching_path_polarization_20261009.check_export")
    actual = checker.audit(repository, root/"completed_HF")
    if actual != json.loads((root/"offline_replay_local.json").read_text()):
        raise ValueError("raw result replay changed")
    result = json.loads((root/"completed_HF/summary.json").read_text())
    if not result["all_longitudinal_gates_passed"] or not result["all_sampled_lifts_unique"] or result["full_G1_certified"]:
        raise ValueError("measured longitudinal gates and conditional unique lifts required")
    rows = []
    for c, (name, observation, *_rest) in zip(result["channels"], flow.CHANNELS):
        if c["name"] != name or c["sampled_branch"]["continuous_path_certified"]:
            raise ValueError("original channel order and conditional interpretation required")
        images, _raw, _digest = read_evaluated_observation(repository/flow.CASE/observation)
        arc = VCNEB(images, pressure=0., climb=False, k=.2).reaction_coordinate()
        arc /= arc[-1]
        for i, audits in enumerate(c["native_Berry_audits_in_source_order"]):
            last = audits[-1]
            rows.append({"channel": name, "image": i, "s_normalized": float(arc[i]),
                "P_modern_SI_C_m2": last["value_modern_SI_C_m2"], "physical_quantum_C_m2": last["physical_quantum_C_m2"],
                "reported_reduced": c["sampled_branch"]["raw_reduced"][i],
                "conditional_lift_reduced": c["sampled_branch"]["lifted_reduced"][i],
                "explicit_branch_integer": c["sampled_branch"]["branch_integers"][i],
                "224_to_228_change_mC_m2": 1000*c["longitudinal_224_to_228_modular_change_C_m2"][i]})
    return {"rows": rows, "new_DFT_calls": 0, "independent_replicates": 0,
            "summary_sha256": sha256(root/"completed_HF/summary.json"),
            "raw_replay_sha256": sha256(root/"offline_replay_local.json"),
            "initial_gauge_integer": 0, "native_reduced_period": 2,
            "continuous_path_full_G1_or_TS_certified": False, "spontaneous_P_selected": False,
            "limitations": flow.LIMITATIONS}


def make_figure(data):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 10, "axes.labelsize": 10.5, "svg.fonttype": "none", "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 3, figsize=(183/25.4, 89/25.4), sharex=True)
    fig.subplots_adjust(left=.08, right=.99, bottom=.20, top=.73, wspace=.60)
    handles = []
    for name, color, marker in (("preserving", "#23877E", "o"), ("reversing", "#C97738", "^")):
        rows = [r for r in data["rows"] if r["channel"] == name]
        x = [r["s_normalized"] for r in rows]
        for j, key in enumerate(("reported_reduced", "conditional_lift_reduced", "224_to_228_change_mC_m2")):
            h, = axes[j].plot(x, [r[key] for r in rows], color=color, marker=marker,
                markersize=4., linewidth=1.35, markerfacecolor="white" if j == 0 else color)
            if j == 0:
                handles.append(h)
    for a, label in zip(axes, ("R3 p (reported class)", "R3 p (conditional lift)", r"224$\to$228 change (mC/m$^2$)")):
        a.set_ylabel(label)
        a.set_xlabel("Normalized arc")
        a.set_xlim(-.04, 1.04)
        a.set_xticks([0, .5, 1.])
        a.tick_params(direction="in", top=True, right=True, length=3.5)
        for s in a.spines.values():
            s.set_visible(True)
    axes[0].set_ylim(-1.45, 1.45)
    axes[1].set_ylim(1., 3.)
    axes[2].set_ylim(0., 1.)
    axes[2].text(.5, 1.10, r"Gate: 10 mC/m$^2$", transform=axes[2].transAxes,
                 ha="center", fontsize=9.5)
    for i, a in enumerate(axes):
        a.text(-.25, 1.055, f"({chr(97+i)})", transform=a.transAxes, fontsize=13,
               fontweight="normal", va="bottom")
    fig.legend(handles, ["Preserving candidate", "Reversing candidate"], ncols=2,
        loc="upper center", bbox_to_anchor=(.51, .985), frameon=True,
        facecolor="white", edgecolor="#9CA3AA", framealpha=.85, fontsize=10)
    fig.canvas.draw()
    return fig


def export(data, output):
    if output.exists():
        raise FileExistsError("refusing existing Berry figure")
    output.mkdir(parents=True)
    fig = make_figure(data)
    stem = output/"hfo2_path_Berry"
    for ext in ("svg", "pdf", "png", "tiff"):
        extra = {"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tiff" else {}
        fig.savefig(stem.with_suffix("."+ext), dpi=600 if ext == "tiff" else 300, facecolor="white", **extra)
    with (output/"source_data.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(data["rows"][0]))
        writer.writeheader()
        writer.writerows(data["rows"])
    qa = {k: v for k, v in data.items() if k != "rows"}
    positions = [a.get_position().bounds for a in fig.axes]
    qa.update({"plot_script_sha256": sha256(Path(__file__)), "source_records": 18,
        "panel_positions_normalized": positions, "all_spines_visible": all(s.get_visible() for a in fig.axes for s in a.spines.values()),
        "panel_titles_empty": all(not a.get_title() for a in fig.axes), "panel_label_weight": "normal",
        "row_alignment_verified": max(p[1] for p in positions)-min(p[1] for p in positions) < 1e-12,
        "figure_dimensions_mm": (fig.get_size_inches()*25.4).tolist(), "legend_frame_alpha": .85,
        "visual_review": "pending; actual PNG inspection required",
        "file_sha256": {p.name: sha256(p) for p in sorted(output.iterdir()) if p.is_file()}})
    flow.save(output/"qa.json", qa)
    import matplotlib.pyplot as plt
    plt.close(fig)
    return qa


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    print(json.dumps({"source_records": export(build_data(a.repository), a.output)["source_records"], "new_DFT_calls": 0}))
