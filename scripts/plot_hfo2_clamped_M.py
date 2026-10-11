"""Draw the frozen E068/E071 training chain; no DFT, fitting or new reference."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

BASE = Path("benchmarks/hfo2_channels/20261008")
MODE_REPORT = BASE / "clamped_reference_E071_20261011/analysis.json"
OBSERVATION = BASE / "clamped_M_terminal_E068_20261011/observations/step_0013/observation.json"
AUDIT = BASE / "clamped_M_terminal_E068_20261011/audit_receipt.json"
MODE_CANONICAL_SHA = "3874505a408f208fba8b58ffd78b3eee9adabf1905fd6232f52b33b55a02536c"
OBSERVATION_SHA = "a8e5d2528a72aa0742d205bc4bdaf76e7d77eb9d688688ae95ed5daf4853ebd5"
AUDIT_SHA = "56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c"
FRAME_IDS = [103, 106, 114, 127]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_sha(value):
    """Pin JSON values independent of checkout newline formatting."""
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode("utf-8")).hexdigest()


def build_data(repository):
    mode_path, obs_path, audit_path = (repository / p for p in (MODE_REPORT, OBSERVATION, AUDIT))
    mode, obs = (json.loads(p.read_text(encoding="utf-8")) for p in (mode_path, obs_path))
    if (canonical_sha(mode) != MODE_CANONICAL_SHA or sha256(obs_path) != OBSERVATION_SHA
            or sha256(audit_path) != AUDIT_SHA):
        raise ValueError("frozen material sources differ")
    if (mode["source_job_id"] != "28722320" or mode["snapshot_step"] != 13
            or mode["terminal_audit_sha256"] != AUDIT_SHA
            or mode["trajectory_sha256"] != obs["evaluated_chain_sha256"]
            or mode["strain_training_condition"] != 0 or mode["pressure_GPa"] != 0
            or not obs["ordinary_residual_pass"] or obs["climb"]
            or mode["replayed_fmax_eV_A"] > .10
            or [f["frame_id"] for f in mode["Cmma_reference_frames"]] != FRAME_IDS):
        raise ValueError("original zero-strain ordinary training chain required")
    energy = np.asarray(mode["relative_energy_meV_fu"], float)
    arc = np.asarray(obs["extended_reaction_coordinate_A"], float)
    green = np.asarray(mode["green_strains_from_original_T"], float)
    triplet = np.asarray(mode["T_triplet_rank3_fraction"], float)
    optical = np.asarray(mode["T_lowest_two_complete_doublets_rank4_fraction"], float)
    cmma = np.asarray([f["common_T_origin_fraction"] for f in mode["Cmma_reference_frames"]])
    if (energy.shape != (9,) or arc.shape != (9,) or green.shape != (9, 3, 3)
            or triplet.shape != (9,) or optical.shape != (9,) or cmma.shape != (4, 9)
            or not all(np.isfinite(a).all() for a in (energy, arc, green, triplet, optical, cmma))
            or not np.allclose(energy, obs["relative_enthalpy_meV_fu"], atol=1e-9, rtol=0)
            or arc[0] != 0 or np.any(np.diff(arc) <= 0)
            or not np.allclose(green, green.transpose(0, 2, 1), atol=1e-12, rtol=0)
            or not np.allclose(green[:, [0, 2, 0], [0, 2, 2]], 0, atol=1e-12, rtol=0)
            or any(np.any((a < 0) | (a > 1)) for a in (triplet, optical, cmma))
            or int(np.argmax(energy)) != mode["highest_image_index"]):
        raise ValueError("frozen path coordinates, fractions or boundary differ")
    rows = []
    for i in range(9):
        rows.append(dict(image_index=i, generalized_arc_A=float(arc[i]),
            normalized_generalized_arc=float(arc[i] / arc[-1]), energy_relative_PO_meV_fu=float(energy[i]),
            T_geometric_triplet_rank3_squared_fraction=float(triplet[i]),
            T_optical_doublets_rank4_squared_fraction=float(optical[i]),
            **{f"Cmma_frame_{frame}_rank4_squared_fraction": float(cmma[j, i])
               for j, frame in enumerate(FRAME_IDS)},
            **{f"Green_E_{a}{b}_original_T": float(green[i, a, b])
               for a, b in ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))}))
    return dict(rows=rows, highest_image_index=mode["highest_image_index"],
        source_job_id="28722320", snapshot_step=13, source_paths=[str(p) for p in (MODE_REPORT, OBSERVATION, AUDIT)],
        mode_report_canonical_sha256=MODE_CANONICAL_SHA,
        source_file_sha256={str(p): sha256(repository / p) for p in (MODE_REPORT, OBSERVATION, AUDIT)},
        new_DFT_calls=0, native_data_reaudit_claimed=False,
        ordinary_fmax_eV_A=mode["replayed_fmax_eV_A"], smooth_MEP_interpolation_used=False,
        TS_certified=False, G3_selection_or_holdout_access=False)


def make_figure(data):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 10, "axes.labelsize": 11, "legend.fontsize": 9.5,
        "axes.linewidth": .9, "svg.fonttype": "none", "pdf.fonttype": 42})
    fig, axes = plt.subplots(3, 1, figsize=(183 / 25.4, 207 / 25.4), sharex=True)
    fig.subplots_adjust(left=.145, right=.975, bottom=.085, top=.91, hspace=.60)
    rows = data["rows"]
    x = np.array([r["normalized_generalized_arc"] for r in rows])
    peak = data["highest_image_index"]
    colors = ("#256B91", "#757D85", "#23877E")
    axes[0].plot(x, [r["energy_relative_PO_meV_fu"] for r in rows], color=colors[0],
                 marker="o", markersize=4.5, lw=1.7, label="9 sampled images (7 interior)")
    axes[0].axhline(0, lw=.7, color="#D4D9DE", zorder=0)
    axes[0].set_ylabel(r"$E-E_{\rm PO^+}$ (meV/f.u.)")
    axes[0].set_ylim(-115, 125)
    axes[0].set_yticks([-100, -50, 0, 50, 100])
    for key, label, color, marker, style in (
        ("T_geometric_triplet_rank3_squared_fraction", "T geometric (rank 3)", colors[0], "o", "-"),
        ("T_optical_doublets_rank4_squared_fraction", "T optical (rank 4)", colors[1], "s", "--")):
        axes[1].plot(x, [100*r[key] for r in rows], color=color, marker=marker,
                     markersize=4, lw=1.4, ls=style, label=label)
    cmma = np.array([[100*r[f"Cmma_frame_{f}_rank4_squared_fraction"] for r in rows] for f in FRAME_IDS])
    for j, y in enumerate(cmma):
        axes[1].plot(x, y, color=colors[2], marker="D" if j == 0 else None,
                     markersize=3.6, lw=1, ls=":", label="Cmma (rank 4; 4 frames)" if j == 0 else None)
    axes[1].fill_between(x, cmma.min(axis=0), cmma.max(axis=0), color=colors[2], alpha=.18)
    axes[1].set_ylabel("Atomic squared-norm coverage (%)")
    axes[1].set_ylim(0, 70)
    axes[1].set_yticks([0, 20, 40, 60])
    for key, label, color, marker, style in (
        ("Green_E_11_original_T", r"$E_{yy}$", colors[0], "o", "-"),
        ("Green_E_01_original_T", r"$E_{xy}$", colors[2], "s", "--"),
        ("Green_E_12_original_T", r"$E_{yz}$", colors[1], "D", ":")):
        axes[2].plot(x, [100*r[key] for r in rows], color=color, marker=marker,
                     markersize=4, lw=1.4, ls=style, label=label)
    axes[2].set_ylabel("Green strain, original T axes (%)")
    axes[2].set_ylim(-8.5, 3)
    axes[2].set_yticks([-8, -6, -4, -2, 0, 2])
    axes[2].set_xlabel("Normalized generalized arc")
    for i, ax in enumerate(axes):
        ax.set_xlim(-.025, 1.025)
        ax.set_xticks([0, .25, .5, .75, 1.])
        ax.tick_params(direction="in", top=True, right=True, length=4)
        for spine in ax.spines.values():
            spine.set_visible(True)
        ax.axvline(x[peak], color="#C97738", lw=.9, ls="--", zorder=0)
        ax.text(-.145, 1.075, f"({chr(97+i)})", transform=ax.transAxes,
                fontsize=13, fontweight="normal", va="bottom")
        ax.legend(loc="lower center", bbox_to_anchor=(.52, 1.075), ncols=3 if i else 1,
                  frameon=True, facecolor="white", edgecolor="#9CA3AA", framealpha=.85,
                  borderpad=.35, columnspacing=1., handlelength=1.7)
    fig.align_ylabels(axes)
    fig.canvas.draw()
    return fig


def export(data, output):
    if output.exists():
        raise FileExistsError("fresh dated figure directory required")
    fig = make_figure(data)
    output.mkdir(parents=True)
    source = output / "source_data.csv"
    with source.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data["rows"][0]))
        writer.writeheader()
        writer.writerows(data["rows"])
    stem = output / "hfo2_clamped_M"
    for extension in ("svg", "pdf", "png", "tiff"):
        extra = {"pil_kwargs": {"compression": "tiff_lzw"}} if extension == "tiff" else {}
        fig.savefig(stem.with_suffix("." + extension), dpi=600 if extension == "tiff" else 300,
                    facecolor="white", **extra)
    positions = [ax.get_position().bounds for ax in fig.axes]
    qa = {k: v for k, v in data.items() if k != "rows"}
    qa.update(plot_script_sha256=sha256(Path(__file__)), source_data_sha256=sha256(source),
        image_records=9, independent_replicates=0, figure_dimensions_mm=(fig.get_size_inches()*25.4).tolist(),
        panel_positions_normalized=positions, axes_columns_aligned=all(
            np.allclose((p[0], p[2]), (positions[0][0], positions[0][2]), atol=1e-12, rtol=0) for p in positions),
        all_spines_visible=all(s.get_visible() for ax in fig.axes for s in ax.spines.values()),
        panel_labels=["(a)", "(b)", "(c)"], panel_label_weight="normal",
        panel_titles_empty=all(not ax.get_title() for ax in fig.axes),
        legends_outside_axes=True, legend_frame_alpha=.85, visual_review="pending",
        Cmma_range="four retained registrations, not numerical/statistical uncertainty",
        limitations=["atomic coverage is not energy or cell-strain decomposition",
                     "connected discrete points are not interpolated or certified continuous MEP",
                     "frozen numeric reports reused; native SCF audit is historical E068/E071 evidence"],
        file_sha256={p.name: sha256(p) for p in sorted(output.glob(stem.name + ".*"))})
    with (output / "qa.json").open("x", encoding="utf-8") as stream:
        json.dump(qa, stream, indent=2, allow_nan=False)
        stream.write("\n")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return qa


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = export(build_data(args.repository), args.output)
    print(json.dumps({"records": result["image_records"], "new_DFT_calls": 0,
                      "columns_aligned": result["axes_columns_aligned"]}))


if __name__ == "__main__":
    main()
