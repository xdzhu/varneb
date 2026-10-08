"""Plot a hash-bound ordinary band, retaining descriptive-mode limitations."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from ase.io import read

from scripts.analyze_hfo2_chain_observations import analyze_observation
from scripts.analyze_hfo2_channel_network import read_evaluated_observation
from scripts.analyze_hfo2_parent_patterns import rotated_t_triplet
from scripts.audit_hfo2_static_replica import sha256


def assert_replay_equal(record, replayed):
    """Compare every descriptor, including strings/hashes, not just energies."""
    if isinstance(record, dict):
        if not isinstance(replayed, dict) or record.keys() != replayed.keys():
            raise ValueError("stale descriptor keys")
        for key in record:
            assert_replay_equal(record[key], replayed[key])
    elif isinstance(record, list):
        if not isinstance(replayed, list) or len(record) != len(replayed):
            raise ValueError("stale descriptor dimensions")
        for left, right in zip(record, replayed):
            assert_replay_equal(left, right)
    elif isinstance(record, (int, float)) and not isinstance(record, bool):
        if (not isinstance(replayed, (int, float)) or not np.isfinite(record)
                or not np.isfinite(replayed)
                or not np.isclose(record, replayed, rtol=1e-10, atol=1e-11)):
            raise ValueError("stale or nonfinite numeric descriptor")
    elif record != replayed:
        raise ValueError("stale descriptor identity")


def optical_group_weights(frequencies, amplitudes, acoustic, peak, tolerance=1e-6):
    """Use complete frequency groups; weights are invariant within a group."""
    frequencies = np.asarray(frequencies, float)
    amplitudes = np.asarray(amplitudes, float)
    acoustic = np.asarray(acoustic)
    if (frequencies.shape != (36,) or amplitudes.ndim != 2
            or amplitudes.shape[1] != 36 or not np.isfinite(frequencies).all()
            or not np.isfinite(amplitudes).all() or acoustic.shape != (3,)
            or not np.issubdtype(acoustic.dtype, np.integer)
            or len(np.unique(acoustic)) != 3 or np.any(acoustic < 0)
            or np.any(acoustic >= 36) or not 0 <= peak < len(amplitudes)):
        raise ValueError("finite full Hf4O8 Gamma descriptors required")
    optical = np.setdiff1d(np.arange(36), acoustic)
    groups = []
    for mode in optical[np.argsort(frequencies[optical])]:
        if not groups or abs(frequencies[mode] - frequencies[groups[-1][0]]) > tolerance:
            groups.append([int(mode)])
        else:
            groups[-1].append(int(mode))
    total = np.sum(amplitudes[:, optical]**2, axis=1)
    valid = total > 1e-20
    if not valid[peak]:
        raise ValueError("undefined optical displacement at sampled maximum")
    squared = np.array([np.sum(amplitudes[:, group]**2, axis=1) for group in groups])
    ranking = sorted(range(len(groups)), key=lambda i: (-squared[i, peak], groups[i][0]))
    selected = ranking[:2]
    fractions = np.full((len(amplitudes), 3), np.nan)
    fractions[valid, :2] = squared[selected][:, valid].T / total[valid, None]
    fractions[valid, 2] = 1 - fractions[valid, :2].sum(axis=1)
    metadata = [{"mode_indices_zero_based": groups[i],
                 "mean_frequency_THz": float(np.mean(frequencies[groups[i]]))}
                for i in selected]
    return metadata, fractions, total


def build_data(dataset, variants, gamma_path):
    """Re-audit frozen physics and replay the archived analysis before drawing."""
    analysis_path = dataset / "analysis.json"
    analysis = json.loads(analysis_path.read_text())
    if analysis["new_DFT_calls"] != 0 or analysis["physical_parameters_changed"]:
        raise ValueError("unaltered frozen analysis required")
    if (sha256(variants / "T.vasp") != analysis["T_reference_sha256"]
            or sha256(gamma_path) != analysis["T_Gamma_source_sha256"]
            or sha256(Path(analyze_observation.__code__.co_filename)) != analysis["analysis_script_sha256"]):
        raise ValueError("reference or analysis source hash changed")
    if len(analysis["observations"]) != 1:
        raise ValueError("one converged band required")
    archived = analysis["observations"][0]
    folder = dataset / archived["observation"]
    images, raw, digest = read_evaluated_observation(folder)
    t = read(variants / "T.vasp", format="vasp")
    parent, patterns, _ = rotated_t_triplet(t)
    with np.load(gamma_path) as gamma:
        replayed = analyze_observation(folder, t, parent, patterns, gamma)
        if not np.array_equal(gamma["frequencies_thz"], analysis["T_Gamma_frequencies_THz"]):
            raise ValueError("reference frequencies changed")
    assert_replay_equal(archived, replayed)
    if (replayed["status"] != "NEB_residual_passed_TS_and_sampling_gates_pending"
            or replayed["replayed_fmax_eV_A"] > .10):
        raise ValueError("ordinary residual convergence required")
    records = replayed["images"]
    arc = np.array([r["extended_reaction_coordinate_A"] for r in records])
    if arc[0] != 0 or np.any(np.diff(arc) <= 0):
        raise ValueError("strictly ordered full generalized arc required")
    groups, fractions, norm = optical_group_weights(
        analysis["T_Gamma_frequencies_THz"],
        [r["T_Gamma_Q_sqrt_amu_A"] for r in records],
        replayed["Gamma_acoustic_mode_indices"], replayed["highest_image_index"])
    rows = []
    for i, r in enumerate(records):
        strain = np.asarray(r["green_strain_from_original_T"])
        row = {"image_index": i, "s_normalized": float(arc[i] / arc[-1]),
               "s_generalized_A": float(arc[i]), "energy_relative_T_meV_fu": r["relative_energy_meV_fu"],
               "Q_pattern_x_A": r["parent_pattern_Q_A"][0],
               "Q_pattern_y_A": r["parent_pattern_Q_A"][1],
               "Q_pattern_z_A": r["parent_pattern_Q_A"][2],
               "triplet_captured_fraction": r["parent_pattern_captured_squared_norm_fraction"],
               "triplet_residual_A": r["parent_pattern_residual_A"],
               "Gamma_optical_squared_norm_amu_A2": float(norm[i])}
        for name, a, b in (("xx", 0, 0), ("yy", 1, 1), ("zz", 2, 2),
                           ("yz", 1, 2), ("xz", 0, 2), ("xy", 0, 1)):
            row[f"green_strain_{name}"] = float(strain[a, b])
        for j, name in enumerate(("Gamma_group1_fraction", "Gamma_group2_fraction", "Gamma_other_fraction")):
            row[name] = float(fractions[i, j]) if np.isfinite(fractions[i, j]) else None
        for name, key in (("NEB_atomic_max_vector_eV_A", "NEB_atomic_block_max_vector_eV_A"),
                          ("NEB_cell_max_vector_eV_A", "NEB_cell_block_max_vector_eV_A")):
            row[name] = r.get(key)  # fixed endpoints have no NEB residual
        rows.append(row)
    peak = replayed["highest_image_index"]
    return {"rows": rows, "Gamma_selected_groups": groups,
            "Gamma_selection_rule": "two largest complete frequency groups at sampled maximum; descriptive",
            "source_analysis_sha256": sha256(analysis_path), "source_observation_sha256": digest,
            "source_T_reference_sha256": analysis["T_reference_sha256"],
            "source_Gamma_sha256": analysis["T_Gamma_source_sha256"],
            "source_job_id": raw["source_job_id"], "snapshot_step": raw["snapshot_step"],
            "replayed_fmax_eV_A": replayed["replayed_fmax_eV_A"], "peak_image_index": peak,
            "largest_residual_image_index": replayed["residual_dominant_image"],
            "peak_triplet_captured_fraction": rows[peak]["triplet_captured_fraction"],
            "PO_triplet_captured_fraction": rows[-1]["triplet_captured_fraction"],
            "new_DFT_calls": 0, "parameters_changed": False, "TS_certified": False,
            "energy_partition_or_prediction_claimed": False, "interpolation": "none; straight connections only"}


def make_figure(data):
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 10, "axes.labelsize": 11, "legend.fontsize": 9,
                         "svg.fonttype": "none", "pdf.fonttype": 42, "axes.linewidth": .9})
    fig, axes = plt.subplots(3, 2, figsize=(7.8, 6.9), sharex=True)
    fig.subplots_adjust(left=.12, right=.98, bottom=.10, top=.96, wspace=.42, hspace=.38)
    axes = axes.ravel()
    rows = data["rows"]
    x = np.array([r["s_normalized"] for r in rows])
    colors = ("#256B91", "#C97738", "#737A83")
    def values(key):
        return np.array([np.nan if r[key] is None else r[key] for r in rows])
    def line(axis, key, color, label=None, scale=1, style="-", marker="o"):
        axis.plot(x, values(key)*scale, style, color=color, marker=marker,
                  markersize=3.7, linewidth=1.5, label=label)
    def legend(axis, **kwargs):
        axis.legend(frameon=True, facecolor="white", edgecolor="#9CA3AA", framealpha=.85, **kwargs)
    for i, axis in enumerate(axes):
        axis.set_xlim(-.025, 1.025)
        axis.set_xticks([0, .25, .5, .75, 1])
        axis.tick_params(direction="in", top=True, right=True, length=4)
        for spine in axis.spines.values():
            spine.set_visible(True)
        axis.text(-.17, 1.04, f"({chr(97+i)})", transform=axis.transAxes,
                  fontsize=13, fontweight="normal", va="bottom")
    line(axes[0], "energy_relative_T_meV_fu", colors[0])
    axes[0].axhline(0, color="#C7CDD1", linewidth=.7, zorder=0)
    axes[0].set_ylabel(r"$E-E_T$ (meV/f.u.)")
    axes[0].set_ylim(-105, 65)
    p = data["peak_image_index"]
    axes[0].annotate("sampled maximum", (x[p], rows[p]["energy_relative_T_meV_fu"]),
                     xytext=(5, 12), textcoords="offset points", fontsize=9)
    axes[0].annotate("T", (x[0], 0), xytext=(6, 7), textcoords="offset points")
    axes[0].annotate(r"PO$^+$", (x[-1], rows[-1]["energy_relative_T_meV_fu"]),
                     xytext=(-29, 8), textcoords="offset points")
    for name, color in zip(("x", "y", "z"), colors):
        line(axes[1], f"Q_pattern_{name}_A", color, f"$Q_{name}$")
    axes[1].set_ylabel(r"Pattern amplitude ($\mathrm{\AA}$)")
    axes[1].set_ylim(-1.12, 1.45)
    legend(axes[1], loc="upper right", ncols=3, columnspacing=.8, handlelength=1.3)
    for name, color in zip(("xx", "yy", "zz"), colors):
        line(axes[2], f"green_strain_{name}", color, rf"$\eta_{{{name}}}$", scale=100)
    axes[2].set_ylabel("Green strain (%)")
    axes[2].set_ylim(-3.8, 7.2)
    legend(axes[2], loc="upper center", ncols=3, columnspacing=.6, handlelength=1.3)
    line(axes[3], "triplet_captured_fraction", colors[0], scale=100)
    axes[3].set_ylabel("Triplet captured norm (%)")
    axes[3].set_ylim(0, 107)
    for i, group in enumerate(data["Gamma_selected_groups"]):
        label = f"{group['mean_frequency_THz']:.2f} THz ({len(group['mode_indices_zero_based'])})"
        line(axes[4], f"Gamma_group{i+1}_fraction", colors[i], label, scale=100)
    line(axes[4], "Gamma_other_fraction", colors[2], "Other optical", scale=100, style="--", marker="s")
    axes[4].set_ylabel(r"T-$\Gamma$ optical norm (%)")
    axes[4].set_ylim(0, 147)
    axes[4].set_yticks([0, 50, 100])
    legend(axes[4], loc="upper center", ncols=1, labelspacing=.18, handlelength=1.7)
    line(axes[5], "NEB_atomic_max_vector_eV_A", colors[0], "Atomic")
    line(axes[5], "NEB_cell_max_vector_eV_A", colors[1], "Scaled cell", marker="s")
    axes[5].axhline(.10, color=colors[2], linestyle="--", linewidth=1)
    axes[5].text(.02, .104, "ordinary target", fontsize=9)
    axes[5].set_ylabel(r"NEB residual (eV/$\mathrm{\AA}$)")
    axes[5].set_ylim(0, .13)
    axes[5].set_yticks([0, .05, .10])
    legend(axes[5], loc="upper right", bbox_to_anchor=(1, .79), labelspacing=.2)
    for axis in axes[-2:]:
        axis.set_xlabel("Normalized generalized arc")
    fig.align_ylabels(axes[::2])
    fig.align_ylabels(axes[1::2])
    fig.canvas.draw()
    return fig


def export(data, output):
    if output.exists():
        raise FileExistsError("refusing existing figure bundle")
    fig = make_figure(data)
    output.mkdir(parents=True)
    stem = output / "hfo2_T_PO_ordinary_modes"
    for extension in ("svg", "pdf", "png", "tiff"):
        options = {"pil_kwargs": {"compression": "tiff_lzw"}} if extension == "tiff" else {}
        fig.savefig(stem.with_suffix(f".{extension}"), dpi=600 if extension == "tiff" else 300,
                    facecolor="white", **options)
    csv_path = output / "source_data.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(data["rows"][0]))
        writer.writeheader()
        writer.writerows(data["rows"])
    axes = fig.axes
    positions = [a.get_position().bounds for a in axes]
    qa = {key: value for key, value in data.items() if key != "rows"}
    qa.update({"source_data_sha256": sha256(csv_path), "plot_script_sha256": sha256(Path(__file__)),
               "panel_positions_normalized": positions, "panel_labels": [f"({chr(97+i)})" for i in range(6)],
               "panel_label_weight": "normal", "all_spines_visible": all(s.get_visible() for a in axes for s in a.spines.values()),
               "panel_titles_empty": all(not a.get_title() for a in axes),
               "row_alignment_verified": all(abs(positions[i][1]-positions[i+1][1]) < 1e-12 for i in (0, 2, 4)),
               "column_alignment_verified": all(abs(positions[i][0]-positions[i+2][0]) < 1e-12 for i in (0, 1, 2, 3)),
               "figure_dimensions_mm": (fig.get_size_inches()*25.4).tolist(),
               "visual_review": "pending; numerical/layout checks do not replace visual review",
               "file_sha256": {p.name: sha256(p) for p in sorted(output.glob("hfo2_T_PO_ordinary_modes.*"))}})
    (output / "qa.json").write_text(json.dumps(qa, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return qa


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008"
    parser.add_argument("--dataset", type=Path, default=root / "converged_gap")
    parser.add_argument("--variants", type=Path, default=root / "reference_variants")
    parser.add_argument("--gamma", type=Path, default=root / "gamma_analysis/T_d0.01.npz")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("refusing existing figure bundle")
    qa = export(build_data(args.dataset, args.variants, args.gamma), args.output)
    print(json.dumps({key: qa[key] for key in ("source_job_id", "replayed_fmax_eV_A",
                                             "Gamma_selected_groups", "TS_certified")}, indent=2))


if __name__ == "__main__":
    main()
