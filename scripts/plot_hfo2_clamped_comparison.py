"""Pinned T/M training observations: geometry coverage is not gradient coverage."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from scripts.plot_hfo2_clamped_M import (
    BASE, FRAME_IDS, AUDIT, AUDIT_SHA, MODE_REPORT, MODE_CANONICAL_SHA,
    OBSERVATION, OBSERVATION_SHA, build_band_data, canonical_sha, sha256,
)

T_REPORT = BASE / "clamped_T_reference_E080_20261011/displacement.json"
T_GRADIENT = BASE / "clamped_T_reference_E080_20261011/gradient.json"
M_GRADIENT = BASE / "clamped_mode_gradient_E078_20261011/analysis.json"
T_OBSERVATION = BASE / "clamped_T_terminal_E079_20261011/observations/step_0010/observation.json"
T_AUDIT = BASE / "clamped_T_terminal_E079_20261011/audit_receipt.json"
T_AUDIT_SHA = "c5653b4759b202199cb75739f99038aaa52b3ccc220ef226efbacc980e5d716e"
T_REPORT_SHA = "2a98aa799bd75f60481acd1725f249937d7d50db5b697b9aad39a2ad434fbd35"
T_OBSERVATION_SHA = "f0628260d2daf3aabac1860a92f17200d00f999609cbdd7dbd93d748bdb5b7d3"
T_GRADIENT_SHA = "8e040045e3c971c6988a19f78401000ea8c73c50c1fb055a1cb3b9a5151bfc48"
M_GRADIENT_SHA = "434332f3f4f267e8f377e0aa4c3339488335f371ea99c4d171d098cc802ae0f4"
SPANS = ["T_geometric_triplet", "T_two_lowest_complete_optical_doublets",
         *[f"Cmma_frame_{i}" for i in FRAME_IDS]]
STYLE = [
    ("T geometric (rank 3)", "#256B91", "o", "-"),
    ("T optical (rank 4)", "#757D85", "s", "--"),
    ("Cmma (rank 4; 4 frames)", "#23877E", "D", ":"),
]


def build_data(repository):
    specs = [
        ("T", T_REPORT, T_OBSERVATION, T_AUDIT, T_REPORT_SHA,
         T_OBSERVATION_SHA, T_AUDIT_SHA, "28692775", 10, T_GRADIENT, T_GRADIENT_SHA),
        ("M", MODE_REPORT, OBSERVATION, AUDIT, MODE_CANONICAL_SHA,
         OBSERVATION_SHA, AUDIT_SHA, "28722320", 13, M_GRADIENT, M_GRADIENT_SHA),
    ]
    bands, observations, initial_geometries = [], [], []
    for channel, mode_path, obs_path, audit_path, mode_pin, obs_pin, audit_pin, job, step, grad_path, grad_pin in specs:
        band = build_band_data(repository, mode_path, obs_path, audit_path,
            mode_digest=mode_pin, observation_digest=obs_pin, audit_digest=audit_pin, job=job, step=step)
        gradient = json.loads((repository / grad_path).read_text(encoding="utf-8"))
        if (canonical_sha(gradient) != grad_pin or gradient["terminal_audit_sha256"] != audit_pin
                or gradient["source_job_id"] != job or gradient["terminal_step"] != step
                or gradient["strain_training_condition"] != 0
                or not gradient["all_E023_Cmma_frames_retained"]
                or gradient["G3_selection_or_holdout_access"] or gradient["new_DFT_calls"]
                or gradient["scheduler_mutations"] or len(gradient["images"]) != 9):
            raise ValueError("frozen gradient or training scope differs")
        for i, (row, grad) in enumerate(zip(band["rows"], gradient["images"])):
            if (grad["image_index"] != i or set(grad["representations"]) != set(SPANS)
                    or not np.isclose(row["energy_relative_PO_meV_fu"], grad["relative_energy_meV_fu"], atol=1e-9, rtol=0)):
                raise ValueError("gradient/image matching differs")
            row.update(channel=channel, source_job_id=job, snapshot_step=step,
                       zero_T_displacement_denominator=band["zero_displacement"][i])
            for span in SPANS:
                result = grad["representations"][span]
                if (result["gradient_unit"] != "eV/(sqrt(amu)*Angstrom)"
                        or result["source_basis_rank"] != result["translation_free_basis_rank"]):
                    raise ValueError("gradient unit or rank differs")
                for key in ("retained_squared_gradient_fraction", "internal_gradient_norm",
                            "retained_gradient_norm", "omitted_gradient_norm",
                            "rigid_translation_gradient_norm", "zero_internal_gradient"):
                    row[f"{span}_{key}"] = result[key]
        band.update(channel=channel, gradient_canonical_sha256=grad_pin)
        band["source_file_sha256"][str(grad_path)] = sha256(repository / grad_path)
        bands.append(band)
        obs = json.loads((repository / obs_path).read_text(encoding="utf-8"))
        snapshot = repository / obs_path.parent / "POSCAR_00"
        if sha256(snapshot) != obs["raw_image_evaluations"][0]["snapshot_POSCAR_sha256"]:
            raise ValueError("pinned initial snapshot changed")
        from ase.io import read
        initial_geometries.append(read(snapshot, format="vasp"))
        observations.append(obs)
    # Equal normalized arcs are not equal absolute positions. Compare the
    # actual native initial record/ensemble, not those normalized coordinates.
    for key in ("pressure_GPa", "cell_scale_A", "k_eV_A2",
                "formula_units", "fmax_target_eV_A", "climb", "n_total_images", "n_active_images"):
        if observations[0][key] != observations[1][key]:
            raise ValueError("same-ensemble T/M comparison required")
    # Native reference filenames/comments differ, not the boundary itself.
    provenance = {"reference_file", "reference_file_sha256"}
    boundaries = [{k: v for k, v in obs["mechanical_boundary"].items() if k not in provenance}
                  for obs in observations]
    if boundaries[0] != boundaries[1]:
        raise ValueError("same mechanical boundary required")
    # POSCAR last-digit serialization differs by ~4e-16 A. Each byte hash is
    # checked separately; common SCF bytes/EFS and ordered geometry must match.
    left, right = initial_geometries
    if (left.get_chemical_symbols() != right.get_chemical_symbols()
            or not np.allclose(left.cell.array, right.cell.array, atol=1e-12, rtol=0)
            or not np.allclose(left.positions, right.positions, atol=1e-12, rtol=0)):
        raise ValueError("same ordered initial geometry required")
    initial_keys = ("input_sha256", "raw_log_sha256",
                    "energy_eV_cell", "forces_eV_A", "stress_ASE_voigt_eV_A3")
    initial = [{k: obs["raw_image_evaluations"][0][k] for k in initial_keys} for obs in observations]
    if initial[0] != initial[1]:
        raise ValueError("same ordered cached PO+ native record required")
    return dict(bands=bands, common_initial_record_canonical_sha256=canonical_sha(initial[0]),
        common_initial_energy_eV_cell=initial[0]["energy_eV_cell"], new_DFT_calls=0,
        initial_snapshot_byte_hashes_checked_separately=True,
        initial_ordered_geometry_match_tolerance_A=1e-12,
        scheduler_mutations=0, G3_selection_or_holdout_access=False, TS_certified=False,
        native_data_reaudit_claimed=False, smooth_MEP_interpolation_used=False,
        source_contract="P0/100Ry/full10auDZP/PBE/Gamma2x2x2; two fixed T-plane rows; open row/tilts",
        normalized_arcs_are_individual_not_common_absolute_coordinates=True)


def make_figure(data):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 10, "axes.labelsize": 11, "legend.fontsize": 10,
        "axes.linewidth": .9, "svg.fonttype": "none", "pdf.fonttype": 42})
    fig, axes = plt.subplots(3, 2, figsize=(183 / 25.4, 202 / 25.4), sharex="col", sharey="row")
    fig.subplots_adjust(left=.115, right=.975, bottom=.085, top=.92, wspace=.18, hspace=.57)
    displacement_keys = ["T_geometric_triplet_rank3_squared_fraction", "T_optical_doublets_rank4_squared_fraction",
                         *[f"Cmma_frame_{f}_rank4_squared_fraction" for f in FRAME_IDS]]
    for col, band in enumerate(data["bands"]):
        rows = band["rows"]
        x = np.array([r["normalized_generalized_arc"] for r in rows])
        peak = band["highest_image_index"]
        axes[0, col].plot(x, [r["energy_relative_PO_meV_fu"] for r in rows],
            color="#256B91", marker="o", markersize=4, lw=1.5,
            label=band.get("plot_label", rf"$\mathrm{{PO^+}}\rightarrow\mathrm{{{band['channel']}}}$"))
        axes[0, col].axhline(0, lw=.6, color="#CDD3D8", zorder=0)
        axes[0, col].set_ylim(*data.get("energy_limits_meV_fu", (-115, 125)))
        axes[0, col].set_yticks([-100, -50, 0, 50, 100])
        axes[0, col].legend(loc="lower center", bbox_to_anchor=(.5, 1.04),
            framealpha=.85, facecolor="white", edgecolor="#9CA3AA", borderpad=.3)
        for j, (disp, span) in enumerate(zip(displacement_keys, SPANS)):
            label, color, marker, linestyle = STYLE[min(j, 2)]
            y_disp = np.array([100*r[disp] for r in rows])
            y_disp[np.array(band["zero_displacement"], dtype=bool)] = np.nan
            y_grad = np.array([100*r[f"{span}_retained_squared_gradient_fraction"] for r in rows])
            y_grad[[bool(r[f"{span}_zero_internal_gradient"]) for r in rows]] = np.nan
            for row_index, values in ((1, y_disp), (2, y_grad)):
                axes[row_index, col].plot(x, values, color=color, ls=linestyle,
                    marker=marker if j < 3 else None, markersize=3.4, lw=1.2 if j < 2 else .9,
                    label=label if j < 3 else None)
        for row_index in (1, 2):
            axes[row_index, col].set_ylim(-4, 104)
            axes[row_index, col].set_yticks([0, 25, 50, 75, 100])
        axes[2, col].set_xlabel("Normalized generalized arc")
        for row_index in range(3):
            ax = axes[row_index, col]
            ax.set_xlim(-.025, 1.025)
            ax.set_xticks([0, .25, .5, .75, 1.])
            ax.tick_params(direction="in", top=True, right=True, length=4)
            ax.axvline(x[peak], color="#C97738", lw=.85, ls="--", zorder=0)
            ax.text(-.16, 1.045, f"({chr(97+2*row_index+col)})", transform=ax.transAxes,
                    fontsize=13, fontweight="normal", va="bottom")
            for spine in ax.spines.values():
                spine.set_visible(True)
    axes[0, 0].set_ylabel(r"$E-E_{\rm PO^+}$ (meV/f.u.)")
    axes[1, 0].set_ylabel("Squared displacement\ncoverage (%)")
    axes[2, 0].set_ylabel("Squared atomic-gradient\ncoverage (%)")
    handles, labels = axes[1, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center", bbox_to_anchor=(.56, .666), ncols=3,
               framealpha=.85, facecolor="white", edgecolor="#9CA3AA", borderpad=.4,
               columnspacing=.8, handlelength=1.5)
    fig.align_ylabels(axes[:, 0])
    fig.canvas.draw()
    return fig


def export(data, output, *, stem_name="hfo2_clamped_TM"):
    if Path(stem_name).name != stem_name or Path(stem_name).suffix:
        raise ValueError("a plain extension-free figure stem is required")
    if output.exists():
        raise FileExistsError("fresh dated comparison figure directory required")
    fig = make_figure(data)
    output.mkdir(parents=True)
    rows = [row for band in data["bands"] for row in band["rows"]]
    source = output / "source_data.csv"
    with source.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    stem = output / stem_name
    for ext in ("pdf", "svg", "png", "tiff"):
        extra = {"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tiff" else {}
        fig.savefig(stem.with_suffix("." + ext), dpi=600 if ext == "tiff" else 300,
                    facecolor="white", **extra)
    positions = np.array([ax.get_position().bounds for ax in fig.axes]).reshape(3, 2, 4)
    qa = {k: v for k, v in data.items() if k != "bands"}
    qa.update(source_bands=[{k: v for k, v in band.items() if k != "rows"} for band in data["bands"]],
        plot_script_sha256=sha256(Path(__file__)), source_data_sha256=sha256(source),
        image_records=18, independent_replicates=0,
        figure_dimensions_mm=(fig.get_size_inches()*25.4).tolist(),
        panel_positions_normalized=positions.tolist(),
        grid_aligned=bool(np.allclose(positions[:, :, [0, 2]], positions[0:1, :, [0, 2]], atol=1e-12, rtol=0)
                          and np.allclose(positions[:, :, [1, 3]], positions[:, 0:1, [1, 3]], atol=1e-12, rtol=0)),
        panel_labels=[f"({chr(97+i)})" for i in range(6)], panel_label_weight="normal",
        all_spines_visible=all(s.get_visible() for ax in fig.axes for s in ax.spines.values()),
        panel_titles_empty=all(not ax.get_title() for ax in fig.axes),
        legends_outside_axes=True, legend_frame_alpha=.85,
        zero_displacement_fraction_rendered_as_missing=True, all_four_Cmma_frames_retained=True,
        Cmma_frame_spread_is_statistical_uncertainty=False, visual_review="pending",
        file_sha256={p.name: sha256(p) for p in sorted(output.glob(stem.name + ".*"))})
    with (output / "qa.json").open("x", encoding="utf-8") as stream:
        json.dump(qa, stream, indent=2, allow_nan=False); stream.write("\n")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return qa


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = export(build_data(args.repository), args.output)
    print(json.dumps({"records": result["image_records"], "grid_aligned": result["grid_aligned"], "new_DFT_calls": 0}))


if __name__ == "__main__":
    main()
