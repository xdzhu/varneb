"""Three pinned zero-strain training profiles, not a forecast or smooth MEP."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from scripts.plot_hfo2_clamped_comparison import (
    BASE, FRAME_IDS, SPANS, STYLE, T_OBSERVATION, build_data as build_TM,
    build_band_data, canonical_sha, sha256,
)

DEPENDENCY_SHA256 = {
    # E105 adds opt-in M strain limits; the legacy renderer default is unchanged.
    # Re-pin this reviewed code upgrade, not the historical figure/data bytes.
    "scripts/plot_hfo2_clamped_comparison.py": "01f9842bd1cd860909dd6847330d732a0bd6bb0c3f4641daa42b83f77a38d664",
    "scripts/plot_hfo2_clamped_M.py": "3aadfdda38c0974acf16908015ad802010823ea06c43bf66f0e34bfff3a86612",
}
R_REPORT = BASE / "clamped_reversing_reference_E085_20261011/displacement.json"
R_GRADIENT = BASE / "clamped_reversing_reference_E085_20261011/gradient.json"
R_OBSERVATION = BASE / "clamped_reversing_terminal_E083_20261011/observations/step_0005/observation.json"
R_AUDIT = BASE / "clamped_reversing_terminal_E083_20261011/audit_receipt.json"
R_REPORT_SHA = "5c513db4766b19ae8bc9950972d469006801dc7ce2382b6f05988f5afb63308b"
R_GRADIENT_SHA = "23ff35ed551bb3fb59a372d72afd9ca98849bb269fd125d9bf04a475304a630c"
R_OBSERVATION_SHA = "669fea17d95a57ba6a2acf5d228806173370a0b1a6a51388aa71d57ee6471b85"
R_AUDIT_SHA = "632e0e9b2b821b3bfc7d66b361c44c7db734988fb5d018ea2f31441fab14ed41"
STEM = "hfo2_clamped_TMR"


def build_data(repository):
    """Reuse the frozen T/M adapter, adding only the audited reversing band."""
    repository = Path(repository)
    for relative, digest in DEPENDENCY_SHA256.items():
        if sha256(repository / relative) != digest:
            raise ValueError("pinned historical plotting dependency changed")
    data = build_TM(repository)
    band = build_band_data(repository, R_REPORT, R_OBSERVATION, R_AUDIT,
        mode_digest=R_REPORT_SHA, observation_digest=R_OBSERVATION_SHA,
        audit_digest=R_AUDIT_SHA, job="28692776", step=5)
    gradient = json.loads((repository / R_GRADIENT).read_text(encoding="utf-8"))
    if (canonical_sha(gradient) != R_GRADIENT_SHA
            or gradient["terminal_audit_sha256"] != R_AUDIT_SHA
            or gradient["source_job_id"] != "28692776" or gradient["terminal_step"] != 5
            or gradient["strain_training_condition"] != 0
            or not gradient["all_E023_Cmma_frames_retained"]
            or gradient["G3_selection_or_holdout_access"] or gradient["new_DFT_calls"]
            or gradient["scheduler_mutations"] or len(gradient["images"]) != 9):
        raise ValueError("frozen reversing gradient or training scope differs")
    for i, (row, grad) in enumerate(zip(band["rows"], gradient["images"])):
        if (grad["image_index"] != i or set(grad["representations"]) != set(SPANS)
                or not np.isclose(row["energy_relative_PO_meV_fu"], grad["relative_energy_meV_fu"], atol=1e-9, rtol=0)):
            raise ValueError("reversing gradient/image matching differs")
        row.update(channel="reversing", source_job_id="28692776", snapshot_step=5,
                   zero_T_displacement_denominator=band["zero_displacement"][i])
        for span in SPANS:
            result = grad["representations"][span]
            if (result["gradient_unit"] != "eV/(sqrt(amu)*Angstrom)"
                    or result["source_basis_rank"] != result["translation_free_basis_rank"]):
                raise ValueError("reversing gradient unit or rank differs")
            for key in ("retained_squared_gradient_fraction", "internal_gradient_norm",
                        "retained_gradient_norm", "omitted_gradient_norm",
                        "rigid_translation_gradient_norm", "zero_internal_gradient"):
                row[f"{span}_{key}"] = result[key]
    band.update(channel="reversing", gradient_canonical_sha256=R_GRADIENT_SHA)
    band["source_file_sha256"][str(R_GRADIENT)] = sha256(repository / R_GRADIENT)
    data["bands"].append(band)

    # Equality of normalized arcs would not establish a shared initial state.
    observations = [json.loads((repository / Path(b["source_paths"][1])).read_text(encoding="utf-8"))
                    for b in data["bands"]]
    first = observations[0]
    keys = ("pressure_GPa", "cell_scale_A", "k_eV_A2", "formula_units",
            "fmax_target_eV_A", "climb", "n_total_images", "n_active_images")
    boundary = {k: v for k, v in first["mechanical_boundary"].items()
                if k not in ("reference_file", "reference_file_sha256")}
    native_keys = ("input_sha256", "raw_log_sha256", "energy_eV_cell",
                   "forces_eV_A", "stress_ASE_voigt_eV_A3")
    initial = {k: first["raw_image_evaluations"][0][k] for k in native_keys}
    from ase.io import read
    geometries = []
    for b, obs in zip(data["bands"], observations):
        if any(obs[key] != first[key] for key in keys):
            raise ValueError("same zero-strain ensemble required")
        if boundary != {k: v for k, v in obs["mechanical_boundary"].items()
                        if k not in ("reference_file", "reference_file_sha256")}:
            raise ValueError("same mechanical boundary required")
        if initial != {k: obs["raw_image_evaluations"][0][k] for k in native_keys}:
            raise ValueError("same ordered cached native PO record required")
        snapshot = repository / Path(b["source_paths"][1]).parent / "POSCAR_00"
        if sha256(snapshot) != obs["raw_image_evaluations"][0]["snapshot_POSCAR_sha256"]:
            raise ValueError("pinned initial snapshot changed")
        geometries.append(read(snapshot, format="vasp"))
    for g in geometries[1:]:
        if (g.get_chemical_symbols() != geometries[0].get_chemical_symbols()
                or not np.allclose(g.cell.array, geometries[0].cell.array, atol=1e-12, rtol=0)
                or not np.allclose(g.positions, geometries[0].positions, atol=1e-12, rtol=0)):
            raise ValueError("same ordered initial geometry required; no repair")
    if canonical_sha(initial) != data["common_initial_record_canonical_sha256"]:
        raise ValueError("historical common native record differs")
    identities = {(r["raw_log_sha256"], r["input_sha256"]["STRU"])
                  for obs in observations for r in obs["raw_image_evaluations"]}
    data.update(image_records=27, distinct_manifest_native_log_STRU_identities=len(identities),
        source_dependency_sha256=DEPENDENCY_SHA256.copy(), independently_audited_new_SCFs=0,
        training_condition=0., reversing_added_without_reference_or_input_change=True)
    return data


def make_figure(data):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 9.5, "axes.labelsize": 10.5, "legend.fontsize": 9.5,
        "axes.linewidth": .9, "svg.fonttype": "none", "pdf.fonttype": 42})
    fig, axes = plt.subplots(3, 3, figsize=(183/25.4, 207/25.4), sharex="col", sharey="row")
    fig.subplots_adjust(left=.13, right=.975, bottom=.085, top=.91, wspace=.17, hspace=.58)
    keys = ["T_geometric_triplet_rank3_squared_fraction", "T_optical_doublets_rank4_squared_fraction",
            *[f"Cmma_frame_{i}_rank4_squared_fraction" for i in FRAME_IDS]]
    labels = [r"$\mathrm{PO^+}\rightarrow\mathrm{T}$", r"$\mathrm{PO^+}\rightarrow\mathrm{M}$",
              r"$\mathrm{PO^+}\rightarrow\mathrm{PO^-}$" + "\nReversing pattern"]
    for col, band in enumerate(data["bands"]):
        rows = band["rows"]
        x = np.array([r["normalized_generalized_arc"] for r in rows])
        axes[0, col].plot(x, [r["energy_relative_PO_meV_fu"] for r in rows],
            color="#256B91", marker="o", markersize=3.8, lw=1.4, label=labels[col])
        axes[0, col].axhline(0., color="#CDD3D8", lw=.6, zorder=0)
        axes[0, col].set_ylim(-115, 370)
        axes[0, col].set_yticks([-100, 0, 100, 200, 300])
        axes[0, col].legend(loc="lower center", bbox_to_anchor=(.5, 1.075), frameon=True,
            facecolor="white", edgecolor="#9CA3AA", framealpha=.85, borderpad=.3)
        for j, (key, span) in enumerate(zip(keys, SPANS)):
            name, color, marker, linestyle = STYLE[min(j, 2)]
            yd = np.array([100*r[key] for r in rows])
            yd[np.array(band["zero_displacement"], dtype=bool)] = np.nan
            yg = np.array([100*r[f"{span}_retained_squared_gradient_fraction"] for r in rows])
            yg[[r[f"{span}_zero_internal_gradient"] for r in rows]] = np.nan
            for row, values in ((1, yd), (2, yg)):
                axes[row, col].plot(x, values, color=color, ls=linestyle,
                    marker=marker if j < 3 else None, markersize=3., lw=1.15 if j < 2 else .9,
                    label=name if j < 3 else None)
        for row in (1, 2):
            axes[row, col].set_ylim(-4, 104)
            axes[row, col].set_yticks([0, 25, 50, 75, 100])
        axes[2, col].set_xlabel("Normalized arc")
        for row in range(3):
            ax = axes[row, col]
            ax.set_xlim(-.025, 1.025)
            ax.set_xticks([0., .5, 1.])
            ax.tick_params(direction="in", top=True, right=True, length=3.6)
            ax.axvline(x[band["highest_image_index"]], color="#C97738", ls="--", lw=.85, zorder=0)
            ax.text(-.22, 1.025, f"({chr(97+3*row+col)})", transform=ax.transAxes,
                    fontsize=13, fontweight="normal", va="bottom")
            for spine in ax.spines.values():
                spine.set_visible(True)
    axes[0, 0].set_ylabel(r"$E-E_{\rm PO^+}$ (meV/f.u.)")
    axes[1, 0].set_ylabel("Squared displacement\ncoverage (%)")
    axes[2, 0].set_ylabel("Squared atomic-gradient\ncoverage (%)")
    handles, _ = axes[1, 0].get_legend_handles_labels()
    fig.legend(handles, ["T geometric (rank 3)", "T optical (rank 4)", "Cmma (rank 4; 4 frames)"],
        loc="center", bbox_to_anchor=(.55, .657), ncols=3, frameon=True, framealpha=.85,
        facecolor="white", edgecolor="#9CA3AA", borderpad=.4, columnspacing=.7, handlelength=1.4)
    fig.align_ylabels(axes[:, 0])
    fig.canvas.draw()
    return fig


def export(data, output):
    output = Path(output)
    if output.exists():
        raise FileExistsError("fresh dated three-band figure directory required")
    fig = make_figure(data)
    output.mkdir(parents=True)
    rows = [r for b in data["bands"] for r in b["rows"]]
    if any(set(r) != set(rows[0]) for r in rows):
        raise ValueError("same complete scalar fields required for all bands")
    source = output / "source_data.csv"
    with source.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    for extension in ("pdf", "svg", "png", "tiff"):
        options = {"pil_kwargs": {"compression": "tiff_lzw"}} if extension == "tiff" else {}
        fig.savefig(output / f"{STEM}.{extension}", dpi=600 if extension == "tiff" else 300,
            facecolor="white", **options)
    positions = np.array([a.get_position().bounds for a in fig.axes]).reshape(3, 3, 4)
    qa = {k: v for k, v in data.items() if k != "bands"}
    qa.update(source_bands=[{k: v for k, v in b.items() if k != "rows"} for b in data["bands"]],
        plot_script_sha256=sha256(Path(__file__)), source_data_sha256=sha256(source),
        independent_replicates=0, figure_dimensions_mm=(fig.get_size_inches()*25.4).tolist(),
        panel_positions_normalized=positions.tolist(),
        grid_aligned=bool(np.allclose(positions[:, :, [0, 2]], positions[0:1, :, [0, 2]], atol=1e-12, rtol=0)
            and np.allclose(positions[:, :, [1, 3]], positions[:, 0:1, [1, 3]], atol=1e-12, rtol=0)),
        panel_labels=[f"({chr(97+i)})" for i in range(9)], panel_label_weight="normal",
        panel_titles_empty=all(not a.get_title() for a in fig.axes),
        all_spines_visible=all(s.get_visible() for a in fig.axes for s in a.spines.values()),
        energy_ylim_meV_fu=[-115, 370], legends_outside_axes=True, legend_frame_alpha=.85,
        all_four_Cmma_frames_retained=True, Cmma_spread_is_statistical_uncertainty=False,
        zero_denominator_fractions_rendered_as_missing=True, visual_review="pending",
        file_sha256={p.name: sha256(p) for p in sorted(output.glob(STEM+".*"))})
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
    qa = export(build_data(args.repository), args.output)
    print(json.dumps({"records": qa["image_records"], "grid_aligned": qa["grid_aligned"], "new_DFT_calls": 0}))


if __name__ == "__main__":
    main()
