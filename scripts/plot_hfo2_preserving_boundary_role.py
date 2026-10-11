"""Compare two already audited preserving bands, without claiming TSs or strain response."""
from __future__ import annotations

import argparse
import csv
import json
from importlib.metadata import version
from pathlib import Path

import numpy as np
from ase.io import read

from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.analyze_hfo2_channel_network import read_evaluated_observation
from scripts.analyze_hfo2_clamped_reference import original_chart_images, pinned_json, terminal_images
from scripts.analyze_hfo2_network_update import structure_audit
from scripts.export_hfo2_chain_observation import replay
from scripts.export_hfo2_clamped_observation import CELL_SCALE
from scripts.audit_hfo2_static_replica import sha256
from vcneb.continuous_projection import continuous_reference_coordinates

BASE = Path("benchmarks/hfo2_channels/20261008")
FREE = BASE / "switching_converged_update_20261009_1255/preserving_step69"
CLAMPED = BASE / "clamped_held_flips_E100_20261011/zero_strain_preserving"
GATE = BASE / "G1_review_clamped_canary_20261009/G1_review_HF.json"
SEEDS = BASE / "clamped_endpoint_seeds/clamped_seed_manifest.json"
PINS = {
    "free": "2058e4045e57ddfa23c95d551593aa861e6ccca169b9e98d3f12bf112bd5b229",
    "clamped": "46145da14f24b2b8bd46aa4656f0982650d6b2dc661282f58d809aab8a6d3202",
    "gate": "6368cbf67e05fc1e61073fd1c5544b5bdccd572d9150865aa7ad4bdd0cd4ed0e",
    "seeds": "46c19e628b71e0a44a595be5f0ec9a2518e0719db0d1b286b294f5160f4af9a4",
}
STEM = "hfo2_preserving_boundary_role"


def build_data(repository):
    """Reuse exact observation/native readers and the unchanged symmetry sweep."""
    repository = Path(repository)
    pinned_json(repository/FREE/"observation.json", PINS["free"])
    gate = pinned_json(repository/GATE, PINS["gate"])
    seeds = pinned_json(repository/SEEDS, PINS["seeds"])
    if (not gate["candidate_G1_gate_passed"] or not gate["raw_on_HF"]
            or gate["physical_contract"] != CONTRACT or gate["fresh_raw_terminal_records"] != 40
            or seeds["physical_contract_sha256"] != CONTRACT):
        raise ValueError("historical native audit or fixed electronic contract differs")
    t_path = repository/BASE/"reference_variants/T.vasp"
    if sha256(t_path) != seeds["sources"]["T"]["sha256"]:
        raise ValueError("original T reference differs")
    tetragonal = read(t_path, format="vasp")
    free, old, _ = read_evaluated_observation(repository/FREE)
    clamped, new, residual, receipt = terminal_images(repository/CLAMPED, PINS["clamped"])
    if (old["source_job_id"] != "28319570" or old["snapshot_step"] != 69
            or receipt["job_id"] != "28727382" or receipt["terminal_step"] != 15
            or old["replayed_fmax_eV_A"] > .10 or residual["fmax_eV_A"] > .10):
        raise ValueError("the two actual ordinary terminal bands are required")
    reference_rows = np.asarray(seeds["unstrained_T_reference_cell_A"])[:2]
    if not np.allclose(new["mechanical_boundary"]["reference_cell_A"][:2],
                       reference_rows, atol=1e-9, rtol=0):
        raise ValueError("only the admitted zero-strain training boundary is used")
    transported = original_chart_images(clamped, seeds["orientation"], tetragonal)
    # The old free-cell metric is replayed; the clamp uses its registered metric.
    free_chain, _, _, _ = replay(free)
    bands = []
    for label, images, original, report, arc in (
            ("Free cell", free, free, old, free_chain.reaction_coordinate()),
            ("Clamped, zero strain", clamped, transported, new,
             np.asarray(new["extended_reaction_coordinate_A"]))):
        if (len(images) != 9 or report["n_active_images"] != 7
                or report["formula_units"] != 4 or report["pressure_GPa"] != 0
                or report["climb"]
                or report["cell_scale_A"] != (free_chain.cell_scale if label == "Free cell" else CELL_SCALE)
                or not np.isfinite(arc).all() or arc[0] != 0 or np.any(np.diff(arc) <= 0)):
            raise ValueError("nine complete ordered P=0 ordinary samples required")
        chart = continuous_reference_coordinates(original, tetragonal)
        rows = []
        for i, (image, point) in enumerate(zip(images, report["raw_image_evaluations"])):
            if any(point["input_sha256"][name] != digest for name, digest in CONTRACT.items()):
                raise ValueError("same six electronic bytes required at every frame")
            audit = structure_audit(image)
            deformation = chart["deformation_gradients"][i]
            diagnostic = report["path_diagnostics"]["images"][i]
            rows.append(dict(boundary=label, image_index=i, s_normalized=float(arc[i]/arc[-1]),
                energy_eV_cell=float(image.get_potential_energy()),
                energy_relative_own_PO_meV_fu=float((image.get_potential_energy()-images[0].get_potential_energy())*250),
                volume_A3=audit["volume_A3"], minimum_periodic_distance_A=audit["minimum_periodic_distance_A"],
                T_F_xx=float(deformation[0, 0]), T_F_yy=float(deformation[1, 1]), T_F_zz=float(deformation[2, 2]),
                raw_atomic_force_max_eV_A=float(np.linalg.norm(image.get_forces(), axis=1).max()),
                true_tangential_force_eV_A=diagnostic.get("true_tangential_force_eV_per_A"),
                ordinary_residual_eV_A=diagnostic.get("neb_residual_generalized_force_eV_per_A"),
                symprec_0p001=audit["symmetry_sweep"][0]["symbol"],
                symprec_0p01=audit["symmetry_sweep"][1]["symbol"],
                symprec_0p05=audit["symmetry_sweep"][2]["symbol"],
                raw_log_sha256=point["raw_log_sha256"], STRU_sha256=point["input_sha256"]["STRU"],
                structure_audit=audit))
        energies = np.array([r["energy_relative_own_PO_meV_fu"] for r in rows])
        bands.append(dict(label=label, source_job_id=report["source_job_id"], step=report["snapshot_step"],
            original_cell_scale_A=report["cell_scale_A"],
            ordinary_fmax_eV_A=report["replayed_fmax_eV_A"],
            sampled_peak_meV_fu=float(energies.max()), central_sample=rows[4],
            centre_minus_mean_adjacent_meV_fu=float(energies[4]-(energies[3]+energies[5])/2), rows=rows))
    return dict(status="descriptive_boundary_comparison_not_stationarity_or_prediction",
        bands=bands, physical_contract=CONTRACT, source_receipt_pins=PINS.copy(),
        T_reference_sha256=sha256(t_path), packages={n: version(n) for n in ("numpy", "ase", "pymatgen", "spglib")},
        free_union_sampling_peak_meV_fu=gate["sampling"]["preserving_union_peak_meV_fu"],
        same_electronic_contract=True, same_mechanical_boundary=False,
        common_initial_well_used=False, each_band_keeps_own_normalized_arc=True,
        new_DFT_calls=0, scheduler_mutations=0, G3_selection_or_holdout_access=False,
        TS_or_stable_intermediate_certified=False, smooth_MEP_interpolation_used=False,
        local_clamped_native_frames_reparsed=9, local_free_original_native_logs_reparsed=0,
        historical_free_native_audit_used=True, undistributed_six_physical_assets_rechecked_locally=False,
        image_records=18, independent_replicates=0)


def make_figure(data):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    with mpl.rc_context({"font.family": "DejaVu Sans", "font.size": 10,
                         "svg.fonttype": "none", "pdf.fonttype": 42}):
        fig, axes = plt.subplots(1, 2, figsize=(183/25.4, 88/25.4))
        fig.subplots_adjust(left=.09, right=.985, bottom=.32, top=.90, wspace=.28)
        handles = []
        for band, color, style in zip(data["bands"], ("#48729A", "#BA784A"), ("-", "--")):
            rows = band["rows"]
            x = [r["s_normalized"] for r in rows]
            for ax, field in zip(axes, ("energy_relative_own_PO_meV_fu", "T_F_xx")):
                y = [r[field] for r in rows]
                line, = ax.plot(x, y, style, color=color, marker="o", markersize=3.3, linewidth=1.5)
                ax.plot(x[4], y[4], marker="s", markersize=6.5, markerfacecolor="white", markeredgecolor=color)
                if ax is axes[0]:
                    handles.append(line)
        for label, ax in zip(("(a)", "(b)"), axes):
            ax.text(-.12, 1.065, label, transform=ax.transAxes, fontsize=13, fontweight="normal")
            ax.set_xlim(-.025, 1.025); ax.set_xticks([0, .25, .5, .75, 1])
            ax.set_xlabel("Normalized path arc", fontsize=10)
            ax.tick_params(direction="in", top=True, right=True, labelsize=10)
            for spine in ax.spines.values():
                spine.set_visible(True)
        axes[0].set_ylabel(r"$E-E_{\mathrm{PO^+}}$ (meV/HfO$_2$)", fontsize=10)
        axes[0].set_ylim(-22, 108); axes[0].axhline(0, color=".75", linewidth=.6, zorder=0)
        axes[1].set_ylabel(r"T-chart stretch $F_{xx}$", fontsize=10)
        axes[1].set_ylim(.99, 1.14); axes[1].set_yticks([1, 1.04, 1.08, 1.12])
        marker = Line2D([], [], color=".4", linestyle="none", marker="s", markerfacecolor="white", markersize=6.5)
        fig.legend(handles+[marker], [b["label"] for b in data["bands"]]+["Pbcn sample"],
                   loc="lower center", bbox_to_anchor=(.53, .025), ncol=3, fontsize=10,
                   frameon=True, facecolor="white", edgecolor=".55", framealpha=.85)
        fig.canvas.draw()
        return fig


def export(data, output):
    output = Path(output)
    if output.exists():
        raise FileExistsError("fresh boundary-comparison output required")
    fig = make_figure(data)
    output.mkdir(parents=True)
    rows = [{k: v for k, v in row.items() if k != "structure_audit"}
            for band in data["bands"] for row in band["rows"]]
    with (output/"source_data.csv").open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    for extension in ("pdf", "svg", "png", "tiff"):
        options = {"pil_kwargs": {"compression": "tiff_lzw"}} if extension == "tiff" else {}
        fig.savefig(output/f"{STEM}.{extension}", dpi=600 if extension == "tiff" else 300, **options)
    data = dict(data, script_sha256=sha256(Path(__file__)), visual_review="pending",
                figure_dimensions_mm=(fig.get_size_inches()*25.4).tolist(),
                file_sha256={p.name: sha256(p) for p in sorted(output.iterdir())})
    with (output/"analysis.json").open("x", encoding="utf-8") as stream:
        json.dump(data, stream, indent=2, allow_nan=False); stream.write("\n")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = export(build_data(args.repository), args.output)
    print(json.dumps({b["label"]: {"centre_meV_fu": b["central_sample"]["energy_relative_own_PO_meV_fu"],
                                "centre_T_F_xx": b["central_sample"]["T_F_xx"]} for b in result["bands"]}))


if __name__ == "__main__":
    main()
