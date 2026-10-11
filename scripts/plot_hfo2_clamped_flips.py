"""Two pinned ordered zero-strain flips; descriptive samples, not TS certificates."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from scripts import plot_hfo2_clamped_comparison as comparison
from scripts import plot_hfo2_clamped_three_bands as historical
from scripts.plot_hfo2_clamped_M import BASE, build_band_data, canonical_sha, sha256

DEPENDENCIES = {
    **historical.DEPENDENCY_SHA256,
    "scripts/plot_hfo2_clamped_three_bands.py": "1bad21a1315baea2835d42ad62325f2ba4b45a29cd706299f78e19c215c19df3",
}
CASE = BASE / "clamped_preserving_reference_E101_20261011"
TERMINAL = BASE / "clamped_held_flips_E100_20261011/zero_strain_preserving"
MODE = CASE / "displacement.json"
GRADIENT = CASE / "gradient.json"
OBSERVATION = TERMINAL / "observations/step_0015/observation.json"
AUDIT = TERMINAL / "audit_receipt.json"
MODE_SHA = "a99ceff3c8a77653fad04b47bff0f4bcf4ec246da81fa877d97c554aaef53966"
GRADIENT_SHA = "15b48f47eccfda196f84f9c8daa9572e9cd0e427f5e4c1ae0f22283c64e5e820"
OBSERVATION_SHA = "babd26a259339ac4a8fb9f0154895b7292294e97bdc8f2ad478fbfc2312e8995"
AUDIT_SHA = "46145da14f24b2b8bd46aa4656f0982650d6b2dc661282f58d809aab8a6d3202"
STEM = "hfo2_clamped_flips"


def build_data(repository):
    """Consume admitted scalar reports; never launch a calculator or read holdouts."""
    repository = Path(repository)
    for relative, digest in DEPENDENCIES.items():
        if sha256(repository / relative) != digest:
            raise ValueError("pinned historical plotting dependency changed")
    # This existing adapter retains all four frames and checks the old records.
    old = historical.build_data(repository)
    reversing = old["bands"][2]
    preserving = build_band_data(repository, MODE, OBSERVATION, AUDIT,
        mode_digest=MODE_SHA, observation_digest=OBSERVATION_SHA,
        audit_digest=AUDIT_SHA, job="28727382", step=15)
    gradient = json.loads((repository / GRADIENT).read_text(encoding="utf-8"))
    if (canonical_sha(gradient) != GRADIENT_SHA
            or gradient["terminal_audit_sha256"] != AUDIT_SHA
            or gradient["source_job_id"] != "28727382" or gradient["terminal_step"] != 15
            or gradient["strain_training_condition"] != 0
            or not gradient["all_E023_Cmma_frames_retained"]
            or gradient["G3_selection_or_holdout_access"] or gradient["new_DFT_calls"]
            or gradient["scheduler_mutations"] or len(gradient["images"]) != 9):
        raise ValueError("frozen preserving gradient or training scope differs")
    for i, (row, grad) in enumerate(zip(preserving["rows"], gradient["images"])):
        if (grad["image_index"] != i or set(grad["representations"]) != set(comparison.SPANS)
                or not np.isclose(row["energy_relative_PO_meV_fu"],
                                  grad["relative_energy_meV_fu"], atol=1e-9, rtol=0)):
            raise ValueError("preserving gradient/image matching differs")
        row.update(channel="preserving", source_job_id="28727382", snapshot_step=15,
                   zero_T_displacement_denominator=preserving["zero_displacement"][i])
        for span in comparison.SPANS:
            result = grad["representations"][span]
            if (result["gradient_unit"] != "eV/(sqrt(amu)*Angstrom)"
                    or result["source_basis_rank"] != result["translation_free_basis_rank"]):
                raise ValueError("gradient unit or rank differs")
            for key in ("retained_squared_gradient_fraction", "internal_gradient_norm",
                        "retained_gradient_norm", "omitted_gradient_norm",
                        "rigid_translation_gradient_norm", "zero_internal_gradient"):
                row[f"{span}_{key}"] = result[key]
    preserving.update(channel="preserving", gradient_canonical_sha256=GRADIENT_SHA)
    preserving["source_file_sha256"][str(GRADIENT)] = sha256(repository / GRADIENT)
    bands = [preserving, reversing]
    observations = [json.loads((repository / Path(b["source_paths"][1])).read_text(encoding="utf-8"))
                    for b in bands]
    first = observations[0]
    from ase.io import read
    geometries = []
    ensemble = ("pressure_GPa", "cell_scale_A", "k_eV_A2", "formula_units",
                "fmax_target_eV_A", "climb", "n_total_images", "n_active_images")
    boundary = lambda obs: {k: v for k, v in obs["mechanical_boundary"].items()
                           if k not in ("reference_file", "reference_file_sha256")}
    native_keys = ("input_sha256", "raw_log_sha256", "energy_eV_cell",
                   "forces_eV_A", "stress_ASE_voigt_eV_A3")
    initial = {k: first["raw_image_evaluations"][0][k] for k in native_keys}
    for band, obs in zip(bands, observations):
        if any(obs[k] != first[k] for k in ensemble) or boundary(obs) != boundary(first):
            raise ValueError("same zero-strain mechanical/electronic ensemble required")
        if initial != {k: obs["raw_image_evaluations"][0][k] for k in native_keys}:
            raise ValueError("same ordered cached native PO+ record required")
        snapshot = repository / Path(band["source_paths"][1]).parent / "POSCAR_00"
        if sha256(snapshot) != obs["raw_image_evaluations"][0]["snapshot_POSCAR_sha256"]:
            raise ValueError("pinned initial snapshot changed")
        geometries.append(read(snapshot, format="vasp"))
    left, right = geometries
    if (left.get_chemical_symbols() != right.get_chemical_symbols()
            or not np.allclose(left.cell.array, right.cell.array, atol=1e-12, rtol=0)
            or not np.allclose(left.positions, right.positions, atol=1e-12, rtol=0)):
        raise ValueError("same ordered initial geometry required; no remapping")
    if canonical_sha(initial) != old["common_initial_record_canonical_sha256"]:
        raise ValueError("historical common initial native record differs")
    identities = {(r["raw_log_sha256"], r["input_sha256"]["STRU"])
                  for obs in observations for r in obs["raw_image_evaluations"]}
    for band, label in zip(bands, ("Preserving pattern", "Reversing pattern")):
        band["plot_label"] = r"$\mathrm{PO^+}\rightarrow\mathrm{PO^-}$" + "\n" + label
    data = {k: v for k, v in old.items() if k not in ("bands", "image_records",
            "distinct_manifest_native_log_STRU_identities", "reversing_added_without_reference_or_input_change")}
    data.update(bands=bands, image_records=18,
        distinct_manifest_native_log_STRU_identities=len(identities),
        source_dependency_sha256=DEPENDENCIES.copy(),
        ordered_minus_endpoints_are_not_merged=True)
    return data


def make_figure(data):
    """Keep the existing six-panel style; only set the preregistered energy scale."""
    figure = comparison.make_figure(data)
    figure.set_size_inches(183 / 25.4, 195 / 25.4)
    for ax in figure.axes[:2]:
        ax.set_ylim(-10, 370)
        ax.set_yticks([0, 100, 200, 300])
    figure.canvas.draw()
    return figure


def export(data, output):
    output = Path(output)
    if output.exists():
        raise FileExistsError("fresh dated two-flip figure directory required")
    rows = [r for band in data["bands"] for r in band["rows"]]
    if any(set(r) != set(rows[0]) for r in rows):
        raise ValueError("same complete scalar fields required")
    fig = make_figure(data)
    output.mkdir(parents=True)
    source = output / "source_data.csv"
    with source.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    for extension in ("pdf", "svg", "png", "tiff"):
        options = {"pil_kwargs": {"compression": "tiff_lzw"}} if extension == "tiff" else {}
        fig.savefig(output / f"{STEM}.{extension}", dpi=600 if extension == "tiff" else 300,
                    facecolor="white", **options)
    positions = np.array([ax.get_position().bounds for ax in fig.axes]).reshape(3, 2, 4)
    qa = {k: v for k, v in data.items() if k != "bands"}
    qa.update(source_bands=[{k: v for k, v in b.items() if k != "rows"} for b in data["bands"]],
        plot_script_sha256=sha256(Path(__file__)), source_data_sha256=sha256(source),
        independent_replicates=0, figure_dimensions_mm=(fig.get_size_inches()*25.4).tolist(),
        panel_positions_normalized=positions.tolist(),
        grid_aligned=bool(np.allclose(positions[:, :, [0, 2]], positions[0:1, :, [0, 2]], atol=1e-12, rtol=0)
            and np.allclose(positions[:, :, [1, 3]], positions[:, 0:1, [1, 3]], atol=1e-12, rtol=0)),
        panel_labels=[f"({chr(97+i)})" for i in range(6)], panel_label_weight="normal",
        all_spines_visible=all(s.get_visible() for ax in fig.axes for s in ax.spines.values()),
        panel_titles_empty=all(not ax.get_title() for ax in fig.axes),
        energy_ylim_meV_fu=[-10, 370], legends_outside_axes=True, legend_frame_alpha=.85,
        all_four_Cmma_frames_retained=True, Cmma_spread_is_statistical_uncertainty=False,
        zero_denominator_fractions_rendered_as_missing=True, visual_review="pending",
        file_sha256={p.name: sha256(p) for p in sorted(output.glob(STEM + ".*"))})
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
    print(json.dumps({k: qa[k] for k in ("image_records", "grid_aligned", "new_DFT_calls")}))


if __name__ == "__main__":
    main()
