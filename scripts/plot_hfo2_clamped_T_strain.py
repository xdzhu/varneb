"""Two pinned seen training strains in unchanged reference spaces; no forecast."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from scripts import plot_hfo2_clamped_comparison as paired
from scripts.plot_hfo2_clamped_M import BASE, FRAME_IDS, build_band_data, canonical_sha, sha256

REPORT = BASE / "clamped_T_strain_reference_E089_20261011/displacement.json"
GRADIENT = REPORT.parent / "gradient.json"
OBSERVATION = BASE / "clamped_plus1_T_terminal_E087_20261011/observations/step_0005/observation.json"
AUDIT = OBSERVATION.parents[2] / "audit_receipt.json"
REPORT_SHA = "357bc8cba76a86acc8a2f9b7aa2d355bfa79ce72d7b43e3f31d4bd6820af594b"
GRADIENT_SHA = "c5346151f69ac43b5910be6f4d77275808d0e7a8b9327c5c125d77328b4f936f"
OBSERVATION_SHA = "0c5e68a923fb3562e380525b3c89ac4f6c0fa3d396d8ef9cbee49851f90db098"
AUDIT_SHA = "32a6e49a618e6e62632ae71a5e0ed093095254f6858221be2ade3d83cc194a28"


def build_data(repository):
    specs = [
        (0., paired.T_REPORT, paired.T_GRADIENT, paired.T_OBSERVATION, paired.T_AUDIT,
         paired.T_REPORT_SHA, paired.T_GRADIENT_SHA, paired.T_OBSERVATION_SHA,
         paired.T_AUDIT_SHA, "28692775", 10),
        (.01, REPORT, GRADIENT, OBSERVATION, AUDIT, REPORT_SHA, GRADIENT_SHA,
         OBSERVATION_SHA, AUDIT_SHA, "28709788", 5),
    ]
    return build_strain_data(repository, specs, channel="T")


def build_strain_data(repository, specs, *, channel):
    """A pinned zero/+1% pair; reference and native-contract guards are shared."""
    if channel not in ("T", "M") or [spec[0] for spec in specs] != [0., .01]:
        raise ValueError("registered T/M zero/+0.01 pair required before label access")
    bands, observations, references = [], [], []
    for strain, mode_path, grad_path, obs_path, audit_path, mode_pin, grad_pin, obs_pin, audit_pin, job, step in specs:
        band = build_band_data(repository, mode_path, obs_path, audit_path,
            mode_digest=mode_pin, observation_digest=obs_pin, audit_digest=audit_pin,
            job=job, step=step, expected_strain=strain)
        mode = json.loads((repository / mode_path).read_text(encoding="utf-8"))
        obs = json.loads((repository / obs_path).read_text(encoding="utf-8"))
        grad = json.loads((repository / grad_path).read_text(encoding="utf-8"))
        if (canonical_sha(grad) != grad_pin or grad["terminal_audit_sha256"] != audit_pin
                or grad["source_job_id"] != job or grad["terminal_step"] != step
                or grad["strain_training_condition"] != strain
                or not grad["all_E023_Cmma_frames_retained"]
                or grad["G3_selection_or_holdout_access"] or grad["new_DFT_calls"]
                or grad["scheduler_mutations"] or len(grad["images"]) != 9):
            raise ValueError("pinned completed training gradient required")
        for i, (row, g, raw) in enumerate(zip(band["rows"], grad["images"], obs["raw_image_evaluations"])):
            if (g["image_index"] != i or raw["image_index"] != i
                    or set(g["representations"]) != set(paired.SPANS)
                    or not np.isclose(g["relative_energy_meV_fu"], row["energy_relative_PO_meV_fu"], atol=1e-9, rtol=0)):
                raise ValueError("ordered training image/gradient identity differs")
            row.update(strain_training_condition=strain, source_job_id=job, snapshot_step=step,
                       raw_energy_eV_cell=raw["energy_eV_cell"],
                       zero_T_displacement_denominator=band["zero_displacement"][i])
            for span in paired.SPANS:
                r = g["representations"][span]
                if (r["gradient_unit"] != "eV/(sqrt(amu)*Angstrom)"
                        or r["source_basis_rank"] != r["translation_free_basis_rank"]
                        or not np.isclose(r["internal_gradient_norm"]**2,
                            r["retained_gradient_norm"]**2 + r["omitted_gradient_norm"]**2,
                            atol=1e-12, rtol=0)):
                    raise ValueError("gradient unit/rank/orthogonal norm identity differs")
                for key in ("retained_squared_gradient_fraction", "internal_gradient_norm",
                            "retained_gradient_norm", "omitted_gradient_norm",
                            "rigid_translation_gradient_norm", "zero_internal_gradient"):
                    row[f"{span}_{key}"] = r[key]
        band.update(channel=channel, strain_training_condition=strain,
            plot_label=rf"$\mathrm{{PO^+}}\rightarrow\mathrm{{{channel}}},\ \epsilon={100*strain:g}\%$",
            gradient_canonical_sha256=grad_pin)
        band["source_file_sha256"][str(grad_path)] = sha256(repository / grad_path)
        bands.append(band)
        observations.append(obs)
        references.append({k: mode[k] for k in ("original_registration_sha256",
            "author_source_sha256", "clamped_seed_manifest_sha256", "fixed_G2_orientation",
            "fixed_T_chart_integer_gauge", "T_lowest_two_complete_doublets_rank4_indices")})
    if references[0] != references[1]:
        raise ValueError("unchanged original reference/ordered-lift contract required")
    for key in ("pressure_GPa", "cell_scale_A", "k_eV_A2", "formula_units",
                "fmax_target_eV_A", "climb", "n_total_images", "n_active_images"):
        if observations[0][key] != observations[1][key]:
            raise ValueError("same physical and numerical training contract required")
    boundaries = [obs["mechanical_boundary"] for obs in observations]
    for key in ("allow_tilt", "cell_dofs", "cell_scale_A", "fixed_cell_rows",
                "kind", "raw_chain_policy"):
        if boundaries[0][key] != boundaries[1][key]:
            raise ValueError("same mechanical boundary family required")
    normal_difference = float(np.max(np.abs(np.array(boundaries[0]["unit_normal"])
                                            - boundaries[1]["unit_normal"])))
    if normal_difference > 1e-12:
        raise ValueError("same mechanical boundary normal required within serialization tolerance")
    left, right = [np.array(b["reference_cell_A"]) for b in boundaries]
    if not np.allclose(right[:2], 1.01*left[:2], atol=1e-9, rtol=0):
        raise ValueError("the prescribed +0.01 substrate increment differs")
    six = {k: v for k, v in observations[0]["raw_image_evaluations"][0]["input_sha256"].items() if k != "STRU"}
    if len(six) != 6 or any(
        {k: v for k, v in raw["input_sha256"].items() if k != "STRU"} != six
        for obs in observations for raw in obs["raw_image_evaluations"]):
        raise ValueError("original six electronic-input hashes differ")
    peaks = [b["rows"][b["highest_image_index"]] for b in bands]
    return dict(bands=bands, training_strains=[0., .01], reference_contract=references[0],
        source_contract="P0/100Ry/full10auDZP/PBE/Gamma2x2x2; fixed T substrate family; open row/tilts",
        prescribed_substrate_strain_only=True, different_strains_reuse_same_initial_SCF=False,
        unit_normal_serialization_max_difference=normal_difference,
        unit_normal_comparison_tolerance=1e-12,
        energies_relative_to_each_matched_PO=True,
        observed_sampled_forward_peak_increment_meV_fu=peaks[1]["energy_relative_PO_meV_fu"]-peaks[0]["energy_relative_PO_meV_fu"],
        observed_PO_energy_increment_eV_cell=bands[1]["rows"][0]["raw_energy_eV_cell"]-bands[0]["rows"][0]["raw_energy_eV_cell"],
        new_DFT_calls=0, scheduler_mutations=0, G3_selection_or_holdout_access=False,
        TS_certified=False, model_accuracy_certified=False, independent_prediction=False,
        native_data_reaudit_claimed=False, smooth_MEP_interpolation_used=False,
        normalized_arcs_are_individual_not_common_absolute_coordinates=True,
        comparison_script_sha256=sha256(Path(__file__)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = paired.export(build_data(args.repository), args.output, stem_name="hfo2_clamped_T_strain")
    print(json.dumps({"records": result["image_records"], "grid_aligned": result["grid_aligned"],
                      "new_DFT_calls": 0, "independent_prediction": False}))


if __name__ == "__main__":
    main()
