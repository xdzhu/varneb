"""Reproduce a provisional four-channel PO+ comparison from frozen SCFs.

No DFT, scheduler mutation, physical parameter change or source overwrite.
The specification names all required candidates; missing coverage, unconverged
paths and absent sampling/error evidence never become a selectivity claim.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read

from examples.hfo2_fixed_input_factory import CONTRACT, same_ordered_geometry
from scripts.audit_hfo2_static_replica import sha256
from scripts.export_hfo2_chain_observation import replay
from scripts.analyze_hfo2_chain_observations import analyze_observation
from scripts.analyze_hfo2_parent_patterns import rotated_t_triplet
from vcneb import ChannelPath, summarize_competing_paths


def read_evaluated_observation(folder):
    """Verify frozen numeric E/F/stress and replay the actual source direction."""
    record = folder / "observation.json"
    report = json.loads(record.read_text())
    trajectory = folder / "evaluated_chain.traj"
    if (report["status"] != "complete_observation_not_final_result" or report["climb"]
            or report["pressure_GPa"] != 0 or report["fmax_target_eV_A"] != .10
            or report["k_eV_A2"] != .2 or report["formula_units"] != 4
            or sha256(trajectory) != report["evaluated_chain_sha256"]):
        raise ValueError("ordinary fixed-contract P=0 frozen Hf4O8 evidence required")
    images = read(trajectory, index=":")
    raw = report["raw_image_evaluations"]
    if len(images) != report["n_total_images"] or len(images) != len(raw):
        raise ValueError("incomplete image evidence")
    for i, (image, point) in enumerate(zip(images, raw)):
        poscar = folder / f"POSCAR_{i:02d}"
        if (point["image_index"] != i or not point["ordered_periodic_geometry_matched"]
                or image.get_chemical_symbols() != ["Hf"] * 4 + ["O"] * 8
                or any(point["input_sha256"][name] != value for name, value in CONTRACT.items())
                or sha256(poscar) != point["snapshot_POSCAR_sha256"]
                or not same_ordered_geometry(image, read(poscar, format="vasp"))
                or not np.isclose(image.get_potential_energy(), point["energy_eV_cell"], rtol=0, atol=1e-10)
                or not np.allclose(image.get_forces(), point["forces_eV_A"], rtol=0, atol=1e-12)
                or not np.allclose(image.get_stress(), point["stress_ASE_voigt_eV_A3"], rtol=0, atol=1e-12)):
            raise ValueError("frozen numeric/geometry/input evidence mismatch")
    chain, forces, _, _ = replay(images)
    residual = float(np.linalg.norm(forces, axis=1).max())
    if (abs(residual - report["replayed_fmax_eV_A"]) > 1e-10
            or abs(residual - report["optimizer_log_row"]["fmax_eV_A"]) > 7e-7
            or abs(float(chain.enthalpies.max()) - report["optimizer_log_row"]["highest_energy_eV_cell"]) > 7e-7):
        raise ValueError("ordinary force/energy replay differs from the production record")
    return images, report, sha256(record)


def analyze(specification, output):
    if output.exists():
        raise FileExistsError("refusing existing network analysis")
    spec = json.loads(specification.read_text())
    root = specification.parent
    if (spec["formula_units"] != 4 or spec["mechanical_family"] != "free_cell_P0_E0"
            or spec["mechanical_parameters"] != {"pressure_GPa": 0., "external_field_V_A": 0.}):
        raise ValueError("this case adapter is limited to the original P=0/E=0 contract")
    t_path = (root / spec["T_reference"]).resolve()
    gamma_path = (root / spec["T_Gamma_reference"]).resolve()
    t = read(t_path)
    parent, patterns, _ = rotated_t_triplet(t)
    paths, sources, mode_records = [], [], []
    with np.load(gamma_path) as gamma:
        for entry in spec["channels"]:
            folder = (root / entry["observation"]).resolve()
            images, record, digest = read_evaluated_observation(folder)
            paths.append(ChannelPath(
                name=entry["name"], role=entry["role"], images=images,
                energies_eV_cell=[p["energy_eV_cell"] for p in record["raw_image_evaluations"]],
                physical_contract=CONTRACT, mechanical_family=spec["mechanical_family"],
                mechanical_parameters=spec["mechanical_parameters"], pressure_eV_A3=0.,
                neb_fmax_eV_A=record["replayed_fmax_eV_A"], source_id=digest,
                reverse=entry["reverse"],
            ))
            sources.append({"name": entry["name"], "relative_observation": entry["observation"],
                            "observation_sha256": digest, "source_job_id": record["source_job_id"],
                            "snapshot_step": record["snapshot_step"],
                            "optimizer_time_CST": record["optimizer_log_row"]["time_CST"]})
            mode = analyze_observation(folder, t, parent, patterns, gamma)
            mode["channel_name"] = entry["name"]
            mode["source_direction_used"] = "reverse" if entry["reverse"] else "forward"
            # Keep original image indices/lift/optimization metric; reversing
            # the thermodynamic view must not silently recompute the NEB metric.
            mode_records.append(mode)
    result = summarize_competing_paths(paths, formula_units=4, required_channels=spec["required_channels"])
    result.update({
        "status": "provisional_unconverged_four_channel_observations",
        "sources": sources, "reference_mode_observations_in_source_direction": mode_records,
        "new_DFT_calls": 0, "physical_parameters_changed": False,
        "existing_image_evaluations": sum(len(p.images) for p in paths),
        "specification_sha256": sha256(specification), "analysis_script_sha256": sha256(Path(__file__)),
        "channel_module_sha256": sha256(Path(__file__).resolve().parents[1] / "vcneb/channel_competition.py"),
        "T_reference_sha256": sha256(t_path), "T_Gamma_reference_sha256": sha256(gamma_path),
        "H1_selectivity_conclusion": "not_evaluated_without_convergence_sampling_and_measured_error",
        "distinct_switching_MEPs_certified": False,
    })
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--specification", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = analyze(args.specification, args.output)
    print(json.dumps({"status": report["status"], "cached_images": report["existing_image_evaluations"],
                      "new_DFT_calls": report["new_DFT_calls"], "ready": report["ready_for_bounded_discrete_comparison"],
                      "channels": [{"name": r["name"], "source_fmax_eV_A": r["source_NEB_fmax_eV_A"],
                                    "provisional_discrete_barrier_meV_fu": r["barrier_from_common_initial_eV_fu"] * 1000}
                                   for r in report["channels"]]}, indent=2))


if __name__ == "__main__":
    main()
