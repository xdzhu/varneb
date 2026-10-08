"""Audit frozen SCFs and describe force/mode evolution, without new DFT.

Use a single initial-image periodic gauge throughout each chain. Parent
patterns and stationary-T Gamma modes are descriptive reference bases, not
local phonons at every image, energy partitions or TS certificates.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read

from scripts.audit_hfo2_static_replica import sha256
from scripts.analyze_hfo2_parent_patterns import rotated_t_triplet
from scripts.export_hfo2_chain_observation import replay
from vcneb.continuous_projection import continuous_reference_coordinates, project_reference_basis
from vcneb.phonons import identify_acoustic_modes


def analyze_observation(folder, tetragonal, parent, patterns, gamma):
    source = folder / "observation.json"
    report = json.loads(source.read_text())
    if (report["status"] != "complete_observation_not_final_result"
            or report["new_DFT_calls"] != 0 or report["pressure_GPa"] != 0
            or report["climb"] or report["formula_units"] != 4):
        raise ValueError("fixed ordinary P=0 Hf4O8 observation required")
    trajectory = folder / "evaluated_chain.traj"
    if sha256(trajectory) != report["evaluated_chain_sha256"]:
        raise ValueError("cached evaluated trajectory hash changed")
    images = read(trajectory, index=":")
    raw = report["raw_image_evaluations"]
    if len(images) != report["n_total_images"] or len(raw) != len(images):
        raise ValueError("incomplete image evidence")
    for i, (a, evaluation) in enumerate(zip(images, raw)):
        if (evaluation["image_index"] != i
                or sha256(folder / f"POSCAR_{i:02d}") != evaluation["snapshot_POSCAR_sha256"]
                or a.get_chemical_symbols() != ["Hf"] * 4 + ["O"] * 8
                or not np.isclose(a.get_potential_energy(), evaluation["energy_eV_cell"], atol=1e-10, rtol=0)
                or not np.allclose(a.get_forces(), evaluation["forces_eV_A"], atol=1e-12, rtol=0)
                or not np.allclose(a.get_stress(), evaluation["stress_ASE_voigt_eV_A3"], atol=1e-12, rtol=0)):
            raise ValueError("frozen numeric E/F/stress/geometry evidence mismatch")
    chain, forces, diagnostics, _ = replay(images)
    fmax = float(np.linalg.norm(forces, axis=1).max())
    if not np.isclose(fmax, report["replayed_fmax_eV_A"], atol=1e-10, rtol=0):
        raise ValueError("local force replay disagrees with HF observation")
    chart = continuous_reference_coordinates(images, parent)
    pattern_projection = project_reference_basis(chart["displacements_A"], patterns.reshape(3, -1).T)
    t_chart = continuous_reference_coordinates(images, tetragonal)
    masses = gamma["masses_amu"]
    # Preserve the original Phonopy masses (O=15.9994), not ASE's default
    # O=15.999. Neither the DFT calculation nor the eigenpairs are retuned.
    if not np.allclose(masses, [178.49] * 4 + [15.9994] * 8, rtol=0, atol=1e-8):
        raise ValueError("T Gamma mass ordering differs from its recorded Phonopy reference")
    eigenvectors = gamma["eigenvectors_mass_weighted"]
    acoustic, translation_overlap = identify_acoustic_modes(eigenvectors, masses)
    optical = np.setdiff1d(np.arange(len(eigenvectors)), acoustic)
    gamma_projection = project_reference_basis(t_chart["displacements_A"], eigenvectors, metric_weights=masses)
    blocks = forces.reshape(len(images) - 2, 15, 3)
    records = []
    for i, a in enumerate(images):
        record = {
            "image_index": i,
            "extended_reaction_coordinate_A": float(chain.reaction_coordinate()[i]),
            "relative_energy_meV_fu": float((chain.enthalpies[i] - chain.enthalpies[0]) * 250),
            "parent_pattern_Q_A": pattern_projection["amplitudes"][i].tolist(),
            "parent_pattern_residual_A": float(pattern_projection["residual_norm"][i]),
            "parent_pattern_captured_squared_norm_fraction": float(pattern_projection["captured_squared_norm_fraction"][i]),
            "parent_chart_max_fractional_offset": float(np.max(np.abs(chart["fractional_offsets"][i]))),
            "inside_local_parent_chart_0p45": bool(np.max(np.abs(chart["fractional_offsets"][i])) <= .45),
            "deformation_from_original_T": chart["deformation_gradients"][i].tolist(),
            "green_strain_from_original_T": chart["green_strains"][i].tolist(),
            "T_Gamma_Q_sqrt_amu_A": gamma_projection["amplitudes"][i].tolist(),
            "T_Gamma_optical_total_squared_norm_amu_A2": float(np.sum(gamma_projection["amplitudes"][i, optical]**2)),
            "T_Gamma_full_basis_reconstruction_residual_sqrt_amu_A": float(gamma_projection["residual_norm"][i]),
        }
        if 0 < i < len(images) - 1:
            norms = np.linalg.norm(blocks[i - 1], axis=1)
            vector = int(np.argmax(norms))
            d = diagnostics["images"][i]
            record.update({
                "NEB_atomic_block_max_vector_eV_A": float(norms[:12].max()),
                "NEB_cell_block_max_vector_eV_A": float(norms[12:].max()),
                "NEB_max_vector_eV_A": float(norms.max()),
                "NEB_dominant_block": "atomic" if vector < 12 else "cell",
                "NEB_dominant_atom_or_cell_row": vector if vector < 12 else vector - 12,
                "true_perpendicular_max_vector_eV_A": d["true_perpendicular_force_max_vector_eV_per_A"],
                "spring_max_vector_eV_A": d["spring_force_max_vector_eV_per_A"],
                "true_tangential_euclidean_eV_A": d["true_tangential_force_eV_per_A"],
            })
        records.append(record)
    dominant = max(records[1:-1], key=lambda r: r["NEB_max_vector_eV_A"])
    peak = int(np.argmax(chain.enthalpies))
    forward = float((chain.enthalpies[peak] - chain.enthalpies[0]) * 250)
    reverse = float((chain.enthalpies[peak] - chain.enthalpies[-1]) * 250)
    return {
        "observation": folder.name, "source_observation_sha256": sha256(source),
        "source_job_id": report["source_job_id"], "snapshot_step": report["snapshot_step"],
        "status": "unconverged_observation" if fmax > .10 else "NEB_residual_passed_TS_and_sampling_gates_pending",
        "replayed_fmax_eV_A": fmax, "highest_image_index": peak,
        "provisional_discrete_forward_barrier_meV_fu": forward,
        "provisional_discrete_reverse_barrier_meV_fu": reverse,
        "reaction_energy_meV_fu": forward - reverse,
        "residual_dominant_image": dominant["image_index"],
        "residual_dominant_block": dominant["NEB_dominant_block"],
        "max_parent_pattern_residual_A": max(r["parent_pattern_residual_A"] for r in records),
        "min_parent_pattern_captured_squared_norm_fraction": min(r["parent_pattern_captured_squared_norm_fraction"] for r in records),
        "all_images_inside_local_parent_chart": all(r["inside_local_parent_chart_0p45"] for r in records),
        "parent_chart_fixed_integer_shifts": chart["fixed_integer_lattice_shifts_by_atom"].tolist(),
        "T_chart_fixed_integer_shifts": t_chart["fixed_integer_lattice_shifts_by_atom"].tolist(),
        "Gamma_acoustic_mode_indices": acoustic.tolist(), "Gamma_acoustic_min_translation_overlap": translation_overlap,
        "images": records,
    }


def analyze(root, variants, gamma_path, output):
    if output.exists():
        raise FileExistsError("refusing existing observation analysis")
    t_path = variants / "T.vasp"
    tetragonal = read(t_path)
    parent, patterns, _ = rotated_t_triplet(tetragonal)
    with np.load(gamma_path) as gamma:
        # pathlib orders case-insensitively on Windows, case-sensitively on
        # Linux. Use an explicit common order for cross-platform evidence.
        folders = sorted((p.parent for p in root.glob("*/observation.json")), key=lambda p: p.name.casefold())
        if not folders:
            raise ValueError("no complete exported observations")
        results = [analyze_observation(folder, tetragonal, parent, patterns, gamma) for folder in folders]
        frequencies = gamma["frequencies_thz"].tolist()
    result = {
        "schema_version": 1, "new_DFT_calls": 0, "physical_parameters_changed": False,
        "analysis_script_sha256": sha256(Path(__file__)), "T_reference_sha256": sha256(t_path),
        "continuous_projection_module_sha256": sha256(Path(__file__).resolve().parents[1] / "vcneb/continuous_projection.py"),
        "T_Gamma_source_sha256": sha256(gamma_path), "T_Gamma_frequencies_THz": frequencies,
        "Gamma_metric_masses_amu": [178.49] * 4 + [15.9994] * 8,
        "Gamma_mass_convention": "original Phonopy masses; not replaced by ASE default atomic weights",
        "pattern_order": ["rotated_T_pattern_x", "rotated_T_pattern_y", "T_pattern_z"],
        "observations": results,
        "limitations": ["descriptive continuous-reference projections, not an energy decomposition",
                        "T Gamma eigenvectors are reference modes, not local phonons or TS unstable modes",
                        "individual degenerate Gamma coordinates are basis-dependent; use complete-subspace weights",
                        "Green strain is separate from atomic pattern amplitude; no silent cell rotation",
                        "provisional discrete peaks need path convergence, sampling and full-variable TS audits",
                        "fixed initial atom gauge preserves winding; no image-wise MIC folding"],
    }
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "variants", "gamma", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.root, args.variants, args.gamma, args.output)
    keys = ("observation", "replayed_fmax_eV_A", "residual_dominant_image", "residual_dominant_block",
            "provisional_discrete_forward_barrier_meV_fu", "max_parent_pattern_residual_A",
            "min_parent_pattern_captured_squared_norm_fraction")
    print(json.dumps([{key: r[key] for key in keys} for r in result["observations"]], indent=2))


if __name__ == "__main__":
    main()
