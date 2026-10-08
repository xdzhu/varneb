"""Describe a frozen G1 network and tolerance-dependent structure identities.

No DFT, remapping, symmetry enforcement or phase relaxation. The cases may
be incompletely optimized; a space group does not certify a stationary phase.
"""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
import warnings

import numpy as np
from ase.io import read

from examples.hfo2_fixed_input_factory import same_ordered_geometry
from scripts.analyze_hfo2_channel_network import read_evaluated_observation
from scripts.analyze_hfo2_chain_observations import analyze_observation
from scripts.analyze_hfo2_parent_patterns import rotated_t_triplet
from scripts.audit_hfo2_static_replica import sha256
from vcneb.continuous_projection import continuous_reference_coordinates


def structure_audit(atoms):
    """Validate ordered periodic Hf4O8 before a fixed, reported tolerance sweep."""
    from pymatgen.core import IStructure, Lattice
    from pymatgen.core.structure import StructureError
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    cell, frac = atoms.cell.array.copy(), atoms.get_scaled_positions(wrap=False).copy()
    if (atoms.get_chemical_symbols() != ["Hf"] * 4 + ["O"] * 8
            or not atoms.pbc.all() or not np.isfinite(cell).all()
            or not np.isfinite(frac).all() or np.linalg.det(cell) <= 1e-12):
        raise ValueError("finite ordered right-handed periodic Hf4O8 required")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            structure = IStructure(Lattice(cell, pbc=(True, True, True)), atoms.get_chemical_symbols(),
                                   frac, coords_are_cartesian=False, to_unit_cell=False,
                                   validate_proximity=True)
        except StructureError as error:
            raise ValueError("invalid periodic structure proximity") from error
        distances = structure.distance_matrix
        minimum = float(distances[np.triu_indices(12, 1)].min())
        if not structure.is_ordered or minimum < 1.6:
            raise ValueError("ordered valid study geometry with minimum distance >=1.6 Angstrom required")
        symmetry = []
        for tolerance in (.001, .01, .05):
            analyzer = SpacegroupAnalyzer(structure, symprec=tolerance, angle_tolerance=1.)
            symmetry.append({"symprec_A": tolerance, "angle_tolerance_deg": 1.,
                             "symbol": analyzer.get_space_group_symbol(),
                             "number": int(analyzer.get_space_group_number()),
                             "point_group": analyzer.get_point_group_symbol()})
    if (not np.array_equal(atoms.cell.array, cell)
            or not np.array_equal(atoms.get_scaled_positions(wrap=False), frac)):
        raise ValueError("structure audit unexpectedly changed the input geometry")
    return {"object_kind": "periodic_IStructure", "coordinate_mode": "fractional_unwrapped",
            "periodic": [True, True, True], "ordered": True, "oxidation_states_guessed": False,
            "volume_A3": float(structure.volume), "minimum_periodic_distance_A": minimum,
            "symmetry_sweep": symmetry, "warnings": [str(x.message) for x in caught],
            "standardized_wrapped_reordered_or_relaxed": False}


def analyze_case(specification: Path) -> dict:
    specification = Path(specification)
    spec = json.loads(specification.read_text())
    if (spec["formula_units"] != 4 or spec["mechanical_family"] != "free_cell_P0_E0"
            or spec["mechanical_parameters"] != {"pressure_GPa": 0., "external_field_V_A": 0.}):
        raise ValueError("fixed P=0/E=0 study ensemble required")
    root = specification.parent
    t_path, gamma_path = root / spec["T_reference"], root / spec["T_Gamma_reference"]
    t = read(t_path, format="vasp")
    parent, patterns, _ = rotated_t_triplet(t)
    reference_audit = structure_audit(t)
    expected_roles = {"PO_to_T": "decay", "PO_to_M": "decay",
                      "PO_flip_T_pattern_preserving": "switching", "PO_flip_T_pattern_reversing": "switching"}
    if (spec["required_channels"] != expected_roles or len(spec["channels"]) != 4
            or {e["name"]: e["role"] for e in spec["channels"]} != expected_roles
            or any(not isinstance(e["reverse"], bool) for e in spec["channels"])):
        raise ValueError("all four registered candidate roles and explicit source directions required")
    paths, reference_energy, common_initial = [], None, None
    with np.load(gamma_path) as gamma:
        for entry in spec["channels"]:
            folder = root / entry["observation"]
            images, observation, digest = read_evaluated_observation(folder)
            mode = analyze_observation(folder, t, parent, patterns, gamma)
            t_chart = continuous_reference_coordinates(images, t)
            u = t_chart["displacements_A"]
            centered = u - u.mean(axis=1, keepdims=True)
            rms = np.sqrt(np.sum(centered**2, axis=(1, 2)) / 12.)
            source_order = list(range(len(images)))
            arc = np.array(observation["extended_reaction_coordinate_A"], dtype=float)
            if arc[0] != 0. or np.any(np.diff(arc) <= 0):
                raise ValueError("strictly increasing source generalized arc required")
            if entry["reverse"]:
                source_order.reverse()
            initial = images[source_order[0]].get_potential_energy()
            if reference_energy is None:
                reference_energy = initial
                common_initial = images[source_order[0]]
            elif (abs(initial - reference_energy) > 1e-8
                  or not same_ordered_geometry(images[source_order[0]], common_initial)):
                raise ValueError("same ordered periodic PO+ geometry and energy required")
            rows = []
            for view_index, index in enumerate(source_order):
                image, m = images[index], mode["images"][index]
                audit = structure_audit(image)
                s = arc[index] / arc[-1]
                rows.append({
                    "view_image_index": view_index, "source_image_index": index,
                    "s_normalized": float(1. - s if entry["reverse"] else s),
                    "energy_relative_PO_meV_fu": float((image.get_potential_energy() - reference_energy) * 1000. / 4.),
                    "source_POSCAR_sha256": observation["raw_image_evaluations"][index]["snapshot_POSCAR_sha256"],
                    "source_log_sha256": observation["raw_image_evaluations"][index]["raw_log_sha256"],
                    "Q_pattern_x_A": m["parent_pattern_Q_A"][0],
                    "Q_pattern_y_A": m["parent_pattern_Q_A"][1],
                    "Q_pattern_z_A": m["parent_pattern_Q_A"][2],
                    "triplet_captured_squared_norm_fraction": m["parent_pattern_captured_squared_norm_fraction"],
                    "T_ordered_chart_translation_free_atomic_RMS_A": float(rms[index]),
                    "T_Green_strain_Frobenius_norm": float(np.linalg.norm(t_chart["green_strains"][index])),
                    "NEB_atomic_max_vector_eV_A": m.get("NEB_atomic_block_max_vector_eV_A"),
                    "NEB_cell_max_vector_eV_A": m.get("NEB_cell_block_max_vector_eV_A"),
                    "NEB_max_vector_eV_A": m.get("NEB_max_vector_eV_A"),
                    "structure_audit": audit,
                })
            peak = max(range(len(rows)), key=lambda i: rows[i]["energy_relative_PO_meV_fu"])
            paths.append({"name": entry["name"], "role": entry["role"], "rows": rows,
                          "source_job_id": observation["source_job_id"], "snapshot_step": observation["snapshot_step"],
                          "source_observation_sha256": digest, "source_direction": "reverse" if entry["reverse"] else "forward",
                          "source_NEB_fmax_eV_A": observation["replayed_fmax_eV_A"],
                          "ordinary_residual_passed": observation["replayed_fmax_eV_A"] <= .10,
                          "peak_view_image_index": peak, "peak_energy_relative_PO_meV_fu": rows[peak]["energy_relative_PO_meV_fu"],
                          "residual_dominant_source_image": mode["residual_dominant_image"],
                          "residual_dominant_block": mode["residual_dominant_block"]})
    return {"format_version": 1, "status": "G1_frozen_structure_observations_not_final_mechanisms",
            "specification_sha256": sha256(specification), "analysis_source_sha256": sha256(Path(__file__)),
            "T_reference_sha256": sha256(t_path), "T_Gamma_sha256": sha256(gamma_path),
            "packages": {name: version(name) for name in ("numpy", "ase", "pymatgen", "spglib")},
            "reference_T_structure_audit": reference_audit,
            "common_PO_energy_eV_cell": reference_energy, "channels": paths,
            "existing_image_records": sum(len(x["rows"]) for x in paths),
            "ordinary_force_target_eV_A": .10, "new_DFT_calls": 0, "physical_parameters_changed": False,
            "source_or_geometry_mutated": False, "full_G1_passed": False, "TS_certified": False,
            "limitations": ["same ordered T chart, one fixed initial gauge, no image-wise MIC or mapping",
                            "normalised arcs use original source metrics, not a common physical distance",
                            "space groups are tolerance-dependent structural observations, not relaxed phases or TS indices",
                            "T geometric patterns are not phonons, polarization or energy contributions",
                            "only the ordinary T-PO residual passed; no final competing barrier or strain prediction verdict"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--specification", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("refusing an existing material-observation report")
    result = analyze_case(args.specification)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"records": result["existing_image_records"], "new_DFT_calls": 0,
                      "peak_groups": {c["name"]: [s["symbol"] for s in c["rows"][c["peak_view_image_index"]]["structure_audit"]["symmetry_sweep"]] for c in result["channels"]}}))


if __name__ == "__main__":
    main()
