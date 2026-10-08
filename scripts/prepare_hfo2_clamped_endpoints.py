"""Prepare the preregistered common-substrate seeds, without calling DFT.

These are affine geometry starters, not relaxed endpoints. In particular,
epsilon=0 references the T substrate, not each phase's own free cell.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read

from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.audit_hfo2_static_replica import sha256
from scripts.prepare_hfo2_reference_variants import write_clean_poscar
from vcneb import clamped_plane_vcneb_boundary, endpoint_structure_record


AXES = (2, 0, 1)  # new Cartesian x/y/z = old z/x/y; det=+1
PHASE_FILES = {
    "T": "reference_variants/T.vasp",
    "PO_plus": "reference_variants/PO.vasp",
    "PO_minus_T_preserving": "reference_variants/PO_minus_T_preserving_inversion.vasp",
    "PO_minus_T_reversing": "reference_variants/PO_minus_T_reversing_inversion.vasp",
    "M": "M_CONTCAR.vasp",
}
STRAINS = (("strain_0000", 0.0), ("strain_p0100", 0.01))


def geometry_guard(images, *, minimum_distance_A=1.6):
    """Reject unsafe geometry before a calculator or optimizer mutation."""
    for image in images:
        cell = image.cell.array
        if (not np.isfinite(cell).all() or not np.isfinite(image.positions).all()
                or np.linalg.det(cell) <= 1e-12 or np.linalg.cond(cell) > 1e12):
            raise ValueError("endpoint must have finite coordinates and a positive, nonsingular cell")
        if not np.all(image.pbc):
            raise ValueError("bulk endpoint requires three periodic directions")
        distances = image.get_all_distances(mic=True)
        np.fill_diagonal(distances, np.inf)
        if len(image) > 1 and distances.min() < minimum_distance_A:
            raise ValueError("endpoint has an unsafe minimum atom distance")


def orient_long_axis_x(atoms):
    """Rotate Cartesian axes and cyclically relabel cell vectors, not atoms."""
    geometry_guard([atoms])
    rotation = np.eye(3)[list(AXES)]
    result = atoms.copy()  # ASE.copy intentionally drops the calculator
    result.set_cell(atoms.cell.array[list(AXES)] @ rotation.T, scale_atoms=False)
    result.set_positions(atoms.positions @ rotation.T)
    np.testing.assert_allclose(result.get_scaled_positions(wrap=False),
                               atoms.get_scaled_positions(wrap=False)[:, AXES], atol=1e-12)
    if result.get_chemical_symbols() != atoms.get_chemical_symbols():
        raise ValueError("orientation must preserve ordered atom identity")
    geometry_guard([result])
    return result


def _json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8", newline="\n")


def prepare(case_root: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing existing endpoint seed namespace: {output}")
    variants_path = case_root / "reference_variants/variant_manifest.json"
    variant_manifest = json.loads(variants_path.read_text(encoding="utf-8"))
    m_audit_path = case_root / "M_endpoint_audit.json"
    m_audit = json.loads(m_audit_path.read_text(encoding="utf-8"))
    if m_audit["status"] != "raw_endpoint_gates_passed":
        raise ValueError("M source must be the raw-audited ABACUS endpoint, not a literature seed")
    sources, oriented = {}, {}
    for phase, relative in PHASE_FILES.items():
        path = case_root / relative
        if phase != "M" and sha256(path) != variant_manifest["generated_structure_sha256"][path.name]:
            raise ValueError(f"registered {phase} source checksum changed")
        atoms = read(path, format="vasp")
        if atoms.get_chemical_symbols() != ["Hf"] * 4 + ["O"] * 8:
            raise ValueError("all phases must retain ordered Hf4O8")
        if phase == "M" and endpoint_structure_record(atoms)["sha256"] != m_audit["endpoint"]["sha256"]:
            raise ValueError("M geometry is not the audited ABACUS endpoint")
        oriented[phase] = orient_long_axis_x(atoms)
        sources[phase] = {"file": relative, "sha256": sha256(path),
                          "ordered_geometry": endpoint_structure_record(atoms)}
    reference = oriented["T"].cell.array.copy()
    if np.linalg.norm(reference[0]) <= max(np.linalg.norm(reference[1]), np.linalg.norm(reference[2])):
        raise ValueError("registered T long axis was not mapped to research x")
    cell_scale = float(np.linalg.det(reference) ** (1 / 3))
    manifest = {
        "format_version": 1, "status": "geometry_prepared_DFT_pending",
        "purpose": "finite_G2_common_substrate_endpoint_seeds",
        "source_audits_sha256": {"variant_manifest.json": sha256(variants_path),
                                 "M_endpoint_audit.json": sha256(m_audit_path)},
        "sources": sources, "preparation_script_sha256": sha256(Path(__file__)),
        "orientation": {"new_axes_from_old": list(AXES),
                        "cartesian_rotation_rows": np.eye(3)[list(AXES)].tolist(),
                        "determinant": 1, "cell_row_permutation": list(AXES),
                        "atom_permutation": list(range(12))},
        "unstrained_T_reference_cell_A": reference.tolist(),
        "physical_contract_sha256": CONTRACT, "pressure_gpa": 0.0,
        "allow_tilt": True, "cell_scale_A": cell_scale,
        "strain_conditions": [strain for _, strain in STRAINS],
        "holdout_strain_reserved_not_generated": 0.005,
        "new_DFT_calls": 0, "n_geometry_seeds": 10, "seeds": [],
        "limitations": [
            "seeds require fresh clamped BFGS and phase/variant audit before any NEB",
            "no free-cell energy, force, stress or Berry result is carried onto a clamped seed",
            "epsilon=0 is T-substrate clamping, not zero strain of each free phase",
            "a phase can disappear on relaxation; no symmetry/mode freezing to retain it",
            "G2 production matrix remains gated on G1, not submitted by this preparation",
        ],
    }
    output.mkdir(parents=True, exist_ok=False)
    for condition, strain in STRAINS:
        substrate = reference.copy()
        substrate[:2] *= 1 + strain
        boundary = clamped_plane_vcneb_boundary(12, substrate, allow_tilt=True)
        for phase, original in oriented.items():
            seed = original.copy()
            cell = seed.cell.array.copy()
            cell[:2] = substrate[:2]
            seed.set_cell(cell, scale_atoms=True)  # explicit affine starting geometry
            geometry_guard([seed])
            boundary.validate_images([seed])
            directory = output / condition / phase
            directory.mkdir(parents=True, exist_ok=False)
            path = directory / "POSCAR.seed"
            write_clean_poscar(path, seed)
            reread = read(path, format="vasp")
            boundary.validate_images([reread])
            distances = seed.get_all_distances(mic=True)
            np.fill_diagonal(distances, np.inf)
            entry = {
                "format_version": 1, "status": "geometry_seed_not_relaxed",
                "phase_label": phase, "strain": strain,
                "seed_file": path.name, "seed_sha256": sha256(path),
                "source_sha256": sources[phase]["sha256"],
                "ordered_species": seed.get_chemical_symbols(), "reference_cell_A": substrate.tolist(),
                "allow_tilt": True, "cell_scale_A": cell_scale, "pressure_gpa": 0.0,
                "minimum_distance_guard_A": 1.6,
                "physical_contract_sha256": CONTRACT,
                "free_phase_in_plane_length_changes": (
                    np.linalg.norm(substrate[:2], axis=1) / np.linalg.norm(original.cell[:2], axis=1) - 1
                ).tolist(),
                "minimum_distance_A": float(distances.min()),
                "volume_A3": float(seed.get_volume()),
                "ordered_geometry": endpoint_structure_record(reread),
                "new_DFT_calls": 0,
            }
            _json(directory / "endpoint_seed.json", entry)
            manifest["seeds"].append({"phase_label": phase, "strain": strain,
                                      "manifest_file": (directory / "endpoint_seed.json").relative_to(output).as_posix(),
                                      "manifest_sha256": sha256(directory / "endpoint_seed.json"),
                                      "seed_sha256": entry["seed_sha256"],
                                      "minimum_distance_A": entry["minimum_distance_A"]})
    _json(output / "clamped_seed_manifest.json", manifest)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-root", type=Path, default=Path("benchmarks/hfo2_channels/20261008"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = prepare(args.case_root, args.output)
    print(json.dumps({key: report[key] for key in
                      ("status", "n_geometry_seeds", "new_DFT_calls", "strain_conditions")}))


if __name__ == "__main__":
    main()
