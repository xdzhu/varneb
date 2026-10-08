"""Inert audit of published Cmma cells and QE Gamma displacement conventions.

Raw author files remain outside Git. This does not map onto a production path,
import LDA force constants into PBE, launch DFT, or claim a prediction baseline
has passed. Only reviewed hashes, declared frames and fixed species mappings
are accepted; no automatic symmetry repair or image-dependent permutation.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

from ase.io import read
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.ase import AseAtomsAdaptor
from scipy.optimize import linear_sum_assignment

from scripts.prepare_hfo2_reference_variants import symmetry_report
from vcneb.qe_modes import read_qe_gamma_modes
from vcneb.reference_variants import apply_parent_operation


AUTHOR_COMMIT = "a438e4ecf63cddfa68d4ae8ea83d4c19bd238969"
SOURCES = {
    "bto.in": ("figure2/bto.in", "8028cd106e5d9eb1c9a27b412b9b1a54fa819701"),
    "bto.out": ("figure2/bto.out", "5a10ac810a65ad83879c47ccc3a2688b7f2d413b"),
    "matdyn.modes": ("figure2/matdyn.modes", "b2c3dddeb3ddc684867b332ee616bc9be9352212"),
    "README": ("figure2/README", "176d734b8fc1b0e311c2cdd43e999b05d40445f9"),
    "Cmma.vasp": ("figure1/Cmma.vasp", "9b337a6b6a2f23785973da76f2585cdefec7d661"),
}


def assign_declared_frame(source, target, basis_change, rotation, shift, *, tolerance_A):
    """One exact-frame/species assignment, rejecting ambiguous near neighbors.

    ASE row convention: H'=B H R.T and f'=f inv(B)+shift. B and R must be
    proper; no uniform scale fitting or reference/path optimization is done.
    A species-preserving permutation only records an equivalent representation.
    """
    basis, rotation = np.asarray(basis_change, float), np.asarray(rotation, float)
    shift = np.asarray(shift, float)
    if (basis.shape != (3, 3) or rotation.shape != (3, 3) or shift.shape != (3,)
            or not np.isfinite(np.r_[basis.ravel(), rotation.ravel(), shift]).all()
            or not np.allclose(basis, np.rint(basis), rtol=0, atol=1e-12)
            or not np.isclose(np.linalg.det(basis), 1, rtol=0, atol=1e-12)
            or not np.allclose(rotation @ rotation.T, np.eye(3), rtol=0, atol=1e-12)
            or not np.isclose(np.linalg.det(rotation), 1, rtol=0, atol=1e-12)
            or not np.isfinite(tolerance_A) or tolerance_A <= 0):
        raise ValueError("proper rotation, integer unimodular basis and finite shift/tolerance required")
    if (len(source) != len(target) or not all(source.pbc) or not all(target.pbc)
            or sorted(source.get_chemical_symbols()) != sorted(target.get_chemical_symbols())):
        raise ValueError("same periodic composition and atom count required")
    for atoms in (source, target):
        if (not len(atoms) or len(atoms) > 128 or not np.isfinite(atoms.positions).all()
                or not np.isfinite(atoms.cell.array).all() or np.linalg.det(atoms.cell.array) <= 0):
            raise ValueError("finite right-handed periodic structure required, at most 128 atoms")
    cell = basis @ source.cell.array @ rotation.T
    cell_error = float(np.max(np.abs(cell - target.cell.array)))
    if cell_error > tolerance_A:
        raise ValueError("declared cells disagree; scaling/strain cannot be silently fitted")
    sites = source.get_scaled_positions(wrap=False) @ np.linalg.inv(basis) + shift
    delta = sites[:, None, :] - target.get_scaled_positions(wrap=False)[None, :, :]
    integers = np.rint(delta)
    distance = np.linalg.norm((delta - integers) @ target.cell.array, axis=2)
    distance[source.numbers[:, None] != target.numbers[None, :]] = np.inf
    row, permutation = linear_sum_assignment(distance)
    matched = distance[row, permutation]
    margins = np.sort(distance, axis=1)[:, 1] - np.sort(distance, axis=1)[:, 0]
    if np.max(matched) > tolerance_A or np.min(margins) <= tolerance_A:
        raise ValueError("declared-frame species assignment fails or is ambiguous")
    return {
        "basis_change_integer": basis.astype(int).tolist(),
        "rotation_Cartesian_proper": rotation.tolist(),
        "translation_fractional": shift.tolist(),
        "source_to_target_indices": permutation.tolist(),
        "integer_shifts_source_minus_target": integers[row, permutation].astype(int).tolist(),
        "maximum_lattice_component_error_A": cell_error,
        "maximum_site_error_A": float(np.max(matched)),
        "minimum_assignment_margin_A": float(np.min(margins)),
        "tolerance_A": tolerance_A,
    }


def reconstruct_table_S2():
    """Literal printed primitive data -> existing 12-atom conventional cell.

    This doubles a 6-atom representation of the same published reference, not
    our production cell or a phonon calculation. Rounding residuals are kept.
    """
    primitive = Structure(
        Lattice([[3.70983, 0, 0], [.12179, 3.70783, 0], [0, 0, 4.84891]]),
        ["Hf", "Hf", "O", "O", "O", "O"],
        [[.25, .75, .78785], [.75, .25, .21215], [0, 0, .5], [.5, .5, .5],
         [.75, .75, 0], [.25, .25, 0]],
        coords_are_cartesian=False, validate_proximity=True,
    )
    before = primitive.as_dict()
    matrix = np.array([[1, -1, 0], [1, 1, 0], [0, 0, 1]])
    conventional = primitive.copy().make_supercell(matrix)
    assert primitive.as_dict() == before
    cell = conventional.lattice.matrix
    x = cell[0] / np.linalg.norm(cell[0])
    z = np.cross(cell[0], cell[1])
    z /= np.linalg.norm(z)
    rotation = np.stack([x, np.cross(z, x), z])
    atoms = AseAtomsAdaptor.get_atoms(conventional)
    atoms.set_cell(cell @ rotation.T, scale_atoms=True)
    return atoms, {
        "source": "arXiv2412.16792v2 ancillary SI p.3 TableS2, literal rounded numbers",
        "primitive_atoms": 6, "conventional_atoms": 12,
        "primitive_to_conventional_integer_matrix": matrix.tolist(),
        "determinant": int(round(np.linalg.det(matrix))),
        "Cartesian_rotation": rotation.tolist(),
        "primitive_input_unchanged": primitive.as_dict() == before,
    }


def audit(source_root: Path):
    sources = {}
    for name, (relative, expected_blob) in SOURCES.items():
        path = source_root / name
        if path.stat().st_size > 12 * 1024**2:
            raise ValueError("author input exceeds reviewed byte bound")
        data = path.read_bytes()
        blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        if blob != expected_blob:
            raise ValueError(f"published {name} Git blob changed")
        sources[name] = {"URL": f"https://raw.githubusercontent.com/yuboqiuab/unstableflatband/{AUTHOR_COMMIT}/{relative}",
                         "bytes": len(data), "git_blob": blob, "sha256": hashlib.sha256(data).hexdigest()}
    qe = read(source_root / "bto.in", format="espresso-in")
    canonical = read(source_root / "Cmma.vasp", format="vasp")
    geometry = {"QE": symmetry_report(qe), "published_POSCAR": symmetry_report(canonical)}
    # README declares the improper y/z exchange. Compose it with canonical
    # Cmma inversion (-I, translation [.5,0,0]) to obtain a proper frame.
    # The translation and permutation are checked, not guessed from labels.
    frame = np.array([[-1, 0, 0], [0, 0, -1], [0, -1, 0]])
    qe_mapping = assign_declared_frame(qe, canonical, frame, frame, [.5, 0, 0], tolerance_A=3e-6)
    reconstructed, table_transform = reconstruct_table_S2()
    geometry["SI_primitive_reconstruction"] = symmetry_report(reconstructed)
    table_mapping = assign_declared_frame(reconstructed, canonical, np.eye(3), np.eye(3),
                                         [.5, .75, 0], tolerance_A=5e-5)
    # ASE espresso-in uses default element masses. Preserve the actual author
    # ATOMIC_SPECIES/QE-output masses instead, verified in the pinned inputs.
    masses = np.array([178.49] * 4 + [15.9994] * 8)
    modes = read_qe_gamma_modes(source_root / "matdyn.modes", masses)
    translations = np.zeros((36, 3))
    for direction in range(3):
        translations[direction::3, direction] = np.sqrt(masses)
    translations /= np.linalg.norm(translations, axis=0)
    scores = np.sum(np.abs(translations.T @ modes.mass_weighted_eigenvectors)**2, axis=0)
    acoustic = np.sort(np.argsort(scores)[-3:])
    acoustic_overlap = float(np.linalg.svd(translations.T @ modes.mass_weighted_eigenvectors[:, acoustic],
                                         compute_uv=False).min())
    if acoustic_overlap < .99:
        raise ValueError("published ordered Gamma rows do not resolve rigid translations")
    centering = apply_parent_operation(qe, qe, np.eye(3), [.5, 0, .5]).source_to_target
    vectors = modes.mass_weighted_eigenvectors.reshape(12, 3, 36)
    translated = vectors[np.argsort(centering)].reshape(36, 36)
    parity = []
    for i in range(4):
        vector = modes.mass_weighted_eigenvectors[:, i]
        even = float(np.linalg.norm(.5 * (translated[:, i] + vector))**2)
        odd = float(np.linalg.norm(.5 * (translated[:, i] - vector))**2)
        parity.append({"source_mode_one_based": i+1, "frequency_THz": float(modes.frequencies_THz[i]),
                       "centering_even_squared_norm": even, "centering_odd_squared_norm": odd})
    return {
        "schema_version": 1, "status": "published_reference_inert_mapping_not_material_prediction",
        "author_commit": AUTHOR_COMMIT, "sources": sources,
        "software": {name: importlib.metadata.version(name) for name in ("pymatgen", "ase", "numpy", "scipy", "spglib")},
        "geometry": geometry, "QE_to_published_POSCAR": qe_mapping,
        "SI_primitive_reconstruction": {**table_transform, "assignment": table_mapping},
        "printed_TableS1_issue": "literal x=0.05000 differs from both pinned author structure and independently reconstructed TableS2; printed table left unchanged",
        "source_QE_coordinate_axes": "source y/z correspond to paper z/y; explicit proper frame/permutation recorded",
        "mode_convention_source": "https://www.quantum-espresso.org/Doc/INPUT_MATDYN.html#flvec",
        "Gamma_modes": {
            "count": 36, "masses_amu": masses.tolist(), "ASE_default_masses_not_used": qe.get_masses().tolist(),
            "printed_normalization_max_defect": modes.printed_normalization_max_defect,
            "mass_weighted_Gram_max_defect": modes.mass_weighted_gram_max_defect,
            "acoustic_indices_zero_based": acoustic.tolist(), "minimum_translation_overlap": acoustic_overlap,
            "four_resolved_negative_reference_modes": parity,
            "centering_translation_in_source_fractional_axes": [.5, 0, .5],
            "centering_source_to_target_indices": list(centering),
            "all_signed_frequencies_THz": modes.frequencies_THz.tolist(),
            "later_q_blocks_not_imported": True, "complex_vectors_preserved_no_QR_or_ASR": True,
        },
        "new_DFT_calls": 0, "physical_parameters_changed": False,
        "limitations": ["author LDA reference eigenvectors are geometrical candidates, not our PBE curvatures",
                        "co-located author files define the association; underlying flfrc/DFPT generation was not audited",
                        "no irrep label inferred from mode number or centering character",
                        "production ordered-atom/path gauge not yet mapped; strong baseline not yet benchmarked",
                        "no whole-spectrum or TS stability certification, prediction, joint elimination or acceleration claim",
                        "no explicit redistribution license seen in repository root; raw data/eigenvectors not redistributed"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("refusing an existing Cmma audit")
    result = audit(args.source_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"QE_mapping": result["QE_to_published_POSCAR"],
                      "SI_mapping": result["SI_primitive_reconstruction"]["assignment"],
                      "Gamma": result["Gamma_modes"]}))


if __name__ == "__main__":
    main()
