"""Project the atomic component of a converged VCNEB path onto Gamma modes.

This postprocessor never launches DFT.  It requires a Gamma force-constant
archive generated at a stationary endpoint and the final complete VCNEB chain.
For a variable-cell path, fractional atomic displacements are expressed in the
reference endpoint cell before projection.  Consequently the report concerns
atomic normal coordinates only; cell degrees of freedom remain a separate
VCNEB coordinate and must not be interpreted as phonon amplitudes.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import tempfile

import numpy as np
from ase import Atoms
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb import (
    diagonalize_gamma_modes,
    load_gamma_force_constants,
    project_displacements_onto_gamma_modes,
    read_chain_trajectory,
    tangent_mode_overlaps,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, help="Stationary endpoint structure readable by ASE")
    parser.add_argument("--trajectory", required=True, help="Flat ASE VCNEB trajectory")
    parser.add_argument("--force-constants", required=True, help="Gamma force_constants.npz archive")
    parser.add_argument("--n-images", type=int, default=7, help="Total images per stored chain")
    parser.add_argument("--trajectory-step", type=int, default=-1, help="Complete chain snapshot to analyze")
    parser.add_argument(
        "--reference-permutation",
        default=None,
        help="Comma-separated reference atom indices in path order; required when NEB mapping reordered atoms",
    )
    parser.add_argument("--output", default="gamma_path_modes.json")
    parser.add_argument("--include-translations", action="store_true")
    return parser.parse_args()


def reference_cell_displacements(images: list[Atoms], reference: Atoms) -> np.ndarray:
    """Return MIC atomic displacements represented in the reference cell.

    The input path may vary its cell.  Each image is converted to fractional
    coordinates, minimum-imaged against the reference endpoint, then expressed
    through the fixed reference lattice.  This is an explicit atomic-only
    convention, appropriate for comparing with endpoint Gamma eigenvectors.
    """

    reference_symbols = reference.get_chemical_symbols()
    reference_pbc = np.asarray(reference.pbc, dtype=bool)
    if reference.cell.rank != 3 or reference.get_volume() <= 0.0:
        raise ValueError("reference must have a non-singular 3D periodic cell")
    q0 = reference.get_scaled_positions(wrap=False)
    values: list[np.ndarray] = []
    for index, image in enumerate(images):
        if image.get_chemical_symbols() != reference_symbols:
            raise ValueError(f"image {index} atom symbols/order differ from reference")
        if not np.array_equal(np.asarray(image.pbc, dtype=bool), reference_pbc):
            raise ValueError(f"image {index} PBC differs from reference")
        delta_q = image.get_scaled_positions(wrap=False) - q0
        for axis, periodic in enumerate(reference_pbc):
            if periodic:
                delta_q[:, axis] -= np.rint(delta_q[:, axis])
        values.append(delta_q @ reference.cell.array)
    return np.asarray(values, dtype=float)


def reorder_reference_force_constants(
    reference: Atoms,
    force_constants: np.ndarray,
    masses_amu: np.ndarray,
    permutation_text: str | None,
) -> tuple[Atoms, np.ndarray, np.ndarray, list[int] | None]:
    """Apply an audited endpoint mapping to both structure and force constants."""

    if permutation_text is None:
        return reference, force_constants, masses_amu, None
    try:
        permutation = [int(value.strip()) for value in permutation_text.split(",")]
    except ValueError as exc:
        raise ValueError("--reference-permutation must be comma-separated integer indices") from exc
    n_atoms = len(reference)
    if len(permutation) != n_atoms or sorted(permutation) != list(range(n_atoms)):
        raise ValueError("--reference-permutation must contain every reference atom index exactly once")
    values = np.asarray(force_constants, dtype=float)
    if values.shape == (n_atoms, 3, n_atoms, 3):
        reordered = values[np.ix_(permutation, range(3), permutation, range(3))]
    elif values.shape == (3 * n_atoms, 3 * n_atoms):
        dofs = [3 * atom + axis for atom in permutation for axis in range(3)]
        reordered = values[np.ix_(dofs, dofs)]
    else:
        raise ValueError("force_constants shape is incompatible with the reference permutation")
    return reference[permutation], reordered, np.asarray(masses_amu)[permutation], permutation


def _write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.replace(path)


def make_report(
    *,
    reference: Atoms,
    images: list[Atoms],
    force_constants: np.ndarray,
    masses_amu: np.ndarray,
    include_translations: bool,
) -> tuple[dict, dict[str, np.ndarray]]:
    if len(reference) != len(masses_amu):
        raise ValueError("force-constant archive mass count differs from reference structure")
    if not np.allclose(reference.get_masses(), masses_amu, rtol=0.0, atol=1e-8):
        raise ValueError("force-constant archive masses differ from reference structure")
    modes = diagonalize_gamma_modes(force_constants, masses_amu, project_translations=True)
    displacements = reference_cell_displacements(images, reference)
    coordinates = project_displacements_onto_gamma_modes(displacements, modes)
    overlaps = tangent_mode_overlaps(coordinates)
    translation_mask = np.abs(modes.frequencies_cm1) < 1e-7
    candidate = np.ones(len(modes.frequencies_cm1), dtype=bool) if include_translations else ~translation_mask
    peak = np.max(np.abs(coordinates), axis=0)
    ranking = [int(index) for index in np.argsort(-peak) if candidate[index]]
    report = {
        "format_version": 1,
        "interpretation": {
            "reference": "stationary endpoint Gamma force constants",
            "atomic_displacement_convention": "fractional image-reference displacement MIC-mapped into the reference cell",
            "cell_degrees_of_freedom": "excluded; analyze separately through VCNEB generalized coordinates",
            "translations_projected_before_diagonalization": True,
        },
        "n_atoms": len(reference),
        "n_images": len(images),
        "frequencies_cm1": [float(value) for value in modes.frequencies_cm1],
        "eigenvalues_eV_per_A2_amu": [float(value) for value in modes.eigenvalues_eV_per_A2_amu],
        "normal_coordinates_sqrt_amu_A": coordinates.tolist(),
        "segment_tangent_mode_overlaps": overlaps.tolist(),
        "dominant_atomic_modes": [
            {
                "mode_index": index,
                "frequency_cm1": float(modes.frequencies_cm1[index]),
                "max_abs_coordinate_sqrt_amu_A": float(peak[index]),
                "max_abs_segment_overlap": float(np.max(np.abs(overlaps[:, index]))),
            }
            for index in ranking[: min(10, len(ranking))]
        ],
    }
    arrays = {
        "reference_cell_displacements_A": displacements,
        "normal_coordinates_sqrt_amu_A": coordinates,
        "segment_tangent_mode_overlaps": overlaps,
        "frequencies_cm1": modes.frequencies_cm1,
        "eigenvectors_mass_weighted": modes.eigenvectors,
    }
    return report, arrays


def main() -> None:
    args = parse_args()
    if args.n_images < 2:
        raise ValueError("--n-images includes endpoints and must be at least 2")
    reference = read(args.reference)
    images = read_chain_trajectory(args.trajectory, n_images=args.n_images, step=args.trajectory_step)
    force_constants, masses_amu = load_gamma_force_constants(args.force_constants)
    reference, force_constants, masses_amu, permutation = reorder_reference_force_constants(
        reference,
        force_constants,
        masses_amu,
        args.reference_permutation,
    )
    report, arrays = make_report(
        reference=reference,
        images=images,
        force_constants=force_constants,
        masses_amu=masses_amu,
        include_translations=args.include_translations,
    )
    output = Path(args.output).resolve()
    report["reference_path"] = str(Path(args.reference).resolve())
    report["trajectory_path"] = str(Path(args.trajectory).resolve())
    report["force_constant_archive"] = str(Path(args.force_constants).resolve())
    report["trajectory_step"] = args.trajectory_step
    report["reference_atom_permutation_in_path_order"] = permutation
    _write_json_atomic(output, report)
    np.savez_compressed(output.with_suffix(".npz"), **arrays)
    print(f"[DONE] projected {len(images)} images onto Gamma modes; report={output}")


if __name__ == "__main__":
    main()
