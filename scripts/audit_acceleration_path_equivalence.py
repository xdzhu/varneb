"""Compare first-threshold VCNEB chains on a common joint arc coordinate.

This is a diagnostic, not a retrospectively chosen path-equivalence pass gate.
The input trajectories must be complete *evaluated* chains from the same
calculator/pressure contract and have identical ordered endpoints.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read
from ase.units import GPa

from vcneb.analysis import path_reaction_coordinate
from vcneb.core import deformation_from_cell


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_chain(path: Path) -> list[Atoms]:
    images = list(read(str(path), index=":"))
    if len(images) < 3:
        raise ValueError(f"{path}: incomplete chain")
    symbols = images[0].get_chemical_symbols()
    for index, image in enumerate(images):
        if image.get_chemical_symbols() != symbols:
            raise ValueError(f"{path}: atom order changes at image {index}")
        if not np.all(np.asarray(image.pbc, dtype=bool) == np.asarray(images[0].pbc, dtype=bool)):
            raise ValueError(f"{path}: periodicity changes at image {index}")
        if not np.isfinite(image.get_potential_energy()) or image.get_volume() <= 0:
            raise ValueError(f"{path}: non-finite energy or invalid cell at image {index}")
    return images


def joint_vectors(images: list[Atoms], reference_cell: np.ndarray, cell_scale: float) -> np.ndarray:
    """Use the same unwrapped fractional/strain metric as path_reaction_coordinate."""

    previous_q = images[0].get_scaled_positions(wrap=False).copy()
    vectors = []
    for index, image in enumerate(images):
        q = image.get_scaled_positions(wrap=False).copy()
        if index:
            delta = q - previous_q
            for axis, periodic in enumerate(np.asarray(images[0].pbc, dtype=bool)):
                if periodic:
                    delta[:, axis] -= np.rint(delta[:, axis])
            q = previous_q + delta
        previous_q = q
        deform = deformation_from_cell(image.cell.array, reference_cell)
        vectors.append(np.r_[(q @ reference_cell).reshape(-1),
                             (cell_scale * (deform - np.eye(3))).reshape(-1)])
    return np.asarray(vectors)


def audit_pair(
    baseline: list[Atoms], variant: list[Atoms], *, n_formula: int, pressure_gpa: float,
) -> dict:
    if len(baseline) != len(variant):
        raise ValueError("chain image counts differ")
    if n_formula < 1 or not np.isfinite(pressure_gpa):
        raise ValueError("invalid formula count or pressure")
    if any(a.get_chemical_symbols() != b.get_chemical_symbols()
           for a, b in zip(baseline, variant)):
        raise ValueError("atom mapping differs between chains")
    reference_cell = baseline[0].cell.array
    cell_scale = float(np.linalg.det(reference_cell) ** (1.0 / 3.0))
    base_x = joint_vectors(baseline, reference_cell, cell_scale)
    var_x = joint_vectors(variant, reference_cell, cell_scale)
    endpoint_errors = np.linalg.norm(base_x[[0, -1]] - var_x[[0, -1]], axis=1)
    if np.max(endpoint_errors) > 1e-8:
        raise ValueError(f"ordered endpoints differ in joint metric: {endpoint_errors}")

    base_s, _ = path_reaction_coordinate(baseline, cell_scale=cell_scale)
    var_s, _ = path_reaction_coordinate(variant, cell_scale=cell_scale)
    for label, s, x in (("baseline", base_s, base_x), ("variant", var_s, var_x)):
        calculated = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(x, axis=0), axis=1))]
        np.testing.assert_allclose(calculated, s, rtol=0, atol=1e-8,
                                   err_msg=f"{label} metric disagrees with vcneb.analysis")
        if not np.all(np.diff(s) > 0):
            raise ValueError(f"{label} has a collapsed/reversed arc segment")
    base_t = base_s / base_s[-1]
    var_t = var_s / var_s[-1]

    def enthalpies(images: list[Atoms]) -> np.ndarray:
        return np.asarray([image.get_potential_energy() + pressure_gpa * GPa * image.get_volume()
                           for image in images], dtype=float)

    base_h = enthalpies(baseline)
    var_h = enthalpies(variant)
    if abs(base_h[0] - var_h[0]) > 1e-8 or abs(base_h[-1] - var_h[-1]) > 1e-8:
        raise ValueError("endpoint enthalpy references differ")
    base_h -= base_h[0]
    var_h -= var_h[0]

    t = np.unique(np.r_[np.linspace(0.0, 1.0, 401), base_t, var_t])
    base_h_fit = np.interp(t, base_t, base_h)
    var_h_fit = np.interp(t, var_t, var_h)
    energy_delta = (var_h_fit - base_h_fit) * (1000.0 / n_formula)
    base_x_fit = np.column_stack([np.interp(t, base_t, base_x[:, k])
                                  for k in range(base_x.shape[1])])
    var_x_fit = np.column_stack([np.interp(t, var_t, var_x[:, k])
                                 for k in range(var_x.shape[1])])
    geometry_delta = np.linalg.norm(var_x_fit - base_x_fit, axis=1)
    direct_image_delta = np.linalg.norm(var_x - base_x, axis=1)
    return {
        "n_total_images": len(baseline),
        "n_formula_units": n_formula,
        "pressure_GPa": pressure_gpa,
        "reference_cell_scale_A": cell_scale,
        "baseline_total_joint_arc_A": float(base_s[-1]),
        "variant_total_joint_arc_A": float(var_s[-1]),
        "baseline_max_relative_enthalpy_eV_per_formula": float(np.max(base_h) / n_formula),
        "variant_max_relative_enthalpy_eV_per_formula": float(np.max(var_h) / n_formula),
        "maximum_relative_enthalpy_difference_meV_per_formula": float(
            (np.max(var_h) - np.max(base_h)) * 1000.0 / n_formula
        ),
        "max_abs_common_arc_energy_difference_meV_per_formula": float(np.max(np.abs(energy_delta))),
        "rms_common_arc_energy_difference_meV_per_formula": float(np.sqrt(np.mean(energy_delta**2))),
        "max_common_arc_joint_geometry_difference_A": float(np.max(geometry_delta)),
        "rms_common_arc_joint_geometry_difference_A": float(np.sqrt(np.mean(geometry_delta**2))),
        "max_same_index_joint_geometry_difference_A": float(np.max(direct_image_delta)),
        "limitation": "Piecewise-linear common normalized arc comparison; not a predeclared equivalence pass gate or proof of the same saddle basin.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--variant", nargs=2, action="append", metavar=("LABEL", "TRAJ"), required=True)
    parser.add_argument("--n-formula", type=int, required=True)
    parser.add_argument("--pressure-gpa", type=float, default=0.0)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true",
                        help="replace only the explicitly selected generated report")
    args = parser.parse_args()
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(args.output)
    baseline = load_chain(args.baseline)
    variants = {}
    hashes = {str(args.baseline): sha256(args.baseline)}
    for label, name in args.variant:
        if label in variants:
            raise ValueError(f"duplicate variant label: {label}")
        path = Path(name)
        variants[label] = audit_pair(baseline, load_chain(path),
                                     n_formula=args.n_formula,
                                     pressure_gpa=args.pressure_gpa)
        hashes[str(path)] = sha256(path)
    report = {
        "status": "first_threshold_common_arc_diagnostic_no_predeclared_equivalence_gate",
        "baseline": str(args.baseline), "variants": variants,
        "source_sha256": hashes,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({label: {
        "max_energy_difference_meV_per_formula": result["max_abs_common_arc_energy_difference_meV_per_formula"],
        "max_joint_geometry_difference_A": result["max_common_arc_joint_geometry_difference_A"],
    } for label, result in variants.items()}))


if __name__ == "__main__":
    main()
