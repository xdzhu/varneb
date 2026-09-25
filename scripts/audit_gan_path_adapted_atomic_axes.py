"""Audited geometric low-rank analysis of the converged GaN VCNEB path.

The returned atomic axes are empirical path singular vectors, *not* Gamma
phonons or transition-state unstable modes. Cell strain is reported separately.
No calculator or DFT job is launched.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from examples.analyze_path_gamma_modes import reference_cell_displacements
from scripts.analyze_gan_qian_path import latest_evaluated_chain
from scripts.audit_gan_ts_candidate import audit as audit_ts_candidate
from vcneb.analysis import path_reaction_coordinate
from vcneb.path_modes import fit_path_atomic_svd, leave_one_image_out_residuals
from vcneb.provenance import validate_static_endpoint_identity


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _cell_report(images) -> tuple[list[list[float]], float, float]:
    reference_cell = np.asarray(images[0].cell.array, dtype=float)
    strains = []
    max_asymmetry = 0.0
    max_cell_reconstruction_error = 0.0
    for image in images:
        cell = np.asarray(image.cell.array, dtype=float)
        deformation = np.linalg.solve(reference_cell, cell)
        asymmetry = float(np.max(np.abs(deformation - deformation.T)))
        symmetric = 0.5 * (deformation + deformation.T)
        if np.linalg.det(cell) <= 0.0 or np.min(np.linalg.eigvalsh(symmetric)) <= 0.0:
            raise ValueError("GaN path has an invalid or inverted cell")
        max_asymmetry = max(max_asymmetry, asymmetry)
        max_cell_reconstruction_error = max(
            max_cell_reconstruction_error,
            float(np.max(np.abs(reference_cell @ symmetric - cell))),
        )
        eta = symmetric - np.eye(3)
        strains.append([eta[0, 0], eta[1, 1], eta[2, 2],
                        eta[1, 2], eta[0, 2], eta[0, 1]])
    return strains, max_asymmetry, max_cell_reconstruction_error


def _reference_fractional_margin(images, reference) -> float:
    q0 = reference.get_scaled_positions(wrap=False)
    margins = []
    for image in images:
        delta = image.get_scaled_positions(wrap=False) - q0
        delta -= np.rint(delta)
        margins.append(float(np.min(0.5 - np.abs(delta))))
    return min(margins)


def audit(
    summary_path: Path, trajectory_path: Path, generic_audit_path: Path,
    route_analysis_path: Path, ts_audit_path: Path,
) -> tuple[dict, list[dict]]:
    ts_report = json.loads(ts_audit_path.read_text(encoding="utf-8"))
    if ts_report != audit_ts_candidate(
        summary_path, trajectory_path, generic_audit_path, route_analysis_path
    ):
        raise ValueError("stored TS-candidate audit differs from fresh source audit")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    images = latest_evaluated_chain(trajectory_path, 29)
    identity = validate_static_endpoint_identity(summary, images)
    if identity["status"] != "passed":
        raise ValueError("GaN chain endpoints differ from archived static endpoints")
    if (len(images) != 29 or any(len(image) != 4 for image in images)
            or any(image.get_chemical_symbols() != ["Ga", "Ga", "N", "N"] for image in images)):
        raise ValueError("expected aligned four-atom GaN chain with 29 images")
    masses = images[0].get_masses()
    initial_d, initial_translations = reference_cell_displacements(
        images, images[0], masses_amu=masses, remove_translations=True
    )
    final_d, final_translations = reference_cell_displacements(
        images, images[-1], masses_amu=masses, remove_translations=True
    )
    if (np.max(np.abs(initial_d[0])) > 1e-10
            or np.max(np.abs(final_d[-1])) > 1e-10
            or _reference_fractional_margin(images, images[0]) <= 1e-5
            or _reference_fractional_margin(images, images[-1]) <= 1e-5):
        raise ValueError("reference/gauge is discontinuous or near an ambiguous half cell")
    initial_fits = {rank: fit_path_atomic_svd(initial_d, masses, n_components=rank)
                    for rank in (1, 2, 3)}
    final_fit = fit_path_atomic_svd(final_d, masses, n_components=3)
    held_out = leave_one_image_out_residuals(initial_d, masses, n_components=3)
    strain, max_asym, max_cell_error = _cell_report(images)
    if max_asym > 1e-5 or max_cell_error > 1e-5:
        raise ValueError("cell rotation/asymmetry exceeds geometric reporting tolerance")
    cell_scale = float(summary["path_diagnostics"]["cell_scale_A"])
    coordinate, segments = path_reaction_coordinate(images, cell_scale=cell_scale)
    if len(coordinate) != 29 or not np.all(np.diff(coordinate) > 0.0):
        raise ValueError("GaN path reaction coordinate is not strictly increasing")
    normalized_s = coordinate / coordinate[-1]
    enthalpy = np.asarray(summary["image_enthalpies_eV"], dtype=float)
    if enthalpy.shape != (29,) or int(np.argmax(enthalpy)) != 15:
        raise ValueError("stored enthalpy does not identify the audited peak")
    relative_per_formula = (enthalpy - enthalpy[0]) / 2.0
    if abs(relative_per_formula[15] - ts_report["barrier_eV_per_GaN"]) > 1e-7:
        raise ValueError("peak barrier conversion disagrees with TS audit")
    fit3 = initial_fits[3]
    rows = []
    for index, image in enumerate(images):
        row = {
            "image_index": index,
            "reaction_coordinate_normalized": float(normalized_s[index]),
            "relative_enthalpy_eV_per_GaN": float(relative_per_formula[index]),
            "volume_A3_per_GaN": float(image.get_volume() / 2.0),
            "empirical_Q1_sqrt_amu_A": float(fit3.coefficients_sqrt_amu_A[index, 0]),
            "empirical_Q2_sqrt_amu_A": float(fit3.coefficients_sqrt_amu_A[index, 1]),
            "empirical_Q3_sqrt_amu_A": float(fit3.coefficients_sqrt_amu_A[index, 2]),
            "rank3_atomic_residual_sqrt_amu_A": float(fit3.residuals_sqrt_amu_A[index]),
            "rank3_held_out_residual_sqrt_amu_A": float(held_out[index]),
            "removed_translation_x_A": float(initial_translations[index, 0]),
            "removed_translation_y_A": float(initial_translations[index, 1]),
            "removed_translation_z_A": float(initial_translations[index, 2]),
            **{f"symmetric_strain_{label}": float(value)
               for label, value in zip(("xx", "yy", "zz", "yz", "xz", "xy"), strain[index])},
        }
        rows.append(row)
    report = {
        "status": "audited_GaN_path_geometric_empirical_axes_not_phonons",
        "route": "GaN B4-to-B1 tetragonal",
        "pressure_GPa": 45.7,
        "n_images_total": 29,
        "n_formula_units": 2,
        "peak_image_index": 15,
        "peak_enthalpy_eV_per_GaN": float(relative_per_formula[15]),
        "endpoint_identity": identity,
        "reference": "initial B4; endpoint atom order and periodic gauge as archived",
        "atomic_axis_kind": "mass_weighted_path_SVD_empirical_not_phonon",
        "atomic_singular_values_sqrt_amu_A": fit3.singular_values_sqrt_amu_A.tolist(),
        "captured_squared_atomic_path_norm_fraction_by_rank": {
            str(rank): fit.captured_squared_norm_fraction for rank, fit in initial_fits.items()
        },
        "max_atomic_reconstruction_residual_sqrt_amu_A_by_rank": {
            str(rank): float(np.max(fit.residuals_sqrt_amu_A))
            for rank, fit in initial_fits.items()
        },
        "max_rank3_held_out_residual_sqrt_amu_A": float(np.max(held_out)),
        "max_rank3_final_reference_residual_sqrt_amu_A": float(
            np.max(final_fit.residuals_sqrt_amu_A)
        ),
        "max_successive_internal_atomic_displacement_A": float(
            np.max(np.linalg.norm(np.diff(initial_d, axis=0), axis=2))
        ),
        "initial_reference_min_fractional_half_cell_margin": _reference_fractional_margin(images, images[0]),
        "final_reference_min_fractional_half_cell_margin": _reference_fractional_margin(images, images[-1]),
        "max_cell_deformation_antisymmetry": max_asym,
        "max_symmetric_cell_reconstruction_error_A": max_cell_error,
        "initial_cell_A": images[0].cell.array.tolist(),
        "final_cell_A": images[-1].cell.array.tolist(),
        "initial_reference_basis_mass_weighted": fit3.basis_mass_weighted.tolist(),
        "final_reference_rank3_captured_fraction": final_fit.captured_squared_norm_fraction,
        "max_removed_initial_reference_translation_A": float(
            np.max(np.linalg.norm(initial_translations, axis=1))
        ),
        "max_removed_final_reference_translation_A": float(
            np.max(np.linalg.norm(final_translations, axis=1))
        ),
        "reaction_coordinate_cell_scale_A": cell_scale,
        "reaction_coordinate_total_metric_A": float(coordinate[-1]),
        "smallest_segment_metric_A": float(np.min(segments)),
        "rows": rows,
        "full_variable_cell_TS_certified": False,
        "phonon_force_constants_available_in_this_analysis": False,
        "interpretation": (
            "Three data-driven atomic directions compactly describe the sampled chain; "
            "this is descriptive and not a phonon basis or an independent mechanism proof. "
            "Cell strain and volume are separate. Endpoint Gamma force constants and "
            "local full-variable-cell TS refinement/Hessian remain required."
        ),
        "source_sha256": {
            "summary": _sha256(summary_path),
            "trajectory": _sha256(trajectory_path),
            "generic_audit": _sha256(generic_audit_path),
            "route_analysis": _sha256(route_analysis_path),
            "ts_audit": _sha256(ts_audit_path),
            "auditor": _sha256(Path(__file__)),
        },
    }
    return report, rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("summary", "trajectory", "generic_audit", "route_analysis",
                 "ts_audit", "output", "csv"):
        parser.add_argument("--" + name.replace("_", "-"), type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.csv.exists():
        raise FileExistsError("refusing to overwrite GaN geometry-mode audit")
    report, rows = audit(
        args.summary, args.trajectory, args.generic_audit,
        args.route_analysis, args.ts_audit,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    with args.csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({
        "status": report["status"],
        "rank3_captured_fraction": report["captured_squared_atomic_path_norm_fraction_by_rank"]["3"],
        "max_rank3_held_out_residual_sqrt_amu_A": report["max_rank3_held_out_residual_sqrt_amu_A"],
        "max_cell_deformation_antisymmetry": report["max_cell_deformation_antisymmetry"],
    }))


if __name__ == "__main__":
    main()
