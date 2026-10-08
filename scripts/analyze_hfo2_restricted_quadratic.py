"""Retrospective two-direction G0 quadratic check; no new DFT or barrier claim."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from scripts.analyze_hfo2_joint_probe_reuse import EVIDENCE, ROOT, analyze, audit_raw
from scripts.audit_hfo2_static_replica import sha256
from vcneb.quadratic_reduction import condition_quadratic_energy


def analyze_restricted(manifest: dict, results: dict[int, dict]) -> dict:
    """Use only the same-reference measured 2x2 block, with its complement fixed."""
    joint = analyze(manifest, results)
    basis = np.array(manifest["directions"], dtype=float).T
    # This conversion is specific to the registered equal xx+yy direction.
    expected_cell = np.array([1., 1., 0., 0., 0., 0.]) / np.sqrt(2.)
    if (not np.allclose(basis[36:, 0], 0., rtol=0, atol=1e-12)
            or not np.allclose(basis[:36, 1], 0., rtol=0, atol=1e-12)
            or not np.allclose(basis[36:, 1], expected_cell, rtol=0, atol=1e-12)
            or manifest["formula_units"] != 4
            or not np.isfinite(manifest["center_energy_eV_cell"])):
        raise ValueError("registered atomic secant, equal xx+yy strain and four formula units required")
    gradient0 = np.array(manifest["center_gradient_eV_A"], dtype=float)
    g2 = basis.T @ gradient0
    jacobian = float(np.sqrt(2.) * manifest["cell_scale_A"])
    # A measured operational screening floor, not a statistical uncertainty bound.
    floor = joint["two_step_full_action_operator_spread_eV_A2"]
    models, release_records = [], []
    for item in joint["step_records"]:
        h2 = np.array(item["symmetric_projected_eV_A2"])
        model = condition_quadratic_energy(
            h2, g2, manifest["center_energy_eV_cell"], np.eye(2)[:, 1:],
            np.eye(2)[:, :1], stability_floor=floor,
        )
        models.append(model)
        release_records.append({
            "step_A": item["step_A"],
            "frozen_cell_curvature_eV_A2": float(model.curvature.frozen[0, 0]),
            "released_cell_curvature_eV_A2": float(model.curvature.relaxed[0, 0]),
            "softening_eV_A2": float(model.curvature.softening[0, 0]),
            "softening_fraction": float(model.curvature.softening[0, 0] / model.curvature.frozen[0, 0]),
            "released_atomic_offset_A": float(model.eliminated_offset[0]),
            "atomic_response_per_cell_coordinate": float(model.curvature.orthogonal_response[0, 0]),
            "release_energy_change_at_zero_cell_coordinate_eV_cell": model.reference_energy_change,
            "retained_gradient_at_zero_eV_A": float(model.retained_gradient_at_zero[0]),
            "frozen_biaxial_strain_curvature_eV_cell": float(model.curvature.frozen[0, 0] * jacobian**2),
            "restricted_released_biaxial_strain_curvature_eV_cell": float(model.curvature.relaxed[0, 0] * jacobian**2),
            "restricted_released_biaxial_strain_gradient_eV_cell": float(model.retained_gradient_at_zero[0] * jacobian),
        })
    hshort = np.array(joint["step_records"][0]["symmetric_projected_eV_A2"])
    action = np.array(joint["step_records"][0]["full_hessian_action_eV_A2"])
    checks = []
    for point in manifest["points"]:
        if point["step_A"] != .02:
            continue
        z = np.zeros(2)
        z[point["direction"]] = point["sign"] * point["step_A"]
        observed = results[point["index"]]
        actual_delta = float(observed["energy_eV_cell"] - manifest["center_energy_eV_cell"])
        predicted_delta = float(g2 @ z + .5 * z @ hshort @ z)
        gradient_error = g2 + hshort @ z - basis.T @ observed["gradient_eV_A"]
        full_gradient_error = gradient0 + action @ z - observed["gradient_eV_A"]
        checks.append({
            "index": point["index"], "direction": point["direction"],
            "coordinate_A": z.tolist(), "observed_energy_change_eV_cell": actual_delta,
            "predicted_energy_change_eV_cell": predicted_delta,
            "energy_residual_eV_cell": predicted_delta - actual_delta,
            "projected_gradient_residual_eV_A": gradient_error.tolist(),
            "projected_gradient_residual_norm_eV_A": float(np.linalg.norm(gradient_error)),
            "full_gradient_action_residual_norm_eV_A": float(np.linalg.norm(full_gradient_error)),
        })
    # The axis probes' convex hull is |z_atomic|+|z_cell| <= .02 A, not a box.
    # Minimise this piecewise linear norm on the requested retained interval
    # by evaluating its endpoints and its two possible kinks.
    halfspan = .02
    model = models[0]
    offset = float(model.eliminated_offset[0])
    response = float(model.curvature.orthogonal_response[0, 0])
    candidates = [-halfspan, 0., halfspan]
    if response != 0. and abs(-offset / response) <= halfspan:
        candidates.append(-offset / response)
    l1 = [float(np.linalg.norm(model.full_displacement(np.array([x])), ord=1)) for x in candidates]
    minimum = min(l1)
    responses = []
    for x in (-halfspan, -.01, 0., .01, halfspan):
        z = model.full_displacement(np.array([x]))
        responses.append({
            "cell_coordinate_A": x, "biaxial_strain": x / jacobian,
            "atomic_coordinate_A": float(z[0]), "model_energy_change_eV_cell": model.energy_change(np.array([x])),
            "axis_probe_convex_hull_L1_A": float(np.linalg.norm(z, ord=1)),
            "within_axis_probe_convex_hull": bool(np.linalg.norm(z, ord=1) <= halfspan + 1e-12),
            "conditional_DFT_evaluation_available": False,
        })
    return {
        "format_version": 1, "status": "same_center_restricted_quadratic_pilot_only",
        "scope": "post_hoc_G0_free_cell_slice_not_G2_clamped_training_or_prospective_prediction",
        "direction_semantics": manifest["direction_semantics"],
        "center_energy_eV_cell": manifest["center_energy_eV_cell"],
        "center_projected_gradient_eV_A": g2.tolist(),
        "center_full_gradient_norm_eV_A": joint["center_full_gradient_norm_eV_A"],
        "pressure_eV_A3": 0., "formula_units": 4,
        "reference_source_directory": manifest["source_directory"],
        "reference_STRU_sha256": manifest["source_STRU_sha256"],
        "reference_log_sha256": manifest["source_log_sha256"],
        "probe_log_sha256": joint["probe_log_sha256"],
        "units": {"coordinates": "Angstrom", "energy": "eV/cell", "gradient": "eV/Angstrom",
                  "hessian": "eV/Angstrom**2", "biaxial_strain": "dimensionless",
                  "strain_curvature": "eV/cell per dimensionless strain squared"},
        "cell_coordinate_per_biaxial_strain_A": jacobian,
        "unit_contract": "z_cell=sqrt(2)*cell_scale_A*epsilon_xx=yy; gradients multiply and curvatures multiply by its square",
        "measured_subspace_dimension": 2, "retained_cell_dimension": 1,
        "released_atomic_dimension": 1, "clamped_unmeasured_atomic_dimension": 35,
        "clamped_unmeasured_cell_dimension": 5,
        "other_reference_Hessians_merged": False,
        "joint_curvature_audit": joint,
        "operational_eliminated_stability_floor_eV_A2": floor,
        "floor_is_total_DFT_error_bound": False,
        "release_records": release_records,
        "retrospective_axis_crosschecks": checks,
        "maximum_axis_energy_residual_eV_cell": max(abs(x["energy_residual_eV_cell"]) for x in checks),
        "maximum_axis_projected_gradient_residual_eV_A": max(x["projected_gradient_residual_norm_eV_A"] for x in checks),
        "long_step_data_seen_before_model": True,
        "short_step_fit_inputs": "center_energy_and_gradient_plus_short_step_full_gradient_columns; no long_step_energy_fit",
        "stability_screen_uses_two_step_data": True,
        "response_samples": responses,
        "sampling_coverage": {"axis_convex_hull_radius_A": halfspan,
                              "requested_cell_interval_A": [-halfspan, halfspan],
                              "minimum_release_L1_A_over_interval": minimum,
                              "all_released_coordinates_outside_probe_convex_hull": bool(minimum > halfspan + 1e-12)},
        "new_DFT_calls": 0, "parameters_changed": False,
        "full_joint_Hessian_measured": False, "unsampled_stability_known": False,
        "saddle_certified": False, "conditional_DFT_surface_established": False,
        "independent_barrier_prediction_validated": False,
        "interpretation": "One measured atomic direction softens the retained strain curvature, but nonstationary offset and sampling coverage preclude a validated released DFT branch. All other directions remain clamped; Schur condensation is established mathematics, not a novelty claim.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("refusing an existing analysis receipt")
    if args.raw_root is None:
        manifest = json.loads((EVIDENCE / "work_probe_manifest.json").read_text())
        results = {i: json.loads((EVIDENCE / "point_audits/calculations" / f"{i:02d}" / "point_audit.json").read_text()) for i in range(8)}
    else:
        manifest, results = audit_raw(args.raw_root)
    report = analyze_restricted(manifest, results)
    report.update({
        "raw_source_audit_performed": args.raw_root is not None,
        "raw_root": str(args.raw_root) if args.raw_root is not None else None,
        "archived_manifest_sha256": sha256(EVIDENCE / "work_probe_manifest.json"),
        "raw_manifest_sha256": sha256(args.raw_root / "probe_manifest.json") if args.raw_root is not None else None,
        "analysis_source_sha256": sha256(Path(__file__)),
        "quadratic_model_source_sha256": sha256(ROOT / "vcneb/quadratic_reduction.py"),
        "stable_reduction_source_sha256": sha256(ROOT / "vcneb/relaxed_curvature.py"),
        "joint_probe_analysis_source_sha256": sha256(ROOT / "scripts/analyze_hfo2_joint_probe_reuse.py"),
    })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), "new_DFT_calls": 0,
                      "softening_fraction": report["release_records"][0]["softening_fraction"],
                      "release_outside_coverage": report["sampling_coverage"]["all_released_coordinates_outside_probe_convex_hull"]}))


if __name__ == "__main__":
    main()
