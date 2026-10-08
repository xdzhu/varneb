"""Zero-DFT implementation checks, runnable on HF without pytest."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import ase
import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT

from scripts.relax_clamped_ase_endpoint import load_seed
from vcneb import (ActiveJointCurvatureCoordinates, assemble_joint_directional_curvature,
                   clamped_plane_vcneb_boundary, joint_curvature_probes)
from vcneb.mode_surface import audit_directional_curvature_consistency


ROOT = Path(__file__).resolve().parents[1]


def verify(seed_root: Path) -> dict:
    records = json.loads((seed_root / "clamped_seed_manifest.json").read_text())["seeds"]
    if len(records) != 10:
        raise ValueError("expected exactly ten registered unrelaxed geometry starters")
    geometries = []
    for item in records:
        atoms, boundary, record = load_seed(seed_root / item["manifest_file"])
        chart = ActiveJointCurvatureCoordinates.for_clamped_plane(
            atoms, cell_scale_A=record["cell_scale_A"], boundary=boundary,
        )
        probes = joint_curvature_probes(chart, chart.translation_free_basis(), step_A=.005)
        if len(probes) != 72 or any(probe.atoms.calc is not None for probe in probes):
            raise AssertionError("expected 72 inert translation-free probes per starter")
        drift = 0.
        for probe in probes:
            boundary.validate_images([probe.atoms])
            if probe.atoms.get_chemical_symbols() != atoms.get_chemical_symbols():
                raise AssertionError("probe changed atom order")
            drift = max(drift, float(np.max(np.abs(probe.atoms.cell[:2] - atoms.cell[:2]))))
        geometries.append({"phase_label": record["phase_label"], "strain": record["strain"],
                           "seed_sha256": record["seed_sha256"], "status": record["status"],
                           "probe_count": len(probes), "max_substrate_drift_A": drift})

    # Two measured columns cannot identify the missing diagonal, regardless
    # of whether their own projected Hessian is positive or indefinite.
    matrix = np.array([[2., .3, .7], [.3, -1., .5], [.7, .5, 4.]])
    basis = np.eye(3)[:, :2]
    h = .01
    gradients = (h * matrix @ basis).T
    selected = assemble_joint_directional_curvature(basis, gradients, -gradients, step_A=h)
    action_error = float(np.max(np.abs(selected.full_hessian_action - matrix @ basis)))
    if action_error > 1e-12:
        raise AssertionError("harmonic Hessian action disagrees")
    if not np.allclose(selected.transverse_action[2], [.7, .5], rtol=0, atol=1e-12):
        raise AssertionError("unsampled coupling was discarded")

    calls, work = 0, []
    for tilt in (False, True):
        for rotated in (False, True):
            atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
            cell = atoms.cell.array.copy()
            cell[0] += [0., .10, .15]
            cell[1] += [.05, 0., .08]
            cell[2] += [.11, -.07, 0.]
            if rotated:
                rotation, _ = np.linalg.qr(np.random.default_rng(81).normal(size=(3, 3)))
                cell = cell @ rotation.T
            atoms.set_cell(cell, scale_atoms=True)
            boundary = clamped_plane_vcneb_boundary(len(atoms), cell, allow_tilt=tilt)
            chart = ActiveJointCurvatureCoordinates.for_clamped_plane(atoms, cell_scale_A=3.6, boundary=boundary)
            directions = chart.translation_free_basis()
            center = np.zeros(chart.size)
            center[0], center[-1] = .003, -.004
            origin = chart.displaced(center)
            origin.calc = EMT()
            pressure = .03  # Cu analytic check only; HfO2 remains P=0
            energy0 = origin.get_potential_energy() + pressure * origin.get_volume()
            calls += 1
            pairs, curvatures = [], []
            steps = (.001, .002)
            for step in steps:
                energies, gradients = [], []
                for probe in joint_curvature_probes(chart, directions, step_A=step, center_delta_A=center):
                    point = probe.atoms
                    boundary.validate_images([point])
                    point.calc = EMT()
                    energies.append(point.get_potential_energy() + pressure * point.get_volume())
                    gradients.append(chart.enthalpy_gradient(point, point.get_forces(), point.get_stress(voigt=False), pressure))
                    calls += 1
                values = np.array(gradients).reshape(directions.shape[1], 2, chart.size)
                pairs.append((np.array(energies).reshape(-1, 2), values))
                curvatures.append(assemble_joint_directional_curvature(directions, values[:, 0], values[:, 1], step_A=step))
            maximum_difference = 0.
            for index in range(directions.shape[1]):
                audit = audit_directional_curvature_consistency(
                    energy0, np.array([pair[0][index] for pair in pairs]),
                    np.array([pair[1][index] @ directions[:, index] for pair in pairs]), steps,
                )
                maximum_difference = max(maximum_difference, float(audit.absolute_energy_gradient_disagreement.max()))
            reciprocity = max(curvature.reciprocity_relative_defect for curvature in curvatures)
            spread = float(np.linalg.norm(curvatures[0].symmetric_projected - curvatures[1].symmetric_projected, ord=2))
            if maximum_difference >= .002 or reciprocity >= 1e-4 or spread >= .005:
                raise AssertionError("Cu/EMT implementation curvature check failed")
            work.append({"allow_tilt": tilt, "rotated": rotated,
                         "translation_free_directions": directions.shape[1],
                         "max_energy_gradient_curvature_difference_eV_A2": maximum_difference,
                         "max_projected_reciprocity_relative_defect": reciprocity,
                         "two_step_projected_operator_spread_eV_A2": spread})
    if max(record["max_substrate_drift_A"] for record in geometries) > 1e-12:
        raise AssertionError("fixed substrate drifted")
    return {
        "format_version": 1, "ASE_version": ase.__version__, "new_DFT_calls": 0,
        "geometry_starters": geometries, "total_inert_HfO2_geometry_probes": 720,
        "HfO2_probe_step_A_is_geometry_check_only": .005,
        "EMT_evaluations": calls, "EMT_two_step_checks": work,
        "synthetic_action_max_error_eV_A2": action_error,
        "synthetic_unsampled_direction_self_curvature_not_measured": True,
        "material_joint_Hessian_measured": False, "saddle_certified": False,
        "material_barrier_prediction_validated": False,
        "module_sha256": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                          for name in ("vcneb/__init__.py", "vcneb/joint_curvature.py", "vcneb/active_curvature.py",
                                       "vcneb/joint_stencil.py", "vcneb/mode_surface.py",
                                       "scripts/verify_joint_curvature_stencil.py")},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=Path, default=ROOT / "benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    if args.report is not None and args.report.exists():
        raise FileExistsError("refusing an existing verification receipt")
    result = verify(args.seeds)
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.report is None:
        print(payload)
    else:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(payload, encoding="utf-8")
        print(json.dumps({"report": str(args.report), "new_DFT_calls": 0,
                          "geometry_probes": result["total_inert_HfO2_geometry_probes"],
                          "EMT_evaluations": result["EMT_evaluations"]}))


if __name__ == "__main__":
    main()
