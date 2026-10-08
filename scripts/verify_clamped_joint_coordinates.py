"""Reproducible zero-DFT checks for the joint probe/path mechanical subspace."""

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
from vcneb import clamped_plane_vcneb_boundary
from vcneb.active_curvature import ActiveJointCurvatureCoordinates as JointCurvatureCoordinates


ROOT = Path(__file__).resolve().parents[1]
WORK_TOLERANCE = 3e-5  # eV/A; prespecified implementation check, not a DFT bound


def verify(seed_root: Path) -> dict:
    records = json.loads((seed_root / "clamped_seed_manifest.json").read_text())["seeds"]
    if len(records) != 10:
        raise ValueError("expected the ten registered, unrelaxed G2 geometry starters")
    geometry = []
    for item in records:
        atoms, boundary, record = load_seed(seed_root / item["manifest_file"])
        chart = JointCurvatureCoordinates.for_clamped_plane(
            atoms, cell_scale_A=record["cell_scale_A"], boundary=boundary,
        )
        block = boundary.mode_basis[3 * len(atoms):, 3 * len(atoms):]
        difference = float(np.max(np.abs(chart.deformation_basis.reshape(3, 9).T - block)))
        drift = 0.0
        for trial in range(5):
            delta = np.random.default_rng(trial).normal(size=chart.size) * 0.005
            probe = chart.displaced(delta)
            boundary.validate_images([probe])
            drift = max(drift, float(np.max(np.abs(probe.cell[:2] - atoms.cell[:2]))))
            if probe.calc is not None or probe.get_chemical_symbols() != atoms.get_chemical_symbols():
                raise AssertionError("probe inherited a calculator or changed atom identity")
        geometry.append({
            "phase_label": record["phase_label"], "strain": record["strain"],
            "seed_sha256": record["seed_sha256"], "status": record["status"],
            "cell_dofs": chart.cell_dofs, "joint_size": chart.size,
            "translation_free_size": chart.translation_free_basis().shape[1],
            "probe_path_basis_max_difference": difference, "max_substrate_drift_A": drift,
        })

    work = []
    calls = 0
    for rotated in (False, True):
        reference = bulk("Cu", "fcc", a=3.6, cubic=True)
        cell = reference.cell.array.copy()
        cell[0] += [0.0, 0.10, 0.15]
        cell[1] += [0.05, 0.0, 0.08]
        cell[2] += [0.11, -0.07, 0.0]
        if rotated:
            rotation, _ = np.linalg.qr(np.random.default_rng(81).normal(size=(3, 3)))
            cell = cell @ rotation.T
        reference.set_cell(cell, scale_atoms=True)
        for tilt in (False, True):
            boundary = clamped_plane_vcneb_boundary(len(reference), cell, allow_tilt=tilt)
            chart = JointCurvatureCoordinates.for_clamped_plane(
                reference, cell_scale_A=3.6, boundary=boundary,
            )
            for deformed in (False, True):
                center = np.zeros(chart.size)
                if deformed:
                    center = np.random.default_rng(314).normal(size=chart.size) * 0.015
                atoms = chart.displaced(center)
                atoms.calc = EMT()
                pressure = 0.03  # analytic Cu check only, not an HfO2 experiment
                gradient = chart.enthalpy_gradient(
                    atoms, atoms.get_forces(), atoms.get_stress(voigt=False), pressure,
                )
                calls += 1
                numeric = []
                h = 1e-5
                for index in range(chart.size):
                    energies = []
                    for sign in (1, -1):
                        delta = center.copy()
                        delta[index] += sign * h
                        probe = chart.displaced(delta)
                        boundary.validate_images([probe])
                        probe.calc = EMT()
                        energies.append(probe.get_potential_energy() + pressure * probe.get_volume())
                        calls += 1
                    numeric.append((energies[0] - energies[1]) / (2 * h))
                error = float(np.max(np.abs(gradient - numeric)))
                if error > WORK_TOLERANCE:
                    raise AssertionError(f"joint stress/force work error {error} exceeds {WORK_TOLERANCE}")
                work.append({"rotated": rotated, "allow_tilt": tilt, "deformed": deformed,
                             "max_work_gradient_error_eV_A": error})

    if max(r["max_substrate_drift_A"] for r in geometry) > 1e-12:
        raise AssertionError("joint perturbations drifted the prescribed substrate")
    return {
        "format_version": 1, "ASE_version": ase.__version__, "new_DFT_calls": 0,
        "EMT_evaluations": calls, "geometry_starters_checked": len(geometry),
        "material_joint_Hessian_measured": False, "material_barrier_prediction_validated": False,
        "saddle_certified": False, "work_tolerance_eV_A": WORK_TOLERANCE,
        "geometry": geometry, "EMT_work": work,
        "module_sha256": {
            name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in ("vcneb/__init__.py", "vcneb/joint_curvature.py",
                         "vcneb/active_curvature.py", "vcneb/epitaxial_boundary.py",
                         "scripts/verify_clamped_joint_coordinates.py")
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=Path, default=ROOT / "benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    if args.report is not None and args.report.exists():
        raise FileExistsError(f"refusing existing verification report: {args.report}")
    result = verify(args.seeds)
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.report is None:
        print(payload, end="")
    else:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(payload, encoding="utf-8")
        print(json.dumps({"report": str(args.report), "new_DFT_calls": 0,
                          "geometry_starters_checked": result["geometry_starters_checked"],
                          "max_work_gradient_error_eV_A": max(
                              row["max_work_gradient_error_eV_A"] for row in result["EMT_work"])}))


if __name__ == "__main__":
    main()
