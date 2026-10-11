"""Bounded EMT work and actual frozen-frame checks; never launches HfO2 DFT."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT

from scripts.analyze_hfo2_clamped_reference import terminal_images
from vcneb import CalculatorModeEvaluator, ClampedModeCoordinates, clamped_plane_vcneb_boundary
from vcneb.active_curvature import ActiveJointCurvatureCoordinates

ROOT = Path(__file__).resolve().parents[1]
BASE = Path("benchmarks/hfo2_channels/20261008")
CASES = [
    ("clamped_T_terminal_E079_20261011", "c5653b4759b202199cb75739f99038aaa52b3ccc220ef226efbacc980e5d716e"),
    ("clamped_M_terminal_E068_20261011", "56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c"),
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(repository, workdir):
    if workdir.exists():
        raise FileExistsError("fresh verification workdir required")
    workdir.mkdir(parents=True)
    work_rows = []
    total_EMT_evaluations = 0
    for rotated in (False, True):
        for tilt in (False, True):
            for pressure in (0., .03):
                reference = bulk("Cu", "fcc", a=3.6, cubic=True)
                cell = reference.cell.array.copy()
                cell += [[0., .12, .09], [.04, 0., .07], [.19, -.13, .06]]
                if rotated:
                    rotation, _ = np.linalg.qr(np.random.default_rng(82).normal(size=(3, 3)))
                    cell = cell @ rotation.T
                reference.set_cell(cell, scale_atoms=True)
                boundary = clamped_plane_vcneb_boundary(len(reference), cell, allow_tilt=tilt)
                chart = ClampedModeCoordinates(reference, boundary=boundary, cell_scale_A=3.6)
                evaluator = CalculatorModeEvaluator(chart, lambda *args: EMT(),
                    workdir=workdir / f"EMT-r{int(rotated)}-tilt{int(tilt)}-p{pressure:g}",
                    calculator_id="Cu-EMT-verification-not-HfO2", pressure_eV_per_A3=pressure)
                y = np.random.default_rng(83).normal(size=chart.coordinate_count)*.01
                energy, gradient = evaluator(y)
                h = 1e-5
                numerical = []
                for i in range(chart.coordinate_count):
                    values = []
                    for sign in (1., -1.):
                        probe = y.copy(); probe[i] += sign*h
                        atoms = chart.to_atoms(probe)
                        boundary.validate_images([atoms])
                        atoms.calc = EMT()
                        values.append(atoms.get_potential_energy()+pressure*atoms.get_volume())
                    numerical.append((values[0]-values[1])/(2*h))
                error = float(np.max(np.abs(gradient-numerical)))
                if error > 3e-5:
                    raise ValueError("EMT independent energy-gradient work check failed")
                _, cached = evaluator(y)
                if evaluator.n_new_evaluations != 1 or not np.array_equal(cached, gradient):
                    raise ValueError("same-coordinate cache reuse failed")
                resumed = CalculatorModeEvaluator(chart, lambda *args: EMT(), workdir=evaluator.workdir,
                    calculator_id="Cu-EMT-verification-not-HfO2", pressure_eV_per_A3=pressure, resume=True)
                if resumed(y)[0] != energy or resumed.n_new_evaluations != 0:
                    raise ValueError("actual serialized-cache recovery failed")
                count = 1+2*chart.coordinate_count
                total_EMT_evaluations += count
                work_rows.append(dict(rotated=rotated, allow_tilt=tilt, pressure_eV_A3=pressure,
                    coordinate_count=chart.coordinate_count, EMT_evaluations=count,
                    maximum_work_gradient_error_eV_A=error, cache_recovery_checked=True))
    datasets = [(case, pin, terminal_images(repository / BASE / case, pin)) for case, pin in CASES]
    reference = datasets[0][2][0][-1]  # Actual same-contract clamped T endpoint.
    native_rows = []
    for case, pin, (images, observation, residual, receipt) in datasets:
        if observation["pressure_GPa"] != 0 or not observation["mechanical_boundary"]["allow_tilt"]:
            raise ValueError("registered P0/open-tilt training band required")
        boundary = clamped_plane_vcneb_boundary(len(reference),
            np.asarray(observation["mechanical_boundary"]["reference_cell_A"]), allow_tilt=True)
        chart = ClampedModeCoordinates(reference, boundary=boundary, cell_scale_A=observation["cell_scale_A"])
        primitive = ActiveJointCurvatureCoordinates.for_clamped_plane(reference,
            boundary=boundary, cell_scale_A=observation["cell_scale_A"])
        for i, atoms in enumerate(images):
            y = chart.from_atoms(atoms)
            recovered = chart.to_atoms(y)
            energy, gradient = chart.energy_and_gradient(y, atoms)
            expected = primitive.enthalpy_gradient(atoms, atoms.get_forces(), atoms.get_stress(voigt=False), 0.)
            geometry_error = float(np.max(np.abs(recovered.positions-atoms.positions)))
            cell_error = float(np.max(np.abs(recovered.cell.array-atoms.cell.array)))
            if (geometry_error > 1e-12 or cell_error > 1e-12 or energy != atoms.get_potential_energy()
                    or not np.array_equal(gradient, expected)):
                raise ValueError("actual ordered native-frame/chart replay differs")
            native_rows.append(dict(job=receipt["job_id"], step=receipt["terminal_step"], image_index=i,
                terminal_receipt_sha256=pin, coordinate_count=chart.coordinate_count,
                ordinary_fmax_eV_A=residual["fmax_eV_A"], energy_eV_cell=energy,
                max_ordered_position_roundtrip_error_A=geometry_error,
                max_cell_roundtrip_error_A=cell_error,
                native_active_gradient_primitive_byte_equal=True))
    sources = ["vcneb/clamped_mode_coordinates.py", "vcneb/mode_evaluator.py", "vcneb/active_curvature.py",
               "vcneb/epitaxial_boundary.py", "scripts/verify_clamped_mode_evaluator.py"]
    return dict(status="clamped_mode_calculator_flow_checked_not_material_conditional_PES",
        source_sha256={name: digest(repository / name) for name in sources},
        EMT_work=work_rows, EMT_evaluations=total_EMT_evaluations,
        native_portable_frames_reparsed=len(native_rows), native_frames=native_rows,
        native_energy_finite_differences_measured=False,
        native_gradient_agreement_is_reuse_of_existing_primitive_not_independent_DFT_validation=True,
        complete_native_runtime_and_six_asset_audits="historical E068/E079, not new local reads",
        new_DFT_calls=0, scheduler_mutations=0, production_runtime_or_inputs_modified=False,
        G3_selection_or_holdout_access=False, material_conditional_surface_computed=False,
        material_Hessian_or_TS_certified=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, default=ROOT)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("fresh verification report required")
    result = verify(args.repository, args.workdir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False); stream.write("\n")
    print(json.dumps({"EMT_evaluations": result["EMT_evaluations"],
        "native_frames_reparsed": result["native_portable_frames_reparsed"], "new_DFT_calls": 0}))


if __name__ == "__main__":
    main()
