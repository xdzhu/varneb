"""Bounded NumPy/ASE check, no pytest, calculator installation or DFT.

Run from an immutable repository source root. This checks implementation,
not the hafnia material model or the complete registered forecast batch.
"""

import json
from pathlib import Path

import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT

from scripts.relax_clamped_ase_endpoint import load_seed
from vcneb import (BiaxialClampedCurvatureCoordinates,
                   assemble_joint_directional_curvature,
                   clamped_plane_vcneb_boundary, joint_curvature_probes)


evaluations = 0


def evaluate(chart, delta):
    global evaluations
    atoms = chart.displaced(delta)
    atoms.calc = EMT()
    energy = atoms.get_potential_energy() + .03*atoms.get_volume()
    gradient = chart.enthalpy_gradient(delta, atoms, atoms.get_forces(),
                                      atoms.get_stress(voigt=False), .03)
    evaluations += 1
    return energy, gradient


checks = []
for tilt in (False, True):
    for rotated in (False, True):
        reference = bulk("Cu", "fcc", a=3.6, cubic=True)
        cell = reference.cell.array.copy()
        cell[0] += [0., .10, .15]
        cell[1] += [.05, 0., .08]
        cell[2] += [.11, -.07, 0.]
        if rotated:
            rotation, _ = np.linalg.qr(np.random.default_rng(81).normal(size=(3, 3)))
            cell = cell @ rotation.T
        reference.set_cell(cell, scale_atoms=True)
        boundary = clamped_plane_vcneb_boundary(len(reference), cell, allow_tilt=tilt)
        chart = BiaxialClampedCurvatureCoordinates(
            reference, boundary=boundary, cell_scale_A=3.6,
            anchor_strain=.01, strain_scale_A=4.2,
        )
        for deformed in (False, True):
            delta = np.zeros(chart.size)
            if deformed:
                delta = np.random.default_rng(17).normal(size=chart.size)*.007
            _, analytic = evaluate(chart, delta)
            numerical = []
            for i in range(chart.size):
                direction = np.eye(chart.size)[i]*1e-5
                plus, _ = evaluate(chart, delta+direction)
                minus, _ = evaluate(chart, delta-direction)
                numerical.append((plus-minus)/2e-5)
            error = float(np.max(np.abs(analytic-numerical)))
            assert error < 3e-5, error
            checks.append({"tilt": tilt, "rotated": rotated, "deformed": deformed,
                           "max_gradient_error_eV_A": error})

directions = np.zeros((chart.size, 2))
directions[:, 0] = chart.controlled_direction()
directions[0, 1], directions[3, 1] = 1/np.sqrt(2), -1/np.sqrt(2)
center = np.random.default_rng(9).normal(size=chart.size)*.005
mixed = []
for h in (1e-4, 2e-4):
    probes = joint_curvature_probes(chart, directions, step_A=h, center_delta_A=center)
    gradients = np.array([evaluate(chart, probe.delta_A)[1] for probe in probes])
    block = assemble_joint_directional_curvature(directions, gradients[::2],
                                                gradients[1::2], step_A=h)
    assert block.reciprocity_relative_defect < 1e-6
    assert np.linalg.norm(block.transverse_action) > 1e-3
    mixed.append(block)
mixed_difference = float(np.max(np.abs(mixed[0].raw_projected-mixed[1].raw_projected)))
assert mixed_difference < 2e-5

seed_root = Path("benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds")
seeds_checked = 0
for condition, epsilon in (("strain_0000", 0.), ("strain_p0100", .01)):
    for phase in ("T", "PO_plus", "PO_minus_T_preserving", "PO_minus_T_reversing", "M"):
        atoms, boundary, seed = load_seed(seed_root/condition/phase/"endpoint_seed.json")
        chart = BiaxialClampedCurvatureCoordinates(
            atoms, boundary=boundary, cell_scale_A=seed["cell_scale_A"],
            anchor_strain=epsilon, strain_scale_A=seed["cell_scale_A"],
        )
        delta = np.zeros(chart.size)
        delta[-1] = .0001*chart.strain_scale_A
        probe = chart.displaced(delta)
        np.testing.assert_allclose(probe.cell[:2], atoms.cell[:2]
                                   * (1+epsilon+.0001)/(1+epsilon), atol=1e-13)
        assert probe.calc is None and chart.size == 40 and chart.internal_size == 39
        assert probe.get_chemical_symbols() == atoms.get_chemical_symbols()
        np.testing.assert_array_equal(chart.internal_relaxation_basis()[-1], 0.)
        seeds_checked += 1

print(json.dumps({"status": "implementation_compatibility_passed",
                  "EMT_gradient_checks": checks, "EMT_evaluations": evaluations,
                  "mixed_reciprocity_relative_defects": [m.reciprocity_relative_defect for m in mixed],
                  "mixed_step_difference_eV_A2": mixed_difference,
                  "actual_HfO2_uncomputed_seed_geometries_checked": seeds_checked,
                  "new_DFT_calls": 0, "material_prediction_or_TS_certified": False}, indent=2))
