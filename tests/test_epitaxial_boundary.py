"""Fixed substrate vectors must survive optimizer proposals and stress work."""

from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read

from vcneb import VCNEB, ClampedPlaneFilter, clamped_plane_vcneb_boundary, cell_work_derivative
from vcneb.core import cell_force_from_arrays


def _cell():
    # No Cartesian component mask can enforce this rotated, oblique plane.
    return np.array([[4.0, 0.2, 0.7], [0.5, 4.2, 1.1], [0.2, -0.3, 5.1]])


@pytest.mark.parametrize("allow_tilt", [False, True])
def test_exact_substrate_subspace_under_arbitrary_proposals(allow_tilt):
    cell = _cell()
    boundary = clamped_plane_vcneb_boundary(2, cell, allow_tilt=allow_tilt)
    assert boundary.cell_dofs == (3 if allow_tilt else 1)
    basis = boundary.mode_basis
    np.testing.assert_allclose(basis.T @ basis, np.eye(basis.shape[1]), atol=1e-14)
    assert np.array_equal(boundary.cell_mask, np.ones((3, 3)))
    images = [Atoms("Ar2", scaled_positions=[[0.1, 0.2, x], [0.4, 0.5, 0.6]],
                    cell=cell, pbc=True) for x in (0.2, 0.3, 0.4)]
    neb = VCNEB(images, cell_scale=5.0, **boundary.vcneb_kwargs(images))
    for seed in range(5):
        proposal = neb.get_x() + np.random.default_rng(seed).normal(size=15) * 0.2
        candidate = neb.candidate_images_from_x(proposal)[1]
        np.testing.assert_allclose(candidate.cell[:2], cell[:2], atol=1e-13)
        if not allow_tilt:
            delta = candidate.cell[2] - cell[2]
            np.testing.assert_allclose(delta, np.dot(delta, boundary.normal) * boundary.normal,
                                       atol=1e-13)
        boundary.validate_images([candidate])
    assert not any(a.flags.writeable for a in
                   (boundary.reference_cell, boundary.normal, basis, boundary.cell_mask))


def test_all_three_row2_directions_remain_available_without_rotation():
    cell = _cell()
    boundary = clamped_plane_vcneb_boundary(1, cell, allow_tilt=True)
    changes = [cell @ basis.reshape(3, 3).T for basis in boundary.mode_basis[3:].T[3:]]
    assert np.linalg.matrix_rank(np.array(changes)[:, 2, :]) == 3
    for change in changes:
        np.testing.assert_allclose(change[:2], 0.0, atol=1e-14)
    rotation = np.array([[0., -1., 2.], [1., 0., -3.], [-2., 3., 0.]])
    assert np.linalg.norm(cell[:2] @ rotation.T) > 0


@pytest.mark.parametrize("bad_index", [0, 1, 2])
def test_reject_incompatible_plane_before_neb_can_project_it(bad_index):
    images = [Atoms("Ar", positions=[[1, 1, 1]], cell=_cell(), pbc=True) for _ in range(3)]
    before = [im.cell.array.copy() for im in images]
    images[bad_index].cell[0, 0] += 0.001
    boundary = clamped_plane_vcneb_boundary(1, _cell(), allow_tilt=True)
    with pytest.raises(ValueError, match=f"image {bad_index}: substrate"):
        boundary.vcneb_kwargs(images)
    assert images[bad_index].cell[0, 0] == pytest.approx(before[bad_index][0, 0] + 0.001)


def test_normal_only_rejects_tilt_in_seed_including_interior():
    boundary = clamped_plane_vcneb_boundary(1, _cell(), allow_tilt=False)
    image = Atoms("Ar", cell=_cell(), pbc=True)
    image.cell[2] += 0.01 * _cell()[0]
    with pytest.raises(ValueError, match="tilt"):
        boundary.vcneb_kwargs([image])


def test_actual_released_hafnia_endpoints_are_not_silently_clamped():
    root = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008"
    initial = read(root / "reference_variants/T.vasp")
    final = read(root / "reference_variants/PO.vasp")
    boundary = clamped_plane_vcneb_boundary(len(initial), initial.cell.array, allow_tilt=True)
    with pytest.raises(ValueError, match="substrate"):
        boundary.vcneb_kwargs([initial, final])


@pytest.mark.parametrize("allow_tilt", [False, True])
def test_open_cell_gradient_matches_finite_difference_in_oblique_cell(allow_tilt):
    reference = bulk("Cu", "fcc", a=3.6, cubic=True)
    cell = reference.cell.array.copy()
    cell[0] += [0.0, 0.10, 0.15]
    cell[1] += [0.05, 0.0, 0.08]
    reference.set_cell(cell, scale_atoms=True)
    boundary = clamped_plane_vcneb_boundary(len(reference), cell, allow_tilt=allow_tilt)
    basis = boundary.mode_basis[3 * len(reference):, 3 * len(reference):]
    cell_scale, pressure = 3.6, 0.03
    for direction in basis.T:
        delta_F = direction.reshape(3, 3) / cell_scale
        h = 1e-5
        pair = []
        for sign in (1, -1):
            atoms = reference.copy()
            atoms.set_cell(cell @ (np.eye(3) + sign * h * delta_F).T, scale_atoms=True)
            atoms.calc = EMT()
            pair.append(atoms.get_potential_energy() + pressure * atoms.get_volume())
        reference.calc = EMT()
        stress = reference.get_stress(voigt=False)
        generalized_force = cell_force_from_arrays(stress, reference, cell, pressure=pressure)
        expected = -float(np.sum(generalized_force * delta_F))
        assert expected == pytest.approx((pair[0] - pair[1]) / (2 * h), abs=3e-6)
        assert cell_work_derivative(stress, cell, cell @ delta_F.T, pressure=pressure) == pytest.approx(
            expected, abs=1e-12)


def test_clamp_reaction_stress_does_not_enter_neb_residual_but_open_stress_does():
    cell = np.diag([4., 4., 5.])
    images = [Atoms("Ar", positions=[[1., 1., x]], cell=cell, pbc=True) for x in (1., 2., 3.)]
    for image in images:
        image.calc = SinglePointCalculator(image, energy=0., forces=np.zeros((1, 3)),
                                           stress=np.diag([10., -20., 0.]))
    boundary = clamped_plane_vcneb_boundary(1, cell, allow_tilt=True)
    neb = VCNEB(images, cell_scale=1., **boundary.vcneb_kwargs(images))
    np.testing.assert_allclose(neb.get_forces(), 0., atol=1e-12)
    images[1].calc = SinglePointCalculator(images[1], energy=0., forces=np.zeros((1, 3)),
                                           stress=np.diag([10., -20., 0.2]))
    neb = VCNEB(images, cell_scale=1., **boundary.vcneb_kwargs(images))
    assert np.linalg.norm(neb.get_forces()) > 1.


def test_biaxial_work_derivative_finite_difference_at_nonzero_strain():
    base = bulk("Cu", "fcc", a=3.6, cubic=True)
    epsilon = 0.012
    cell = base.cell.array.copy()
    cell[:2] *= 1 + epsilon
    cell[2] += [0.05, 0.06, -0.01]
    base.set_cell(cell, scale_atoms=True)
    base.calc = EMT()
    derivative = np.zeros((3, 3))
    derivative[:2] = cell[:2] / (1 + epsilon)
    pressure, h = 0.03, 1e-5
    expected = cell_work_derivative(base.get_stress(voigt=False), cell, derivative,
                                    pressure=pressure)
    pair = []
    for sign in (1, -1):
        atoms = base.copy()
        atoms.set_cell(cell + sign * h * derivative, scale_atoms=True)
        atoms.calc = EMT()
        pair.append(atoms.get_potential_energy() + pressure * atoms.get_volume())
    assert expected == pytest.approx((pair[0] - pair[1]) / (2 * h), abs=3e-6)


@pytest.mark.parametrize("n_atoms,cell,tilt,match", [
    (0, _cell(), True, "positive integer"),
    (True, _cell(), True, "positive integer"),
    (1, np.zeros((3, 3)), True, "positive volume"),
    (1, np.full((3, 3), np.nan), True, "finite"),
    (1, _cell(), 1, "boolean"),
])
def test_invalid_boundary_inputs(n_atoms, cell, tilt, match):
    with pytest.raises(ValueError, match=match):
        clamped_plane_vcneb_boundary(n_atoms, cell, allow_tilt=tilt)


def test_invalid_work_inputs():
    with pytest.raises(ValueError, match="symmetric"):
        cell_work_derivative(np.arange(9).reshape(3, 3), _cell(), np.eye(3))
    with pytest.raises(ValueError, match="finite"):
        cell_work_derivative(np.eye(3), _cell(), np.full((3, 3), np.nan))


@pytest.mark.parametrize("allow_tilt", [False, True])
def test_endpoint_filter_energy_gradient_at_deformed_point(allow_tilt):
    atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
    cell = atoms.cell.array.copy()
    cell[0] += [0.0, 0.10, 0.15]
    cell[1] += [0.05, 0.0, 0.08]
    atoms.set_cell(cell, scale_atoms=True)
    atoms.calc = EMT()
    boundary = clamped_plane_vcneb_boundary(len(atoms), cell, allow_tilt=allow_tilt)
    target = ClampedPlaneFilter(atoms, boundary, cell_scale_A=3.6, pressure_eV_per_A3=0.03)
    center = target.get_positions()
    center[0, 0] += 0.01
    center[-3:] += 0.01
    target.set_positions(center)
    center = target.get_positions()
    gradient = -target.get_forces()
    for index in range(center.size):
        pair = []
        for sign in (1, -1):
            trial = center.copy().reshape(-1)
            trial[index] += sign * 1e-5
            target.set_positions(trial.reshape(center.shape))
            pair.append(target.get_potential_energy())
        assert gradient.reshape(-1)[index] == pytest.approx((pair[0] - pair[1]) / 2e-5, abs=3e-5)
    target.set_positions(center)
    np.testing.assert_allclose(atoms.cell[:2], cell[:2], atol=1e-13)


def test_bfgs_endpoint_keeps_substrate_at_every_step():
    from ase.optimize import BFGS

    atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
    cell = atoms.cell.array.copy()
    cell[:2] *= 1.01
    cell[2, 2] *= 1.03
    atoms.set_cell(cell, scale_atoms=True)
    atoms.positions[0, 2] += 0.02
    atoms.calc = EMT()
    boundary = clamped_plane_vcneb_boundary(len(atoms), cell, allow_tilt=True)
    target = ClampedPlaneFilter(atoms, boundary, cell_scale_A=3.6)
    optimizer = BFGS(target, logfile=None, maxstep=0.03)
    observations = []

    def check():
        boundary.validate_images([atoms])
        observations.append(atoms.cell.array.copy())

    optimizer.attach(check)
    assert optimizer.run(fmax=0.02, steps=30)
    assert len(observations) > 2
    for observed in observations:
        np.testing.assert_allclose(observed[:2], cell[:2], atol=1e-13)


def test_filter_rejects_invalid_step_without_mutating_atoms():
    atoms = Atoms("Ar", positions=[[1, 1, 1]], cell=np.diag([4., 4., 5.]), pbc=True)
    boundary = clamped_plane_vcneb_boundary(1, atoms.cell.array, allow_tilt=True)
    target = ClampedPlaneFilter(atoms, boundary, cell_scale_A=5.)
    positions, cell = atoms.positions.copy(), atoms.cell.array.copy()
    proposal = target.get_positions()
    proposal[-1, -1] = -5.
    with pytest.raises(ValueError, match="positive volume"):
        target.set_positions(proposal)
    np.testing.assert_array_equal(atoms.positions, positions)
    np.testing.assert_array_equal(atoms.cell.array, cell)
