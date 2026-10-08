import numpy as np
import pytest
from ase import Atoms

from vcneb.continuous_projection import continuous_reference_coordinates, project_reference_basis


def toy_chain():
    ref = Atoms("NaCl", scaled_positions=[[.1, .1, .1], [.6, .6, .6]], cell=np.diag([5., 6., 7.]), pbc=True)
    images = []
    for amplitude in [0, .2, .4, .6, .8, 1.]:
        a = ref.copy()
        q = ref.get_scaled_positions(wrap=False)
        q[0, 0] += amplitude
        a.set_scaled_positions(q)
        images.append(a)
    return ref, images


def test_full_winding_is_not_folded_image_by_image():
    ref, images = toy_chain()
    before = [(a.cell.array.copy(), a.positions.copy()) for a in images]
    result = continuous_reference_coordinates(images, ref)
    assert result["displacements_A"][-1, 0, 0] == pytest.approx(5.)
    assert np.allclose(result["displacements_A"][:, 0, 0], np.arange(6))
    assert np.array_equal(result["fixed_integer_lattice_shifts_by_atom"], np.zeros((2, 3)))
    for a, (cell, positions) in zip(images, before):
        assert np.array_equal(a.cell.array, cell) and np.array_equal(a.positions, positions)


def test_single_integer_representative_change_is_gauge_invariant():
    ref, images = toy_chain()
    expected = continuous_reference_coordinates(images, ref)
    for a in images:
        q = a.get_scaled_positions(wrap=False)
        q += [[1, -2, 0], [-1, 0, 3]]
        a.set_scaled_positions(q)
    actual = continuous_reference_coordinates(images, ref)
    assert np.allclose(actual["displacements_A"], expected["displacements_A"])
    assert np.array_equal(actual["fixed_integer_lattice_shifts_by_atom"], [[-1, 2, 0], [1, 0, -3]])


@pytest.mark.parametrize("fault", ["wrap", "half", "order", "pbc"])
def test_invalid_gauge_is_rejected_not_repaired(fault):
    ref, images = toy_chain()
    if fault == "wrap":
        images[-1].wrap()
    elif fault == "half":
        q = ref.get_scaled_positions()
        q[0, 1] += .5
        ref.set_scaled_positions(q)
    elif fault == "order":
        ref = ref[[1, 0]]
    else:
        ref.pbc[2] = False
    with pytest.raises(ValueError):
        continuous_reference_coordinates(images, ref)


def test_oblique_finite_strain_and_non_affine_displacement_are_separate():
    ref, _ = toy_chain()
    ref.set_cell([[5, .3, .2], [.1, 6, .4], [.5, .2, 7]], scale_atoms=True)
    f = np.array([[1.01, .03, -.01], [-.02, .99, .04], [.02, 0, 1.02]])
    deformed = ref.copy()
    deformed.set_cell(ref.cell.array @ f.T, scale_atoms=True)
    result = continuous_reference_coordinates([ref, deformed], ref)
    assert np.allclose(result["deformation_gradients"][1], f)
    assert np.allclose(result["green_strains"][1], .5 * (f.T @ f - np.eye(3)))
    assert np.allclose(result["displacements_A"], 0, atol=1e-14)


def test_mass_metric_translation_removal_and_parseval():
    masses = np.array([2., 5.])
    u = np.random.default_rng(2).normal(size=(3, 2, 3))
    basis = np.eye(6)
    result = project_reference_basis(u, basis, metric_weights=masses)
    translated = project_reference_basis(u + [1.3, -2., .7], basis, metric_weights=masses)
    assert np.allclose(result["amplitudes"], translated["amplitudes"])
    assert np.allclose(np.sum(result["amplitudes"]**2, axis=1), result["total_norm"]**2)
    assert np.allclose(result["residual_norm"], 0)
    assert np.allclose(result["captured_squared_norm_fraction"], 1)


def test_degenerate_subspace_weight_not_individual_coordinate_is_invariant():
    u = np.random.default_rng(3).normal(size=(2, 2, 3))
    basis = np.eye(6)[:, :2]
    rotation = np.array([[.6, -.8], [.8, .6]])
    first = project_reference_basis(u, basis)
    second = project_reference_basis(u, basis @ rotation)
    assert not np.allclose(first["amplitudes"], second["amplitudes"])
    assert np.allclose(np.sum(first["amplitudes"]**2, axis=1), np.sum(second["amplitudes"]**2, axis=1))
    assert np.allclose(first["residual_norm"], second["residual_norm"])


def test_zero_displacement_is_explicit_and_invalid_metric_fails():
    zero = project_reference_basis(np.zeros((1, 2, 3)), np.eye(6))
    assert zero["zero_displacement"].all() and zero["captured_squared_norm_fraction"][0] == 0
    for weights in ([1, 0], [1], [1, float("nan")]):
        with pytest.raises(ValueError, match="metric_weights"):
            project_reference_basis(np.zeros((1, 2, 3)), np.eye(6), metric_weights=weights)
    with pytest.raises(ValueError, match="orthonormal"):
        project_reference_basis(np.zeros((1, 2, 3)), np.ones((6, 2)))
