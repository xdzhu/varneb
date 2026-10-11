"""Independent virtual-work, gauge, frame and subspace tests; no DFT."""
import numpy as np
import pytest

from vcneb.continuous_projection import project_reference_atomic_gradient


def fixture():
    rng = np.random.default_rng(78)
    masses = np.array([4., 16., 25.])
    translation = np.zeros((9, 3))
    for k in range(3):
        translation[k::3, k] = np.sqrt(masses/masses.sum())
    basis = rng.normal(size=(9, 3))
    basis -= translation @ (translation.T @ basis)
    basis = np.linalg.qr(basis)[0]
    reference = np.array([[3., .2, .1], [0., 4., .3], [0., 0., 5.]])
    current = np.array([[3.2, -.1, .5], [.2, 3.7, .4], [-.3, .6, 4.8]])
    forces = rng.normal(size=(3, 3)); forces -= forces.mean(axis=0)
    return forces, current, reference, basis, masses, translation


def call(values, basis=None):
    forces, current, reference, original, masses, _ = values
    return project_reference_atomic_gradient(forces, current, reference,
                    original if basis is None else basis, masses_amu=masses)


@pytest.mark.parametrize("same_cell", [False, True])
def test_single_direction_matches_independent_pair_energy_difference(same_cell):
    forces, current, reference, basis, masses, _ = fixture()
    if same_cell:
        current = reference.copy()
    positions = np.array([[.2, .4, .1], [1.1, -.3, .7], [-.5, .8, 1.4]])
    forces = -3*positions + positions.sum(axis=0)
    direction = (basis[:, 0]/np.repeat(np.sqrt(masses), 3)).reshape(3, 3)
    # Independently move fractional atoms in H_ref, then form current positions.
    delta = np.linalg.solve(reference.T, direction.T).T @ current
    def energy(p):
        return .5*sum(np.dot(p[i]-p[j], p[i]-p[j]) for i in range(3) for j in range(i))
    h = 1e-5
    derivative = (energy(positions+h*delta)-energy(positions-h*delta))/(2*h)
    report = project_reference_atomic_gradient(forces, current, reference, basis[:, :1], masses_amu=masses)
    assert report["retained_gradient_norm"] == pytest.approx(abs(derivative), abs=1e-9, rel=0)
    if not same_cell:
        # A displacement-vector pullback on F is not a dual-covector pullback.
        wrong = (-forces @ np.linalg.solve(current, reference)/np.sqrt(masses)[:, None]).ravel()
        assert abs(abs(basis[:, 0] @ wrong)-abs(derivative)) > .01


def test_norm_partition_and_internal_basis_rotation_invariance():
    values = fixture(); a = call(values)
    rotation = np.linalg.qr(np.random.default_rng(79).normal(size=(3, 3)))[0]
    b = call(values, values[3] @ rotation)
    for key in ("internal_gradient_norm", "retained_gradient_norm", "omitted_gradient_norm", "retained_squared_gradient_fraction"):
        assert a[key] == pytest.approx(b[key], abs=1e-13, rel=0)
    assert abs(a["norm_identity_defect"]) < 1e-14
    assert a["source_basis_rank"] == a["translation_free_basis_rank"] == 3
    assert a["gradient_unit"] == "eV/(sqrt(amu)*Angstrom)"


def test_oblique_cartesian_rotation_and_ordered_permutation_covariance():
    forces, current, reference, basis, masses, tr = fixture()
    rotation = np.linalg.qr(np.random.default_rng(80).normal(size=(3, 3)))[0]
    rotated_basis = np.einsum("nak,ba->nbk", basis.reshape(3, 3, 3), rotation).reshape(9, 3)
    rotated = (forces @ rotation.T, current @ rotation.T, reference @ rotation.T, rotated_basis, masses, tr)
    order = [2, 0, 1]
    permuted = (forces[order], current, reference, basis.reshape(3, 3, 3)[order].reshape(9, 3), masses[order], tr)
    a, b, c = call(fixture()), call(rotated), call(permuted)
    for key in ("internal_gradient_norm", "retained_gradient_norm", "omitted_gradient_norm", "retained_squared_gradient_fraction"):
        assert b[key] == pytest.approx(a[key], abs=1e-13, rel=0)
        assert c[key] == pytest.approx(a[key], abs=1e-13, rel=0)


def test_explicit_translation_admixture_is_restricted_without_source_repair():
    values = fixture()
    mixture = values[3][:, 0]*np.cos(.05) + values[5][:, 0]*np.sin(.05)
    before = mixture.copy()
    a, b = call(values, mixture[:, None]), call(values, values[3][:, :1])
    assert a["source_translation_overlap_Frobenius_norm"] == pytest.approx(np.sin(.05), abs=1e-14)
    assert a["restricted_basis_min_singular_value"] == pytest.approx(np.cos(.05), abs=1e-14)
    assert a["retained_gradient_norm"] == pytest.approx(b["retained_gradient_norm"], abs=1e-14)
    assert a["source_modes_or_forces_ASR_corrected"] is False
    assert np.array_equal(mixture, before)


def test_mass_scale_changes_gradient_units_not_fraction():
    values = fixture(); a = call(values)
    b = call((*values[:4], values[4]*9, values[5]))
    assert b["internal_gradient_norm"] == pytest.approx(a["internal_gradient_norm"]/3)
    assert b["retained_squared_gradient_fraction"] == pytest.approx(a["retained_squared_gradient_fraction"])


def test_zero_and_pure_rigid_gradient_have_no_false_coverage_fraction():
    forces, current, reference, basis, masses, tr = fixture()
    for force in (np.zeros_like(forces), masses[:, None]*np.array([.2, -.3, .1])):
        report = call((force, reference, reference, basis, masses, tr))
        assert report["zero_internal_gradient"] is True
        assert report["retained_squared_gradient_fraction"] is None
    assert report["rigid_translation_gradient_norm"] > 1


@pytest.mark.parametrize("bad", ["mass", "force", "left_cell", "singular_cell", "basis", "translation"])
def test_invalid_contract_or_lost_internal_rank_rejected(bad):
    forces, current, reference, basis, masses, tr = fixture()
    if bad == "mass": masses[0] = 0
    if bad == "force": forces[0, 0] = np.nan
    if bad == "left_cell": current[0] *= -1
    if bad == "singular_cell": current[0] = current[1]
    if bad == "basis": basis[:, 0] *= 2
    if bad == "translation": basis = tr[:, :1]
    with pytest.raises(ValueError):
        call((forces, current, reference, basis, masses, tr))
