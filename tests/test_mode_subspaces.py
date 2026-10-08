import numpy as np
import pytest

from vcneb.mode_subspaces import real_mode_subspace, mass_translation_overlaps


def test_complex_phases_and_unitary_mixing_preserve_real_projector():
    rng = np.random.default_rng(14)
    q = np.linalg.qr(rng.normal(size=(12, 4)))[0]
    complex_rotation = np.linalg.qr(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))[0]
    values = q @ complex_rotation
    before = values.copy()
    result = real_mode_subspace(values, expected_rank=4)
    assert np.array_equal(values, before)
    assert np.allclose(result.basis_columns @ result.basis_columns.T, q @ q.T, atol=1e-14)
    assert result.relative_column_projection_errors.max() < 1e-14
    assert result.orthogonality_max_defect < 1e-14
    assert not result.basis_columns.flags.writeable


def test_printed_noise_is_reported_not_claimed_as_eigenpairs():
    values = np.eye(9, 3, dtype=complex)
    values[7, 0] = 2e-6j
    result = real_mode_subspace(values, expected_rank=3)
    assert result.real_singular_values[3] > 0
    assert result.relative_column_projection_errors.max() == pytest.approx(2e-6)
    with pytest.raises(ValueError, match="projection"):
        real_mode_subspace(values, expected_rank=3, max_relative_projection_error=1e-7)
    values[7, 0] = .01j
    with pytest.raises(ValueError, match="real rank"):
        real_mode_subspace(values, expected_rank=3)


def test_source_translation_admixture_is_preserved_not_ASR_repaired():
    masses = np.array([2., 7.])
    vector = np.array([1., 0, 0, -.3, 0, 0])[:, None]
    result = real_mode_subspace(vector, expected_rank=1)
    expected = vector / np.linalg.norm(vector)
    assert np.allclose(result.basis_columns @ result.basis_columns.T, expected @ expected.T)
    assert mass_translation_overlaps(result.basis_columns, masses)[0] > .1


@pytest.mark.parametrize("fault", ["nan", "empty", "zero", "wide", "rank", "bool", "fraction", "tolerance"])
def test_bad_or_ambiguous_inputs_fail(fault):
    values, rank, kwargs = np.eye(6, 2), 2, {}
    if fault == "nan": values[0, 0] = np.nan
    elif fault == "empty": values = np.empty((0, 2))
    elif fault == "zero": values[:, 1] = 0
    elif fault == "wide": values = np.eye(3, 7)
    elif fault == "rank": rank = 1
    elif fault == "bool": rank = True
    elif fault == "fraction": rank = 1.5
    elif fault == "tolerance": kwargs["relative_rank_tolerance"] = 0
    with pytest.raises(ValueError):
        real_mode_subspace(values, expected_rank=rank, **kwargs)


def test_nonfinite_projection_tolerance_and_incompatible_masses_fail():
    with pytest.raises(ValueError):
        real_mode_subspace(np.eye(6, 2), expected_rank=2, max_relative_projection_error=np.nan)
    for masses in ([2], [2, 0], [2, np.nan]):
        with pytest.raises(ValueError):
            mass_translation_overlaps(np.eye(6, 2), masses)


def test_independent_large_column_scales_do_not_change_projector():
    values = np.eye(8, 3) * [1e200, 1e-200, 3.]
    result = real_mode_subspace(values, expected_rank=3)
    assert np.allclose(result.basis_columns @ result.basis_columns.T, np.diag([1., 1., 1., 0, 0, 0, 0, 0]))
