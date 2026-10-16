"""Real-data soft-direction diagnosis, not an instability or TS certificate."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / 'benchmarks/hfo2_channels/20261008/G3_released_complete_E254_20261016'


@pytest.fixture(scope='module')
def evidence():
    spec = importlib.util.spec_from_file_location('released_soft_diagnosis', CASE / 'soft_direction_analysis.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    result = module.analyze()
    stored = json.loads((CASE / 'soft_direction.json').read_bytes())
    assert result.keys() == stored.keys()
    for key in result:
        if key == 'coordinate_weights':
            for field in result[key]:
                np.testing.assert_allclose(result[key][field], stored[key][field], atol=2e-12, rtol=0)
        elif isinstance(result[key], (dict, str, bool)):
            assert result[key] == stored[key]
        else:
            np.testing.assert_allclose(result[key], stored[key], atol=2e-12, rtol=0)
    return result, module


def test_actual_two_steps_track_the_same_negative_direction(evidence):
    r, _ = evidence
    np.testing.assert_allclose(r['minimum_eigenvalues_eV_A2'],
                               [-.022077212863866463, -.024089732924225806], atol=2e-12, rtol=0)
    assert r['two_step_soft_squared_overlap'] == pytest.approx(.9999978887323483, abs=2e-12, rel=0)
    v = np.asarray(r['soft_direction_full40'])
    assert v @ v == pytest.approx(1., abs=2e-12, rel=0)
    pivot = np.flatnonzero(abs(v) >= np.max(abs(v))-1e-12)[0]
    assert v[pivot] > 0 and v[-1] == 0
    np.testing.assert_allclose(r['original_Q_overlap'], 0., atol=2e-12, rtol=0)


def test_species_weights_are_coordinate_weights_not_energy_partition(evidence):
    r, _ = evidence
    v = np.asarray(r['soft_direction_full40'])
    w = r['coordinate_weights']
    assert w['O'] == pytest.approx(np.sum(v[12:36]**2), abs=2e-12, rel=0)
    assert w['Hf'] == pytest.approx(np.sum(v[:12]**2), abs=2e-12, rel=0)
    assert w['open_cell'] == pytest.approx(np.sum(v[36:39]**2), abs=2e-12, rel=0)
    assert w['O'] > .9999 and w['open_cell'] < .00001
    assert sum(w[k] for k in ('Hf', 'O', 'open_cell')) == pytest.approx(1., abs=2e-12, rel=0)
    assert not r['coordinate_weights_are_energy_fractions_or_phonon_labels']


def test_negative_curvature_remains_in_the_fixed_cell_restriction(evidence):
    r, _ = evidence
    assert r['atomic_restriction_dimension'] == 31
    np.testing.assert_allclose(r['atomic_restriction_minima_eV_A2'],
                               [-.022058208479049088, -.02407000861844597], atol=2e-12, rtol=0)
    assert all(x < 0 for x in r['atomic_restriction_minima_eV_A2'])
    assert abs(r['soft_gradient_eV_A']) < 4e-7
    assert r['full_released_gradient_norm_eV_A'] < .03


def test_formal_promotion_is_not_a_new_stable_surface_or_full_index(evidence):
    r, _ = evidence
    assert r['after_formal_third_coordinate_remaining_dimension'] == 33
    np.testing.assert_allclose(r['remaining_raw_minima_eV_A2'],
                               [.3848778969656874, .37509315973880214], atol=2e-12, rtol=0)
    assert not r['full_new_retained_Hessian_stationarity_or_TS_certified']
    assert not r['independent_accuracy_or_stable_nonlinear_branch_certified']
    assert not r['original_independent_sign_admission_passed']


def test_analysis_does_not_reserve_new_points_or_unlock_old_tests(evidence):
    r, _ = evidence
    assert r['new_DFT_inputs_or_scheduler_writes'] == r['new_training_reservations'] == 0
    assert r['original_independent_slots_computed'] == 0
    assert not r['physical_plus005_generated_or_read']
    assert not r['radius_interval_is_complete_geometry_or_DFT_admission']
