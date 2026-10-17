"""Synthetic geometry decomposition plus read-only fixed native evidence."""
from pathlib import Path
import importlib.util
import numpy as np
import pytest

PATH = Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/G3_preserving_branches_E260_20261017/diagnose_prefix_E262.py'
spec = importlib.util.spec_from_file_location('E262_prefix_diagnosis_test', PATH)
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)


def test_frozen_component_excluded_from_radial_and_tangent_gradient():
    free = np.eye(3)[:, 1:]
    value = module.radial_gradient(np.array([.3, .4, 0.]), np.array([99., -2., 3.]), free)
    assert value['frozen_radius_A'] == .3
    assert value['released_radius_A'] == .4
    assert value['radial_gradient_eV_A'] == -2.
    assert value['tangent_gradient_norm_eV_A'] == 3.
    assert value['released_gradient_norm_eV_A'] == pytest.approx(np.sqrt(13.))


def test_zero_released_radius_retains_full_gradient_without_division():
    value = module.radial_gradient(np.array([.3, 0., 0.]), np.array([99., -2., 3.]), np.eye(3)[:, 1:])
    assert value['radial_gradient_eV_A'] == 0.
    assert value['tangent_gradient_norm_eV_A'] == pytest.approx(np.sqrt(13.))


def test_all_first_four_native_calls_accounted_not_only_successful_candidates():
    result = module.diagnose()
    assert result['native_calls'] == 36
    assert [r['attempted_native_SCFs'] for r in result['rows']] == [12, 12, 6, 6]
    assert all(r['total_radius_margin_A'] >= 0 for r in result['rows'])
    assert all(r['energy_rejected_native_SCFs'] == 0 for r in result['rows'])
    assert result['new_DFT_or_scheduler_writes'] == 0
    assert not result['complete_cohort_or_stability_or_predictive_claim']
