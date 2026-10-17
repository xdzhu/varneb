"""Synthetic target/selection checks, not DFT labels or material admission."""
from copy import deepcopy
from pathlib import Path
import importlib.util
import numpy as np
import pytest

P = Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/G3_P_branch_independent_E263_20261018/target_policy.py'
s = importlib.util.spec_from_file_location('E263_training_only_targets', P)
module = importlib.util.module_from_spec(s); s.loader.exec_module(module)


def fixture():
    axes = np.eye(36)[:, :2]
    rows = []
    for i in range(10):
        q = list(module.TRAINING_Q[i//2]); x = np.zeros(36); x[:2] = q; x[2] = .03
        row = dict(id=i, Q_id=i//2, q_A=q, start_sign=1 if i % 2 == 0 else -1,
                   status='attempt_budget_exhausted_no_retry', label=None)
        if i >= 2:
            row.update(status='gradient_candidate_stability_unchecked', label=dict(
                coordinates_internal_A=x.tolist(), energy_eV_cell=-10.+i,
                released_gradient_norm_eV_A=.001, zeta_gradient_eV_A=1e-5, gradient_candidate=True))
        rows.append(row)
    return rows, axes


def test_exact_eight_original_targets_no_training_grid_overlap():
    values = module.targets()
    assert [r['id'] for r in values] == list(range(1, 9))
    assert len({tuple(r['q_A']) for r in values}) == 8
    assert all(tuple(r['q_A']) not in module.TRAINING_Q for r in values)


def test_training_parent_selection_uses_distance_and_id_not_energy():
    rows, axes = fixture()
    first = module.training_only_starts(rows, axes, lambda x: np.linalg.norm(x) <= .5)
    assert first[0]['training_parent_id'] == 2
    changed = deepcopy(rows)
    for row in changed[2:]: row['label']['energy_eV_cell'] *= -10000
    assert module.training_only_starts(changed, axes, lambda x: True) == first
    for row in first:
        np.testing.assert_allclose(axes.T@row['starting_internal_coordinates_A'], row['q_A'], atol=2e-12, rtol=0)


def test_all_outcomes_required_and_directional_stop_not_replaced_by_norm_only():
    rows, axes = fixture()
    with pytest.raises(ValueError, match='all ten'): module.training_only_starts(rows[:-1], axes, lambda x: True)
    rows[2]['label']['zeta_gradient_eV_A'] = .001
    with pytest.raises(ValueError, match='both original'): module.training_only_starts(rows, axes, lambda x: True)


def test_guard_rejection_is_explicit_and_never_repairs_geometry():
    rows, axes = fixture()
    result = module.training_only_starts(rows, axes, lambda x: False)
    assert len(result) == 8
    assert all(r['training_parent_id'] is None and r['starting_internal_coordinates_A'] is None for r in result)
    assert all(len(r['rejected_parent_ids']) == 8 for r in result)
