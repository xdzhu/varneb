"""Terminal handle gating and unavailable arms are not scientific results."""
import importlib.util
import json
from pathlib import Path

import pytest

path = Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/matched_pattern_seeded_E271_20261018/collect_outputs.py'
spec = importlib.util.spec_from_file_location('release_collection', path)
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)


@pytest.mark.parametrize('state,code', [('RUNNING', '0:0'), ('TIMEOUT', '0:0'), ('FAILED', '1:0'), ('COMPLETED', '1:0')])
def test_only_original_completed_handle_can_be_collected(state, code):
    with pytest.raises(ValueError):
        module.require_terminal(f'28947338|{state}|{code}|32|hfacnormal01|iai806|1|3600|start|end|', '28947338')


def test_complete_handle_accounting_is_kept_exactly():
    row = module.require_terminal('28947338|COMPLETED|0:0|32|hfacnormal01|iai806|1|3600|start|end|', '28947338')
    assert row[7] == '3600'
    with pytest.raises(ValueError): module.require_terminal('|'.join(row), 'other-job')


@pytest.mark.parametrize('status', ['attempt_budget_exhausted', 'released_gradient_target_not_reached', 'not_attempted_atomic_candidate_unavailable'])
def test_missing_candidate_energy_and_residual_stay_unavailable(status):
    row = dict(treatment='A', status=status, candidate=None, new_attempts=12)
    r = module.candidate_summary(None, None, None, 'A', row, [])
    assert r['candidate_energy_eV_cell'] is None and r['released_gradient_eV_A'] is None
    assert r['energy_change_from_F_meV_fu'] is None


def test_actual_completed_A_candidate_reproduces_energy_and_held_variables():
    case = path.parent.with_name('matched_pattern_release_E270_20261018')
    data = json.loads((case/'prepared/seeds/1/seed.json').read_bytes())
    row = json.loads((case/'midpoint_A_prefix/treatment_A.json').read_bytes())
    native_path = next((case/'midpoint_A_prefix/evaluations').glob('*/result.json'))
    native = json.loads(native_path.read_bytes()); native['export_source'] = str(native_path)
    data['frozen_energy_eV_cell'] = -9782.60723549071
    r = module.candidate_summary(None, None, data, 'A', row, [native])
    assert r['energy_change_from_F_meV_fu'] == pytest.approx(-.0088583010438)
    assert r['released_gradient_eV_A'] == pytest.approx(.023923315666125)
    assert r['physical_retained_defect_A'] < 1e-15 and r['cell_coordinate_norm_A'] == 0
    native['gradient'][0] += .1
    with pytest.raises(ValueError): module.candidate_summary(None, None, data, 'A', row, [native])
