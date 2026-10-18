"""One actual completed midpoint arm, with failed arms kept unavailable."""
import importlib.util
import json
from pathlib import Path
import shutil

import pytest

CASE = Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/matched_pattern_release_E270_20261018'
spec = importlib.util.spec_from_file_location('E270_arm_review', CASE/'review_arm.py')
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)


def test_actual_midpoint_atomic_candidate_native_replay():
    r = module.review(CASE, CASE/'midpoint_A_prefix', 1, 'A')
    assert r['strain'] == .005 and r['total_arm_attempts_reported'] == 5
    assert r['energy_change_meV_fu'] == pytest.approx(-.008858301007)
    assert r['released_gradient_norm_eV_A'] == pytest.approx(.023923315666125003)
    assert r['cell_coordinate_norm_A'] == 0 and r['physical_retained_defect_A'] < 1e-15
    assert r['local_native_SCFs_reparsed'] == 1 and r['new_DFT_calls'] == 0
    assert not r['full_orbital_pseudopotential_bytes_rechecked_locally']
    assert not r['energy_partition_or_stationary_TS_or_minimum_certified']
    assert not r['changed_path_label_or_forecast_score']


def test_unavailable_zero_strain_arms_have_no_imputed_energy():
    for arm in ('A', 'J'):
        r = module.review(CASE, CASE/'partial_outcomes_2139/0', 0, arm)
        assert r['energy_change_meV_fu'] is None and r['released_gradient_norm_eV_A'] is None
        assert r['local_native_SCFs_reparsed'] == 0


@pytest.mark.parametrize('fault', ['energy', 'force', 'coordinate'])
def test_unreplayed_or_false_candidate_refused(tmp_path, fault):
    destination = tmp_path/'export'
    shutil.copytree(CASE/'midpoint_A_prefix', destination)
    path = destination/'treatment_A.json'
    row = json.loads(path.read_bytes())
    if fault == 'energy': row['candidate']['energy_eV_cell'] += .001
    elif fault == 'force': row['candidate']['released_gradient_norm_eV_A'] = .001
    else: row['candidate']['coordinates'][0] += .001
    path.write_text(json.dumps(row), encoding='utf-8')
    with pytest.raises(ValueError): module.review(CASE, destination, 1, 'A')
