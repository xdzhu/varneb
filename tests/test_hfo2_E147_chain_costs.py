"""Real cumulative production costs plus accounting/source fault fixtures."""
import json
from pathlib import Path

import pytest

from scripts import audit_hfo2_clamped_chain_costs as audit
from vcneb import ResponseEvaluationCost


@pytest.fixture(scope='module')
def result():
    return audit.audit()


def test_all_seven_actual_histories_count_every_segment_once(result):
    assert len(result['bands']) == 7 and result['cumulative_production_segments'] == 16
    assert result['total']['unique_source_DFT_evaluations'] == 1918
    assert result['total']['completed_evaluations'] == 1918
    assert result['total']['failed_evaluations'] == 0
    assert result['total']['recorded_transport_cpu_core_hours_sum'] == pytest.approx(1966.956847628798, abs=1e-9)
    assert result['allocation_core_hours'] == pytest.approx(1976.3644444444444, abs=1e-9)
    assert [r['fresh_interior_SCFs'] for r in result['bands']] == [308, 77, 42, 42, 742, 231, 476]
    assert [r['optimizer_updates'] for r in result['bands']] == [43, 10, 5, 5, 105, 32, 67]
    assert all(r['fresh_initial_interiors'] == 7 for r in result['bands'])
    for r in result['bands']:
        assert r['fresh_interior_SCFs'] == 7*(r['optimizer_updates']+1)
        assert r['SCF_core_hour_proxy'] <= r['allocation_core_hours']
    assert result['bands'][4]['jobs'] == ['28574708', '28661019', '28722810', '28723655', '28727382']
    p = result['bands'][-1]
    assert p['last_segment_fresh_SCFs'] == 119 and p['fresh_interior_SCFs'] == 476
    assert p['SCF_core_hour_proxy'] == pytest.approx(415.4959169591508)
    assert p['last_segment_SCF_core_hour_proxy'] == pytest.approx(99.55591450291375)


def test_complete_recorded_subset_is_not_whole_study_or_method_advantage(result):
    assert result['total']['cost_complete']
    assert not result['all_preparation_or_full_study_cost_complete']
    assert not result['all_native_assets_or_raw_EFS_reaudited_now']
    assert result['pinned_prior_native_audits_reused']
    assert not result['acceleration_TS_physical_error_or_prediction_accuracy_certified']
    assert not result['G3_or_holdout_access'] and not result['physical_parameters_changed']
    assert result['new_DFT_calls'] == result['scheduler_mutations'] == 0
    assert 'unfinished +1% reversing chain' in result['excluded_costs']
    assert set(result['native_accounting']) == set(result['sources']) == audit.JOBS
    assert '28852667' not in audit.JOBS and '28764704' not in audit.JOBS


def native_body():
    return json.loads((audit.CASE/'native_accounting.json').read_text())['stdout']


@pytest.mark.parametrize('index,value', [(1,'RUNNING'), (2,'1:0'), (3,'16'), (4,'2'),
    (5,'other'), (6,'other'), (7,'0'), (7,'-1'), (10,'0'), (9,'2026-10-10T14:12:03')])
def test_native_identity_and_allocation_mismatches_refused(index, value):
    lines = native_body().splitlines(); fields = lines[0].split('|'); fields[index] = value
    lines[0] = '|'.join(fields)
    with pytest.raises(ValueError): audit.parse_accounting('\n'.join(lines))


@pytest.mark.parametrize('mode', ['missing', 'duplicate', 'unknown', 'malformed'])
def test_missing_extra_duplicate_accounting_never_looks_complete(mode):
    lines = native_body().splitlines()
    if mode == 'missing': lines.pop()
    elif mode == 'duplicate': lines.append(lines[0])
    elif mode == 'unknown': lines.append(lines[0].replace('28574708', '99999999'))
    else: lines[0] += '|extra'
    with pytest.raises(ValueError): audit.parse_accounting('\n'.join(lines))


def record(name='scf_000000', **changes):
    r = dict(path=f'/actual/band/image_0001/{name}/call_audit.json',
             audit_sha256='a'*64, input_sha256=dict(audit.CONTRACT), elapsed_seconds=120.)
    r.update(changes)
    return r


def test_fresh_identical_geometries_are_not_deduplicated_by_log_or_structure():
    # Different actual calls remain distinct even if the diagnostic hashes
    # and geometry are identical. Cache accesses would keep the SAME source.
    a, b = record(), record('scf_000001')
    rows = audit.call_records([a, b])
    assert len(rows) == 2 and rows[0].raw_audit_sha256 == rows[1].raw_audit_sha256
    assert audit.segment_cost(rows, 3)['unique_source_DFT_evaluations'] == 2
    with pytest.raises(ValueError, match='duplicate'):
        audit.call_records([a, a])


@pytest.mark.parametrize('value', [None, True, -1., float('nan')])
def test_missing_or_invalid_call_time_is_not_free(value):
    with pytest.raises(ValueError): audit.segment_cost(audit.call_records([record(elapsed_seconds=value)]), 3)


def test_serial_call_sum_must_fit_native_allocation():
    with pytest.raises(ValueError): audit.segment_cost(audit.call_records([record()]), .5)
    r = record(); r['input_sha256']['INPUT'] = 'b'*64
    with pytest.raises(ValueError): audit.call_records([r])


def test_grid_cannot_hide_a_missing_image_with_extra_other_image_calls():
    rows = [ResponseEvaluationCost(f'hf:/band/image_{i:04d}/scf_{j:06d}', 'a'*64,
             'completed', 1., 32) for i in range(1,8) for j in range(2)]
    audit.validate_call_grid(rows, 2)
    with pytest.raises(ValueError): audit.validate_call_grid(rows[:-1], 2)
    with pytest.raises(ValueError): audit.validate_call_grid(rows+[rows[-1]], 2)


def test_original_receipt_changes_refused_before_interpreting_cost(tmp_path):
    path = tmp_path/'source.json'; path.write_text('{"synthetic":true}')
    with pytest.raises(ValueError): audit.pinned(path, 'a'*64)
