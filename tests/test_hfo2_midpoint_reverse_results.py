"""Postprocessing contracts and arithmetic, not synthetic material evidence."""
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'benchmarks/hfo2_channels/20261008/midpoint_reverse_E269_20261018'


def module(name):
    spec = importlib.util.spec_from_file_location('E269_'+name, CASE/f'{name}.py')
    result = importlib.util.module_from_spec(spec); spec.loader.exec_module(result)
    return result


def test_paired_difference_uses_both_paths_and_keeps_old_labels():
    compare = module('compare_paths')
    # Deliberately simple arithmetic fixtures, not DFT or a predicted outcome.
    rows = [dict(channel='preserving', from_zero_peak_meV_fu=100., from_upper_peak_meV_fu=110.),
            dict(channel='T', from_zero_peak_meV_fu=90., from_upper_peak_meV_fu=95.)]
    result = compare.paired_summary(rows)
    assert result['from_zero_T_minus_preserving_meV_fu'] == -10.
    assert result['from_upper_T_minus_preserving_meV_fu'] == -15.
    assert result['seed_gap_difference_meV_fu'] == -5.
    assert result['per_channel_upper_minus_zero_meV_fu'] == {'preserving': 10., 'T': 5.}
    assert result['physical_barrier_error_bound'] is None
    assert not result['original_prediction_scores_replaced']
    assert not result['complete_four_channel_ordering_claimed']


@pytest.mark.parametrize('kind', ['missing', 'duplicate', 'nonfinite', 'negative'])
def test_incomplete_or_invalid_pair_is_not_scored(kind):
    compare = module('compare_paths')
    rows = [dict(channel='preserving', from_zero_peak_meV_fu=100., from_upper_peak_meV_fu=110.),
            dict(channel='T', from_zero_peak_meV_fu=90., from_upper_peak_meV_fu=95.)]
    if kind == 'missing': rows.pop()
    if kind == 'duplicate': rows[1]['channel'] = 'preserving'
    if kind == 'nonfinite': rows[0]['from_upper_peak_meV_fu'] = float('nan')
    if kind == 'negative': rows[0]['from_upper_peak_meV_fu'] = -1.
    with pytest.raises(ValueError): compare.paired_summary(rows)


@pytest.mark.parametrize('state', ['RUNNING', 'PENDING', 'TIMEOUT', 'FAILED'])
def test_live_or_failed_allocation_never_becomes_a_completed_result(state):
    collect = module('collect_results')
    good = ['28946180', 'COMPLETED', '0:0', '32', 'hfacnormal01', 'iai806', '1',
            '1200', '2026-10-18T19:27:47', '2026-10-18T19:47:47']
    collect.require_completed(good)
    bad = good.copy(); bad[1] = state
    with pytest.raises(ValueError): collect.require_completed(bad)


def test_accounting_requires_original_handle_not_a_nearby_job(monkeypatch):
    collect = module('collect_results')
    monkeypatch.setattr(collect.subprocess, 'check_output', lambda *a, **k:
        '999|COMPLETED|0:0|32|hfacnormal01|iai806|1|1200|start|end|\n')
    with pytest.raises(ValueError): collect.accounting('28946180')


def test_actual_initial_snapshots_are_not_promoted_to_terminal_paths():
    compare = module('compare_paths')
    exported = CASE/'initial_complete_chains'
    with pytest.raises(ValueError, match='initial snapshots'):
        compare.compare(ROOT, exported)
    result = compare.compare(ROOT, exported, initial=True)
    assert result['status'] == 'initial_seed_comparison_not_terminal_paths'
    assert result['upper_stage'] == 'initial_evaluation'
    assert not any(r['upper_terminal_path'] or r['upper_ordinary_converged'] for r in result['outcomes'])
    assert result['paired']['per_channel_upper_minus_zero_meV_fu'] == pytest.approx({
        'preserving': 40.3582262370, 'T': .0851851300})
    assert result['outcomes'][0]['from_upper_Qx_A'] == pytest.approx(-1.0393884673)
    assert result['outcomes'][0]['from_upper_open_y_percent'] == pytest.approx(.6547077522)


def test_actual_initial_figure_keeps_all_samples_and_labels():
    compare = module('compare_paths')
    data = compare.compare(ROOT, CASE/'initial_complete_chains', initial=True)
    fig = compare.figure(data)
    try:
        assert len(fig.axes) == 6 and len(fig.legends) == 1
        assert any('initial' in t.get_text() for t in fig.legends[0].get_texts())
        for i, ax in enumerate(fig.axes):
            assert all(len(line.get_xdata()) == 9 for line in ax.lines)
            assert ax.spines['top'].get_visible() and ax.spines['right'].get_visible()
            assert ax.texts[0].get_text() == f'({chr(97+i)})'
            assert ax.texts[0].get_fontweight() == 'normal'
    finally:
        import matplotlib.pyplot as plt
        plt.close(fig)


def test_actual_completed_T_can_be_processed_without_imputing_preserving(tmp_path):
    compare = module('compare_paths')
    exported = CASE/'completed_channels'
    result = compare.compare(ROOT, exported, channels=('T',))
    assert result['status'] == 'one_completed_midpoint_channel_compared'
    assert not result['complete_selected_pair_available'] and result['paired'] is None
    assert [r['channel'] for r in result['outcomes']] == ['T']
    row = result['outcomes'][0]
    assert row['upper_ordinary_converged'] and row['upper_terminal_path']
    assert row['from_upper_peak_meV_fu'] == pytest.approx(90.28624292841414)
    assert row['from_zero_peak_meV_fu'] == pytest.approx(90.41817833576715)
    assert row['upper_SCFs'] == 21 and row['from_upper_updates'] == 2
    with pytest.raises(ValueError): compare.figure(result)
    with pytest.raises(ValueError): compare.write_result(result, tmp_path/'refused_figure', draw=True)
    assert not (tmp_path/'refused_figure').exists()
    with pytest.raises(FileNotFoundError): compare.compare(ROOT, exported)


def test_every_completed_T_native_call_is_retained_and_replays():
    import json
    import numpy as np
    from scripts.audit_hfo2_static_replica import audited_results, sha256
    from examples.hfo2_fixed_input_factory import CONTRACT
    folder = CASE/'completed_channels/T'
    audit = json.loads((folder/'audit.json').read_bytes())
    assert len(audit['calls']) == audit['new_SCF_calls'] == 21
    assert audit['job_id'] == '28946181' and audit['ordinary_converged']
    assert audit['scheduler_row'][1:7] == ['COMPLETED', '0:0', '32', 'hfacnormal01', 'iai806', '1']
    for i, record in enumerate(audit['calls']):
        raw = folder/'all_calls'/f'{i:04d}'
        assert sha256(raw/'call_audit.json') == record['call_audit_sha256']
        call = json.loads((raw/'call_audit.json').read_bytes())
        assert all(call['input_sha256'][n] == h for n, h in CONTRACT.items())
        assert all(sha256(raw/n) == call['input_sha256'][n] for n in ('INPUT', 'KPT', 'STRU'))
        assert sha256(raw/'OUT.ABACUS/running_scf.log') == call['raw_log_sha256']
        parsed = audited_results(raw)
        for key in parsed:
            np.testing.assert_allclose(parsed[key], call['results'][key], atol=1e-12, rtol=0)
