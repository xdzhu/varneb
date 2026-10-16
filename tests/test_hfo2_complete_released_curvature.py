"""Replay both actual E181 stencils, including the refused preserving sign gate."""
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CASES = ROOT / 'benchmarks/hfo2_channels/20261008'
CASE = CASES / 'G3_released_complete_E254_20261016'
T_CASE = CASES / 'G3_T_released_complete_E252_20261015'
CONFIG = CASES / 'G3_released_curvature_E181_20261015/config.json'
READER = CASES / 'G3_released_audit_E190_20261015/audit.py'
ATOL = 2e-12


@pytest.fixture(scope='module')
def evidence():
    pins = {
        CASE / 'native_audit.json': 'a389154749a07109a3a91499343b86f06a2c185fd310611570637d71193e3381',
        CASE / 'preserving_curvature.json': '588e6527936e2ceb8e94d8785293ae39fd36f561c0583004d0709a9fe1a561f2',
        CASE / 'completed.json': '70d3e665c86ed3c73fb2f560ec5ba6902b594b40db47bcfd763e64fa21cfa9bb',
        T_CASE / 'T_curvature.json': '66f7cdc9648afbfc3b61fc63a3835ec5868350c61aa921e56feea578a8498767',
        T_CASE / 'native_audit.json': '7f2dc871917734f3c2f3826c58d70239ef02545a29363b8b9fac8d7902bb68d4',
        CONFIG: '71a2e1b8a46b2ef903e4023d329109c5e74975aaf5bab54f4d4638d15c1add39',
        READER: '51d0452968394ff03d26efface7705ee31c7679215cb2204f6c69e13a5391d4f',
    }
    for path, digest in pins.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
    audit = json.loads((CASE / 'native_audit.json').read_bytes())
    summaries = {
        'T': json.loads((T_CASE / 'T_curvature.json').read_bytes()),
        'preserving': json.loads((CASE / 'preserving_curvature.json').read_bytes()),
    }
    config = json.loads(CONFIG.read_bytes())
    previous = json.loads((T_CASE / 'native_audit.json').read_bytes())
    terminal = json.loads((CASE / 'completed.json').read_bytes())
    spec = importlib.util.spec_from_file_location('complete_released_reader', READER)
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    return audit, summaries, config, previous, terminal, reader


def channel_data(evidence, channel):
    audit, summaries, config, _, _, reader = evidence
    offset = 0 if channel == 'T' else 136
    sources = audit['native_sources'][offset:offset + 136]
    summary = summaries[channel]
    records = [dict(gradient=g, enthalpy_eV=e, validation=s['validation'])
               for g, e, s in zip(summary['full_gradients_eV_A'],
                                  summary['energies_eV_cell'], sources)]
    row = next(r for r in audit['channels'] if r['channel'] == channel)
    return summary, config['channels'][channel], row, records, sources, reader


def test_normal_whole_job_requires_both_complete_native_cohorts(evidence):
    audit, _, _, _, terminal, reader = evidence
    assert audit['mode'] == 'complete'
    assert audit['status'] == 'complete_native_released_stencils_not_stability_certificate'
    assert audit['native_new_stencil_SCFs_reparsed'] == len(audit['native_sources']) == 272
    assert audit['job_id'] == terminal['job_id'] == '28868827'
    account = audit['accounting']
    assert account == reader.accounting(account['raw'])
    assert account['normal_completion'] and account['state'] == 'COMPLETED'
    assert (account['exit'], account['CPUs'], account['nodes']) == ('0:0', 32, 1)
    assert account['node_list'] == 'node9' and account['elapsed_seconds'] == 31120
    assert terminal['status'] == 'complete_both_released_stencils'
    assert terminal['new_native_SCFs'] == 272 and terminal['new_center_SCFs'] == 0
    for row in audit['channels']:
        assert row['completed_native_records'] == row['reserved_attempts'] == 136
        assert row['pending'] is None and row['failures'] == []


def test_complete_read_preserves_every_previously_audited_source(evidence):
    audit, _, _, previous, _, _ = evidence
    assert audit['native_sources'][:138] == previous['native_sources']
    assert len({s['result_path'] for s in audit['native_sources']}) == 272
    assert len({s['result_sha256'] for s in audit['native_sources']}) == 272
    fixed = {k: v for k, v in previous['native_sources'][0]['validation']['input_sha256'].items()
             if k != 'STRU'}
    assert len(fixed) == 6
    for source in audit['native_sources']:
        v = source['validation']
        assert v['raw_SCFs'] == 1 and v['genuine_MPI_ranks'] == 32
        assert v['elapsed_seconds'] > 0
        assert {k: h for k, h in v['input_sha256'].items() if k != 'STRU'} == fixed


@pytest.mark.parametrize('channel', ['T', 'preserving'])
def test_exact_original_136_probe_order_and_literal_anchors(evidence, channel):
    summary, spec, _, _, sources, _ = channel_data(evidence, channel)
    directions = np.asarray(spec['literal_directions_40_by_34'])
    anchor = np.asarray(spec['literal_native_anchor_40_A'])
    assert directions.shape == (40, 34) and anchor.shape == (40,)
    expected = [(h, i, sign) for h in (.01, .02) for i in range(34) for sign in (1, -1)]
    assert len(spec['probes']) == len(expected) == 136
    for index, (probe, (h, i, sign), source) in enumerate(zip(spec['probes'], expected, sources)):
        assert (probe['id'], probe['step_A'], probe['direction'], probe['sign']) == (index, h, i, sign)
        np.testing.assert_allclose(np.asarray(probe['coordinates_full_A']) - anchor,
                                   sign * h * directions[:, i], atol=ATOL, rtol=0)
        assert f'/released_curvature/{channel}/evaluations/eval-{index:06d}-' in source['result_path']
    assert summary['new_native_SCFs'] == 136 and summary['new_center_SCFs'] == 0


@pytest.mark.parametrize('channel', ['T', 'preserving'])
def test_full_native_gradients_reconstruct_both_unrepaired_actions(evidence, channel):
    summary, spec, row, records, _, reader = channel_data(evidence, channel)
    g = np.asarray(summary['full_gradients_eV_A'])
    directions = np.asarray(spec['literal_directions_40_by_34'])
    assert g.shape == (136, 40) and np.isfinite(g).all()
    replayed = reader.reconstruct(spec, records)
    for j, h in enumerate((.01, .02)):
        action = np.column_stack([(g[j*68+2*i] - g[j*68+2*i+1]) / (2*h) for i in range(34)])
        raw = directions.T @ action
        for stored in (summary['blocks'][j], row['blocks'][j], replayed[j]):
            np.testing.assert_allclose(stored['full_hessian_action_eV_A2'], action, atol=ATOL, rtol=0)
            np.testing.assert_allclose(stored['raw_projected_eV_A2'], raw, atol=ATOL, rtol=0)
        eig = np.linalg.eigvalsh((raw + raw.T)/2)
        block = row['blocks'][j]
        np.testing.assert_allclose(block['symmetric_projected_eigenvalues_eV_A2'], eig, atol=ATOL, rtol=0)
        assert block['raw_negative_count'] == np.count_nonzero(eig < 0)
        assert block['reciprocity_relative_defect'] == pytest.approx(
            np.linalg.norm(raw - raw.T)/np.linalg.norm(raw), abs=ATOL, rel=0)
        change = np.linalg.norm(raw - (raw + raw.T)/2, ord=2)
        assert change > 0
        assert block['symmetrization_operator_change_eV_A2'] == pytest.approx(change, abs=ATOL, rel=0)
    reader.check_summary(summary, channel, spec, records, replayed)


@pytest.mark.parametrize('channel,minima,negative,difference', [
    ('T', [.1640645284775998, .1575959364082007], 0, .054468926313563405),
    ('preserving', [-.022077212863866463, -.024089732924225806], 1, .31310777342198165),
])
def test_both_steps_keep_observed_signs_without_selection(evidence, channel, minima, negative, difference):
    summary, spec, row, _, _, _ = channel_data(evidence, channel)
    for j, block in enumerate(row['blocks']):
        assert min(block['symmetric_projected_eigenvalues_eV_A2']) == pytest.approx(minima[j], abs=ATOL, rel=0)
        assert block['raw_negative_count'] == negative
    raw = [np.asarray(b['raw_projected_eV_A2']) for b in row['blocks']]
    assert np.linalg.norm(raw[0]-raw[1], ord=2) == pytest.approx(difference, abs=ATOL, rel=0)
    assert summary['two_step_raw_block_operator_difference_eV_A2'] == pytest.approx(difference, abs=ATOL, rel=0)
    assert summary['base_released_gradient_norm_eV_A'] == spec['base_released_gradient_norm_eV_A']
    assert 0 < summary['base_released_gradient_norm_eV_A'] < .03


def test_negative_preserving_spectrum_cannot_be_clipped_in_summary(evidence):
    summary, spec, _, records, _, reader = channel_data(evidence, 'preserving')
    blocks = reader.reconstruct(spec, records)
    for j in (0, 1):
        broken = deepcopy(summary)
        broken['blocks'][j]['symmetric_projected_eigenvalues_eV_A2'][0] = 0.
        with pytest.raises(ValueError):
            reader.check_summary(broken, 'preserving', spec, records, blocks)


def test_failed_necessary_sign_condition_is_not_repaired_by_gradient_target(evidence):
    audit, summaries, _, _, _, _ = evidence
    passed = {r['channel']: all(min(b['symmetric_projected_eigenvalues_eV_A2']) > 0
                               for b in r['blocks']) for r in audit['channels']}
    assert passed == {'T': True, 'preserving': False}
    assert not all(passed.values())
    assert not audit['stable_branch_TS_or_independent_prediction_certified']
    assert not audit['physical_reserved_geometry_generated_or_read']
    assert audit['new_DFT_or_scheduler_writes'] == 0
    assert all(not s['physical_stability_TS_or_independent_prediction_certified'] for s in summaries.values())


def test_completed_cost_separates_two_cohorts_allocation_and_prior_preparation(evidence):
    audit, summaries, _, _, terminal, _ = evidence
    proxies = []
    for channel in ('T', 'preserving'):
        _, _, row, _, sources, _ = channel_data(evidence, channel)
        seconds = sum(s['validation']['elapsed_seconds'] for s in sources)
        assert summaries[channel]['sum_native_call_seconds'] == seconds
        proxy = seconds*32/3600
        assert summaries[channel]['call_time_proxy_core_hours'] == row['completed_call_proxy_core_hours'] == proxy
        proxies.append(proxy)
    assert sum(proxies) == pytest.approx(275.26456540390644, abs=ATOL, rel=0)
    assert audit['accounting']['allocation_core_hours'] == pytest.approx(31120*32/3600, abs=ATOL, rel=0)
    assert terminal['allocation_core_hours'] is None  # original producer did not know final sacct
    assert audit['all_other_curvature_training_preparation_costs_additional']
    assert terminal['all_other_training_and_preparation_costs_additional']
