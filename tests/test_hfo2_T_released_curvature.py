"""Replay the first complete released stencil, not a whole-job or TS claim."""
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CASES = ROOT / 'benchmarks/hfo2_channels/20261008'
CASE = CASES / 'G3_T_released_complete_E252_20261015'
READER = CASES / 'G3_released_audit_E190_20261015/audit.py'
CONFIG = CASES / 'G3_released_curvature_E181_20261015/config.json'
PREVIOUS = CASES / 'G3_released_audit_E190_20261015/native_prefix.json'
ATOL = 2e-12


@pytest.fixture(scope='module')
def evidence():
    expected = {
        CASE / 'native_audit.json': '7f2dc871917734f3c2f3826c58d70239ef02545a29363b8b9fac8d7902bb68d4',
        CASE / 'T_curvature.json': '66f7cdc9648afbfc3b61fc63a3835ec5868350c61aa921e56feea578a8498767',
        CONFIG: '71a2e1b8a46b2ef903e4023d329109c5e74975aaf5bab54f4d4638d15c1add39',
        PREVIOUS: 'b80d7571848f151f326f2c19576b8f9d2446800fcea776808bcc06aa97255c52',
        READER: '51d0452968394ff03d26efface7705ee31c7679215cb2204f6c69e13a5391d4f',
    }
    for path, digest in expected.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
    audit = json.loads((CASE / 'native_audit.json').read_bytes())
    summary = json.loads((CASE / 'T_curvature.json').read_bytes())
    config = json.loads(CONFIG.read_bytes())
    previous = json.loads(PREVIOUS.read_bytes())
    spec = importlib.util.spec_from_file_location('actual_T_released_reader', READER)
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    return audit, summary, config['channels']['T'], previous, reader


def T_row(audit):
    return next(row for row in audit['channels'] if row['channel'] == 'T')


def records(audit, summary):
    sources = audit['native_sources'][:136]
    return [dict(gradient=g, enthalpy_eV=e, validation=s['validation'])
            for g, e, s in zip(summary['full_gradients_eV_A'],
                               summary['energies_eV_cell'], sources)]


def test_complete_T_is_still_a_prefix_of_the_running_two_center_job(evidence):
    audit, summary, _, _, _ = evidence
    assert audit['mode'] == 'prefix'
    assert audit['status'] == 'native_released_stencil_prefix_not_completion'
    assert audit['job_id'] == '28868827'
    assert audit['accounting']['state'] == 'RUNNING'
    assert audit['accounting']['CPUs'] == 32 and audit['accounting']['nodes'] == 1
    assert not audit['accounting']['normal_completion']
    assert T_row(audit)['completed_native_records'] == T_row(audit)['reserved_attempts'] == 136
    assert T_row(audit)['pending'] is None and T_row(audit)['failures'] == []
    other = next(row for row in audit['channels'] if row['channel'] == 'preserving')
    assert other['completed_native_records'] == 2 and other['reserved_attempts'] == 3
    assert other['blocks'] is None and other['pending'] == 'reserved_probe_2_output_pending'
    assert audit['native_new_stencil_SCFs_reparsed'] == len(audit['native_sources']) == 138
    assert summary['new_native_SCFs'] == 136 and summary['new_center_SCFs'] == 0


def test_all136_literal_probes_keep_original_two_steps_directions_and_order(evidence):
    _, _, spec, _, _ = evidence
    directions = np.asarray(spec['literal_directions_40_by_34'])
    anchor = np.asarray(spec['literal_native_anchor_40_A'])
    assert directions.shape == (40, 34) and anchor.shape == (40,)
    expected = [(h, i, sign) for h in (.01, .02) for i in range(34) for sign in (1, -1)]
    assert len(spec['probes']) == len(expected) == 136
    for index, (probe, (h, i, sign)) in enumerate(zip(spec['probes'], expected)):
        assert (probe['id'], probe['step_A'], probe['direction'], probe['sign']) == (index, h, i, sign)
        np.testing.assert_allclose(np.asarray(probe['coordinates_full_A']) - anchor,
                                   sign * h * directions[:, i], atol=ATOL, rtol=0)


def test_all_native_sources_keep_six_inputs_MPI_and_the_older13_identities(evidence):
    audit, _, _, previous, _ = evidence
    assert len(previous['native_sources']) == 13
    assert audit['native_sources'][:13] == previous['native_sources']
    sources = audit['native_sources'][:136]
    assert len({s['result_path'] for s in sources}) == 136
    assert len({s['result_sha256'] for s in sources}) == 136
    fixed = {k: v for k, v in previous['native_sources'][0]['validation']['input_sha256'].items()
             if k != 'STRU'}
    assert len(fixed) == 6
    for index, source in enumerate(sources):
        assert f'/released_curvature/T/evaluations/eval-{index:06d}-' in source['result_path']
        validation = source['validation']
        assert validation['genuine_MPI_ranks'] == 32 and validation['raw_SCFs'] == 1
        assert validation['elapsed_seconds'] > 0
        assert {k: v for k, v in validation['input_sha256'].items() if k != 'STRU'} == fixed


def test_all136_full_gradients_independently_reconstruct_both_raw_blocks(evidence):
    audit, summary, spec, _, reader = evidence
    gradients = np.asarray(summary['full_gradients_eV_A'])
    directions = np.asarray(spec['literal_directions_40_by_34'])
    assert gradients.shape == (136, 40) and np.isfinite(gradients).all()
    replayed = reader.reconstruct(spec, records(audit, summary))
    for j, h in enumerate((.01, .02)):
        action = np.column_stack([(gradients[j*68+2*i] - gradients[j*68+2*i+1]) / (2*h)
                                  for i in range(34)])
        raw = directions.T @ action
        for stored in (T_row(audit)['blocks'][j], summary['blocks'][j], replayed[j]):
            assert stored['step_A'] == h and stored['complete_pairs'] == 34
            np.testing.assert_allclose(stored['full_hessian_action_eV_A2'], action, atol=ATOL, rtol=0)
            np.testing.assert_allclose(stored['raw_projected_eV_A2'], raw, atol=ATOL, rtol=0)


def test_both_spectra_and_unrepaired_reciprocity_diagnostics_are_retained(evidence):
    audit, summary, _, _, _ = evidence
    minima = [0.1640645284775998, 0.1575959364082007]
    for j, block in enumerate(T_row(audit)['blocks']):
        raw = np.asarray(block['raw_projected_eV_A2'])
        symmetric = (raw + raw.T) / 2
        eigenvalues = np.linalg.eigvalsh(symmetric)
        np.testing.assert_allclose(block['symmetric_projected_eigenvalues_eV_A2'],
                                   eigenvalues, atol=ATOL, rtol=0)
        assert block['raw_negative_count'] == np.count_nonzero(eigenvalues < 0) == 0
        assert eigenvalues[0] == pytest.approx(minima[j], abs=ATOL, rel=0)
        assert block['reciprocity_relative_defect'] == pytest.approx(
            np.linalg.norm(raw-raw.T) / np.linalg.norm(raw), abs=ATOL, rel=0)
        change = np.linalg.norm(raw-symmetric, ord=2)
        assert change > 0  # raw nonsymmetric entries are not repaired away
        assert block['symmetrization_operator_change_eV_A2'] == pytest.approx(change, abs=ATOL, rel=0)
    raw = [np.asarray(b['raw_projected_eV_A2']) for b in T_row(audit)['blocks']]
    difference = np.linalg.norm(raw[0]-raw[1], ord=2)
    assert difference == pytest.approx(0.054468926313563405, abs=ATOL, rel=0)
    assert summary['two_step_raw_block_operator_difference_eV_A2'] == pytest.approx(difference, abs=ATOL, rel=0)


def test_original_summary_checker_accepts_actual_data_and_refuses_mutations(evidence):
    audit, summary, spec, _, reader = evidence
    rows = records(audit, summary)
    blocks = reader.reconstruct(spec, rows)
    reader.check_summary(summary, 'T', spec, rows, blocks)
    for fault in ('count', 'gradient', 'block', 'missing_step'):
        broken = deepcopy(summary)
        if fault == 'count': broken['new_native_SCFs'] = 135
        elif fault == 'gradient': broken['full_gradients_eV_A'][0][0] += .001
        elif fault == 'block': broken['blocks'][0]['raw_projected_eV_A2'][0][0] += .001
        else: broken['blocks'] = broken['blocks'][:1]
        with pytest.raises(ValueError):
            reader.check_summary(broken, 'T', spec, rows, blocks)


def test_actual_cost_keeps_parent_allocation_and_earlier_work_separate(evidence):
    audit, summary, _, _, _ = evidence
    seconds = sum(s['validation']['elapsed_seconds'] for s in audit['native_sources'][:136])
    assert summary['sum_native_call_seconds'] == seconds
    proxy = seconds * 32 / 3600
    assert proxy == pytest.approx(137.83936245814795, abs=ATOL, rel=0)
    assert summary['call_time_proxy_core_hours'] == T_row(audit)['completed_call_proxy_core_hours'] == proxy
    assert audit['all_other_curvature_training_preparation_costs_additional']
    assert audit['accounting']['allocation_core_hours'] == 140.24
    assert not audit['accounting']['normal_completion']  # not a final whole-job cost


def test_positive_local_curvature_does_not_promote_stationarity_TS_or_predictions(evidence):
    audit, summary, spec, _, _ = evidence
    assert summary['base_released_gradient_norm_eV_A'] == spec['base_released_gradient_norm_eV_A']
    assert summary['base_released_gradient_norm_eV_A'] == pytest.approx(0.02629995445855942, abs=ATOL, rel=0)
    assert 0 < summary['base_released_gradient_norm_eV_A'] < .03
    assert not summary['physical_stability_TS_or_independent_prediction_certified']
    assert not audit['stable_branch_TS_or_independent_prediction_certified']
    assert not audit['physical_reserved_geometry_generated_or_read']
    assert audit['new_DFT_or_scheduler_writes'] == 0
