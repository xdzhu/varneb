"""Actual first-amplitude replay and arithmetic guards, not TS certification."""
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from scripts import summarize_hfo2_G3_curvature as summary

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'benchmarks/hfo2_channels/20261008'
CASE = BASE/'G3_first_amplitude_E155_20261014'
PINS = dict(T='1cc4269a8f1437acddb94c9a7c264f52bae888669db2bcc36521982d280666ea',
            preserving='9bfeb3c8a7eca7533246ffd61e6e43df17262fb8ef274dfa87965c1301a10efa')


def assert_derived_arithmetic_same(actual, expected):
    """BLAS-roundoff policy for derived scalars only, never native file bytes."""
    assert type(actual) is type(expected)
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys()
        for key in expected: assert_derived_arithmetic_same(actual[key], expected[key])
    elif isinstance(expected, list):
        assert len(actual) == len(expected)
        for a, e in zip(actual, expected): assert_derived_arithmetic_same(a, e)
    elif isinstance(expected, float):
        assert actual == pytest.approx(expected, abs=2e-12, rel=0)
    else:
        assert actual == expected


@pytest.fixture(scope='module')
def material():
    mp = BASE/'G3_curvature_E150_20261014/manifest.json'
    assert hashlib.sha256(mp.read_bytes()).hexdigest() == summary.MANIFEST_SHA256
    reports = {}
    for channel, pin in PINS.items():
        data = (CASE/(channel+'_native.json')).read_bytes()
        assert hashlib.sha256(data).hexdigest() == pin
        reports[channel] = json.loads(data)
    return json.loads(mp.read_bytes()), reports


def test_actual_native_prefix_replay_keeps_incomplete_second_amplitude(material):
    m, reports = material
    result = summary.summarize(m, reports)
    saved = json.loads((CASE/'first_amplitude_summary.json').read_bytes())
    assert_derived_arithmetic_same(result, {k:v for k,v in saved.items() if k != 'input_sha256'})
    assert [r['audited_prefix_records'] for r in result['rows']] == [76, 77]
    for row in result['rows']:
        assert row['complete_amplitudes'] == [.01]
        assert not row['two_step_comparison_available']
        block = row['blocks'][0]
        assert block['raw_negative_count'] == 2
        assert block['tangent_complement_raw_negative_count'] == 1
        assert block['minimum_retained_dimension_for_raw_positive_complement'] == 2
        assert block['positive_spectral_complement_dimension'] == 34
        assert block['positive_spectral_complement_gradient_norm_eV_A'] > .1
        assert not block['physical_error_bound_or_stationary_index_or_TS_certificate']
    assert not result['new_raw_SCF_reparse']
    assert not result['conditional_axes_or_forecasts_selected']


@pytest.mark.parametrize('channel', ['T', 'preserving'])
def test_internal_rotation_and_external_control_do_not_change_raw_signs(material, channel):
    m, reports = material
    center = m['centers'][channel]
    r = reports[channel]
    a = np.asarray(r['blocks'][0]['raw_projected_eV_A2'])
    b = np.asarray(center['directions'])
    g = np.asarray(r['raw_center_gradient_eV_A'])
    t = np.asarray(center['ordered_local_tangent'])
    baseline = summary.block_diagnostics(a, b, g, t)
    rotation = np.eye(37)
    rotation[:36, :36] = np.linalg.qr(np.random.default_rng(155).normal(size=(36, 36)))[0]
    changed_a = rotation.T @ a @ rotation
    changed_a[-1, -1] = -1e6  # prescribed control is never an internal instability
    changed_g = g.copy(); changed_g[-1] += 1e6
    changed = summary.block_diagnostics(changed_a, b @ rotation, changed_g, t)
    for key in ('symmetrized_internal_eigenvalues_eV_A2', 'tangent_complement_eigenvalues_eV_A2'):
        np.testing.assert_allclose(changed[key], baseline[key], atol=3e-12, rtol=0)
    for key in ('raw_negative_count', 'tangent_complement_raw_negative_count',
                'minimum_retained_dimension_for_raw_positive_complement'):
        assert changed[key] == baseline[key]
    assert changed['measured_internal_gradient_norm_eV_A'] == pytest.approx(
        baseline['measured_internal_gradient_norm_eV_A'], abs=1e-14)


def test_tangent_translation_is_explicit_not_a_force_or_path_repair(material):
    m, reports = material
    c = m['centers']['T']; r = reports['T']
    b = np.asarray(c['directions']); t = np.asarray(c['ordered_local_tangent'])
    g = np.asarray(r['raw_center_gradient_eV_A']); original_g = g.copy()
    baseline = summary.block_diagnostics(r['blocks'][0]['raw_projected_eV_A2'], b, g, t)
    excluded = t-b[:, :36] @ (b[:, :36].T @ t)
    assert baseline['native_tangent_excluded_translation_norm'] == pytest.approx(.10168809680225134)
    changed = summary.block_diagnostics(r['blocks'][0]['raw_projected_eV_A2'], b, g, t+excluded)
    np.testing.assert_allclose(changed['tangent_complement_eigenvalues_eV_A2'],
                               baseline['tangent_complement_eigenvalues_eV_A2'], atol=1e-12, rtol=0)
    assert changed['native_tangent_excluded_translation_norm'] == pytest.approx(
        2*baseline['native_tangent_excluded_translation_norm'])
    np.testing.assert_array_equal(g, original_g)


@pytest.mark.parametrize('fault', ['partial_matrix', 'nonfinite', 'nonorthogonal', 'mixed_control', 'control_tangent', 'zero_tangent'])
def test_malformed_blocks_are_refused(material, fault):
    m, reports = material
    c = m['centers']['T']; r = reports['T']
    a = np.array(r['blocks'][0]['raw_projected_eV_A2'])
    b = np.array(c['directions']); g = np.array(r['raw_center_gradient_eV_A'])
    t = np.array(c['ordered_local_tangent'])
    if fault == 'partial_matrix': a = a[:-1, :-1]
    elif fault == 'nonfinite': a[0, 0] = np.nan
    elif fault == 'nonorthogonal': b[:, 0] *= 2
    elif fault == 'mixed_control': b[0, -1] = .1
    elif fault == 'control_tangent': t[-1] = .01
    elif fault == 'zero_tangent': t *= 0
    with pytest.raises(ValueError): summary.block_diagnostics(a, b, g, t)


@pytest.mark.parametrize('fault', ['failure', 'cached_spectrum', 'partial_spectrum', 'wrong_count',
                                 'wrong_status', 'wrong_center', 'missing_block', 'false_full_flag'])
def test_report_faults_are_not_repaired(material, fault):
    m, reports = material
    changed = copy.deepcopy(reports)
    r = changed['T']
    if fault == 'failure': r['failures_present'] = ['explicit test failure, not native']
    elif fault == 'cached_spectrum': r['blocks'][0]['fixed_strain_internal_eigenvalues_eV_A2'][0] += .1
    elif fault == 'partial_spectrum': r['blocks'][1]['fixed_strain_internal_eigenvalues_eV_A2'] = [1.]
    elif fault == 'wrong_count': r['raw_SCFs_reparsed'] -= 1
    elif fault == 'wrong_status': r['status'] = 'complete_local_curvature_diagnostic'
    elif fault == 'wrong_center': r['center_source'] = reports['preserving']['center_source']
    elif fault == 'missing_block': r['blocks'].pop()
    elif fault == 'false_full_flag': r['blocks'][0]['full_stencil_at_this_step'] = False
    with pytest.raises(ValueError): summary.summarize(m, changed)


def test_actual_CLI_arithmetic_replay_and_no_overwrite_or_wrong_input_bytes(tmp_path):
    output = tmp_path/'summary.json'
    command = [sys.executable, '-m', 'scripts.summarize_hfo2_G3_curvature', '--manifest',
        str(BASE/'G3_curvature_E150_20261014/manifest.json'), '--output', str(output)]
    for channel, pin in PINS.items():
        command += ['--'+channel+'-report', str(CASE/(channel+'_native.json')), '--'+channel+'-sha256', pin]
    first = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    assert first.returncode == 0, first.stderr
    assert_derived_arithmetic_same(json.loads(output.read_bytes()),
                                  json.loads((CASE/'first_amplitude_summary.json').read_bytes()))
    before = output.read_bytes()
    assert subprocess.run(command, cwd=ROOT, capture_output=True).returncode != 0
    assert output.read_bytes() == before
    wrong = command.copy(); wrong[wrong.index('--output')+1] = str(tmp_path/'refused.json')
    wrong[wrong.index('--T-sha256')+1] = '0'*64
    assert subprocess.run(wrong, cwd=ROOT, capture_output=True).returncode != 0
    assert not (tmp_path/'refused.json').exists()
