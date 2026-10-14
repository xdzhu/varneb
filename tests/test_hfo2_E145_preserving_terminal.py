"""Actual +1% preserving terminal: portable replay, not a saddle certificate."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
CASE = REPO/'benchmarks/hfo2_channels/20261008/clamped_plus1_preserving_terminal_E145_20261012'
ROOT = CASE/'completed_HF'
AUDIT_SHA = '0c8d91532bc9051decdfd600eae12ce409fd84a86c677158077948fa352c9d8a'
FMAX = 0.09928337557669162


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def receipt():
    assert sha(ROOT/'audit_receipt.json') == AUDIT_SHA
    return json.loads((ROOT/'audit_receipt.json').read_text())


def test_native_receipt_covers_all_fresh_calls_and_original_runtime():
    from examples.hfo2_fixed_input_factory import CONTRACT
    from scripts.audit_hfo2_clamped_terminal import validate_terminal

    r = receipt()
    assert r['job_id'] == '28812952' and r['terminal_step'] == 17
    assert r['ordinary_converged']
    summary = json.loads((ROOT/'vcneb_summary.json').read_text())
    assert sha(ROOT/'vcneb_summary.json') == r['terminal_summary_sha256']
    validate_terminal('28812952', r['scheduler_row'], summary, False)
    assert r['fresh_interior_SCFs'] == r['expected_fresh_interior_SCFs'] == 119
    assert r['initial_fresh_interior_SCFs'] == 0
    assert len(r['calls']) == 119 and len(r['runtime_code_sha256']) == 338
    assert {(Path(row['path']).parts[-3], Path(row['path']).parts[-2]) for row in r['calls']} == {
        (f'image_{i:04d}', f'scf_{j:06d}') for i in range(1, 8) for j in range(17)}
    assert all({n: row['input_sha256'][n] for n in CONTRACT} == CONTRACT for row in r['calls'])
    assert r['all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked']
    assert r['terminal_observation']['n_total_images'] == 9
    assert r['new_DFT_calls'] == r['scheduler_mutations'] == 0
    assert not r['physical_inputs_or_running_source_changed']
    assert not r['holdout_generated_or_read']
    assert not r['continuous_TS_or_barrier_error_bound_certified']


def test_portable_nine_frames_reparse_exact_geometry_energy_force_stress():
    from ase.io import read
    from examples.hfo2_fixed_input_factory import read_fixed_hfo2_stru, same_ordered_geometry
    from scripts.analyze_hfo2_clamped_residual import analyze
    from scripts.audit_hfo2_static_replica import audited_results

    r = receipt()
    observation = ROOT/'observations/step_0017'
    report = json.loads((observation/'observation.json').read_text())
    assert report == r['terminal_observation']
    residual = analyze(observation)
    assert residual['ordinary_residual_pass']
    assert residual['fmax_eV_A'] == pytest.approx(FMAX, abs=1e-12)
    images = read(observation/'evaluated_chain.traj', index=':')
    assert len(images) == 9
    for image, row in zip(images, report['raw_image_evaluations']):
        raw = observation/'raw'/f"image_{row['image_index']:04d}"
        assert sha(raw/'OUT.ABACUS/running_scf.log') == row['raw_log_sha256']
        assert sha(raw/'source_audit.json') == row['audit_sha256']
        assert all(sha(raw/n) == row['input_sha256'][n] for n in ('INPUT', 'KPT', 'STRU'))
        assert same_ordered_geometry(image, read_fixed_hfo2_stru(raw/'STRU'))
        actual = audited_results(raw)
        for key, field in [('energy', 'energy_eV_cell'), ('forces', 'forces_eV_A'),
                           ('stress', 'stress_ASE_voigt_eV_A3')]:
            np.testing.assert_allclose(actual[key], row[field], atol=1e-12, rtol=0)
    assert not (ROOT/'seed').exists()  # terminal, not an unconverged rolling seed


def test_actual_segment_cost_and_driver_ancestry_are_not_whole_path_cost():
    r = receipt()
    assert r['fresh_SCF_wall_seconds_sum'] == pytest.approx(11200.040381577797, abs=1e-8)
    assert r['fresh_SCF_wall_seconds_sum'] == pytest.approx(sum(x['elapsed_seconds'] for x in r['calls']), abs=1e-8)
    assert r['fresh_SCF_core_hour_proxy_32'] == pytest.approx(99.55591450291375, abs=1e-9)
    driver = json.loads((ROOT/'driver_receipt.json').read_text())
    assert driver['job_id'] == '28812952' and driver['parent_job'] == '28737465'
    assert driver['terminal_audit_receipt_sha256'] == AUDIT_SHA
    helper = REPO/'benchmarks/hfo2_channels/20261008/clamped_preserving_continue_E141_20261012/continue_preserving_HF.py'
    assert driver['executed_driver_sha256'] == sha(helper)
    assert driver['new_DFT_calls'] == driver['scheduler_writes'] == driver['new_independent_chains'] == 0


def test_all_scalar_reports_share_training_receipt_and_no_prediction_claim():
    names = ('displacement', 'gradient', 'path_work', 'external_work')
    reports = [json.loads((CASE/'analysis'/f'{n}.json').read_text()) for n in names]
    for r in reports:
        assert r['terminal_audit_sha256'] == AUDIT_SHA
        assert r['source_job_id'] == '28812952' and r['strain_training_condition'] == .01
        assert r['new_DFT_calls'] == 0
    d, g, w, e = reports
    p = d['highest_image_index']
    assert p == g['sampled_peak_image'] == w['sampled_peak_image'] == e['sampled_peak_image'] == 4
    assert d['relative_energy_meV_fu'][p] == pytest.approx(125.42443115808055, abs=1e-9)
    assert d['T_triplet_rank3_fraction'][p] == pytest.approx(.5217508561554873, abs=1e-12)
    spans = g['sampled_peak']['representations']
    assert spans['T_geometric_triplet']['retained_squared_gradient_fraction'] == pytest.approx(.32878978855005697)
    assert spans['T_two_lowest_complete_optical_doublets']['source_basis_rank'] == 4
    for frame in ('103', '106', '114', '127'):
        span = spans['Cmma_frame_'+frame]
        assert span['source_basis_rank'] == 4 and not span['source_modes_or_forces_ASR_corrected']
        assert .49466696 < span['retained_squared_gradient_fraction'] < .49466700
    assert not w['G3_selection_or_holdout_access']
    assert not d['physical_parameters_changed']
