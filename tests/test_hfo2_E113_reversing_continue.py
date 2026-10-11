"""Real capped-chain evidence; a successful cap is not a converged label."""
import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
CASE = REPO / 'benchmarks/hfo2_channels/20261008/clamped_reversing_continue_E113_20261012'
ROOT = CASE / 'plus1_reversing_starter'
AUDIT_SHA = '9536903dc519e3abfb597ca3072821ef16786272cb2bef517333681803979c76'
SEED_SHA = '571070ea12afeab882b82420ec02d71e7afab3addac3c150761f4d9edeeda043'
HELPER_SHA = '6c64c0155800082e28cbbb280511ac761929bdb6ccc5757fa8d68d156b32eac7'
WRAPPER_SHA = 'ab914119aaa7918bd0cb466ed8dd86daafecb20609a077a5d249d89a3228eed0'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_actual_cap_and_all_call_receipt_are_not_convergence():
    from examples.hfo2_fixed_input_factory import CONTRACT
    from scripts.audit_hfo2_clamped_terminal import validate_terminal

    assert sha(ROOT / 'audit_receipt.json') == AUDIT_SHA
    r = json.loads((ROOT / 'audit_receipt.json').read_text())
    s = json.loads((ROOT / 'vcneb_summary.json').read_text())
    assert sha(ROOT / 'vcneb_summary.json') == r['terminal_summary_sha256']
    validate_terminal('28709791', r['scheduler_row'], s, False, outcome='step-cap')
    with pytest.raises(ValueError):
        validate_terminal('28709791', r['scheduler_row'], s, False)
    assert r['terminal_step'] == 10 and not r['ordinary_converged']
    assert r['initial_fresh_interior_SCFs'] == 7
    assert r['fresh_interior_SCFs'] == r['expected_fresh_interior_SCFs'] == 77
    assert len(r['runtime_code_sha256']) == 338
    assert len(r['calls']) == 77
    assert {(Path(c['path']).parts[-3], Path(c['path']).parts[-2]) for c in r['calls']} == {
        (f'image_{i:04d}', f'scf_{j:06d}') for i in range(1, 8) for j in range(11)}
    assert all({k: c['input_sha256'][k] for k in CONTRACT} == CONTRACT for c in r['calls'])
    assert r['all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked'] is True
    assert r['current_frame_exact_caches_checked'] == 9
    assert r['new_DFT_calls'] == r['scheduler_mutations'] == r['new_independent_chains'] == 0
    assert not r['physical_inputs_or_running_source_changed']
    assert not r['holdout_generated_or_read'] and not r['FIRE_state_restored']


def test_nine_retained_raw_scfs_and_same_boundary_residual():
    from scripts.analyze_hfo2_clamped_residual import analyze
    from scripts.audit_hfo2_static_replica import audited_results

    observation = ROOT / 'observations/step_0010'
    report = json.loads((observation / 'observation.json').read_text())
    r = json.loads((ROOT / 'audit_receipt.json').read_text())
    assert report == r['terminal_observation']
    replay = analyze(observation)
    assert replay['fmax_eV_A'] == pytest.approx(0.12692249290222424, abs=1e-12)
    assert not replay['ordinary_residual_pass']
    assert report['mechanical_boundary']['fixed_cell_rows'] == [0, 1]
    assert report['mechanical_boundary']['allow_tilt'] is True
    assert report['pressure_GPa'] == 0 and report['climb'] is False
    assert report['formula_units'] == 4
    assert report['n_total_images'] == 9 and report['n_active_images'] == 7
    assert len(report['raw_image_evaluations']) == 9
    for row in report['raw_image_evaluations']:
        raw = observation / 'raw' / f"image_{row['image_index']:04d}"
        assert sha(raw / 'OUT.ABACUS/running_scf.log') == row['raw_log_sha256']
        assert sha(raw / 'source_audit.json') == row['audit_sha256']
        assert all(sha(raw / n) == row['input_sha256'][n] for n in ('INPUT', 'KPT', 'STRU'))
        parsed = audited_results(raw)
        for key, expected in (('energy', row['energy_eV_cell']),
                              ('forces', row['forces_eV_A']),
                              ('stress', row['stress_ASE_voigt_eV_A3'])):
            assert np.allclose(parsed[key], expected, atol=1e-12, rtol=0)


def test_latest_seed_geometry_and_exact_cache_manifest():
    from ase.io import read
    from examples.hfo2_fixed_input_factory import same_ordered_geometry

    seed = ROOT / 'seed'
    assert sha(seed / 'manifest.json') == SEED_SHA
    m = json.loads((seed / 'manifest.json').read_text())
    assert m['source_job_id'] == '28709791' and m['source_step'] == 10
    assert m['current_frame_exact_caches'] == 9
    assert not m['optimizer_state_restored'] and not m['physical_inputs_changed']
    assert not m['holdout_generated_or_read']
    assert all(sha(seed / n) == h for n, h in m['files_sha256'].items())
    observation = ROOT / 'observations/step_0010'
    assert sha(observation / 'observation.json') == m['source_observation_sha256']
    old = read(observation / 'evaluated_chain.traj', index=':')
    current = read(seed / 'seed.traj', index=':')
    assert len(old) == len(current) == 9
    for a, b in zip(old, current):
        assert same_ordered_geometry(a, b)
        for key in ('energy', 'forces', 'stress'):
            assert np.allclose(a.calc.results[key], b.calc.results[key], atol=1e-12, rtol=0)


def test_reused_production_helper_and_terminal_driver_bytes():
    helper = REPO / 'benchmarks/hfo2_channels/20261008/clamped_flip_continue_E103_20261012/continue_flips_HF.py'
    assert sha(helper) == HELPER_SHA
    assert sha(REPO / 'cluster/hf_hfo2_clamped_continue_E103_20261012.slurm') == WRAPPER_SHA
    driver = json.loads((ROOT / 'driver_receipt.json').read_text())
    assert driver['job_id'] == driver['original_parent'] == '28709791'
    assert driver['executed_driver_sha256'] == HELPER_SHA
    assert driver['terminal_audit_receipt_sha256'] == AUDIT_SHA
    assert driver['new_DFT_calls'] == driver['scheduler_writes'] == driver['new_independent_chains'] == 0


def capture_fixture(tmp_path, *, change=None):
    path = CASE / 'capture_start_HF.py'
    spec = importlib.util.spec_from_file_location('E113_startup_fixture', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    work = tmp_path / 'band'
    work.mkdir()
    caches = [dict(directory='fixture') for _ in range(9)]
    if change == 'cache':
        caches[4] = None
    (work / 'vcneb_preflight.json').write_text(json.dumps(dict(
        factory='wrong' if change == 'factory' else
        'examples.hfo2_fixed_input_factory:make_clamped_resume_cached_factory',
        calculator_parameters=dict(seed_cache_records=caches))))
    (work / 'vcneb.opt.log').write_text('' if change == 'log' else
        'CheckedFIRE: 0 04:00:00 -9781.593155 0.126922\n')
    raw = work / 'image_0001/scf_000000/OUT.ABACUS/running_scf.log'
    raw.parent.mkdir(parents=True)
    raw.write_text('DSIZE = 1\n' if change == 'DSIZE' else 'DSIZE = 32\n')
    root = tmp_path / '28709791'
    root.mkdir()
    for name in ('dispatch_receipt.json', 'dispatch_journal.jsonl'):
        (root / name).write_text('fixture\n')
    (root / 'audit_receipt.json').write_text('fixture\n')
    stderr = tmp_path / 'stderr'
    stderr.write_text('')
    calls = []
    def run(argv):
        calls.append(argv)
        if argv[0] == 'sacct':
            return '99900001|' + ('PENDING' if change == 'state' else 'RUNNING') + '|0:0|32|hfacnormal01|iai806|start|Unknown'
        if argv[0] == 'squeue':
            return 'fixture queue'
        assert argv[:3] == ['scontrol', 'show', 'job']
        return 'fixture control'
    m = SimpleNamespace(ROOT=tmp_path, accepted=lambda parent: {'job_id': '99900001'},
        paths=lambda parent: dict(work=work, stderr=stderr, terminal=root),
        run=run, slots=lambda *args, **kwargs: ['28737465', '99900001'])
    return module, m, calls


@pytest.mark.parametrize('change', ['state', 'cache', 'factory', 'log', 'DSIZE'])
def test_snapshot_rejects_unstarted_or_inconsistent_native_evidence(tmp_path, change):
    module, m, calls = capture_fixture(tmp_path, change=change)
    with pytest.raises(ValueError):
        module.capture(m)
    assert all(c[0] in {'sacct', 'squeue', 'scontrol'} for c in calls)
    assert not any(c[:2] == ['scontrol', 'release'] for c in calls)


def test_snapshot_synthetic_success_never_launches_dft_or_changes_scheduler(tmp_path):
    module, m, calls = capture_fixture(tmp_path)
    report = module.capture(m)
    assert report['seed_cache_count'] == 9 and report['none_seed_cache_count'] == 0
    assert report['native_first_DSIZE'] == ['32']
    assert report['fresh_call_audits_at_capture'] == 0
    assert not report['first_fresh_call_audit_exists']
    assert report['new_readonly_DFT_calls'] == report['scheduler_writes_for_capture'] == 0
    assert report['startup_is_not_ordinary_convergence_or_full_fresh_call_audit'] is True
    assert [c[0] for c in calls] == ['sacct', 'squeue', 'scontrol']


def test_actual_single_admission_preserves_two_slot_limit():
    r = json.loads((CASE / 'dispatch_receipt.json').read_text())
    journal = CASE / 'dispatch_journal.jsonl'
    rows = [json.loads(line) for line in journal.read_text().splitlines()]
    assert r['job_id'] == '28764704' and r['parent_job'] == '28709791'
    assert r['source_step'] == 10 and r['max_new_updates'] == 40 and r['wall_hours'] == 12
    assert r['total_submissions'] == r['release_writes'] == 1
    assert r['dependency_writes'] == r['new_independent_chains'] == 0
    assert r['audit_receipt_sha256'] == AUDIT_SHA
    assert r['seed_manifest_sha256'] == SEED_SHA
    assert r['executed_helper_sha256'] == HELPER_SHA and r['script_sha256'] == WRAPPER_SHA
    assert sha(journal) == r['journal_sha256']
    assert [row['event'] for row in rows] == [
        'one_held_submission_intent', 'one_accepted_handle', 'one_release_intent',
        'one_existing_handle_released', 'released_slot_independently_corroborated']
    assert rows[0]['committed_slots'] == ['28737465']
    assert rows[-1]['committed_slots'] == ['28737465', '28764704']
    assert all(row['job_id'] == '28764704' for row in rows[1:])
    assert rows[2]['accepted_native_state']['Reason'] == 'JobHeldUser'
    native = rows[-1]['native_control']
    assert native['Reason'] != 'JobHeldUser' and native['Priority'] != '0'
    assert native['NumCPUs'] == '32' and native['TimeLimit'] == '12:00:00'
    assert not r['FIRE_optimizer_state_restored'] and not r['execution_started_claimed']
    assert not r['physical_inputs_or_running_source_changed'] and not r['holdout_generated_or_read']


def test_actual_native_startup_cache_replay_is_not_final_convergence():
    assert sha(CASE / 'startup_observation.json') == '47c33b0eadff46fa5ec5dd2800cc55a8091ac3a7a5ba269e4d27db6894869883'
    r = json.loads((CASE / 'startup_observation.json').read_text())
    assert r['source_job_id'] == '28764704' and r['original_parent'] == '28709791'
    assert '28764704|RUNNING|0:0|32|hfacnormal01|iai806|2026-10-12T04:27:20' in r['native_accounting']
    assert r['committed_study_jobs'] == ['28737465', '28764704']
    assert r['native_first_DSIZE'] == ['32']
    assert r['seed_cache_count'] == 9 and r['none_seed_cache_count'] == 0
    assert r['executed_capture_sha256'] == sha(CASE / 'capture_start_HF.py')
    assert r['runtime_preflight_sha256'] == sha(CASE / 'child_runtime_preflight.json')
    info = json.loads((CASE / 'child_runtime_preflight.json').read_text())
    assert info['factory'] == r['factory']
    assert info['calculator_parameters'] == json.loads((ROOT / 'seed/factory_parameters.json').read_text())
    assert r['optimizer_log_raw'] == (CASE / 'child_optimizer_log_snapshot.log').read_text()
    row = r['optimizer_log_raw'].splitlines()[-1].split()
    assert row[:2] == ['CheckedFIRE:', '0']
    assert float(row[-1]) == pytest.approx(0.12692249290222424, abs=5e-7)
    assert r['admission_receipt_sha256'] == sha(CASE / 'dispatch_receipt.json')
    assert r['admission_journal_sha256'] == sha(CASE / 'dispatch_journal.jsonl')
    assert r['parent_audit_sha256'] == AUDIT_SHA
    assert r['stderr_bytes'] == 0 and not r['failure_exists']
    assert not r['terminal_summary_exists']
    assert r['startup_is_not_ordinary_convergence_or_full_fresh_call_audit'] is True
    assert r['new_readonly_DFT_calls'] == r['scheduler_writes_for_capture'] == 0
    assert not r['running_source_or_physical_inputs_changed'] and not r['holdout_generated_or_read']
