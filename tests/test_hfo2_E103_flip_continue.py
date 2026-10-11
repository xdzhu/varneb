"""Synthetic admission/startup tests; none constitutes a physical DFT result."""
import importlib.util
import json
from pathlib import Path
import shlex
import subprocess

import pytest

REPO = Path(__file__).resolve().parents[1]
CASE = REPO/'benchmarks/hfo2_channels/20261008/clamped_flip_continue_E103_20261012'
SCRIPT = REPO/'cluster/hf_hfo2_clamped_continue_E103_20261012.slurm'


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def manager():
    return load(CASE/'continue_flips_HF.py', 'E103_fixture_manager')


def recipe():
    m = load(REPO/'tests/test_hfo2_rolling_continue_recipe.py', 'E103_reused_recipe_fixture')
    m.SCRIPT = SCRIPT
    m.PARENTS = {p: (10, 7, 338, m.RUNS+'20261010-varneb-nested-controls-E055-r1/source',
        m.RUNS+'20261010-varneb-clamped-G2-E060-r1/'+c+'/band') for p, c in manager().CHANNELS.items()}
    return m


@pytest.mark.parametrize('parent', ['28709790', '28709791'])
def test_startup_exact_program_with_original_parent_geometry(tmp_path, monkeypatch, parent, capsys):
    m = recipe()
    m.run_validator(*m.fixture(tmp_path, parent), parent, monkeypatch)
    assert float(capsys.readouterr().out) == 5.12968067458423


@pytest.mark.parametrize('key,value', [
    ('ordinary_converged', True), ('terminal_step', 9), ('fresh_interior_SCFs', 76),
    ('current_frame_exact_caches_checked', 8), ('initial_fresh_interior_SCFs', 0),
    ('physical_inputs_or_running_source_changed', True), ('holdout_generated_or_read', True),
    ('FIRE_state_restored', True), ('scheduler_mutations', 1), ('new_DFT_calls', 1),
    ('scheduler_row', '28709790|RUNNING|0:0|32|hfacnormal01|iai806')])
def test_unfinished_changed_or_incomplete_parent_rejected(tmp_path, monkeypatch, key, value):
    m = recipe()
    seed, audit, runtime, source = m.fixture(tmp_path, '28709790')
    receipt = json.loads(audit.read_text())
    receipt[key] = value
    m.write_json(audit, receipt)
    with pytest.raises(ValueError):
        m.run_validator(seed, audit, runtime, source, '28709790', monkeypatch)


def test_hash_and_production_driver_byte_unchanged():
    m = manager()
    assert m.sha(SCRIPT) == m.SCRIPT_SHA
    old = (REPO/'cluster/hf_hfo2_clamped_continue_E090_20261011.slurm').read_text()
    assert SCRIPT.read_text().split('module purge', 1)[1] == old.split('module purge', 1)[1]
    assert m.sha(REPO/'benchmarks/hfo2_channels/20261008/clamped_child_audit_E093_20261011/audit_child_HF.py') == m.LEGACY_SHA


def test_future_native_accounting_keeps_channel_named_logs():
    from scripts.audit_hfo2_clamped_terminal import validate_accounting_source
    m = manager()
    p = m.paths('28709790')
    script = m.SCRIPT
    env = '--export=ALL,RUN_DFT=1,WORKDIR='+str(p['work'])
    submit = 'sbatch '+shlex.quote(env)+' '+shlex.quote(str(script))
    body = '|'.join(['99900001', 'COMPLETED', 'iai806', '32', 'hfacnormal01',
        'hfo2-G2-E103-28709790', str(p['stdout']), str(p['stderr']), submit])
    assert validate_accounting_source(body, '99900001', p['work'], script)['NumCPUs'] == '32'


def test_audit_parent_rejects_live_before_loading_support(tmp_path, monkeypatch):
    m = manager()
    monkeypatch.setattr(m, 'ROOT', tmp_path/'audit-only')
    monkeypatch.setattr(m, 'run', lambda *a: native_parent(m, state='RUNNING'))
    monkeypatch.setattr(m, 'load_auditor', lambda: pytest.fail('no imports for a live parent'))
    with pytest.raises(ValueError): m.audit('28709790')
    assert not m.ROOT.exists()


def native_parent(m, parent='28709790', state='COMPLETED', change=None):
    channel = m.CHANNELS[parent]
    env = dict(RUN_DFT='1', SOURCE_ROOT=str(m.SOURCE),
        SEED_ROOT=str(m.RUNS/'20261010-varneb-clamped-G2-remaining-E056-r1/seeds/strain_p0100'/channel),
        WORKDIR=str(m.paths(parent)['old_work']))
    if change == 'override': env['ECUTWFC'] = '120'
    if change in env: env[change] = 'different'
    submit = 'sbatch '+shlex.quote('--export=ALL,'+','.join(k+'='+v for k, v in env.items()))+' '+shlex.quote(str(m.PILOT))
    if change == 'duplicate': submit += ' --export=ALL'
    fields = [parent, state, 'iai806', '32', 'hfacnormal01', 'hfo2-G2-E060-'+channel,
        str(m.ORIGINAL/(channel+'.slurm.out')), str(m.ORIGINAL/(channel+'.slurm.err')), submit]
    if isinstance(change, int): fields[change] = 'different'
    return '|'.join(fields)


@pytest.mark.parametrize('change', [None, 0, 2, 3, 4, 5, 6, 7, 'SOURCE_ROOT', 'WORKDIR', 'SEED_ROOT', 'RUN_DFT', 'override', 'duplicate'])
def test_native_original_submission_contract(change):
    m = manager()
    if change is None:
        m.identity(native_parent(m), '28709790', '28709790')
    else:
        with pytest.raises(ValueError): m.identity(native_parent(m, change=change), '28709790', '28709790')


@pytest.mark.parametrize('state', ['RUNNING', 'PENDING', 'TIMEOUT', 'FAILED', 'CANCELLED'])
def test_live_terminal_or_timeout_is_not_a_successful_cap(state):
    m = manager()
    with pytest.raises(ValueError): m.identity(native_parent(m, state=state), '28709790', '28709790')


def queue(m, *, child=None, released=False):
    body = '28709791|hfo2-G2-E060-PO_flip_T_pattern_reversing|RUNNING|hfacnormal01|32|node1'
    if child:
        body += '\n'+child+'|hfo2-G2-E103-28709790|PENDING|hfacnormal01|32|'+('(Priority)' if released else '(JobHeldUser)')
    return body


def test_runnable_pending_and_unknown_jobs_guard_the_ceiling():
    m = manager()
    extra = {'99900001': 'hfo2-G2-E103-28709790'}
    assert m.slots(queue(m, child='99900001'), '28709790', spare=True, extra=extra) == ['28709791']
    with pytest.raises(ValueError): m.slots(queue(m, child='99900001', released=True), '28709790', spare=True, extra=extra)
    with pytest.raises(ValueError): m.slots(queue(m)+'\nunknown|hfo2-foreign|RUNNING|hfacnormal01|32|node', '28709790', spare=True)
    with pytest.raises(ValueError): m.slots(queue(m)+'\n'+queue(m), '28709790', spare=False)


def dispatch_fixture(tmp_path, monkeypatch, mode='success'):
    m = manager()
    monkeypatch.setattr(m, 'ROOT', tmp_path/'study')
    monkeypatch.setattr(m, 'LOCK', tmp_path/'shared_admission.lock')
    monkeypatch.setattr(m, 'SCRIPT', tmp_path/'wrapper.slurm')
    m.SCRIPT.write_text('# fixture, not production\n')
    root = m.ROOT/'28709790'
    root.mkdir(parents=True)
    monkeypatch.setattr(m, 'verify_audit', lambda *a: {'seed_manifest_sha256': 'a'*64})
    calls, released = [], []
    def response(argv):
        calls.append(argv)
        if argv[0] == 'sacct': return '28709790|'+('RUNNING' if mode == 'live' else 'COMPLETED')+'|0:0|32|hfacnormal01|iai806'
        if argv[0] == 'bash': return ''
        if argv[0] == 'squeue':
            if mode == 'capacity' and sum(c[0] == 'squeue' for c in calls) >= 2:
                return queue(m)+'\n28727368|hfo2-G2-E091-28709789|RUNNING|hfacnormal01|32|node6'
            return queue(m, child='99900001' if any(c[0] == 'sbatch' for c in calls) else None, released=bool(released))
        if argv[0] == 'sbatch':
            records = [json.loads(l) for l in (root/'dispatch_journal.jsonl').read_text().splitlines()]
            assert records[-1]['event'] == m.EVENTS[0]
            return '99900001;ambiguous' if mode == 'ambiguous' else '99900001'
        if argv[:2] == ['scontrol', 'release']:
            if mode == 'timeout': raise subprocess.TimeoutExpired(argv, 45)
            if mode != 'stale': released.append(True)
            return ''
        p = m.paths('28709790')
        row = dict(JobId='99900001', JobName='hfo2-G2-E103-28709790',
            JobState='PENDING', Reason='Priority' if released else 'JobHeldUser',
            Priority='1' if released else '0', Partition='hfacnormal01',
            NumCPUs='16' if mode == 'wrong_cpu' else '32', NumTasks='1',
            Command=str(m.SCRIPT), Dependency='(null)', StdOut=str(p['stdout']), StdErr=str(p['stderr']),
            TimeLimit='12:00:00', Restarts='0', NumNodes='1', UserId='iai806(1000)', **{'CPUs/Task': '32'})
        return ' '.join(k+'='+shlex.quote(v) for k, v in row.items())
    monkeypatch.setattr(m, 'run', response)
    return m, root, calls


def test_exact_one_submission_release_and_native_slot_journal(tmp_path, monkeypatch):
    m, root, calls = dispatch_fixture(tmp_path, monkeypatch)
    result = m.dispatch('28709790', 'b'*64)
    assert result['total_submissions'] == result['release_writes'] == 1
    assert result['new_independent_chains'] == result['dependency_writes'] == 0
    assert not result['FIRE_optimizer_state_restored'] and not result['execution_started_claimed']
    assert not m.LOCK.exists()
    assert m.accepted('28709790')['job_id'] == '99900001'
    before = len(calls)
    with pytest.raises(FileExistsError): m.dispatch('28709790', 'b'*64)
    assert len(calls) == before
    submit = next(c for c in calls if c[0] == 'sbatch')
    assert '--hold' in submit and '--parsable' in submit
    exported = next(c for c in submit if c.startswith('--export='))
    assert 'SOURCE_ROOT='+str(m.SOURCE) in exported and 'PARENT_JOB_ID=28709790' in exported
    assert len(m.slots(queue(m, child='99900001', released=True), '28709790', spare=False)) == 2


@pytest.mark.parametrize('mode,submits,releases', [('live', 0, 0), ('capacity', 0, 0),
    ('ambiguous', 1, 0), ('wrong_cpu', 1, 0), ('timeout', 1, 1), ('stale', 1, 1)])
def test_failed_or_ambiguous_intent_never_repeated(tmp_path, monkeypatch, mode, submits, releases):
    m, root, calls = dispatch_fixture(tmp_path, monkeypatch, mode)
    with pytest.raises((ValueError, subprocess.TimeoutExpired)): m.dispatch('28709790', 'b'*64)
    assert sum(c[0] == 'sbatch' for c in calls) == submits
    assert sum(c[:2] == ['scontrol', 'release'] for c in calls) == releases
    assert not (root/'dispatch_receipt.json').exists()
    if mode != 'live':
        assert m.LOCK.exists() and (root/'dispatch_journal.jsonl').exists()
        before = len(calls)
        with pytest.raises(FileExistsError): m.dispatch('28709791', 'b'*64)
        assert len(calls) == before


def test_unknown_parent_rejects_before_io(monkeypatch):
    m = manager()
    monkeypatch.setattr(m, 'run', lambda *a: pytest.fail('no IO for unknown parent'))
    with pytest.raises(ValueError): m.audit('28727368')
    with pytest.raises(ValueError): m.dispatch('28727368', 'a'*64)


@pytest.fixture(scope='module')
def actual_starter():
    m = manager()
    root = CASE/'plus1_preserving_starter'
    assert m.sha(root/'audit_receipt.json') == '8eaac3f12f513caaa24fc6dac72ba124e3b42f00646b2d31bad8c2c5acde1d15'
    receipt = json.loads((root/'audit_receipt.json').read_text())
    return m, root, receipt


def test_actual_successful_cap_is_not_convergence_and_has_77_native_calls(actual_starter):
    from examples.hfo2_fixed_input_factory import CONTRACT
    from scripts.audit_hfo2_clamped_terminal import validate_terminal
    m, root, receipt = actual_starter
    summary = json.loads((root/'vcneb_summary.json').read_text())
    assert m.sha(root/'vcneb_summary.json') == receipt['terminal_summary_sha256']
    validate_terminal('28709790', receipt['scheduler_row'], summary, False, outcome='step-cap')
    with pytest.raises(ValueError): validate_terminal('28709790', receipt['scheduler_row'], summary, False)
    assert receipt['terminal_step'] == 10 and not receipt['ordinary_converged']
    assert receipt['fresh_interior_SCFs'] == receipt['expected_fresh_interior_SCFs'] == 77
    assert receipt['initial_fresh_interior_SCFs'] == 7 and receipt['current_frame_exact_caches_checked'] == 9
    assert len(receipt['runtime_code_sha256']) == 338
    assert receipt['new_DFT_calls'] == receipt['scheduler_mutations'] == receipt['new_independent_chains'] == 0
    assert receipt['all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked'] is True
    assert not receipt['physical_inputs_or_running_source_changed'] and not receipt['FIRE_state_restored']
    assert {(Path(c['path']).parts[-3], Path(c['path']).parts[-2]) for c in receipt['calls']} == {
        (f'image_{i:04d}', f'scf_{j:06d}') for i in range(1, 8) for j in range(11)}
    assert all({k: c['input_sha256'][k] for k in CONTRACT} == CONTRACT for c in receipt['calls'])


def test_actual_portable_nine_raw_SCFS_and_same_boundary_force_replay(actual_starter):
    import numpy as np
    from scripts.analyze_hfo2_clamped_residual import analyze
    from scripts.audit_hfo2_static_replica import audited_results
    m, root, receipt = actual_starter
    observation = root/'observations/step_0010'
    report = json.loads((observation/'observation.json').read_text())
    assert report == receipt['terminal_observation']
    result = analyze(observation)
    assert not result['ordinary_residual_pass']
    assert result['fmax_eV_A'] == pytest.approx(1.325704069665602, abs=1e-12)
    assert report['mechanical_boundary']['fixed_cell_rows'] == [0, 1]
    assert report['mechanical_boundary']['allow_tilt'] is True
    assert report['pressure_GPa'] == 0 and report['climb'] is False
    for row in report['raw_image_evaluations']:
        i = row['image_index']
        raw = observation/'raw'/f'image_{i:04d}'
        assert m.sha(raw/'OUT.ABACUS/running_scf.log') == row['raw_log_sha256']
        assert m.sha(raw/'source_audit.json') == row['audit_sha256']
        assert all(m.sha(raw/n) == row['input_sha256'][n] for n in ('INPUT', 'KPT', 'STRU'))
        parsed = audited_results(raw)
        for key, expected in (('energy', row['energy_eV_cell']), ('forces', row['forces_eV_A']),
                              ('stress', row['stress_ASE_voigt_eV_A3'])):
            assert np.allclose(parsed[key], expected, atol=1e-12, rtol=0)


def test_actual_seed_is_latest_nine_cached_geometry_only_not_historical_best(actual_starter):
    from ase.io import read
    from examples.hfo2_fixed_input_factory import same_ordered_geometry
    m, root, receipt = actual_starter
    seed = root/'seed'
    assert m.sha(seed/'manifest.json') == receipt['seed_manifest_sha256']
    manifest = json.loads((seed/'manifest.json').read_text())
    assert manifest['source_job_id'] == '28709790' and manifest['source_step'] == 10
    assert manifest['current_frame_exact_caches'] == 9 and not manifest['optimizer_state_restored']
    assert not manifest['physical_inputs_changed'] and not manifest['holdout_generated_or_read']
    assert all(m.sha(seed/n) == h for n, h in manifest['files_sha256'].items())
    observation = root/'observations/step_0010'
    assert m.sha(observation/'observation.json') == manifest['source_observation_sha256']
    original = read(observation/'evaluated_chain.traj', index=':')
    current = read(seed/'seed.traj', index=':')
    assert len(original) == len(current) == 9
    assert all(same_ordered_geometry(a, b) for a, b in zip(original, current))


def test_actual_single_admission_matches_exact_helper_and_durable_native_journal():
    m = manager()
    result = json.loads((CASE/'dispatch_receipt.json').read_text())
    journal = CASE/'dispatch_journal.jsonl'
    records = [json.loads(s) for s in journal.read_text().splitlines()]
    assert result['job_id'] == '28737465' and result['parent_job'] == '28709790'
    assert result['total_submissions'] == result['release_writes'] == 1
    assert result['new_independent_chains'] == result['dependency_writes'] == 0
    assert result['executed_helper_sha256'] == m.sha(CASE/'continue_flips_HF.py')
    assert result['script_sha256'] == m.sha(SCRIPT)
    assert result['journal_sha256'] == m.sha(journal)
    assert [r['event'] for r in records] == m.EVENTS
    assert records[0]['committed_slots'] == ['28727368']
    assert records[-1]['committed_slots'] == ['28727368', '28737465']
    assert records[2]['accepted_native_state']['Reason'] == 'JobHeldUser'
    assert records[-1]['native_control']['Reason'] != 'JobHeldUser'


def test_actual_start_proof_is_not_an_electronic_or_force_convergence_claim():
    m = manager()
    report = json.loads((CASE/'startup_observation.json').read_text())
    assert report['source_job_id'] == '28737465'
    assert '28737465|RUNNING|0:0|32|hfacnormal01|iai806|2026-10-12T00:45:12' in report['native_accounting']
    assert '28709791|PENDING|0:0|0|hfacnormal01|iai806' in report['native_accounting']
    assert report['native_first_DSIZE'] == ['32']
    assert report['seed_cache_count'] == 9 and report['none_seed_cache_count'] == 0
    assert report['factory'].endswith(':make_clamped_resume_cached_factory')
    assert 'CheckedFIRE:    0 00:46:02    -9782.063713        1.325704' in report['optimizer_log_raw']
    assert report['fresh_call_audits_at_capture'] == 0 and not report['first_fresh_call_audit_exists']
    assert not report['terminal_summary_exists'] and not report['failure_exists']
    assert report['stderr_bytes'] == 0
    assert report['admission_receipt_sha256'] == m.sha(CASE/'dispatch_receipt.json')
