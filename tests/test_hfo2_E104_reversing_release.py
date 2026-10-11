"""Reviewed slot-policy adapter; scheduler fixtures are not material evidence."""
import importlib.util
import json
from pathlib import Path
import shlex
import shutil

import pytest

REPO = Path(__file__).resolve().parents[1]
CASE = REPO/'benchmarks/hfo2_channels/20261008/clamped_reversing_release_E104_20261012'
E100 = REPO/'benchmarks/hfo2_channels/20261008/clamped_held_flips_E100_20261011'
E103 = REPO/'benchmarks/hfo2_channels/20261008/clamped_flip_continue_E103_20261012'


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def adapter():
    return load(CASE/'release_reversing_HF.py', 'E104_fixture_adapter')


def components(tmp_path, monkeypatch):
    m = adapter()
    legacy = m.load_pinned(E100/'release_existing_once.py', m.RELEASER_SHA, 'E104_fixture_E100')
    manager = m.load_pinned(E103/'continue_flips_HF.py', m.MANAGER_SHA, 'E104_fixture_E103')
    monkeypatch.setattr(legacy, 'ROOT', tmp_path/'E100')
    legacy.ROOT.mkdir()
    monkeypatch.setattr(manager, 'LOCK', legacy.ROOT/'admission.lock')
    monkeypatch.setattr(manager, 'ROOT', tmp_path/'E103')
    monkeypatch.setattr(manager, 'SCRIPT', REPO/'cluster/hf_hfo2_clamped_continue_E103_20261012.slurm')
    root = manager.ROOT/'28709790'
    (root/'terminal').mkdir(parents=True)
    for n in ('dispatch_receipt.json', 'dispatch_journal.jsonl'):
        shutil.copyfile(E103/n, root/n)
    shutil.copyfile(E103/'plus1_preserving_starter/audit_receipt.json', root/'terminal/audit_receipt.json')
    monkeypatch.setattr(m, 'ROOT', tmp_path/'E104')
    m.ROOT.mkdir()
    loaded = []
    def loader(path, digest, name):
        loaded.append((path, digest, name))
        if path == m.RELEASER and digest == m.RELEASER_SHA: return legacy
        if path == m.MANAGER and digest == m.MANAGER_SHA: return manager
        pytest.fail('only the two pinned originals may be imported')
    monkeypatch.setattr(m, 'load_pinned', loader)
    return m, legacy, manager, loaded


def native_child(m, manager, *, state='RUNNING', change=None):
    admitted = manager.accepted(m.PARENT)
    p = manager.paths(m.PARENT)
    receipt = json.loads((p['terminal']/'audit_receipt.json').read_text())
    env = manager.environment(m.PARENT, receipt, admitted['audit_receipt_sha256'])
    if change == 'override': env['ECUTWFC'] = '120'
    if change in env: env[change] = 'different'
    submit = 'sbatch '+shlex.quote('--export=ALL,'+','.join(k+'='+v for k,v in env.items()))+' '+shlex.quote(str(manager.SCRIPT))
    if change == 'duplicate': submit += ' --export=ALL'
    values = [admitted['job_id'], state, 'iai806', '32', 'hfacnormal01', 'hfo2-G2-E103-'+m.PARENT,
        str(p['stdout']), str(p['stderr']), submit]
    if isinstance(change, int): values[change] = 'different'
    return '|'.join(values)


def queue(*, include_M=False, released=False, child_pending=False):
    body = '28737465|hfo2-G2-E103-28709790|'+('PENDING' if child_pending else 'RUNNING')+'|hfacnormal01|32|'+('(Resources)' if child_pending else 'node25')
    body += '\n28709791|hfo2-G2-E060-PO_flip_T_pattern_reversing|PENDING|hfacnormal01|32|'+('(Priority)' if released else '(JobHeldUser)')
    if include_M: body += '\n28727368|hfo2-G2-E091-28709789|RUNNING|hfacnormal01|32|node6'
    return body


@pytest.mark.parametrize('change', [None, 0, 2, 3, 4, 5, 6, 7, 'SOURCE_ROOT', 'WORKDIR', 'SEED_ROOT',
    'PARENT_JOB_ID', 'PARENT_AUDIT_PATH', 'PARENT_AUDIT_SHA256', 'SEED_MANIFEST_SHA256', 'override', 'duplicate'])
def test_exact_accepted_native_child_ancestry(tmp_path, monkeypatch, change):
    m, _, manager, _ = components(tmp_path, monkeypatch)
    body = native_child(m, manager, change=change)
    if change is None: m.validate_child(body, manager, manager.accepted(m.PARENT))
    else:
        with pytest.raises(ValueError): m.validate_child(body, manager, manager.accepted(m.PARENT))


@pytest.mark.parametrize('state', ['FAILED', 'TIMEOUT', 'CANCELLED', 'PENDING'])
def test_unexpected_native_child_state_rejects(tmp_path, monkeypatch, state):
    m, _, manager, _ = components(tmp_path, monkeypatch)
    with pytest.raises(ValueError): m.validate_child(native_child(m, manager, state=state), manager, manager.accepted(m.PARENT))


def test_original_policy_rejects_new_child_and_adapter_preserves_ceiling(tmp_path, monkeypatch):
    m, legacy, manager, _ = components(tmp_path, monkeypatch)
    with pytest.raises(ValueError): legacy.committed_slots(queue(), require_spare=True)
    monkeypatch.setattr(legacy, 'run', lambda *a: native_child(m, manager))
    new, admitted, raw = m.reviewed_control()
    assert new is legacy and admitted['job_id'] == '28737465' and '|RUNNING|' in raw
    assert new.committed_slots(queue(), require_spare=True) == ['28737465']
    assert new.committed_slots(queue(child_pending=True), require_spare=True) == ['28737465']
    with pytest.raises(ValueError): new.committed_slots(queue(include_M=True), require_spare=True)
    assert len(new.committed_slots(queue(released=True), require_spare=False)) == 2
    with pytest.raises(ValueError): new.committed_slots(queue(released=True, include_M=True), require_spare=False)


def test_changed_pinned_helper_rejects_before_import(tmp_path):
    m = adapter()
    path = tmp_path/'changed.py'
    path.write_text('raise AssertionError("must never import")\n')
    with pytest.raises(ValueError): m.load_pinned(path, 'a'*64, 'forbidden')


def fixture(tmp_path, monkeypatch, mode='success'):
    m, legacy, manager, _ = components(tmp_path, monkeypatch)
    tools = load(REPO/'tests/test_hfo2_E100_held_flips.py', 'E104_original_control_fixture')
    calls, released = [], []
    def audit(path, digest, job):
        calls.append(('audit', str(path), digest, job))
        if job == '28727368' and mode == 'cap': raise ValueError('M must be ordinary-converged')
    monkeypatch.setattr(legacy, 'verified_terminal', audit)
    def preflight(job, require_spare=False):
        assert job == m.TARGET
        committed = legacy.committed_slots(queue(include_M=mode=='two_live', released=bool(released)), require_spare=require_spare)
        return {'synthetic_not_material_seed_proof': True, 'committed_study_slots': committed}
    monkeypatch.setattr(legacy, 'preflight', preflight)
    def run(argv):
        calls.append(tuple(argv))
        if argv[0] == 'sacct': return native_child(m, manager)
        if argv[0] == 'squeue': return queue(include_M=mode=='two_live', released=bool(released))
        if argv[:3] == ['scontrol', 'show', 'job']:
            return tools.control(legacy, job=m.TARGET, Priority='1' if released else '0', Reason='Priority' if released else 'JobHeldUser')
        if argv == ['scontrol', 'release', m.TARGET]:
            if mode == 'timeout': raise TimeoutError('ambiguous release transport')
            if mode != 'stale': released.append(True)
            return ''
        pytest.fail(str(argv))
    monkeypatch.setattr(legacy, 'run', run)
    return m, legacy, manager, calls


def test_readiness_two_live_not_release_or_convergence(tmp_path, monkeypatch):
    m, legacy, _, calls = fixture(tmp_path, monkeypatch, 'two_live')
    proof = m.preflight()
    assert len(proof['original_held_seed_proof']['committed_study_slots']) == 2
    assert proof['new_DFT_calls'] == proof['scheduler_writes'] == 0
    assert not (legacy.ROOT/'28709791.journal.jsonl').exists()
    assert not any(c[:2] == ('scontrol', 'release') for c in calls)


def test_success_reuses_original_single_release_and_explains_policy_substitution(tmp_path, monkeypatch):
    m, legacy, _, calls = fixture(tmp_path, monkeypatch)
    result = m.release_once('M_actual_audit_fixture_digest')
    assert result['target_job'] == '28709791' and result['slot_policy_process_substitution_explicit']
    assert result['new_submissions'] == result['dependency_writes'] == result['new_independent_chains'] == 0
    assert result['release_writes'] == 1 and not result['old_files_or_live_source_changed']
    assert calls[1] == ('audit', str(m.M_AUDIT), 'M_actual_audit_fixture_digest', '28727368')
    assert sum(c[:2] == ('scontrol', 'release') for c in calls) == 1
    assert not (legacy.ROOT/'admission.lock').exists()
    records = [json.loads(s) for s in (legacy.ROOT/'28709791.journal.jsonl').read_text().splitlines()]
    assert records[-1]['committed_slots'] == ['28709791', '28737465']
    before = len(calls)
    with pytest.raises(FileExistsError): m.release_once('same')
    assert len(calls) == before


@pytest.mark.parametrize('mode', ['cap', 'two_live', 'timeout', 'stale'])
def test_nonordinary_full_capacity_or_ambiguous_reply_never_released_twice(tmp_path, monkeypatch, mode):
    m, legacy, _, calls = fixture(tmp_path, monkeypatch, mode)
    with pytest.raises((ValueError, TimeoutError)): m.release_once('fixture')
    writes = sum(c[:2] == ('scontrol', 'release') for c in calls)
    assert writes == (1 if mode in ('timeout', 'stale') else 0)
    assert not (m.ROOT/'release_driver_receipt.json').exists()
    if writes:
        assert (legacy.ROOT/'admission.lock').exists()
        with pytest.raises(FileExistsError): m.release_once('fixture')
        assert sum(c[:2] == ('scontrol', 'release') for c in calls) == 1


def test_actual_read_only_HF_seed_and_accepted_child_preflight():
    m = adapter()
    path = CASE/'preflight.json'
    assert m.sha(path) == 'ba231ea7d47adaf1512f13cec8fa613351ca2ee274a84fbb8a8387dfdbf07f44'
    r = json.loads(path.read_text())
    assert r['target_job'] == '28709791' and r['scheduler_writes'] == r['new_DFT_calls'] == 0
    assert r['executed_adapter_sha256'] == m.sha(CASE/'release_reversing_HF.py')
    assert r['original_release_helper_sha256'] == m.RELEASER_SHA
    assert r['accepted_slot_policy_sha256'] == m.MANAGER_SHA
    proof = r['original_held_seed_proof']
    assert proof['original_runtime_files_checked'] == 338
    assert proof['committed_study_slots'] == ['28727368', '28737465']
    assert proof['seed_manifest_sha256'] == 'cbbb66a64a1a3f98503929152b35800d42b04cfc9a523cab63b83627f0251f08'
    assert proof['native_accounting'] == '28709791|PENDING|0:0|0|hfacnormal01|iai806'
    assert proof['native_control']['Reason'] == 'JobHeldUser' and proof['new_submissions'] == 0
    assert '|RUNNING|iai806|32|hfacnormal01|hfo2-G2-E103-28709790|' in r['accepted_child_native_source_row']
    assert r['accepted_preserving_child']['job_id'] == '28737465'


@pytest.fixture(scope='module')
def actual_M_terminal():
    m = adapter()
    root = CASE/'plus1_M_terminal'
    receipt = json.loads((root/'audit_receipt.json').read_text())
    assert m.sha(root/'audit_receipt.json') == 'b86cdb72eab9e3d8ff2125e16e6fbb2d0954b34d54c0f05099a581c066108389'
    return m, root, receipt


def test_actual_M_ordinary_completion_all_native_calls_and_fixed_contract(actual_M_terminal):
    from examples.hfo2_fixed_input_factory import CONTRACT
    from scripts.audit_hfo2_clamped_terminal import validate_terminal
    m, root, r = actual_M_terminal
    summary = json.loads((root/'vcneb_summary.json').read_text())
    assert m.sha(root/'vcneb_summary.json') == r['terminal_summary_sha256']
    validate_terminal('28727368', r['scheduler_row'], summary, False)
    assert r['ordinary_converged'] and r['terminal_step'] == 22
    assert r['fresh_interior_SCFs'] == r['expected_fresh_interior_SCFs'] == 154
    assert r['initial_fresh_interior_SCFs'] == 0
    assert len(r['runtime_code_sha256']) == 338
    assert r['all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked']
    assert r['new_DFT_calls'] == r['scheduler_mutations'] == 0
    assert not r['physical_inputs_or_running_source_changed']
    assert not r['continuous_TS_or_barrier_error_bound_certified']
    assert all({k: c['input_sha256'][k] for k in CONTRACT} == CONTRACT for c in r['calls'])


def test_actual_M_portable_nine_raw_SCFS_and_same_boundary_force_replay(actual_M_terminal):
    import numpy as np
    from scripts.analyze_hfo2_clamped_residual import analyze
    from scripts.audit_hfo2_static_replica import audited_results
    m, root, r = actual_M_terminal
    observation = root/'observations/step_0022'
    report = json.loads((observation/'observation.json').read_text())
    assert report == r['terminal_observation']
    result = analyze(observation)
    assert result['ordinary_residual_pass']
    assert result['fmax_eV_A'] <= .10
    assert result['fmax_eV_A'] == pytest.approx(report['replayed_fmax_eV_A'], abs=1e-12)
    assert report['mechanical_boundary']['fixed_cell_rows'] == [0, 1]
    assert report['mechanical_boundary']['allow_tilt'] is True
    assert report['pressure_GPa'] == 0 and report['climb'] is False
    assert len(report['raw_image_evaluations']) == 9
    for row in report['raw_image_evaluations']:
        raw = observation/'raw'/f"image_{row['image_index']:04d}"
        assert m.sha(raw/'OUT.ABACUS/running_scf.log') == row['raw_log_sha256']
        assert m.sha(raw/'source_audit.json') == row['audit_sha256']
        assert all(m.sha(raw/n) == row['input_sha256'][n] for n in ('INPUT', 'KPT', 'STRU'))
        parsed = audited_results(raw)
        for key, expected in (('energy', row['energy_eV_cell']), ('forces', row['forces_eV_A']),
                              ('stress', row['stress_ASE_voigt_eV_A3'])):
            assert np.allclose(parsed[key], expected, atol=1e-12, rtol=0)


def test_actual_one_original_reversing_release_preserves_shared_journal_and_policy():
    m = adapter()
    r = json.loads((CASE/'release_driver_receipt.json').read_text())
    old = json.loads((CASE/'28709791.receipt.json').read_text())
    journal = CASE/'28709791.journal.jsonl'
    records = [json.loads(s) for s in journal.read_text().splitlines()]
    assert r['original_release_result'] == old
    assert r['M_terminal_audit_sha256'] == m.sha(CASE/'plus1_M_terminal/audit_receipt.json')
    assert r['target_job'] == old['job_id'] == '28709791'
    assert r['executed_adapter_sha256'] == m.sha(CASE/'release_reversing_HF.py')
    assert r['original_release_helper_sha256'] == old['executed_helper_sha256'] == m.RELEASER_SHA
    assert r['accepted_slot_policy_sha256'] == m.MANAGER_SHA
    assert old['journal_sha256'] == m.sha(journal)
    assert r['slot_policy_process_substitution_explicit']
    assert r['release_writes'] == old['release_writes'] == 1
    assert r['new_submissions'] == r['dependency_writes'] == r['new_independent_chains'] == 0
    assert [v['event'] for v in records] == ['verified_original_held_target',
        'one_existing_release_intent', 'one_existing_release_accepted',
        'released_slot_independently_corroborated']
    assert records[0]['proof']['committed_study_slots'] == ['28737465']
    assert records[-1]['committed_slots'] == ['28709791', '28737465']
    assert records[-1]['native_control']['Reason'] != 'JobHeldUser'


def test_actual_reversing_start_uses_two_cached_endpoints_and_32_native_ranks():
    from examples.hfo2_fixed_input_factory import CONTRACT
    m = adapter()
    r = json.loads((CASE/'startup_observation.json').read_text())
    assert '28709791|RUNNING|0:0|32|hfacnormal01|iai806|2026-10-12T01:27:26' in r['native_accounting']
    assert '28727368|COMPLETED|0:0|32|hfacnormal01|iai806' in r['native_accounting']
    assert '28737465|RUNNING|0:0|32|hfacnormal01|iai806' in r['native_accounting']
    assert r['seed_cache_count'] == 2 and r['none_seed_cache_count'] == 7
    assert r['native_first_DSIZE'] == ['32'] and r['fresh_call_audits_at_capture'] == 1
    assert r['first_fresh_call_elapsed_seconds'] == pytest.approx(125.57138762343675)
    assert {k:r['first_fresh_call_input_sha256'][k] for k in CONTRACT} == CONTRACT
    assert r['adapter_release_receipt_sha256'] == m.sha(CASE/'release_driver_receipt.json')
    assert r['original_release_receipt_sha256'] == m.sha(CASE/'28709791.receipt.json')
    assert r['original_journal_sha256'] == m.sha(CASE/'28709791.journal.jsonl')
    assert not r['original_shared_admission_lock_exists']
    assert not r['failure_exists'] and not r['terminal_summary_exists']
    assert r['stderr_bytes'] == r['new_readonly_DFT_calls'] == r['scheduler_writes_for_capture'] == 0
    assert r['startup_or_release_is_not_SCF_completion_or_convergence']
