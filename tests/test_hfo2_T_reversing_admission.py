import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import pytest

PATH = Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/clamped_T_terminal_E075_20261011/release_existing_reversing_once.py'
spec = importlib.util.spec_from_file_location('E075_test', PATH)
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)


def report():
    return dict(job_id=helper.DONE, ordinary_converged=True,
        status='audited_ordinary_converged_not_TS_or_sampling_certificate',
        new_DFT_calls=0, scheduler_mutations=0, physical_inputs_or_running_source_changed=False,
        holdout_generated_or_read=False,
        all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked=True,
        initial_fresh_interior_SCFs=7, terminal_step=6, fresh_interior_SCFs=49,
        terminal_observation=dict(ordinary_residual_pass=True, replayed_fmax_eV_A=.09))


def queue():
    return '\n'.join([helper.LIVE+'|hfo2-G2-flip-E069|RUNNING|hfacnormal01|32',
        helper.TARGET+'|hfo2-G2-E058-PO_flip_T_pattern_reversing|PENDING|hfacnormal01|32',
        *[j+'|hfo2-G2-E060-'+c+'|PENDING|hfacnormal01|32' for j,c in helper.DOWNSTREAM.items()]])


def test_exact_positive_audit_and_spare_slot():
    helper.verified_audit(report())
    assert helper.one_used_slot(queue()+'\n999|pc-unrelated|RUNNING|other|128') == [helper.LIVE]


@pytest.mark.parametrize('field,value', [('ordinary_converged',False), ('job_id','123'),
    ('new_DFT_calls',1), ('scheduler_mutations',1), ('holdout_generated_or_read',True),
    ('initial_fresh_interior_SCFs',0), ('fresh_interior_SCFs',42), ('terminal_step',5)])
def test_missing_or_miscounted_starter_evidence_rejected(field, value):
    with pytest.raises(ValueError):
        helper.verified_audit({**report(),field:value})


@pytest.mark.parametrize('force', [.103118, True, None, float('nan'), -1])
def test_near_threshold_not_a_converged_release(force):
    value = report()
    value['terminal_observation']['replayed_fmax_eV_A'] = force
    with pytest.raises(ValueError):
        helper.verified_audit(value)


@pytest.mark.parametrize('change', ['missing_live','extra_active','unknown','missing_waiter','duplicate'])
def test_capacity_and_waiters_must_be_native_and_exact(change):
    body = queue()
    if change == 'missing_live': body = '\n'.join(body.splitlines()[1:])
    if change == 'extra_active': body = body.replace('PO_to_T|PENDING','PO_to_T|RUNNING')
    if change == 'unknown': body += '\n999|hfo2-unknown|RUNNING|hfacnormal01|32'
    if change == 'missing_waiter': body = '\n'.join(body.splitlines()[:-1])
    if change == 'duplicate': body += '\n'+body.splitlines()[0]
    with pytest.raises(ValueError):
        helper.one_used_slot(body)


def test_exclusive_intent_never_replayed(tmp_path, monkeypatch):
    monkeypatch.setattr(helper, 'ROOT', tmp_path)
    (tmp_path/'admission_journal.jsonl').write_text('{}\n')
    with pytest.raises(FileExistsError):
        helper.release_once('a'*64)


@pytest.fixture
def transaction(tmp_path, monkeypatch):
    root, wait, source, seed = [tmp_path/name for name in ('audit','wait','source','seed')]
    for path in (root, wait, source, seed): path.mkdir()
    monkeypatch.setattr(helper, 'ROOT', root)
    monkeypatch.setattr(helper, 'WAIT', wait)
    common_file = tmp_path/'common.py'
    common_file.write_text('# synthetic nonexecuted dependency\n')
    monkeypatch.setattr(helper, 'COMMON', common_file)
    monkeypatch.setattr(helper, 'COMMON_SHA', helper.sha(common_file))
    terminal = tmp_path/'terminal'
    terminal.mkdir()
    (terminal/'vcneb_summary.json').write_text('{}\n')
    value = {**report(), 'workdir':str(terminal),
        'terminal_summary_sha256':helper.sha(terminal/'vcneb_summary.json')}
    (root/'audit_receipt.json').write_text(json.dumps(value))
    script, binary = source/'production.slurm', source/'binary'
    script.write_text('# synthetic only\n')
    binary.write_text('not executable\n')
    (seed/'manifest.json').write_text('{}\n')
    (seed/'seed.traj').write_bytes(b'synthetic fixture, no physical calculation')
    proof = dict(source_root=str(source), pilot_sha256=helper.sha(script),
        ABACUS_binary_sha256=helper.sha(binary), source_code_sha256={}, reports=[dict(
        channel=helper.CHANNEL, seed_root=str(seed), manifest_sha256=helper.sha(seed/'manifest.json'),
        seed_files_sha256={'seed.traj':helper.sha(seed/'seed.traj')})])
    (wait/'preflight.json').write_text(json.dumps(proof))
    calls = []
    state = dict(ambiguous=False, bad_downstream=False)
    def fields(body):
        return dict(token.split('=',1) for token in body.split())
    def run(argv):
        calls.append(argv)
        if argv[0] == 'sacct':
            job = argv[argv.index('-j')+1]
            return job+'|'+('COMPLETED' if job == helper.DONE else 'RUNNING')+'|0:0|32|hfacnormal01|iai806'
        if argv[0] == 'squeue': return queue()
        if argv[:2] == ['scontrol','release']:
            if state['ambiguous']: raise TimeoutError('possibly accepted native write')
            return ''
        job = argv[-1]
        channel = helper.DOWNSTREAM.get(job,helper.CHANNEL)
        row = dict(JobId=job, JobName=('hfo2-G2-E060-' if job in helper.DOWNSTREAM else 'hfo2-G2-E058-')+channel,
            UserId='iai806(16284)', JobState='PENDING', Priority='0', Reason='JobHeldUser',
            NumCPUs='32', Partition='hfacnormal01', Command=str(script),
            StdOut=str(helper.RUNS/'20261010-varneb-clamped-G2-E060-r1'/(channel+'.slurm.out')))
        if state['bad_downstream'] and job == '28709788': row['Priority'] = '1163'
        return ' '.join(k+'='+v for k,v in row.items())
    def held(body, job, *, required, allowed):
        row = fields(body)
        assert job == helper.TARGET and not required and not allowed
        assert row['JobState'] == 'PENDING' and row['Priority'] == '0'
        return row
    rules = SimpleNamespace(run=run, PILOT=str(script), BINARY=binary, fields=fields, held_waiter=held)
    common = SimpleNamespace(PREFLIGHT_SHA=helper.sha(wait/'preflight.json'), load_rules=lambda:rules,
        pending_seed_preflight=lambda path:dict(mock_only=True,new_DFT_calls=0))
    monkeypatch.setattr(helper, 'load_common', lambda:common)
    return root, calls, state, helper.sha(root/'audit_receipt.json')


def test_full_mock_admission_writes_only_one_release(transaction):
    root, calls, state, digest = transaction
    receipt = helper.release_once(digest)
    assert receipt['release_writes'] == 1 and receipt['dependency_writes'] == 0
    assert [call for call in calls if call[:2] == ['scontrol','release']] == [
        ['scontrol','release',helper.TARGET]]
    assert not any(call[0] == 'sbatch' or call[:2] == ['scontrol','update'] for call in calls)
    with pytest.raises(FileExistsError): helper.release_once(digest)


def test_unprotected_downstream_blocks_actual_release(transaction):
    root, calls, state, digest = transaction
    state['bad_downstream'] = True
    with pytest.raises(ValueError): helper.release_once(digest)
    assert not any(call[:2] == ['scontrol','release'] for call in calls)


def test_ambiguous_actual_release_is_journalled_without_retry(transaction):
    root, calls, state, digest = transaction
    state['ambiguous'] = True
    with pytest.raises(TimeoutError): helper.release_once(digest)
    assert 'reconcile_do_not_repeat' in (root/'admission_journal.jsonl').read_text()
    with pytest.raises(FileExistsError): helper.release_once(digest)
    assert sum(call[:2] == ['scontrol','release'] for call in calls) == 1


def test_real_native_starter_preflight_requires_seven_initial_SCFs():
    from scripts.audit_hfo2_clamped_terminal import fresh_interior_call_count
    from scripts.export_hfo2_clamped_observation import PILOT_SCRIPT_SHA256
    native = json.loads((PATH.parent/'starter_runtime_preflight.json').read_text())
    assert fresh_interior_call_count(6, native, PILOT_SCRIPT_SHA256) == (7,49)
