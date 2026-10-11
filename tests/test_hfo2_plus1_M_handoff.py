"""Mock dispatch and registered-source checks; no DFT/native scheduler calls."""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO = Path(__file__).resolve().parents[1]
CASE = REPO/'benchmarks/hfo2_channels/20261008/clamped_plus1_T_handoff_E086_20261011'
spec = importlib.util.spec_from_file_location('E086_M_handoff', CASE/'release_existing_M_once.py')
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)
spec = importlib.util.spec_from_file_location('E067_fields_only_M',
    REPO/'benchmarks/hfo2_channels/20261008/clamped_flip_continue_E067_20261011/submit_held_once.py')
parser = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parser)


def report():
    return dict(job_id=helper.DONE, workdir=str(helper.WAIT/'PO_to_T/band'),
        ordinary_converged=True, status='audited_ordinary_converged_not_TS_or_sampling_certificate',
        new_DFT_calls=0, scheduler_mutations=0, physical_inputs_or_running_source_changed=False,
        holdout_generated_or_read=False,
        all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked=True,
        initial_fresh_interior_SCFs=7, terminal_step=5, fresh_interior_SCFs=42,
        terminal_observation=dict(source_job_id=helper.DONE, ordinary_residual_pass=True,
            replayed_fmax_eV_A=.09, mechanical_boundary=dict(reference_file_sha256=helper.SUBSTRATE_SHA,
                kind='clamped_plane', allow_tilt=True, fixed_cell_rows=[0, 1], cell_dofs=3)))


def queue(terminal=False):
    rows = [] if terminal else [helper.LIVE+'|hfo2-G2-flip-E069|RUNNING|hfacnormal01|32']
    return '\n'.join(rows+[j+'|hfo2-G2-E060-'+c+'|PENDING|hfacnormal01|32'
                           for j, c in helper.WAITERS.items()])


def held_row(job, pilot):
    channel = helper.WAITERS[job]
    return dict(JobId=job, JobName='hfo2-G2-E060-'+channel, JobState='PENDING',
        Priority='0', Reason='JobHeldUser', UserId='iai806(16284)', Partition='hfacnormal01',
        NumCPUs='32', NumTasks='1', **{'CPUs/Task': '32'}, Command=pilot,
        TimeLimit='04:00:00', StdOut=str(helper.WAIT/(channel+'.slurm.out')),
        StdErr=str(helper.WAIT/(channel+'.slurm.err')),
        Dependency='(null)' if job == helper.TARGET else 'afterok:'+helper.TARGET+'(unfulfilled)')


def body(row):
    return ' '.join(k+'='+v for k, v in row.items())


def test_completed_T_dependency_can_disappear_but_M_is_still_protected():
    rules = SimpleNamespace(fields=parser.fields, dependencies=parser.dependencies, PILOT='registered.slurm')
    helper.verified_parent(report())
    for job in helper.WAITERS:
        row = held_row(job, rules.PILOT)
        helper.held_waiter(body(row), job, rules)
        if job != helper.TARGET:
            row['Dependency'] += ',afterok:'+helper.DONE+'(fulfilled)'
            helper.held_waiter(body(row), job, rules)
    assert helper.spare_slot(queue()+'\n999|pc-unrelated|RUNNING|other|128', live_terminal=False) == [helper.LIVE]
    assert helper.spare_slot(queue(True), live_terminal=True) == []


@pytest.mark.parametrize('field,value', [('job_id','123'), ('ordinary_converged',False),
    ('fresh_interior_SCFs',35), ('terminal_step',True), ('initial_fresh_interior_SCFs',0),
    ('holdout_generated_or_read',True), ('new_DFT_calls',1), ('scheduler_mutations',1),
    ('workdir','another/band')])
def test_terminal_receipt_cannot_be_relabelled(field, value):
    with pytest.raises(ValueError): helper.verified_parent({**report(), field:value})


@pytest.mark.parametrize('force', [.1001, -1, True, None, float('nan'), float('inf')])
def test_force_criterion_is_not_rounded_or_inferred_from_exit_zero(force):
    r = report()
    r['terminal_observation']['replayed_fmax_eV_A'] = force
    with pytest.raises(ValueError): helper.verified_parent(r)


@pytest.mark.parametrize('field,value', [('reference_file_sha256','a'*64), ('allow_tilt',False),
    ('fixed_cell_rows',[0]), ('cell_dofs',6), ('kind','free_cell')])
def test_same_registered_plus1_boundary_is_required(field, value):
    r = report()
    r['terminal_observation']['mechanical_boundary'][field] = value
    with pytest.raises(ValueError): helper.verified_parent(r)


@pytest.mark.parametrize('field,value', [('Priority','16'), ('NumTasks','32'),
    ('NumCPUs','64'), ('CPUs/Task','1'), ('StdErr','wrong'), ('TimeLimit','12:00:00'),
    ('Dependency','afterok:28723655'), ('UserId','other(1)')])
def test_held_identity_cannot_be_widened(field, value):
    rules = SimpleNamespace(fields=parser.fields, dependencies=parser.dependencies, PILOT='registered.slurm')
    with pytest.raises(ValueError):
        helper.held_waiter(body({**held_row(helper.TARGET, rules.PILOT),field:value}),helper.TARGET,rules)


def test_pending_M_dependency_cannot_disappear_from_flips():
    rules = SimpleNamespace(fields=parser.fields, dependencies=parser.dependencies, PILOT='registered.slurm')
    with pytest.raises(ValueError):
        helper.held_waiter(body({**held_row('28709790',rules.PILOT),'Dependency':'(null)'}),'28709790',rules)


@pytest.mark.parametrize('change', ['parent_still_active','extra_active','missing_waiter','unknown','duplicate'])
def test_two_slot_and_registered_queue_gate(change):
    q = queue()
    if change == 'parent_still_active': q += '\n'+helper.DONE+'|hfo2-G2-T|RUNNING|hfacnormal01|32'
    if change == 'extra_active': q = q.replace('PO_to_M|PENDING','PO_to_M|RUNNING')
    if change == 'missing_waiter': q = '\n'.join(q.splitlines()[:-1])
    if change == 'unknown': q += '\n999|hfo2-unknown|RUNNING|hfacnormal01|32'
    if change == 'duplicate': q += '\n'+q.splitlines()[0]
    with pytest.raises(ValueError): helper.spare_slot(q,live_terminal=False)


@pytest.fixture
def transaction(tmp_path, monkeypatch):
    root, wait, source, seed = [tmp_path/n for n in ('audit','wait','source','seed')]
    for p in (root,wait,source,seed): p.mkdir()
    monkeypatch.setattr(helper,'ROOT',root)
    monkeypatch.setattr(helper,'WAIT',wait)
    terminal = wait/'PO_to_T/band'
    terminal.mkdir(parents=True)
    (terminal/'vcneb_summary.json').write_text('{}\n')
    r = {**report(), 'terminal_summary_sha256':helper.sha(terminal/'vcneb_summary.json')}
    (root/'audit_receipt.json').write_text(json.dumps(r))
    script, binary, code = [source/n for n in ('production.slurm','binary','helper.py')]
    for p in (script,binary,code): p.write_text('# synthetic non-executable fixture\n')
    (seed/'manifest.json').write_text('{}\n')
    (seed/'substrate.vasp').write_text('synthetic, not a physical structure\n')
    monkeypatch.setattr(helper,'SUBSTRATE_SHA',helper.sha(seed/'substrate.vasp'))
    r['terminal_observation']['mechanical_boundary']['reference_file_sha256'] = helper.SUBSTRATE_SHA
    (root/'audit_receipt.json').write_text(json.dumps(r))
    proof = dict(condition='strain_p0100',source_root=str(source),pilot_sha256=helper.sha(script),
        ABACUS_binary_sha256=helper.sha(binary),source_code_sha256={'helper.py':helper.sha(code)},
        reports=[dict(channel=helper.WAITERS[helper.TARGET],condition='strain_p0100',seed_root=str(seed),
            manifest_sha256=helper.sha(seed/'manifest.json'),
            seed_files_sha256={'substrate.vasp':helper.sha(seed/'substrate.vasp')})])
    (wait/'preflight.json').write_text(json.dumps(proof))
    monkeypatch.setattr(helper,'PREFLIGHT_SHA',helper.sha(wait/'preflight.json'))
    calls, state = [], dict(ambiguous=False,race=False,reads=0,companion_terminal=False,parent_running=False)
    def run(argv):
        calls.append(argv)
        if argv[0] == 'sacct':
            j = argv[argv.index('-j')+1]
            s = ('COMPLETED' if state['companion_terminal'] else 'RUNNING') if j == helper.LIVE else (
                'RUNNING' if state['parent_running'] else 'COMPLETED')
            return j+'|'+s+'|0:0|32|hfacnormal01|iai806'
        if argv[0] == 'squeue':
            state['reads'] += 1
            q = queue(state['companion_terminal'])
            return q.replace('PO_to_M|PENDING','PO_to_M|RUNNING') if state['race'] and state['reads']>1 else q
        if argv[:2] == ['scontrol','release']:
            if state['ambiguous']: raise TimeoutError('possibly accepted native write')
            return ''
        return body(held_row(argv[-1],str(script)))
    rules = SimpleNamespace(run=run,PILOT=str(script),BINARY=binary,
        fields=parser.fields,dependencies=parser.dependencies)
    monkeypatch.setattr(helper,'load_common',lambda:SimpleNamespace(load_rules=lambda:rules,
        pending_seed_preflight=lambda p:dict(mock_only=True,new_DFT_calls=0)))
    return root,wait,source,calls,state,helper.sha(root/'audit_receipt.json')


def test_read_only_preflight_never_checks_or_claims_parent_convergence(transaction):
    root,wait,source,calls,state,digest = transaction
    state['parent_running'] = True
    _, result = helper.pending_inputs()
    assert result['admission_or_parent_convergence_claimed'] is False
    assert result['new_DFT_calls'] == result['scheduler_mutations'] == 0
    assert not (root/'admission_journal.jsonl').exists()
    assert all(c[:3] == ['scontrol','show','job'] for c in calls)


@pytest.mark.parametrize('terminal', [False,True])
def test_only_existing_M_is_released_once_with_zero_dependency_writes(transaction,terminal):
    root,wait,source,calls,state,digest = transaction
    state['companion_terminal'] = terminal
    result = helper.release_once(digest)
    assert result['release_writes'] == 1 and result['dependency_writes'] == result['new_submissions'] == 0
    assert not result['startup_claimed']
    assert [c for c in calls if c[:2] == ['scontrol','release']] == [['scontrol','release',helper.TARGET]]
    assert not any(c[0] == 'sbatch' or c[:2] == ['scontrol','update'] for c in calls)
    events = [json.loads(l)['event'] for l in (root/'admission_journal.jsonl').read_text().splitlines()]
    assert events == ['admission_verified','release_intent','release_accepted']
    with pytest.raises(FileExistsError): helper.release_once(digest)


@pytest.mark.parametrize('failure', ['race','parent_running','ambiguous'])
def test_real_state_or_ambiguous_reply_cannot_trigger_a_second_write(transaction,failure):
    root,wait,source,calls,state,digest = transaction
    state[failure] = True
    with pytest.raises(TimeoutError if failure == 'ambiguous' else ValueError): helper.release_once(digest)
    releases = sum(c[:2] == ['scontrol','release'] for c in calls)
    assert releases == (1 if failure == 'ambiguous' else 0)
    if (root/'admission_journal.jsonl').exists():
        with pytest.raises(FileExistsError): helper.release_once(digest)
        assert 'reconcile_do_not_repeat' in (root/'admission_journal.jsonl').read_text()


@pytest.mark.parametrize('change', ['runtime','production_exists','seed','preflight'])
def test_changed_original_pending_inputs_block_any_dispatch(transaction,change):
    root,wait,source,calls,state,digest = transaction
    if change == 'runtime': (source/'helper.py').write_text('# changed\n')
    if change == 'production_exists': (wait/helper.WAITERS[helper.TARGET]/'band').mkdir(parents=True)
    if change == 'seed': (source.parent/'seed/substrate.vasp').write_text('changed\n')
    if change == 'preflight': (wait/'preflight.json').write_text('{}\n')
    with pytest.raises(ValueError): helper.release_once(digest)
    assert not any(c[:2] == ['scontrol','release'] for c in calls)
    assert not (root/'admission_journal.jsonl').exists()


def test_registered_plus1_reference_is_not_zero_strain_or_holdout():
    p = REPO/'benchmarks/hfo2_channels/20261008/clamped_G2_training_queue_E060_20261010/preflight.json'
    assert helper.sha(p) == helper.PREFLIGHT_SHA
    proof = json.loads(p.read_text(encoding='utf-8'))
    assert proof['condition'] == 'strain_p0100' and len(proof['source_code_sha256']) == 338
    assert all(r['seed_files_sha256']['substrate.vasp'] == helper.SUBSTRATE_SHA for r in proof['reports'])
    assert [r['channel'] for r in proof['reports'][1:]] == list(helper.WAITERS.values())
