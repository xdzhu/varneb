"""Validate the exact finite same-chain recipe, without Slurm or DFT calls."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys

import pytest
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read
from examples.hfo2_fixed_input_factory import CONTRACT, read_fixed_hfo2_stru, same_ordered_geometry
from scripts.audit_hfo2_static_replica import audited_results
from scripts import export_hfo2_clamped_observation as observer
from scripts.analyze_hfo2_clamped_residual import analyze
from vcneb import VCNEB, clamped_plane_vcneb_boundary

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('flip_terminal_audit',
    ROOT/'benchmarks/hfo2_channels/20261008/clamped_flip_continue_E067_20261011/audit_HF.py')
terminal = importlib.util.module_from_spec(spec)
spec.loader.exec_module(terminal)
dispatch_spec = importlib.util.spec_from_file_location('submit_held_once',
    ROOT/'benchmarks/hfo2_channels/20261008/clamped_flip_continue_E067_20261011/submit_held_once.py')
dispatch = importlib.util.module_from_spec(dispatch_spec)
dispatch_spec.loader.exec_module(dispatch)
sys.modules['submit_held_once'] = dispatch
handoff_spec = importlib.util.spec_from_file_location('flip_handoff',
    ROOT/'benchmarks/hfo2_channels/20261008/clamped_flip_continue_E067_20261011/handoff_release.py')
handoff = importlib.util.module_from_spec(handoff_spec)
handoff_spec.loader.exec_module(handoff)
SCRIPT = ROOT/'cluster/hf_hfo2_clamped_flip_continue_E067_20261011.slurm'
SOURCE = SCRIPT.read_text()
CODE = re.search(r"python -c \\\n+ '(.+?)' \"\$\{seed_root\}/manifest.json\"", SOURCE, flags=re.DOTALL).group(1)


@pytest.fixture
def seed(tmp_path):
    names = ('initial.vasp', 'final.vasp', 'seed.traj', 'substrate.vasp', 'factory_parameters.json')
    for name in names:
        (tmp_path/name).write_text('synthetic '+name)
    manifest = {'status': 'audited_G2_geometry_continuation_prepared',
        'source_job_id': '28661019', 'source_step': 20, 'current_frame_exact_caches': 9,
        'new_independent_chains': 0, 'new_DFT_calls': 0, 'physical_inputs_changed': False,
        'climb': False, 'holdout_generated_or_read': False, 'optimizer_state_restored': False,
        'pressure_GPa': 0, 'ordinary_fmax_eV_A': .10, 'n_total_images': 9,
        'n_fixed_endpoints': 2, 'n_active_images': 7, 'cell_scale_A': observer.CELL_SCALE,
        'files_sha256': {name: observer.sha256(tmp_path/name) for name in names}}
    return tmp_path, manifest


def validate(seed, *, digest=None, optimize=False):
    root, manifest = seed
    path = root/'manifest.json'
    path.write_text(json.dumps(manifest))
    env = dict(os.environ)
    if optimize:
        env['PYTHONOPTIMIZE'] = '2'
    return subprocess.run([sys.executable, '-c', CODE, str(path),
                           observer.sha256(path) if digest is None else digest],
                          capture_output=True, text=True, timeout=10, env=env)


def test_registered_script_and_runner_flags_match_original_unchanged_M_recipe():
    assert observer.sha256(SCRIPT) == observer.FLIP_CONTINUE_SCRIPT_SHA256
    assert observer.REGISTERED_FACTORIES[observer.sha256(SCRIPT)].endswith(':make_clamped_resume_cached_factory')
    old = (ROOT/'cluster/hf_hfo2_clamped_M_continue_E062_20261011.slurm').read_text()
    assert SOURCE.split('-m vcneb.material_runner',1)[1] == old.split('-m vcneb.material_runner',1)[1]
    assert '#SBATCH -p hfacnormal01' in SOURCE and '#SBATCH --cpus-per-task=32' in SOURCE
    assert '#SBATCH --time=12:00:00' in SOURCE and 'source-fixed' in SOURCE
    assert 'sbatch' not in SOURCE and 'srun' not in SOURCE and 'scontrol' not in SOURCE


def test_valid_seed_replays_scale_without_any_external_calculation(seed):
    reply = validate(seed)
    assert reply.returncode == 0 and float(reply.stdout) == observer.CELL_SCALE


@pytest.mark.parametrize('key,value', [
    ('source_job_id','28661020'), ('source_step',10), ('current_frame_exact_caches',2),
    ('new_independent_chains',1), ('new_DFT_calls',1), ('physical_inputs_changed',True),
    ('climb',True), ('holdout_generated_or_read',True), ('optimizer_state_restored',True),
    ('pressure_GPa',1), ('ordinary_fmax_eV_A',.2), ('n_total_images',7),
    ('n_active_images',5), ('cell_scale_A',6), ('status','unaudited'),
])
def test_changed_seed_contract_refused_even_with_assertions_disabled(seed, key, value):
    seed[1][key] = value
    reply = validate(seed, optimize=True)
    assert reply.returncode != 0 and 'seed contract differs' in reply.stderr


@pytest.mark.parametrize('name', ['initial.vasp','final.vasp','seed.traj','substrate.vasp','factory_parameters.json'])
def test_each_pinned_seed_file_change_refused(seed, name):
    (seed[0]/name).write_text('changed')
    reply = validate(seed)
    assert reply.returncode != 0 and 'seed file changed' in reply.stderr


def test_manifest_and_file_inventory_tampering_refused(seed):
    reply = validate(seed, digest='0'*64)
    assert reply.returncode != 0 and 'seed manifest changed' in reply.stderr
    seed[1]['files_sha256']['../outside'] = hashlib.sha256(b'anything').hexdigest()
    reply = validate(seed)
    assert reply.returncode != 0 and 'seed file set changed' in reply.stderr


@pytest.mark.parametrize('row,status,converged,failure', [
    ('28661019|RUNNING|0:0|32|hfacnormal01|iai806','max_steps_reached',False,False),
    ('28661019|TIMEOUT|0:0|32|hfacnormal01|iai806','max_steps_reached',False,False),
    ('28661019|COMPLETED|1:0|32|hfacnormal01|iai806','max_steps_reached',False,False),
    ('28661019|COMPLETED|0:0|64|hfacnormal01|iai806','max_steps_reached',False,False),
    ('28661019|COMPLETED|0:0|32|hfacnormal01|iai806','converged',True,False),
    ('28661019|COMPLETED|0:0|32|hfacnormal01|iai806','max_steps_reached',False,True),
    ('28661020|COMPLETED|0:0|32|hfacnormal01|iai806','max_steps_reached',False,False),
])
def test_terminal_audit_cannot_resume_live_failed_converged_or_wrong_parent(row,status,converged,failure):
    with pytest.raises(ValueError):
        terminal.validate_terminal(row, {'status':status,'converged':converged}, failure)


def test_completed_but_unconverged_segment_is_not_mislabeled_as_convergence():
    terminal.validate_terminal('28661019|COMPLETED|0:0|32|hfacnormal01|iai806',
                               {'status':'max_steps_reached','converged':False}, False)


def waiter_body(job, dep='afterok:28722320(unfulfilled)', state='PENDING'):
    channel = dispatch.WAITERS[job]
    return (f'JobId={job} JobName=hfo2-G2-E058-{channel} UserId=iai806(16284) '
            f'JobState={state} Partition=hfacnormal01 NumCPUs=32 Priority=0 Reason=JobHeldUser '
            f'Dependency={dep} Command={dispatch.PILOT} '
            f'StdOut={dispatch.WAIT_ROOT}/{channel}.slurm.out StdErr={dispatch.WAIT_ROOT}/{channel}.slurm.err')


@pytest.mark.parametrize('dep', ['afterany:28722320', 'afterok:28722320?afterok:999',
    'afterok:999', 'afterok:28722320,afterok:123', '(null)'])
def test_held_successor_requires_original_live_parent_without_unrelated_dependencies(dep):
    with pytest.raises(ValueError):
        dispatch.held_waiter(waiter_body('28692775',dep), '28692775',
                             required={dispatch.OTHER},allowed={dispatch.OTHER,dispatch.PARENT})


def test_handoff_acknowledges_stale_read_without_repeated_write():
    target = {dispatch.OTHER,'999'}
    replies = iter([waiter_body('28692775'), waiter_body('28692775','afterok:28722320(unfulfilled),afterok:999(unfulfilled)')])
    calls, pauses = [], []
    def lookup(argv):
        calls.append(argv)
        return next(replies)
    row = handoff.acknowledge('28692775',target,lookup=lookup,pause=pauses.append)
    assert dispatch.dependencies(row['Dependency']) == target and pauses == [1]
    assert all(a[:3] == ['scontrol','show','job'] for a in calls)
    with pytest.raises(TimeoutError):
        handoff.acknowledge('28692775',target,lookup=lambda a:waiter_body('28692775'),pause=lambda n:None)
    with pytest.raises(ValueError):
        handoff.acknowledge('28692775',target,lookup=lambda a:waiter_body('28692775',state='RUNNING'),pause=lambda n:None)


def test_one_held_submission_and_later_handoff_never_release_successors(tmp_path, monkeypatch):
    monkeypatch.setattr(dispatch,'ROOT',tmp_path)
    monkeypatch.setattr(dispatch,'SCRIPT',tmp_path/'recipe.slurm')
    (tmp_path/'audit_receipt.json').write_text('synthetic audited receipt')
    (tmp_path/'recipe.slurm').write_bytes(SCRIPT.read_bytes())
    (tmp_path/'submit_held_once.py').write_bytes(Path(dispatch.__file__).read_bytes())
    audit = {'seed_manifest_sha256':'a'*64}
    monkeypatch.setattr(dispatch,'verify_audit_seed',lambda:audit)
    dependencies = {job:'afterok:28722320(unfulfilled)' for job in dispatch.WAITERS}
    calls = []
    def lookup(argv):
        calls.append(argv)
        if argv[0]=='sacct':
            return '28661019|COMPLETED|0:0|32|hfacnormal01|iai806'
        if argv[0]=='squeue':
            return '28722320|hfo2-G2-M|RUNNING|hfacnormal01|32'
        if argv[0]=='bash':
            return ''
        if argv[0]=='sbatch':
            assert '--hold' in argv
            return '999'
        if argv[:2]==['scontrol','show']:
            if argv[-1]=='999':
                return ('JobId=999 JobName=hfo2-G2-flip-E067 UserId=iai806(16284) '
                        'JobState=PENDING Reason=JobHeldUser Priority=0 Partition=hfacnormal01 '
                        f'NumCPUs=32 Command={dispatch.SCRIPT} StdOut={tmp_path / "flip.slurm.out"} '
                        f'StdErr={tmp_path / "flip.slurm.err"} Dependency=(null)')
            return waiter_body(argv[-1],dependencies[argv[-1]])
        if argv[:2]==['scontrol','update']:
            dependencies[argv[2].split('=',1)[1]] = argv[3].split('=',1)[1]
            return ''
        if argv[:2]==['scontrol','release']:
            assert argv[-1]=='999'
            return ''
        raise AssertionError(argv)
    monkeypatch.setattr(dispatch,'run',lookup)
    dispatch.submit_once()
    assert len([a for a in calls if a[0]=='sbatch'])==1
    assert not any(a[:2] in (['scontrol','release'],['scontrol','update']) for a in calls)
    with pytest.raises(FileExistsError):
        dispatch.submit_once()
    assert len([a for a in calls if a[0]=='sbatch'])==1
    # Acknowledge defaults bind at definition; replace with explicit mock lookup.
    original_ack = handoff.acknowledge
    monkeypatch.setattr(handoff,'acknowledge',lambda j,e:original_ack(j,e,lookup=lookup,pause=lambda n:None))
    handoff.handoff_once()
    assert [a for a in calls if a[:2]==['scontrol','release']] == [['scontrol','release','999']]
    assert len([a for a in calls if a[:2]==['scontrol','update']])==2
    receipt = json.loads((tmp_path/'handoff_receipt.json').read_text())
    assert receipt['successors_still_held']==list(dispatch.WAITERS)
    assert receipt['new_submissions_during_handoff']==0 and not receipt['execution_started_claimed']
    assert all(dispatch.dependencies(d)=={dispatch.OTHER,'999'} for d in dependencies.values())
    with pytest.raises(ValueError):
        handoff.handoff_once()
    assert len([a for a in calls if a[:2]==['scontrol','update']])==2


def test_ambiguous_submission_keeps_intent_and_cannot_be_repeated(tmp_path,monkeypatch):
    monkeypatch.setattr(dispatch,'ROOT',tmp_path)
    monkeypatch.setattr(dispatch,'verify_audit_seed',lambda:{'seed_manifest_sha256':'a'*64})
    (tmp_path/'audit_receipt.json').write_text('synthetic audit')
    calls=[]
    def lookup(argv):
        calls.append(argv)
        if argv[0]=='sacct': return '28661019|COMPLETED|0:0|32|hfacnormal01|iai806'
        if argv[0]=='squeue': return '28722320|hfo2-G2-M|RUNNING|hfacnormal01|32'
        if argv[0]=='bash': return ''
        if argv[0]=='sbatch': raise TimeoutError('accepted but response lost')
        return waiter_body(argv[-1])
    monkeypatch.setattr(dispatch,'run',lookup)
    with pytest.raises(TimeoutError): dispatch.submit_once()
    journal=[json.loads(r) for r in (tmp_path/'submission_journal.jsonl').read_text().splitlines()]
    assert journal[-1]['event']=='reconcile_existing_intent_do_not_resubmit'
    with pytest.raises(FileExistsError): dispatch.submit_once()
    assert len([a for a in calls if a[0]=='sbatch'])==1


def test_separate_verified_release_never_repeats_submission_or_dependency_writes(tmp_path,monkeypatch):
    spec = importlib.util.spec_from_file_location('flip_verified_release',
        ROOT/'benchmarks/hfo2_channels/20261008/clamped_flip_continue_E067_20261011/release_verified.py')
    release = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(release)
    monkeypatch.setattr(dispatch,'ROOT',tmp_path)
    monkeypatch.setattr(dispatch,'verify_audit_seed',lambda:None)
    (tmp_path/'audit_receipt.json').write_text('synthetic audit')
    (tmp_path/'submit_held_once.py').write_bytes(Path(dispatch.__file__).read_bytes())
    accepted={'job_id':'999','total_submissions':1,'parent_job_id':dispatch.PARENT,
        'audit_receipt_sha256':dispatch.sha(tmp_path/'audit_receipt.json'),
        'executed_helper_sha256':dispatch.sha(tmp_path/'submit_held_once.py')}
    (tmp_path/'accepted_receipt.json').write_text(json.dumps(accepted))
    journal=[{'event':'dependency_update_intent','job_id':j,'dependency':'afterok:28722320:999'}
             for j in dispatch.WAITERS]+[{'event':'reconcile_do_not_repeat'}]
    (tmp_path/'handoff_journal.jsonl').write_text('\n'.join(json.dumps(r) for r in journal))
    calls=[]
    def lookup(argv):
        calls.append(argv)
        if argv[0]=='squeue': return '28722320|hfo2-G2-M|RUNNING|hfacnormal01|32'
        if argv[:2]==['scontrol','release']: return ''
        if argv[-1]=='999':
            return ('JobId=999 JobName=hfo2-G2-flip-E067 UserId=iai806(16284) '
                    'JobState=PENDING Reason=JobHeldUser Priority=0 Partition=hfacnormal01 '
                    f'NumCPUs=32 Command={dispatch.SCRIPT} StdOut={tmp_path / "flip.slurm.out"} '
                    f'StdErr={tmp_path / "flip.slurm.err"} Dependency=(null)')
        return waiter_body(argv[-1],'afterok:28722320,afterok:999')
    monkeypatch.setattr(dispatch,'run',lookup)
    release.release_once()
    assert [a for a in calls if a[:2]==['scontrol','release']]==[['scontrol','release','999']]
    assert not any(a[0]=='sbatch' or a[:2]==['scontrol','update'] for a in calls)
    actual=json.loads((tmp_path/'handoff_receipt.json').read_text())
    assert actual['dependency_writes_during_reconciliation']==actual['new_submissions_during_reconciliation']==0
    assert actual['successors_still_held']==list(dispatch.WAITERS) and not actual['execution_started_claimed']
    with pytest.raises(FileExistsError): release.release_once()
    assert len([a for a in calls if a[:2]==['scontrol','release']])==1


@pytest.mark.parametrize('step',[0,15,20])
def test_actual_HF_terminal_segment_raw_EFS_and_residual_replay(step):
    case=ROOT/'benchmarks/hfo2_channels/20261008/clamped_flip_continue_E067_20261011'
    directory=case/'observations'/f'step_{step:04d}'
    report=json.loads((directory/'observation.json').read_text())
    images=read(directory/'evaluated_chain.traj',index=':')
    assert len(images)==9 and report['source_job_id']=='28661019'
    assert observer.sha256(directory/'evaluated_chain.traj')==report['evaluated_chain_sha256']
    for i,image in enumerate(images):
        raw=directory/'raw'/f'image_{i:04d}'
        pinned=report['raw_image_evaluations'][i]
        assert set(pinned['input_sha256'])=={*CONTRACT,'STRU'}
        assert all(pinned['input_sha256'][n]==h for n,h in CONTRACT.items())
        # Full six-file verification occurred on HF. Only these portable
        # actual bytes, native EFS and declared digests are replayed locally.
        for name in ('INPUT','KPT','STRU'):
            assert observer.sha256(raw/name)==pinned['input_sha256'][name]
        assert observer.sha256(raw/'OUT.ABACUS/running_scf.log')==pinned['raw_log_sha256']
        assert same_ordered_geometry(image,read_fixed_hfo2_stru(raw/'STRU'))
        parsed=audited_results(raw)
        for key,expected in (('energy',pinned['energy_eV_cell']),('forces',pinned['forces_eV_A']),
                             ('stress',pinned['stress_ASE_voigt_eV_A3'])):
            np.testing.assert_allclose(parsed[key],expected,atol=1e-12,rtol=0)
        image.calc=SinglePointCalculator(image,**parsed)
    reference=read(case/'seed/substrate.vasp').cell.array
    boundary=clamped_plane_vcneb_boundary(12,reference,allow_tilt=True)
    chain=VCNEB(images,cell_scale=report['cell_scale_A'],k=.2,climb=False,pressure=0,
                **boundary.vcneb_kwargs(images))
    assert np.linalg.norm(chain.get_forces(),axis=1).max()==pytest.approx(report['replayed_fmax_eV_A'],abs=1e-10,rel=0)
    np.testing.assert_allclose((chain.enthalpies-chain.enthalpies[0])*1000/4,
                              report['relative_enthalpy_meV_fu'],atol=1e-9,rtol=0)
    actual=analyze(directory)
    historical=json.loads((directory/'residual.json').read_text())
    for key in ('fmax_eV_A','perpendicular_only_fmax_same_geometry_eV_A',
                'spring_only_fmax_same_geometry_eV_A'):
        assert actual[key]==pytest.approx(historical[key],abs=1e-10,rel=0)
    assert not report['ordinary_residual_pass'] and report['new_DFT_calls']==0


def test_actual_HF_single_dispatch_delayed_acknowledgement_and_separate_release():
    case=ROOT/'benchmarks/hfo2_channels/20261008/clamped_flip_continue_E067_20261011'
    audit=json.loads((case/'audit_receipt.json').read_text())
    accepted=json.loads((case/'accepted_receipt.json').read_text())
    receipt=json.loads((case/'handoff_receipt.json').read_text())
    startup=json.loads((case/'startup_verification.json').read_text())
    journal=[json.loads(r) for r in (case/'submission_journal.jsonl').read_text().splitlines()]
    handoff_rows=[json.loads(r) for r in (case/'handoff_journal.jsonl').read_text().splitlines()]
    assert [r['job_id'] for r in journal if r['event']=='accepted_held_handle']==['28722810']
    assert [r['job_id'] for r in handoff_rows if r['event']=='dependency_update_intent']==list(dispatch.WAITERS)
    assert handoff_rows[-1]['event']=='reconcile_do_not_repeat'
    assert receipt['job_id']==accepted['job_id']=='28722810' and receipt['total_E067_submissions']==1
    assert receipt['new_submissions_during_reconciliation']==receipt['dependency_writes_during_reconciliation']==0
    assert receipt['accepted_receipt_sha256']==observer.sha256(case/'accepted_receipt.json')
    assert receipt['original_handoff_journal_sha256']==observer.sha256(case/'handoff_journal.jsonl')
    assert accepted['audit_receipt_sha256']==observer.sha256(case/'audit_receipt.json')
    assert audit['fresh_parent_interior_SCFs']==140 and audit['current_frame_exact_caches_checked']==9
    assert len(audit['runtime_code_sha256'])==334
    assert audit['executed_audit_source_sha256']==observer.sha256(case/'audit_HF.py')
    assert audit['analysis_source_sha256']['export_hfo2_clamped_observation.py']==observer.sha256(case/'executed_exporter.py')
    assert receipt['executed_release_sha256']==observer.sha256(case/'release_verified.py')
    assert startup['first_fresh_SCF_native_DSIZE_lines']==['                                    DSIZE = 32']
    assert '0.270056' in startup['new_flip_optimizer_log'] and startup['stderr_tail']==''
    assert all('|RUNNING|0:0|' in row and '|32|' in row for row in startup['raw_accounting'].splitlines())
    assert not startup['first_fresh_SCF_completion_claimed'] and not receipt['execution_started_claimed']
    assert not any(audit[k] for k in ('new_DFT_calls','physical_inputs_changed','holdout_generated_or_read','ordinary_converged'))
