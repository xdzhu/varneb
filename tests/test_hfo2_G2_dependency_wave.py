"""Offline guards for the bounded second G2 wave; never contact HF or submit."""
import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'benchmarks/hfo2_channels/20261008/clamped_G2_queue_E058_20261010'
spec = importlib.util.spec_from_file_location('G2_wave_E058', CASE/'queue_next_pair.py')
wave = importlib.util.module_from_spec(spec)
spec.loader.exec_module(wave)

ACCT = '\n'.join(f'{j}|RUNNING|0:0|32|hfacnormal01|iai806' for j in wave.PARENTS)
QUEUE = '\n'.join(f'{j}|hfo2-G2-cont-{i}|RUNNING|hfacnormal01|32' for i,j in enumerate(wave.PARENTS))


def test_dependency_waits_for_both_parents_and_only_remaining_zero_strain_starts():
    for channel in wave.CHANNELS:
        argv = wave.sbatch_argv(channel)
        assert '--dependency=afterok:28661019:28661020' in argv
        exports = next(a for a in argv if a.startswith('--export='))
        assert 'strain_0000/'+channel in exports and 'RUN_DFT=1' in exports
        assert 'SOURCE_ROOT='+wave.SOURCE.as_posix() in exports
        assert (wave.ROOT/channel/'band').as_posix() in exports
        assert argv[-1] == (wave.SOURCE/wave.PILOT).as_posix()
    for channel in ('PO_to_M','PO_flip_T_pattern_preserving','strain_p0050','new_channel'):
        with pytest.raises(ValueError): wave.sbatch_argv(channel)


def test_healthy_parents_can_be_live_or_healthy_completed_and_unrelated_jobs_are_untouched():
    assert len(wave.reconcile(ACCT,QUEUE+'\n999|pc-other|RUNNING|other|128')['active_study_rows']) == 2
    assert not wave.reconcile(ACCT.replace('RUNNING','COMPLETED'),'')['active_study_rows']


@pytest.mark.parametrize('bad', [
    ACCT.replace('RUNNING','TIMEOUT'), ACCT.replace('RUNNING','FAILED'),
    ACCT.replace('0:0','1:0'), ACCT.replace('|32|','|64|'),
    ACCT.replace('hfacnormal01','other'), ACCT.replace('iai806','someone'),
    ACCT.splitlines()[0], ACCT+'\n'+ACCT.splitlines()[0], ACCT.replace('28661019','1234')])
def test_unknown_failed_or_wrong_resource_parents_prevent_submission(bad):
    with pytest.raises(ValueError): wave.reconcile(bad,QUEUE)


def test_existing_other_study_job_even_pending_prevents_duplicate_wave():
    with pytest.raises(ValueError):
        wave.reconcile(ACCT,QUEUE+'\n123|hfo2-G2-E058-PO_to_T|PENDING|hfacnormal01|32')


def setup_submission(monkeypatch,tmp_path):
    monkeypatch.setattr(wave,'ROOT',tmp_path/'fresh')
    wave.ROOT.mkdir()
    seed=tmp_path/'seed';seed.mkdir()
    (seed/'manifest.json').write_text('{}')
    receipt=dict(status='preflight_passed_NOT_submitted',queue_helper_sha256='helper',
        registered_channels=list(wave.CHANNELS),reports=[dict(channel=c,seed_root=str(seed),
        manifest_sha256='manifest',seed_files_sha256={}) for c in wave.CHANNELS])
    (wave.ROOT/'preflight.json').write_text(json.dumps(receipt))
    def digest(path):
        if str(path).endswith('queue_next_pair.py'): return 'helper'
        if Path(path) == wave.SOURCE/wave.PILOT: return wave.PILOT_SHA
        if str(path).endswith('manifest.json'): return 'manifest'
        return 'receipt'
    monkeypatch.setattr(wave,'sha',digest)
    monkeypatch.setattr(wave,'scheduler_snapshot',lambda: wave.reconcile(ACCT,QUEUE))


def test_partial_accepted_handle_survives_ambiguous_next_reply_and_replay_is_refused(monkeypatch,tmp_path):
    setup_submission(monkeypatch,tmp_path)
    replies=iter(['12345','unclear scheduler response'])
    calls=[]
    def fake_run(argv): calls.append(argv);return next(replies)
    monkeypatch.setattr(wave,'run',fake_run)
    with pytest.raises(ValueError,match='ambiguous'): wave.submit_once()
    lines=[json.loads(s) for s in (wave.ROOT/'submission_journal.jsonl').read_text().splitlines()]
    assert any(r['event']=='accepted_handle' and r['job_id']=='12345' for r in lines)
    assert lines[-1]['event']=='ambiguous_or_failed_reply_RECONCILE_DO_NOT_RETRY'
    assert len(calls)==2 and not (wave.ROOT/'submission_receipt.json').exists()
    with pytest.raises(FileExistsError): wave.submit_once()
    assert len(calls)==2


def test_submission_is_single_use_and_records_both_accepted_jobs(monkeypatch,tmp_path):
    setup_submission(monkeypatch,tmp_path)
    replies=iter(['12345','12346'])
    monkeypatch.setattr(wave,'run',lambda argv:next(replies))
    r=wave.submit_once()
    assert [j['job_id'] for j in r['jobs']]==['12345','12346']
    assert r['G2_unique_chains_with_handles']==4 and r['G2_training_starts_remaining_without_handles']==4
    assert r['segment_steps']==10 and r['max_study_simultaneous_allocations']==2
    assert r['job_completion_is_not_NEB_convergence']
    assert not r['physical_parameters_changed'] and not r['held_out_condition_generated_or_read']
    with pytest.raises(FileExistsError): wave.submit_once()


def test_delivered_actual_HF_preflight_pins_original_source_seeds_and_four_caches():
    receipt=json.loads((CASE/'preflight.json').read_text())
    assert receipt['status']=='preflight_passed_NOT_submitted' and receipt['code_files_byte_checked']==338
    assert len(receipt['source_code_sha256'])==338 and receipt['archive_sha256']==wave.ARCHIVE_SHA
    assert receipt['pilot_sha256']==wave.PILOT_SHA and receipt['ABACUS_binary_sha256']==wave.BINARY_SHA
    assert receipt['queue_helper_sha256']==wave.sha(CASE/'queue_next_pair.py')==wave.sha(CASE/'executed_queue_next_pair.py')
    assert not receipt['new_DFT_calls'] and not receipt['existing_jobs_mutated']
    assert not receipt['physical_parameters_changed'] and not receipt['held_out_condition_generated_or_read']
    seeds=CASE.parent/'clamped_G2_remaining_E056_20261010/seeds/strain_0000'
    for report in receipt['reports']:
        seed=seeds/report['channel']
        assert wave.sha(seed/'manifest.json')==report['manifest_sha256']
        assert all(wave.sha(seed/n)==h for n,h in report['seed_files_sha256'].items())
        assert report['internal_images_fresh']==7 and report['geometry_checked_images']==9
        assert report['physical_bytes_native_32MPI_raw_EFS_rechecked']
        for cache in report['endpoint_caches']:
            path=CASE/'cache_preflight'/report['channel']/f"image_{cache['image_index']:04d}/seed_cache_audit.json"
            assert wave.sha(path)==cache['audit_sha256'] and not cache['new_DFT_calls']
            assert not json.loads(path.read_text())['new_DFT_calls']


def test_delivered_actual_job_handles_have_both_dependencies_and_original_resources():
    receipt=json.loads((CASE/'submission_receipt.json').read_text())
    proof=json.loads((CASE/'queued_jobs_verification.json').read_text())
    assert receipt['preflight_sha256']==wave.sha(CASE/'preflight.json')
    assert receipt['submission_journal_sha256']==wave.sha(CASE/'submission_journal.jsonl')
    assert proof['submission_receipt_sha256']==wave.sha(CASE/'submission_receipt.json')
    assert [(j['channel'],j['job_id']) for j in receipt['jobs']]==[
        ('PO_to_T','28692775'),('PO_flip_T_pattern_reversing','28692776')]
    assert proof['active_study_allocations']==proof['pending_study_starts']==2
    assert not proof['study_mutated'] and not proof['recurring_monitor_created']
    assert receipt['dependency']==wave.DEPENDENCY and receipt['job_completion_is_not_NEB_convergence']
    for job in proof['jobs']:
        f=job['actual_fields']
        assert job['resource_and_dependency_checks_passed'] and f['JobState']=='PENDING' and f['Reason']=='Dependency'
        assert f['Dependency']=='afterok:28661019(unfulfilled),afterok:28661020(unfulfilled)'
        assert f['NumCPUs']==f['CPUs/Task']=='32' and f['NumTasks']=='1' and f['NumNodes']=='1-1'
        assert f['TimeLimit']=='04:00:00' and f['Partition']=='hfacnormal01'
    assert not any(p.suffix in ('.upf','.orb') for p in CASE.rglob('*') if p.is_file())
