"""Pure queue guards: no HF, scheduler or DFT calls."""
import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'benchmarks/hfo2_channels/20261008/clamped_G2_training_queue_E060_20261010'
spec = importlib.util.spec_from_file_location('G2_E060',CASE/'queue_training_waves.py')
wave = importlib.util.module_from_spec(spec)
spec.loader.exec_module(wave)
ACCT = '\n'.join(f'{j}|'+('RUNNING|0:0|32' if i < 2 else 'PENDING|0:0|0')+'|hfacnormal01|iai806'
                 for i,j in enumerate(wave.KNOWN))
QUEUE = '\n'.join(f'{j}|hfo2-G2-existing|'+('RUNNING' if i < 2 else 'PENDING')+'|hfacnormal01|32'
                  for i,j in enumerate(wave.KNOWN))


def test_registered_training_sources_and_two_allocation_guard():
    result = wave.reconcile(ACCT,QUEUE+'\n999|pc-other|RUNNING|other|128')
    assert result['active_study_allocations'] == 2 and result['unrelated_jobs_untouched']
    assert wave.PREDECESSORS == ('28692775','28692776')
    for channel in wave.CHANNELS:
        argv = wave.sbatch_argv(channel,wave.PREDECESSORS)
        assert '--dependency=afterok:28692775:28692776' in argv
        exports = next(a for a in argv if a.startswith('--export='))
        assert 'strain_p0100/'+channel in exports and 'RUN_DFT=1' in exports
        assert argv[-1] == (wave.SOURCE/wave.PILOT).as_posix()
    assert wave.reconcile(ACCT.replace('RUNNING|0:0|32','COMPLETED|0:0|32'),
                          '\n'.join(QUEUE.splitlines()[2:]))['active_study_allocations'] == 0


@pytest.mark.parametrize('bad',[
    ACCT.replace('RUNNING','TIMEOUT'),ACCT.replace('RUNNING','FAILED'),
    ACCT.replace('0:0','1:0'),ACCT.replace('|32|','|64|'),
    ACCT.replace('hfacnormal01','other'),ACCT.replace('iai806','someone'),
    '\n'.join(ACCT.splitlines()[:3]),ACCT+'\n'+ACCT.splitlines()[0],
    ACCT.replace('28692775','9999')])
def test_failed_unknown_missing_or_wrong_resource_handles_prevent_submission(bad):
    with pytest.raises(ValueError):wave.reconcile(bad,QUEUE)


@pytest.mark.parametrize('extra',[
    '999|hfo2-unregistered|PENDING|hfacnormal01|32',
    '28661019|hfo2-existing|RUNNING|hfacnormal01|32'])
def test_unknown_or_duplicate_study_queue_is_rejected(extra):
    with pytest.raises(ValueError):wave.reconcile(ACCT,QUEUE+'\n'+extra)


def test_more_than_two_running_allocations_are_rejected():
    with pytest.raises(ValueError):wave.reconcile(ACCT,QUEUE.replace('PENDING','RUNNING'))


def predecessor_controls(dependency='afterok:28661019(unfulfilled),afterok:28661020(unfulfilled)'):
    return {j:' '.join((f'JobId={j}','JobName=hfo2-G2-E058-'+('PO_to_T','PO_flip_T_pattern_reversing')[i],
        'JobState=PENDING','Partition=hfacnormal01','NumCPUs=32','NumTasks=1','CPUs/Task=32','NumNodes=1-1',
        'Command='+(wave.SOURCE/wave.PILOT).as_posix(),'Dependency='+dependency)) for i,j in enumerate(wave.PREDECESSORS)}


def test_dependency_graph_requires_every_still_active_original_parent():
    rows = wave.reconcile(ACCT,QUEUE)['existing_accounting_rows']
    assert len(wave.validate_predecessor_controls(rows,predecessor_controls())) == 2
    rows[0][1] = 'COMPLETED'
    assert len(wave.validate_predecessor_controls(rows,predecessor_controls('afterok:28661020(unfulfilled)'))) == 2
    rows[1][1] = 'COMPLETED'
    assert len(wave.validate_predecessor_controls(rows,predecessor_controls('(null)'))) == 2


@pytest.mark.parametrize('bad',['(null)','afterany:28661019','afterok:9999','afterok:28661019',
                               'afterok:28661019,afterok:28661019'])
def test_removed_weakened_unknown_or_duplicate_parent_dependencies_block_submission(bad):
    rows = wave.reconcile(ACCT,QUEUE)['existing_accounting_rows']
    with pytest.raises(ValueError):wave.validate_predecessor_controls(rows,predecessor_controls(bad))


def test_missing_changed_or_premature_zero_strain_controls_block_submission():
    rows = wave.reconcile(ACCT,QUEUE)['existing_accounting_rows']
    controls = predecessor_controls()
    with pytest.raises(ValueError):wave.validate_predecessor_controls(rows,{wave.PREDECESSORS[0]:controls[wave.PREDECESSORS[0]]})
    controls[wave.PREDECESSORS[0]] = controls[wave.PREDECESSORS[0]].replace('NumCPUs=32','NumCPUs=64')
    with pytest.raises(ValueError):wave.validate_predecessor_controls(rows,controls)
    controls = {j:t.replace('JobState=PENDING','JobState=RUNNING').replace(
        'Dependency=afterok:28661019(unfulfilled),afterok:28661020(unfulfilled)','Dependency=(null)')
        for j,t in predecessor_controls().items()}
    for r in rows[2:]:r[1] = 'RUNNING';r[3] = '32'
    with pytest.raises(ValueError):wave.validate_predecessor_controls(rows,controls)


@pytest.mark.parametrize('deps',[('1',),('1','1'),('0','2'),('x','2'),['1','2'],(1,'2')])
def test_malformed_or_implicit_dependencies_are_rejected(deps):
    with pytest.raises(ValueError):wave.sbatch_argv(wave.CHANNELS[0],deps)


def test_holdout_and_new_channels_cannot_be_submitted():
    for channel in ('strain_p0050','new_channel'):
        with pytest.raises(ValueError):wave.sbatch_argv(channel,wave.PREDECESSORS)


def setup_submission(monkeypatch,tmp_path):
    monkeypatch.setattr(wave,'ROOT',tmp_path/'fresh');wave.ROOT.mkdir()
    seed = tmp_path/'seed';seed.mkdir();(seed/'manifest.json').write_text('{}')
    receipt = dict(status='four_training_starts_preflight_passed_NOT_submitted',queue_helper_sha256='helper',
        waves=[list(p) for p in wave.WAVES],condition=wave.CONDITION,
        reports=[dict(channel=c,seed_root=str(seed),manifest_sha256='manifest',seed_files_sha256={}) for c in wave.CHANNELS])
    (wave.ROOT/'preflight.json').write_text(json.dumps(receipt))
    def digest(path):
        if str(path).endswith('queue_training_waves.py'):return 'helper'
        if Path(path) == wave.SOURCE/wave.PILOT:return wave.PILOT_SHA
        if str(path).endswith('manifest.json'):return 'manifest'
        return 'receipt'
    monkeypatch.setattr(wave,'sha',digest)
    monkeypatch.setattr(wave,'scheduler_snapshot',lambda:wave.reconcile(ACCT,QUEUE))


def test_actual_returned_handles_chain_second_pair_and_single_use(monkeypatch,tmp_path):
    setup_submission(monkeypatch,tmp_path)
    replies = iter(['12345','12346','12347','12348']);calls = []
    def fake_run(argv):calls.append(argv);return next(replies)
    monkeypatch.setattr(wave,'run',fake_run)
    receipt = wave.submit_once()
    assert [j['job_id'] for j in receipt['jobs']] == ['12345','12346','12347','12348']
    assert all('--dependency=afterok:28692775:28692776' in a for a in calls[:2])
    assert all('--dependency=afterok:12345:12346' in a for a in calls[2:])
    assert receipt['G2_unique_chains_with_handles'] == 8 and receipt['G2_training_starts_remaining_without_handles'] == 0
    assert receipt['job_completion_is_not_NEB_convergence'] and not receipt['parent_or_previous_jobs_modified']
    assert receipt['segment_steps'] == 10 and receipt['segment_walltime_hours'] == 4
    assert not receipt['physical_parameters_changed'] and not receipt['held_out_condition_generated_or_read']
    with pytest.raises(FileExistsError):wave.submit_once()
    assert len(calls) == 4


@pytest.mark.parametrize('bad_reply',['unclear','12345','28692775','0'])
def test_partial_acceptance_survives_failed_reply_without_blind_retry(monkeypatch,tmp_path,bad_reply):
    setup_submission(monkeypatch,tmp_path)
    replies = iter(['12345',bad_reply]);calls = []
    def fake_run(argv):calls.append(argv);return next(replies)
    monkeypatch.setattr(wave,'run',fake_run)
    with pytest.raises(ValueError,match='sbatch reply'):wave.submit_once()
    records = [json.loads(line) for line in (wave.ROOT/'submission_journal.jsonl').read_text().splitlines()]
    assert records[-1]['event'] == 'ambiguous_or_failed_reply_RECONCILE_DO_NOT_RETRY'
    assert any(r['event'] == 'accepted_handle' and r['job_id'] == '12345' for r in records)
    assert not (wave.ROOT/'submission_receipt.json').exists()
    with pytest.raises(FileExistsError):wave.submit_once()
    assert len(calls) == 2


def test_delivered_actual_HF_preflight_checks_four_registered_starts_and_eight_caches():
    receipt = json.loads((CASE/'preflight.json').read_text())
    assert receipt['status'] == 'four_training_starts_preflight_passed_NOT_submitted'
    assert receipt['condition'] == 'strain_p0100' and receipt['code_files_byte_checked'] == 338
    assert receipt['archive_sha256'] == wave.ARCHIVE_SHA and receipt['pilot_sha256'] == wave.PILOT_SHA
    assert receipt['ABACUS_binary_sha256'] == wave.BINARY_SHA
    assert receipt['queue_helper_sha256'] == wave.sha(CASE/'queue_training_waves.py') == wave.sha(CASE/'executed_queue_training_waves.py')
    assert not receipt['new_DFT_calls'] and not receipt['existing_jobs_mutated']
    assert not receipt['physical_parameters_changed'] and not receipt['held_out_condition_generated_or_read']
    assert len(receipt['scheduler']['zero_strain_predecessor_controls']) == 2
    assert all(r['dependency_graph_verified'] for r in receipt['scheduler']['zero_strain_predecessor_controls'])
    seeds = CASE.parent/'clamped_G2_remaining_E056_20261010/seeds/strain_p0100'
    assert [r['channel'] for r in receipt['reports']] == list(wave.CHANNELS)
    for report in receipt['reports']:
        seed = seeds/report['channel']
        assert wave.sha(seed/'manifest.json') == report['manifest_sha256']
        assert all(wave.sha(seed/n) == h for n,h in report['seed_files_sha256'].items())
        assert report['internal_images_fresh'] == 7 and report['geometry_checked_images'] == 9
        assert report['physical_bytes_native_32MPI_raw_EFS_rechecked']
        assert [r['image_index'] for r in report['endpoint_caches']] == [0,8]
        for cache in report['endpoint_caches']:
            path = CASE/'cache_preflight'/report['channel']/f"image_{cache['image_index']:04d}/seed_cache_audit.json"
            assert wave.sha(path) == cache['audit_sha256'] and not cache['new_DFT_calls']
            assert not json.loads(path.read_text())['new_DFT_calls']


def test_delivered_actual_eight_handle_graph_resources_and_bounded_matrix():
    receipt = json.loads((CASE/'submission_receipt.json').read_text())
    proof = json.loads((CASE/'queued_jobs_verification.json').read_text())
    assert receipt['preflight_sha256'] == wave.sha(CASE/'preflight.json')
    assert receipt['submission_journal_sha256'] == wave.sha(CASE/'submission_journal.jsonl')
    assert proof['submission_receipt_sha256'] == wave.sha(CASE/'submission_receipt.json')
    assert [(j['channel'],j['job_id']) for j in receipt['jobs']] == [
        ('PO_to_T','28709788'),('PO_to_M','28709789'),
        ('PO_flip_T_pattern_preserving','28709790'),('PO_flip_T_pattern_reversing','28709791')]
    assert [j['dependency'] for j in receipt['jobs']] == [
        'afterok:28692775:28692776']*2 + ['afterok:28709788:28709789']*2
    assert proof['active_study_allocations'] == 2 and proof['pending_study_starts'] == 6
    assert proof['all_eight_registered_G2_chains_have_handles']
    assert receipt['G2_unique_chains_with_handles'] == receipt['total_registered_G2_matrix_chains'] == 8
    assert not receipt['G2_training_starts_remaining_without_handles']
    assert receipt['job_completion_is_not_NEB_convergence'] and proof['job_completion_is_not_NEB_convergence']
    assert not proof['study_mutated'] and not proof['recurring_monitor_created']
    for job in proof['jobs']:
        f = job['actual_fields']
        assert job['resource_and_dependency_checks_passed'] and f['JobState'] == 'PENDING' and f['Reason'] == 'Dependency'
        assert f['NumCPUs'] == f['CPUs/Task'] == '32' and f['NumTasks'] == '1'
        assert f['NumNodes'] == '1-1' and f['TimeLimit'] == '04:00:00' and f['Partition'] == 'hfacnormal01'
    assert not any(p.suffix in ('.upf','.orb') for p in CASE.rglob('*') if p.is_file())
