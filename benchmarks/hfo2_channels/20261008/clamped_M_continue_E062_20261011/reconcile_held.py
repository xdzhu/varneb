"""Finish the recorded accepted held job; never submit another allocation."""
from datetime import datetime
import json
import os
from pathlib import Path
import reviewed_handoff_rules as rules


def reconcile():
    root = rules.ROOT
    records = [json.loads(line) for line in (root/'submission_journal.jsonl').read_text().splitlines()]
    accepted = [r for r in records if r['event'] == 'accepted_held_handle']
    if (len(accepted) != 1 or accepted[0]['job_id'] != '28722320'
            or records[-1]['event'] != 'reconcile_do_not_retry'
            or records[-1]['waiters_already_updated'] != []
            or rules.sha(root/'queue_continuation.py') != 'fd6ee53777bc8b13b4fad37e6b572ad01824de916295ad9aef797d8d4f30ae5b'
            or rules.sha(rules.SCRIPT) != rules.SCRIPT_SHA):
        raise ValueError('accepted-handle reconciliation evidence differs')
    handle = accepted[0]['job_id']
    actual = rules.run(['sacct','-X','-n','-P','-j',handle,
                        '--format=JobIDRaw,State,ExitCode,Partition,User'])
    if actual != handle+'|PENDING|0:0|hfacnormal01|iai806':
        raise ValueError('accepted held job not corroborated by accounting')
    held = rules.fields(rules.run(['scontrol','show','job','-o',handle]))
    if (held.get('JobState') != 'PENDING' or held.get('Reason') != 'JobHeldUser'
            or held.get('JobName') != 'hfo2-G2-M-E062' or held.get('NumCPUs') != '32'
            or held.get('Command') != str(rules.SCRIPT) or held.get('Partition') != 'hfacnormal01'):
        raise ValueError('same held allocation identity changed')
    receipt = json.loads((root/'audit_receipt.json').read_text())
    seed = root/'seed'
    manifest = json.loads((seed/'manifest.json').read_text())
    if (receipt['seed_manifest_sha256'] != rules.sha(seed/'manifest.json')
            or any(rules.sha(seed/n) != h for n,h in manifest['files_sha256'].items())
            or not receipt['all_parent_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked']
            or receipt['current_frame_exact_caches_checked'] != 9 or receipt['physical_inputs_changed']
            or rules.sha(rules.BINARY) != rules.BINARY_SHA):
        raise ValueError('audited continuation or binary changed')
    runtime = root.parent/'20261010-varneb-clamped-G2-E054-r1/source-fixed'
    if any(rules.sha(runtime/n) != h for n,h in receipt['runtime_code_sha256'].items()):
        raise ValueError('immutable scientific runtime changed')
    for job in rules.WAITING:
        rules.validate_waiter(rules.run(['scontrol','show','job','-o',job]), job, {'28661019'})
    rules.validate_active_queue(rules.run(['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C']), rules.ALLOWED|{handle})
    dependency = 'afterok:28661019:'+handle
    with (root/'reconciliation_journal.jsonl').open('x') as journal:
        def record(event, **values):
            journal.write(json.dumps(dict(event=event, check_CST=datetime.now().isoformat(), **values))+'\n')
            journal.flush(); os.fsync(journal.fileno())
        record('same_accepted_handle_verified', job_id=handle, accounting=actual, held_fields=held)
        for job in rules.WAITING:
            rules.validate_waiter(rules.run(['scontrol','show','job','-o',job]), job, {'28661019'})
            record('dependency_update_intent', job_id=job, dependency=dependency)
            rules.run(['scontrol','update','JobId='+job,'Dependency='+dependency])
            checked = rules.validate_waiter(rules.run(['scontrol','show','job','-o',job]), job, {'28661019',handle})
            record('dependency_update_verified', job_id=job, actual_dependency=checked['Dependency'])
        rules.validate_active_queue(rules.run(['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C']), rules.ALLOWED|{handle})
        record('release_intent', job_id=handle)
        rules.run(['scontrol','release',handle])
        record('released', job_id=handle)
    result = dict(status='same_accepted_continuation_released_after_reconciliation',
        check_CST=datetime.now().isoformat(), job_id=handle, parent_job_id='28661020',
        new_submissions_during_reconciliation=0, total_E062_submissions=1,
        waiting_jobs_updated=list(rules.WAITING), dependency=dependency,
        max_study_simultaneous_allocations=2, new_independent_chains=0,
        new_segment_steps=20, walltime_hours=6, MPI_ranks=32, ordinary_fmax=.10,
        physical_inputs_changed=False, running_source_overwritten=False,
        held_out_condition_generated_or_read=False, old_parent_not_modified=True,
        executed_original_helper_sha256=rules.sha(root/'queue_continuation.py'),
        reviewed_handoff_rules_sha256=rules.sha(Path(rules.__file__)),
        executed_reconciliation_sha256=rules.sha(Path(__file__)),
        audit_receipt_sha256=rules.sha(root/'audit_receipt.json'),
        original_journal_sha256=rules.sha(root/'submission_journal.jsonl'),
        reconciliation_journal_sha256=rules.sha(root/'reconciliation_journal.jsonl'))
    with (root/'submission_receipt.json').open('x') as stream:
        json.dump(result, stream, indent=2); stream.write('\n')
    print(json.dumps(result))


if __name__ == '__main__':
    reconcile()
