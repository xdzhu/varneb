"""Complete the proven partial handoff; no new submission or repeated update."""
from datetime import datetime
import json
import os
from pathlib import Path
import reviewed_handoff_rules_v2 as rules


def finish():
    root, handle = rules.ROOT, '28722320'
    original = [json.loads(r) for r in (root/'submission_journal.jsonl').read_text().splitlines()]
    partial = [json.loads(r) for r in (root/'reconciliation_journal.jsonl').read_text().splitlines()]
    if ([r['job_id'] for r in original if r['event']=='accepted_held_handle'] != [handle]
            or partial[-1]['event'] != 'dependency_update_intent' or partial[-1]['job_id'] != '28692775'
            or rules.sha(root/'queue_continuation.py') != 'fd6ee53777bc8b13b4fad37e6b572ad01824de916295ad9aef797d8d4f30ae5b'
            or rules.sha(rules.SCRIPT) != rules.SCRIPT_SHA):
        raise ValueError('original single accepted handle/partial handoff evidence differs')
    held = rules.fields(rules.run(['scontrol','show','job','-o',handle]))
    if (held.get('JobState') != 'PENDING' or held.get('Reason') != 'JobHeldUser'
            or held.get('JobName') != 'hfo2-G2-M-E062' or held.get('NumCPUs') != '32'
            or held.get('Command') != str(rules.SCRIPT) or held.get('Partition') != 'hfacnormal01'):
        raise ValueError('same held job identity changed')
    receipt = json.loads((root/'audit_receipt.json').read_text())
    seed = root/'seed'
    manifest = json.loads((seed/'manifest.json').read_text())
    runtime = root.parent/'20261010-varneb-clamped-G2-E054-r1/source-fixed'
    if (receipt['seed_manifest_sha256'] != rules.sha(seed/'manifest.json')
            or any(rules.sha(seed/n) != h for n,h in manifest['files_sha256'].items())
            or any(rules.sha(runtime/n) != h for n,h in receipt['runtime_code_sha256'].items())
            or rules.sha(rules.BINARY) != rules.BINARY_SHA):
        raise ValueError('audited data/runtime/binary changed')
    expected = {'28661019',handle}
    first = rules.confirm_waiter('28692775', expected)
    second = rules.validate_waiter(rules.run(['scontrol','show','job','-o','28692776']), '28692776', {'28661019'})
    rules.validate_active_queue(rules.run(['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C']), rules.ALLOWED|{handle})
    dependency = 'afterok:28661019:'+handle
    with (root/'finish_journal.jsonl').open('x') as journal:
        def record(event, **values):
            journal.write(json.dumps(dict(event=event, check_CST=datetime.now().isoformat(), **values))+'\n')
            journal.flush(); os.fsync(journal.fileno())
        record('first_update_already_visible_NOT_repeated', fields=first, second_before=second)
        record('missing_second_update_intent', job_id='28692776', dependency=dependency)
        rules.run(['scontrol','update','JobId=28692776','Dependency='+dependency])
        for job in rules.WAITING:
            checked = rules.confirm_waiter(job, expected)
            record('dependency_update_verified', job_id=job, actual_dependency=checked['Dependency'])
        rules.validate_active_queue(rules.run(['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C']), rules.ALLOWED|{handle})
        record('release_intent', job_id=handle)
        rules.run(['scontrol','release',handle])
        record('released', job_id=handle)
    result = dict(status='same_accepted_continuation_released_after_partial_handoff_reconciliation',
        check_CST=datetime.now().isoformat(), job_id=handle, parent_job_id='28661020',
        total_E062_submissions=1, new_submissions_during_reconciliation=0,
        first_dependency_update_repeated=False, dependency=dependency,
        waiting_jobs_verified=list(rules.WAITING), max_study_simultaneous_allocations=2,
        new_independent_chains=0, new_segment_steps=20, walltime_hours=6, MPI_ranks=32, ordinary_fmax=.10,
        physical_inputs_changed=False, running_source_overwritten=False, holdout_generated_or_read=False,
        executed_original_helper_sha256=rules.sha(root/'queue_continuation.py'),
        reviewed_handoff_rules_sha256=rules.sha(Path(rules.__file__)),
        executed_finish_sha256=rules.sha(Path(__file__)), audit_receipt_sha256=rules.sha(root/'audit_receipt.json'),
        original_journal_sha256=rules.sha(root/'submission_journal.jsonl'),
        partial_journal_sha256=rules.sha(root/'reconciliation_journal.jsonl'),
        finish_journal_sha256=rules.sha(root/'finish_journal.jsonl'))
    with (root/'submission_receipt.json').open('x') as stream:
        json.dump(result, stream, indent=2); stream.write('\n')
    print(json.dumps(result))


if __name__ == '__main__':
    finish()
