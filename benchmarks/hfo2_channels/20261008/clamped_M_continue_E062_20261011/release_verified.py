"""Release the original held job only after both prior writes are visible."""
from datetime import datetime
import json
from pathlib import Path
import reviewed_handoff_rules_v2 as rules


def release():
    root, handle = rules.ROOT, '28722320'
    original = [json.loads(r) for r in (root/'submission_journal.jsonl').read_text().splitlines()]
    partial = [json.loads(r) for r in (root/'finish_journal.jsonl').read_text().splitlines()]
    if ([r['job_id'] for r in original if r['event']=='accepted_held_handle'] != [handle]
            or not any(r['event']=='missing_second_update_intent' for r in partial)
            or any(r['event']=='release_intent' for r in partial)):
        raise ValueError('original single accepted job/partial release history differs')
    held = rules.fields(rules.run(['scontrol','show','job','-o',handle]))
    if (held.get('JobState') != 'PENDING' or held.get('Reason') != 'JobHeldUser'
            or held.get('JobName') != 'hfo2-G2-M-E062' or held.get('NumCPUs') != '32'
            or held.get('Command') != str(rules.SCRIPT) or held.get('Partition') != 'hfacnormal01'
            or rules.sha(rules.SCRIPT) != rules.SCRIPT_SHA or rules.sha(rules.BINARY) != rules.BINARY_SHA):
        raise ValueError('same accepted held job identity/binary changed')
    receipt = json.loads((root/'audit_receipt.json').read_text())
    seed = root/'seed'
    manifest = json.loads((seed/'manifest.json').read_text())
    runtime = root.parent/'20261010-varneb-clamped-G2-E054-r1/source-fixed'
    if (receipt['seed_manifest_sha256'] != rules.sha(seed/'manifest.json')
            or any(rules.sha(seed/n) != h for n,h in manifest['files_sha256'].items())
            or any(rules.sha(runtime/n) != h for n,h in receipt['runtime_code_sha256'].items())):
        raise ValueError('audited data/runtime changed')
    waiters = {job:rules.validate_waiter(rules.run(['scontrol','show','job','-o',job]), job,
                                        {'28661019',handle}) for job in rules.WAITING}
    active = rules.validate_active_queue(rules.run(['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C']), rules.ALLOWED|{handle})
    intent = dict(check_CST=datetime.now().isoformat(), job_id=handle, actual_waiters=waiters,
                  actual_held_fields=held, active_study_jobs=active, audit_receipt_sha256=rules.sha(root/'audit_receipt.json'))
    with (root/'release_intent.json').open('x') as stream:
        json.dump(intent, stream, indent=2); stream.write('\n')
        stream.flush()
        import os
        os.fsync(stream.fileno())
    rules.run(['scontrol','release',handle])
    result = dict(status='release_acknowledged_same_original_handle', check_CST=datetime.now().isoformat(),
        job_id=handle, parent_job_id='28661020', waiting_jobs_verified=list(rules.WAITING),
        total_E062_submissions=1, new_submissions_during_reconciliation=0, dependency_updates_repeated=False,
        max_study_simultaneous_allocations=2, new_independent_chains=0,
        new_segment_steps=20, walltime_hours=6, MPI_ranks=32, ordinary_fmax=.10,
        physical_inputs_changed=False, running_source_overwritten=False, holdout_generated_or_read=False,
        execution_started_claimed=False, executed_release_sha256=rules.sha(Path(__file__)),
        reviewed_handoff_rules_sha256=rules.sha(Path(rules.__file__)),
        audit_receipt_sha256=rules.sha(root/'audit_receipt.json'), release_intent_sha256=rules.sha(root/'release_intent.json'))
    with (root/'submission_receipt.json').open('x') as stream:
        json.dump(result, stream, indent=2); stream.write('\n')
    print(json.dumps(result))


if __name__ == '__main__':
    release()
