"""Read-reconcile two already accepted dependency writes; release existing ID only."""
from datetime import datetime
import json
import os
from pathlib import Path
import submit_held_once as rules


def release_once():
    root = rules.ROOT
    if (root/'verified_release_journal.jsonl').exists() or (root/'handoff_receipt.json').exists():
        raise FileExistsError('existing release evidence; do not repeat a mutation')
    rules.verify_audit_seed()
    accepted = json.loads((root/'accepted_receipt.json').read_text())
    original = [json.loads(r) for r in (root/'handoff_journal.jsonl').read_text().splitlines()]
    handle = accepted['job_id']
    intents = [r for r in original if r['event'] == 'dependency_update_intent']
    if (accepted['total_submissions'] != 1 or accepted['parent_job_id'] != rules.PARENT
            or accepted['audit_receipt_sha256'] != rules.sha(root/'audit_receipt.json')
            or accepted['executed_helper_sha256'] != rules.sha(root/'submit_held_once.py')
            or [r['job_id'] for r in intents] != list(rules.WAITERS)
            or any(r['dependency'] != 'afterok:'+rules.OTHER+':'+handle for r in intents)
            or original[-1]['event'] != 'reconcile_do_not_repeat'
            or any('release' in r['event'] for r in original)):
        raise ValueError('partial original handoff evidence differs; do not release')
    rules.validate_held_continuation(rules.run(['scontrol','show','job','-o',handle]), handle)
    expected = {rules.OTHER,handle}
    verified = {job: rules.held_waiter(rules.run(['scontrol','show','job','-o',job]), job,
                                     required=expected, allowed=expected) for job in rules.WAITERS}
    active = rules.active_queue(rules.run(['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C']),
                               allowed=rules.ALLOWED|{handle})
    with (root/'verified_release_journal.jsonl').open('x') as journal:
        def record(event, **data):
            journal.write(json.dumps(dict(event=event, check_CST=datetime.now().isoformat(), **data))+'\n')
            journal.flush()
            os.fsync(journal.fileno())
        record('independently_verified_release_intent', job_id=handle,
               both_successors_still_held=verified, active_study_jobs=active,
               original_handoff_journal_sha256=rules.sha(root/'handoff_journal.jsonl'))
        try:
            rules.run(['scontrol','release',handle])
            record('existing_continuation_only_released', job_id=handle)
        except Exception as exc:
            record('reconcile_release_do_not_repeat', exception=type(exc).__name__, message=str(exc))
            raise
    receipt = dict(status='existing_continuation_released_after_independent_dependency_reads',
        check_CST=datetime.now().isoformat(), job_id=handle, parent_job_id=rules.PARENT,
        total_E067_submissions=1, new_submissions_during_reconciliation=0,
        dependency_writes_during_reconciliation=0, successors_still_held=list(rules.WAITERS),
        registered_dependency=sorted(expected), max_study_simultaneous_allocations=2,
        new_independent_chains=0, running_source_or_physical_inputs_changed=False,
        holdout_generated_or_read=False, execution_started_claimed=False,
        accepted_receipt_sha256=rules.sha(root/'accepted_receipt.json'),
        executed_release_sha256=rules.sha(Path(__file__)),
        original_handoff_journal_sha256=rules.sha(root/'handoff_journal.jsonl'),
        verified_release_journal_sha256=rules.sha(root/'verified_release_journal.jsonl'))
    with (root/'handoff_receipt.json').open('x') as stream:
        json.dump(receipt, stream, indent=2)
        stream.write('\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    release_once()
