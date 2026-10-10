"""Rewire two held successors once and release only the existing continuation."""
from datetime import datetime
import json
import os
from pathlib import Path
import time
import submit_held_once as rules


def acknowledge(job, expected, *, lookup=rules.run, pause=time.sleep):
    for delay in (0, 1, 2, 4, 8):
        if delay:
            pause(delay)
        body = lookup(['scontrol','show','job','-o',job])
        row = rules.held_waiter(body, job, required={rules.OTHER}, allowed=expected|{rules.PARENT})
        if rules.dependencies(row['Dependency']) == expected:
            return row
        if rules.dependencies(row['Dependency']) - {rules.OTHER,rules.PARENT}:
            raise ValueError('partial/unexpected dependency acknowledgement; reconcile without repeat')
    raise TimeoutError('stale dependency acknowledgement; never repeat accepted update')


def handoff_once():
    root = rules.ROOT
    rules.verify_audit_seed()
    accepted = json.loads((root/'accepted_receipt.json').read_text())
    handle = accepted['job_id']
    if (not handle.isdecimal() or accepted['total_submissions'] != 1
            or accepted['status'] != 'one_accepted_held_same_chain_continuation_not_released'
            or accepted['parent_job_id'] != rules.PARENT or accepted['script_sha256'] != rules.SCRIPT_SHA
            or accepted['dependencies_changed'] or accepted['execution_started_claimed']
            or accepted['audit_receipt_sha256'] != rules.sha(root/'audit_receipt.json')
            or accepted['executed_helper_sha256'] != rules.sha(root/'submit_held_once.py')):
        raise ValueError('accepted same-chain handoff receipt changed')
    rules.validate_held_continuation(rules.run(['scontrol','show','job','-o',handle]), handle)
    for job in rules.WAITERS:
        rules.held_waiter(rules.run(['scontrol','show','job','-o',job]), job,
                          required={rules.OTHER}, allowed={rules.OTHER,rules.PARENT})
    rules.active_queue(rules.run(['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C']), allowed=rules.ALLOWED|{handle})
    expected = {rules.OTHER,handle}
    completed = []
    with (root/'handoff_journal.jsonl').open('x') as journal:
        def record(event, **data):
            journal.write(json.dumps(dict(event=event, check_CST=datetime.now().isoformat(), **data))+'\n')
            journal.flush()
            os.fsync(journal.fileno())
        try:
            for job in rules.WAITERS:
                rules.held_waiter(rules.run(['scontrol','show','job','-o',job]), job,
                                  required={rules.OTHER}, allowed={rules.OTHER,rules.PARENT})
                dependency = 'afterok:'+rules.OTHER+':'+handle
                record('dependency_update_intent', job_id=job, dependency=dependency)
                rules.run(['scontrol','update','JobId='+job,'Dependency='+dependency])
                row = acknowledge(job, expected)
                completed.append(job)
                record('dependency_update_verified', job_id=job, actual=row)
            # Both successors stay HELD after this operation: reaching the
            # next step cap is still not a convergence proof.
            for job in rules.WAITERS:
                rules.held_waiter(rules.run(['scontrol','show','job','-o',job]), job,
                                  required=expected, allowed=expected)
            rules.validate_held_continuation(rules.run(['scontrol','show','job','-o',handle]), handle)
            rules.active_queue(rules.run(['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C']), allowed=rules.ALLOWED|{handle})
            record('continuation_only_release_intent', job_id=handle, both_successors_still_held=True)
            rules.run(['scontrol','release',handle])
            record('continuation_only_released', job_id=handle)
        except Exception as exc:
            record('reconcile_do_not_repeat', exception=type(exc).__name__, message=str(exc),
                   dependency_writes_verified=completed)
            raise
    receipt = dict(status='continuation_released_successors_still_held_for_material_audit',
        check_CST=datetime.now().isoformat(), job_id=handle, parent_job_id=rules.PARENT,
        new_submissions_during_handoff=0, dependency_writes_verified=completed,
        successors_still_held=list(rules.WAITERS), registered_dependency=sorted(expected),
        max_study_simultaneous_allocations=2, new_independent_chains=0,
        running_source_or_physical_inputs_changed=False, holdout_generated_or_read=False,
        execution_started_claimed=False, accepted_receipt_sha256=rules.sha(root/'accepted_receipt.json'),
        executed_handoff_sha256=rules.sha(Path(__file__)), journal_sha256=rules.sha(root/'handoff_journal.jsonl'))
    with (root/'handoff_receipt.json').open('x') as stream:
        json.dump(receipt, stream, indent=2)
        stream.write('\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    handoff_once()
