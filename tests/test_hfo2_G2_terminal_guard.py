"""Offline tests of exact-target scheduler guard; never invoke Slurm or DFT."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'benchmarks/hfo2_channels/20261008/clamped_terminal_guard_E066_20261011'
spec = importlib.util.spec_from_file_location('terminal_guard', CASE/'hold_waiters_HF.py')
guard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(guard)


def body(job, *, held=False):
    channel = guard.WAITERS[job]
    return (f'JobId={job} JobName=hfo2-G2-E058-{channel} UserId=iai806(16284) '
            f'JobState=PENDING Partition=hfacnormal01 NumCPUs=32 Priority={0 if held else 1167} '
            f'Reason={"JobHeldUser" if held else "Dependency"} '
            'Dependency=afterok:28661019(unfulfilled),afterok:28722320(unfulfilled) '
            f'Command={guard.PILOT} StdOut={guard.OUTPUT_ROOT}/{channel}.slurm.out '
            f'StdErr={guard.OUTPUT_ROOT}/{channel}.slurm.err')


QUEUE = ('28661019|hfo2-G2-flip|RUNNING|hfacnormal01|32\n'
         '28722320|hfo2-G2-M|RUNNING|hfacnormal01|32\n'
         '28692775|hfo2-G2-E058-PO_to_T|PENDING|hfacnormal01|32\n'
         '900|pc-other|RUNNING|hfacnormal01|64\n')


def fake_scheduler(*, ambiguous_job=None):
    calls, held = [], set()
    def lookup(argv):
        calls.append(argv)
        if argv[0] == 'squeue':
            return QUEUE
        job = argv[-1]
        if argv[1] == 'show':
            return body(job, held=job in held)
        if argv[1] == 'hold':
            held.add(job)
            if job == ambiguous_job:
                raise TimeoutError('accepted write but lost response')
            return ''
        raise AssertionError('unexpected command')
    return lookup, calls, held


def test_single_use_actual_sequence_never_submits_releases_or_changes_dependencies(tmp_path):
    lookup, calls, held = fake_scheduler()
    receipt = guard.hold_once(tmp_path, lookup=lookup, pause=lambda n: None)
    assert held == set(guard.WAITERS) and receipt['new_submissions'] == 0
    mutations = [argv for argv in calls if argv[:2] == ['scontrol', 'hold']]
    assert mutations == [['scontrol', 'hold', job] for job in guard.WAITERS]
    assert not receipt['dependencies_changed'] and not receipt['running_jobs_changed']
    journal = [json.loads(r) for r in (tmp_path/'hold_journal.jsonl').read_text().splitlines()]
    assert [r['event'] for r in journal] == ['guard_intent', 'hold_intent', 'hold_verified', 'hold_intent', 'hold_verified']
    assert receipt['journal_sha256'] == guard.sha(tmp_path/'hold_journal.jsonl')
    n = len(mutations)
    with pytest.raises(ValueError):
        guard.hold_once(tmp_path, lookup=lookup, pause=lambda n: None)
    assert len([a for a in calls if a[:2] == ['scontrol', 'hold']]) == n


@pytest.mark.parametrize('before,after', [
    ('JobState=PENDING', 'JobState=RUNNING'), ('UserId=iai806(', 'UserId=someone('),
    ('NumCPUs=32', 'NumCPUs=64'), ('Partition=hfacnormal01', 'Partition=other'),
    ('Priority=1167', 'Priority=0'), ('Reason=Dependency', 'Reason=JobHeldUser'),
    ('afterok:28661019', 'afterany:28661019'),
    (',afterok:28722320', '?afterok:28722320'),
    ('28722320(unfulfilled)', '999(unfulfilled)'),
    ('Command='+guard.PILOT, 'Command=/other'),
    ('JobName=hfo2-G2-E058-PO_to_T', 'JobName=pc-other'),
    ('PO_to_T.slurm.out', 'other.slurm.out'),
])
def test_changed_identity_or_dependency_refuses_mutation(tmp_path, before, after):
    calls = []
    def lookup(argv):
        calls.append(argv)
        return QUEUE if argv[0] == 'squeue' else body(argv[-1]).replace(before, after)
    with pytest.raises(ValueError):
        guard.hold_once(tmp_path, lookup=lookup, pause=lambda n: None)
    assert not any(a[:2] == ['scontrol', 'hold'] for a in calls)
    assert not (tmp_path/'hold_journal.jsonl').exists()


@pytest.mark.parametrize('mutated', [
    QUEUE.replace('28661019|hfo2-G2-flip|RUNNING', '28661019|hfo2-G2-flip|PENDING'),
    QUEUE.replace('28722320|hfo2-G2-M|RUNNING', '28722320|hfo2-G2-M|PENDING'),
    QUEUE+'999|hfo2-unknown|RUNNING|hfacnormal01|32\n',
    QUEUE.replace('hfacnormal01|32', 'hfacnormal01|64', 1),
    QUEUE.replace('28692775|hfo2-G2-E058-PO_to_T|PENDING', '28692775|hfo2-G2-E058-PO_to_T|RUNNING'),
])
def test_live_parent_or_allocation_guard(mutated):
    with pytest.raises(ValueError):
        guard.queue(mutated)


@pytest.mark.parametrize('job', list(guard.WAITERS))
def test_ambiguous_write_is_never_repeated_and_partial_evidence_survives(tmp_path, job):
    lookup, calls, held = fake_scheduler(ambiguous_job=job)
    with pytest.raises(TimeoutError):
        guard.hold_once(tmp_path, lookup=lookup, pause=lambda n: None)
    assert job in held and calls.count(['scontrol', 'hold', job]) == 1
    rows = [json.loads(r) for r in (tmp_path/'hold_journal.jsonl').read_text().splitlines()]
    assert rows[-1]['event'] == 'reconcile_do_not_rerun'
    assert not (tmp_path/'hold_receipt.json').exists()
    count = len([a for a in calls if a[:2] == ['scontrol', 'hold']])
    with pytest.raises(ValueError):
        guard.hold_once(tmp_path, lookup=lookup, pause=lambda n: None)
    assert len([a for a in calls if a[:2] == ['scontrol', 'hold']]) == count


def test_bounded_acknowledgement_is_read_only_and_rejects_other_changes():
    replies = iter([body('28692775'), body('28692775', held=True)])
    calls, pauses = [], []
    def lookup(argv):
        calls.append(argv)
        return next(replies)
    row = guard.acknowledge('28692775', lookup=lookup, pause=pauses.append)
    assert row['Priority'] == '0' and pauses == [1]
    assert all(a[:3] == ['scontrol', 'show', 'job'] for a in calls)
    with pytest.raises(TimeoutError):
        guard.acknowledge('28692775', lookup=lambda a: body('28692775'), pause=lambda n: None)
    with pytest.raises(ValueError):
        guard.acknowledge('28692775', lookup=lambda a: body('28692775', held=True).replace('NumCPUs=32','NumCPUs=64'), pause=lambda n: None)


def test_existing_journal_cannot_be_overwritten(tmp_path):
    (tmp_path/'hold_journal.jsonl').write_text('preserve accepted operation\n')
    lookup, calls, _ = fake_scheduler()
    with pytest.raises(FileExistsError):
        guard.hold_once(tmp_path, lookup=lookup, pause=lambda n: None)
    assert not any(a[:2] == ['scontrol', 'hold'] for a in calls)
    assert (tmp_path/'hold_journal.jsonl').read_text() == 'preserve accepted operation\n'


def test_actual_HF_guard_receipt_preserves_two_live_parents_and_exact_dependencies():
    receipt = json.loads((CASE/'hold_receipt.json').read_text())
    verification = json.loads((CASE/'independent_verification.json').read_text())
    rows = [json.loads(r) for r in (CASE/'hold_journal.jsonl').read_text().splitlines()]
    assert receipt['helper_sha256'] == guard.sha(CASE/'executed_hold_waiters_HF.py')
    assert receipt['helper_sha256'] == 'f231b82c3826af4c0b457d5913541df2c65947a21e4938114b0507753d7e7832'
    assert receipt['journal_sha256'] == guard.sha(CASE/'hold_journal.jsonl')
    assert receipt['held_jobs'] == list(guard.WAITERS)
    assert [r['job_id'] for r in rows if r['event'] == 'hold_intent'] == list(guard.WAITERS)
    assert [r['job_id'] for r in rows if r['event'] == 'hold_verified'] == list(guard.WAITERS)
    assert verification['check_CST'] > receipt['check_CST']
    assert guard.queue(verification['raw_registered_queue']) == sorted(guard.PARENTS)
    for job, raw in verification['raw_waiter_rows'].items():
        guard.waiter(raw, job, held=True)
    assert len(verification['raw_parent_accounting'].splitlines()) == 2
    assert all('|RUNNING|0:0|' in r and r.endswith('|32|hfacnormal01')
               for r in verification['raw_parent_accounting'].splitlines())
    assert receipt['new_submissions'] == receipt['new_DFT_calls'] == 0
    assert verification['new_scheduler_mutations'] == verification['new_DFT_calls'] == 0
    assert not any(receipt[k] for k in ('dependencies_changed', 'running_jobs_changed',
                                      'physical_inputs_changed', 'holdout_generated_or_read'))
