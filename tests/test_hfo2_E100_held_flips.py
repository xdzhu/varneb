"""E100 control fixtures are not native convergence or SCF evidence."""
import importlib.util
import json
from pathlib import Path
import shlex

import pytest

REPO = Path(__file__).resolve().parents[1]
CASE = REPO/'benchmarks/hfo2_channels/20261008/clamped_held_flips_E100_20261011'


def module():
    spec = importlib.util.spec_from_file_location('E100', CASE/'release_existing_once.py')
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def control(m, job='28709790', **changes):
    ch = m.CHANNELS[job]
    row = dict(JobId=job, JobName='hfo2-G2-E060-'+ch, JobState='PENDING',
        Priority='0', Reason='JobHeldUser', Partition='hfacnormal01', NumCPUs='32',
        NumTasks='1', Command=str(m.PILOT), StdOut=str(m.WAIT/(ch+'.slurm.out')),
        StdErr=str(m.WAIT/(ch+'.slurm.err')), TimeLimit='04:00:00', Dependency='(null)',
        UserId='iai806(1000)', NumNodes='1', Restarts='0', **{'CPUs/Task': '32'})
    row.update(changes)
    return ' '.join(k+'='+shlex.quote(v) for k, v in row.items())


def queue(m, jobs=('28727368', '28727382')):
    rows = [j+'|hfo2-G2-E091-'+m.LIVE[j]+'|RUNNING|hfacnormal01|32|node6' for j in jobs]
    rows += [j+'|hfo2-G2-E060-'+c+'|PENDING|hfacnormal01|32|(JobHeldUser)' for j, c in m.CHANNELS.items()]
    return '\n'.join(rows)


@pytest.mark.parametrize('job', ['28709790', '28709791'])
def test_exact_original_held_identity(job):
    m = module()
    assert m.held_identity(control(m, job), job)['JobId'] == job


@pytest.mark.parametrize('key', ['JobName', 'JobState', 'Reason', 'Priority', 'Partition',
    'NumCPUs', 'NumTasks', 'CPUs/Task', 'Command', 'StdOut', 'StdErr', 'TimeLimit',
    'Dependency', 'UserId', 'NumNodes', 'Restarts'])
def test_changed_held_contract_rejected(key):
    m = module()
    with pytest.raises(ValueError):
        m.held_identity(control(m, **{key: 'different'}), '28709790')


def test_readiness_accepts_two_used_slots_but_release_cannot():
    m = module()
    assert len(m.committed_slots(queue(m), require_spare=False)) == 2
    with pytest.raises(ValueError, match='no free committed slot'):
        m.committed_slots(queue(m), require_spare=True)
    assert m.committed_slots(queue(m, ('28727368',)), require_spare=True) == ['28727368']


def test_unheld_pending_reserves_slot_and_unrelated_job_untouched():
    m = module()
    body = queue(m, ('28727368',)).replace('(JobHeldUser)', '(Resources)', 1)
    with pytest.raises(ValueError, match='no free committed slot'):
        m.committed_slots(body, require_spare=True)
    assert m.committed_slots(queue(m, ())+'\npc123|pc-other|RUNNING|other|128|node', require_spare=True) == []


@pytest.mark.parametrize('extra', [
    'unknown|hfo2-foreign|RUNNING|hfacnormal01|32|node',
    '28727368|hfo2-G2-E091-28709789|RUNNING|hfacnormal01|32|node',
    'bad row'])
def test_unknown_duplicate_malformed_queue_rejected(extra):
    m = module()
    with pytest.raises(ValueError):
        m.committed_slots(queue(m)+'\n'+extra, require_spare=False)


def test_unknown_handles_reject_before_io(monkeypatch):
    m = module()
    def forbidden(*a, **kw):
        raise AssertionError('no native IO for unregistered handle')
    monkeypatch.setattr(m, 'run', forbidden)
    with pytest.raises(ValueError):
        m.preflight('28709788')
    with pytest.raises(ValueError):
        m.release_once('new', 'digest')


def fixture(tmp_path, monkeypatch):
    m = module()
    monkeypatch.setattr(m, 'ROOT', tmp_path)
    calls = []
    monkeypatch.setattr(m, 'verified_terminal', lambda *a: calls.append(('audit', a)))
    monkeypatch.setattr(m, 'preflight', lambda *a, **kw: {'fixture_not_material': True})
    released = []
    def run(argv):
        calls.append(tuple(argv))
        if argv[0] == 'squeue':
            body = queue(m, ('28727368',))
            return body.replace('(JobHeldUser)', '(Priority)', 1) if released else body
        if argv[:3] == ['scontrol', 'show', 'job']:
            return control(m, Priority='1', Reason='Priority') if released else control(m)
        if argv == ['scontrol', 'release', '28709790']:
            released.append(True)
            return ''
        raise AssertionError(argv)
    monkeypatch.setattr(m, 'run', run)
    return m, calls


def test_one_release_no_submit_update_and_no_replay(tmp_path, monkeypatch):
    m, calls = fixture(tmp_path, monkeypatch)
    result = m.release_once('28709790', 'fixture_digest')
    assert result['new_submissions'] == result['dependency_writes'] == 0
    assert result['release_writes'] == 1
    assert not (tmp_path/'admission.lock').exists()
    before = list(calls)
    with pytest.raises(FileExistsError): m.release_once('28709790', 'fixture_digest')
    assert calls == before
    events = [json.loads(s)['event'] for s in (tmp_path/'28709790.journal.jsonl').read_text().splitlines()]
    assert events == ['verified_original_held_target', 'one_existing_release_intent',
        'one_existing_release_accepted', 'released_slot_independently_corroborated']
    assert not any(c[0] == 'sbatch' or c[:2] == ('scontrol', 'update') for c in calls)


def test_ambiguous_write_keeps_lock_and_journal_and_never_retries(tmp_path, monkeypatch):
    m, calls = fixture(tmp_path, monkeypatch)
    original = m.run
    def ambiguous(argv):
        if argv[:2] == ['scontrol', 'release']:
            raise TimeoutError('reply missing')
        return original(argv)
    monkeypatch.setattr(m, 'run', ambiguous)
    with pytest.raises(TimeoutError): m.release_once('28709790', 'fixture_digest')
    assert (tmp_path/'admission.lock').exists()
    assert 'one_existing_release_intent' in (tmp_path/'28709790.journal.jsonl').read_text()
    with pytest.raises(FileExistsError): m.release_once('28709791', 'fixture_digest')


def test_capacity_change_before_write_stops_release(tmp_path, monkeypatch):
    m, calls = fixture(tmp_path, monkeypatch)
    monkeypatch.setattr(m, 'run', lambda argv: queue(m))
    with pytest.raises(ValueError, match='no free committed slot'):
        m.release_once('28709790', 'fixture_digest')
    assert 'one_existing_release_intent' not in (tmp_path/'28709790.journal.jsonl').read_text()
    assert (tmp_path/'admission.lock').exists()


def test_empty_release_reply_without_native_change_keeps_shared_lock(tmp_path, monkeypatch):
    m, calls = fixture(tmp_path, monkeypatch)
    original = m.run
    def stale(argv):
        if argv[:2] == ['scontrol', 'release']: return ''
        return original(argv)
    monkeypatch.setattr(m, 'run', stale)
    with pytest.raises(ValueError, match='independently corroborated'):
        m.release_once('28709790', 'fixture_digest')
    assert (tmp_path/'admission.lock').exists()
    with pytest.raises(FileExistsError): m.release_once('28709791', 'fixture_digest')


def test_terminal_cap_is_not_convergence(tmp_path, monkeypatch):
    m = module()
    receipt = tmp_path/'audit.json'
    receipt.write_text(json.dumps({'job_id': '28727382', 'ordinary_converged': False,
        'status': 'audited_step_cap_not_converged_same_chain_seed'}))
    monkeypatch.setattr(m, 'run', lambda *a: pytest.fail('reject before native terminal query'))
    with pytest.raises(ValueError, match='ordinary-converged'):
        m.verified_terminal(receipt, m.sha(receipt), '28727382')


@pytest.mark.parametrize('change', [None, 'SOURCE_ROOT', 'SEED_ROOT', 'WORKDIR', 'RUN_DFT',
    'physical_override', 'duplicate_export', 'wrapper'])
def test_original_submit_environment_checked_before_seed_io(change, monkeypatch):
    m = module()
    ch = m.CHANNELS['28709790']
    env = dict(RUN_DFT='1', SOURCE_ROOT=str(m.SOURCE),
        SEED_ROOT=str(m.RUNS/'20261010-varneb-clamped-G2-remaining-E056-r1/seeds/strain_p0100'/ch),
        WORKDIR=str(m.WAIT/ch/'band'))
    if change in env: env[change] = 'different'
    if change == 'physical_override': env['ECUTWFC'] = '120'
    submit = 'sbatch '+shlex.quote('--export=ALL,'+','.join(k+'='+v for k, v in env.items()))+' '
    submit += shlex.quote(str(m.PILOT) if change != 'wrapper' else '/different/wrapper')
    if change == 'duplicate_export': submit = submit.replace('sbatch ', 'sbatch --export=ALL ', 1)
    def run(argv):
        if argv[0] == 'scontrol': return control(m)
        if argv[0] == 'squeue': return queue(m)
        if argv[-1] == '--format=JobIDRaw,SubmitLine%4000': return '28709790|'+submit
        return '28709790|PENDING|0:0|0|hfacnormal01|iai806'
    monkeypatch.setattr(m, 'run', run)
    monkeypatch.setattr(m, 'pending_readiness', lambda job: {'fixture_seed': True})
    if change is None:
        assert m.preflight('28709790')['scheduler_writes'] == 0
    else:
        with pytest.raises(ValueError, match='submitted source/seed/environment'):
            m.preflight('28709790')


@pytest.fixture(scope='module')
def actual_preserving_terminal():
    from scripts.analyze_hfo2_clamped_reference import terminal_images
    return terminal_images(CASE/'zero_strain_preserving',
        '46145da14f24b2b8bd46aa4656f0982650d6b2dc661282f58d809aab8a6d3202')


def test_actual_terminal_full_native_audit_and_segment_cost(actual_preserving_terminal):
    from examples.hfo2_fixed_input_factory import CONTRACT
    _, report, residual, audit = actual_preserving_terminal
    assert audit['job_id'] == '28727382' and audit['terminal_step'] == 15
    assert audit['scheduler_row'] == '28727382|COMPLETED|0:0|32|hfacnormal01|iai806'
    assert audit['initial_fresh_interior_SCFs'] == 0
    assert audit['fresh_interior_SCFs'] == audit['expected_fresh_interior_SCFs'] == 105
    assert len(audit['runtime_code_sha256']) == 334
    calls = audit['calls']
    assert len(calls) == 105
    assert {(Path(r['path']).parts[-3], Path(r['path']).parts[-2]) for r in calls} == {
        (f'image_{i:04d}', f'scf_{j:06d}') for i in range(1, 8) for j in range(15)}
    assert all({k: r['input_sha256'][k] for k in CONTRACT} == CONTRACT for r in calls)
    assert sum(r['elapsed_seconds'] for r in calls) == pytest.approx(11800.9455248853, abs=1e-9)
    assert audit['fresh_SCF_core_hour_proxy_32'] == pytest.approx(104.897293554536, abs=1e-9)
    assert audit['new_DFT_calls'] == audit['scheduler_mutations'] == 0
    assert audit['continuous_TS_or_barrier_error_bound_certified'] is False
    assert residual['ordinary_residual_pass'] is True


def test_actual_energy_units_boundary_and_shared_initial_well(actual_preserving_terminal):
    import numpy as np
    from examples.hfo2_fixed_input_factory import same_ordered_geometry
    from scripts.analyze_hfo2_clamped_reference import terminal_images
    images, report, residual, _ = actual_preserving_terminal
    assert len(images) == report['n_total_images'] == 9 and report['n_active_images'] == 7
    assert report['formula_units'] == 4 and report['pressure_GPa'] == 0
    assert report['climb'] is False and report['fmax_target_eV_A'] == .10
    assert report['mechanical_boundary']['fixed_cell_rows'] == [0, 1]
    assert report['mechanical_boundary']['allow_tilt'] is True
    assert residual['fmax_eV_A'] == pytest.approx(.0988168085874911, abs=1e-12)
    energy = np.array([a.get_potential_energy() for a in images])
    profile = (energy-energy[0])*1000/4
    assert np.argmax(profile) == 4
    assert np.allclose(profile, report['relative_enthalpy_meV_fu'], atol=1e-9, rtol=0)
    assert max(profile) == pytest.approx(96.82818300279905, abs=1e-9)
    assert max(profile)-profile[-1] == pytest.approx(96.86344180590822, abs=1e-9)
    t_images, t_report, _, _ = terminal_images(CASE.parent/'clamped_T_terminal_E079_20261011',
        'c5653b4759b202199cb75739f99038aaa52b3ccc220ef226efbacc980e5d716e')
    assert same_ordered_geometry(images[0], t_images[0])
    assert energy[0] == t_images[0].get_potential_energy()
    for key in ('reference_cell_A', 'fixed_cell_rows', 'allow_tilt', 'raw_chain_policy'):
        assert report['mechanical_boundary'][key] == t_report['mechanical_boundary'][key]


def test_actual_single_release_journal_and_source_bytes():
    m = module()
    result = json.loads((CASE/'28709790.receipt.json').read_text())
    journal = CASE/'28709790.journal.jsonl'
    rows = [json.loads(s) for s in journal.read_text().splitlines()]
    assert result['job_id'] == '28709790' and result['release_writes'] == 1
    assert result['new_submissions'] == result['dependency_writes'] == result['new_independent_chains'] == 0
    assert result['executed_helper_sha256'] == m.sha(CASE/'release_existing_once.py')
    assert result['journal_sha256'] == m.sha(journal)
    assert rows[0]['proof']['committed_study_slots'] == ['28727368']
    assert rows[0]['preserving_audit_sha256'] == '46145da14f24b2b8bd46aa4656f0982650d6b2dc661282f58d809aab8a6d3202'
    assert [r['event'] for r in rows] == ['verified_original_held_target',
        'one_existing_release_intent', 'one_existing_release_accepted', 'released_slot_independently_corroborated']
    assert rows[-1]['committed_slots'] == ['28709790', '28727368']


def test_actual_start_is_separate_from_release_and_scf_completion():
    report = json.loads((CASE/'startup_observation.json').read_text())
    assert '28709790|RUNNING|32|node1|node1' in report['squeue_raw']
    assert '28709791|PENDING|32||(JobHeldUser)' in report['squeue_raw']
    assert report['native_first_DSIZE'] == ['32']
    assert report['stderr_bytes'] == 0 and report['failure_exists'] is False
    assert report['first_SCF_complete_at_capture'] is False
