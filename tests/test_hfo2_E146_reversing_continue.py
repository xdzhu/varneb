"""Review one finite reversal continuation; fixtures never call a scheduler."""
import hashlib
import importlib.util
import json
from pathlib import Path
import shlex
import subprocess

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
CASE = REPO/'benchmarks/hfo2_channels/20261008/clamped_reversing_cap_continue_E146_20261014'
SCRIPT = REPO/'cluster/hf_hfo2_clamped_continue_E146_20261014.slurm'


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def manager():
    return load(CASE/'continue_reversing_HF.py', 'E146_review')


def recipe():
    m = load(REPO/'tests/test_hfo2_rolling_continue_recipe.py', 'E146_startup')
    m.SCRIPT = SCRIPT
    m.PARENTS = {'28764704': (40, 0, 338, m.RUNS+'20261010-varneb-nested-controls-E055-r1/source',
        m.RUNS+'20261012-varneb-clamped-flip-continue-E103-r1/28709791/PO_flip_T_pattern_reversing/band')}
    return m


def test_exact_startup_and_unchanged_production_command(tmp_path, monkeypatch, capsys):
    f = recipe()
    f.run_validator(*f.fixture(tmp_path, '28764704'), '28764704', monkeypatch)
    assert float(capsys.readouterr().out) == 5.12968067458423
    old = REPO/'cluster/hf_hfo2_clamped_continue_E103_20261012.slurm'
    assert SCRIPT.read_bytes().split(b'module purge', 1)[1] == old.read_bytes().split(b'module purge', 1)[1]
    assert manager().sha(SCRIPT) == manager().SCRIPT_SHA


@pytest.mark.parametrize('field,value', [('terminal_step', 10), ('ordinary_converged', True),
    ('fresh_interior_SCFs', 77), ('initial_fresh_interior_SCFs', 7),
    ('current_frame_exact_caches_checked', 8), ('holdout_generated_or_read', True),
    ('physical_inputs_or_running_source_changed', True), ('FIRE_state_restored', True),
    ('new_DFT_calls', 1), ('scheduler_mutations', 1),
    ('scheduler_row', '28764704|RUNNING|0:0|32|hfacnormal01|iai806')])
def test_startup_refuses_stale_cap_or_changed_contract(tmp_path, monkeypatch, field, value):
    f = recipe()
    seed, audit, runtime, source = f.fixture(tmp_path, '28764704')
    record = json.loads(audit.read_text()); record[field] = value; f.write_json(audit, record)
    with pytest.raises(ValueError):
        f.run_validator(seed, audit, runtime, source, '28764704', monkeypatch)


def test_actual_cap_raw_frames_latest_seed_and_scope():
    from ase.io import read
    from examples.hfo2_fixed_input_factory import CONTRACT, same_ordered_geometry
    from scripts.analyze_hfo2_clamped_residual import analyze
    from scripts.audit_hfo2_clamped_terminal import validate_terminal
    from scripts.audit_hfo2_static_replica import audited_results
    m = manager(); root = CASE/'completed_HF'
    assert m.sha(root/'audit_receipt.json') == m.AUDIT_SHA
    assert m.sha(root/'seed/manifest.json') == m.SEED_SHA
    r = json.loads((root/'audit_receipt.json').read_text())
    assert r['job_id'] == m.PARENT == '28764704' and r['terminal_step'] == 40
    assert not r['ordinary_converged']
    summary = json.loads((root/'vcneb_summary.json').read_text())
    validate_terminal(m.PARENT, r['scheduler_row'], summary, False, outcome='step-cap')
    with pytest.raises(ValueError):
        validate_terminal(m.PARENT, r['scheduler_row'], summary, False)
    assert len(r['calls']) == r['fresh_interior_SCFs'] == r['expected_fresh_interior_SCFs'] == 280
    assert len(r['runtime_code_sha256']) == 338 and r['current_frame_exact_caches_checked'] == 9
    assert r['all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked']
    assert all({n: x['input_sha256'][n] for n in CONTRACT} == CONTRACT for x in r['calls'])
    obs = root/'observations/step_0040'
    assert analyze(obs)['fmax_eV_A'] == pytest.approx(.1099599030476814, abs=1e-12)
    old = read(obs/'evaluated_chain.traj', index=':')
    new = read(root/'seed/seed.traj', index=':')
    assert len(old) == len(new) == 9
    for index, (a, b) in enumerate(zip(old, new)):
        assert same_ordered_geometry(a, b)
        parsed = audited_results(obs/'raw'/f'image_{index:04d}')
        for key in ('energy', 'forces', 'stress'):
            np.testing.assert_allclose(a.calc.results[key], b.calc.results[key], atol=1e-12, rtol=0)
            np.testing.assert_allclose(a.calc.results[key], parsed[key], atol=1e-12, rtol=0)
    manifest = json.loads((root/'seed/manifest.json').read_text())
    assert all(m.sha(root/'seed'/n) == h for n, h in manifest['files_sha256'].items())
    assert manifest['source_job_id'] == m.PARENT and not manifest['optimizer_state_restored']
    assert r['new_DFT_calls'] == r['scheduler_mutations'] == r['new_independent_chains'] == 0
    assert not r['holdout_generated_or_read'] and not r['physical_inputs_or_running_source_changed']
    assert r['fresh_SCF_core_hour_proxy_32'] == pytest.approx(316.20010658119287)


def fixture(tmp_path, monkeypatch, mode='success'):
    m = manager(); m.ROOT = tmp_path/'study'; m.ROOT.mkdir()
    m.LOCK = tmp_path/'shared.lock'; m.SCRIPT = tmp_path/'wrapper'
    monkeypatch.setattr(m, 'verify_parent', lambda: None)
    old = load(REPO/'benchmarks/hfo2_channels/20261008/clamped_flip_continue_E103_20261012/continue_flips_HF.py', 'E146_slots')
    old.ROOT = tmp_path/'empty-original'
    monkeypatch.setattr(m, 'legacy', lambda: old)
    calls, released = [], []
    def run(argv):
        calls.append(argv)
        if argv[0] == 'sacct':
            return m.PARENT+'|'+('RUNNING' if mode == 'live' else 'COMPLETED')+'|0:0|32|hfacnormal01|iai806'
        if argv[0] == 'bash': return ''
        if argv[0] == 'squeue':
            if mode == 'unknown': return '99900002|hfo2-unknown|RUNNING|hfacnormal01|32|node1'
            if any(c[0] == 'sbatch' for c in calls):
                return '99900001|hfo2-G2-E146-28764704|PENDING|hfacnormal01|32|'+('(Priority)' if released else '(JobHeldUser)')
            return ''
        if argv[0] == 'sbatch':
            assert json.loads((m.ROOT/'dispatch_journal.jsonl').read_text().splitlines()[-1])['event'] == m.EVENTS[0]
            return '99900001;ambiguous' if mode == 'ambiguous' else '99900001'
        if argv[:2] == ['scontrol', 'release']:
            if mode == 'timeout': raise subprocess.TimeoutExpired(argv, 45)
            if mode != 'stale': released.append(True)
            return ''
        assert argv[:3] == ['scontrol', 'show', 'job']
        p = m.paths()
        row = dict(JobId='99900001', JobName='hfo2-G2-E146-28764704', JobState='PENDING',
            Reason='Priority' if released else 'JobHeldUser', Priority='1' if released else '0',
            Partition='hfacnormal01', NumCPUs='16' if mode == 'cpu' else '32', NumTasks='1',
            Command=str(m.SCRIPT), Dependency='(null)', StdOut=str(p['stdout']), StdErr=str(p['stderr']),
            TimeLimit='12:00:00', Restarts='0', NumNodes='1', UserId='iai806(1000)', **{'CPUs/Task': '32'})
        return ' '.join(k+'='+shlex.quote(v) for k, v in row.items())
    monkeypatch.setattr(m, 'run', run)
    return m, calls


def test_one_durable_submission_release_and_no_replay(tmp_path, monkeypatch):
    m, calls = fixture(tmp_path, monkeypatch)
    r = m.dispatch()
    assert m.accepted() == r and not m.LOCK.exists()
    assert r['total_submissions'] == r['release_writes'] == 1
    assert r['source_step'] == 40 and r['new_independent_chains'] == 0
    before = len(calls)
    with pytest.raises(FileExistsError): m.dispatch()
    assert len(calls) == before


@pytest.mark.parametrize('mode,submits,releases', [('live', 0, 0), ('unknown', 0, 0),
    ('ambiguous', 1, 0), ('cpu', 1, 0), ('timeout', 1, 1), ('stale', 1, 1)])
def test_ambiguous_or_unsafe_admission_never_retries(tmp_path, monkeypatch, mode, submits, releases):
    m, calls = fixture(tmp_path, monkeypatch, mode)
    with pytest.raises((ValueError, subprocess.TimeoutExpired)): m.dispatch()
    assert sum(c[0] == 'sbatch' for c in calls) == submits
    assert sum(c[:2] == ['scontrol', 'release'] for c in calls) == releases
    assert not (m.ROOT/'dispatch_receipt.json').exists()
    if submits:
        assert m.LOCK.exists()
        before = len(calls)
        with pytest.raises(FileExistsError): m.dispatch()
        assert len(calls) == before


def test_actual_one_accepted_handle_keeps_original_ancestry_and_no_start_claim():
    m = manager()
    r = json.loads((CASE/'dispatch_receipt.json').read_text())
    journal = CASE/'dispatch_journal.jsonl'
    assert r['job_id'] == '28852667' and r['parent_job'] == '28764704'
    assert r['source_step'] == 40 and r['max_new_updates'] == 40 and r['wall_hours'] == 12
    assert r['audit_receipt_sha256'] == m.AUDIT_SHA and r['seed_manifest_sha256'] == m.SEED_SHA
    assert r['script_sha256'] == m.sha(SCRIPT)
    assert r['executed_helper_sha256'] == m.sha(CASE/'continue_reversing_HF.py')
    assert r['total_submissions'] == r['release_writes'] == 1
    assert r['dependency_writes'] == r['new_independent_chains'] == 0
    assert not r['execution_started_claimed'] and not r['FIRE_optimizer_state_restored']
    assert not r['holdout_generated_or_read'] and not r['physical_inputs_or_running_source_changed']
    assert m.sha(journal) == r['journal_sha256'] == 'd54189e9e0eb271fdbdfed5020df1acaf62ba49f9e1ef6610fec4a6bf52dd5b5'
    records = [json.loads(line) for line in journal.read_text().splitlines()]
    assert [x['event'] for x in records] == m.EVENTS
    assert records[0]['committed_slots'] == []
    assert records[2]['accepted_native_state']['Reason'] == 'JobHeldUser'
    assert records[-1]['native_control']['Reason'] != 'JobHeldUser'
    assert records[-1]['committed_slots'] == ['28852667']
