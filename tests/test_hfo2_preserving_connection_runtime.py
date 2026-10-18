"""Finite runtime admission, not evidence that a real job has run."""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'benchmarks/hfo2_channels/20261008/preserving_connection_E273_20261019'


def load(name):
    spec = importlib.util.spec_from_file_location('E273_test_'+name, CASE/(name+'.py'))
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


control = load('control')


@pytest.mark.parametrize('state,code', [('RUNNING','0:0'), ('TIMEOUT','0:0'), ('FAILED','1:0'), ('COMPLETED','1:0')])
def test_no_parent_live_timeout_or_failure_promotion(state, code):
    with pytest.raises(ValueError, match='original E272 whole job'):
        control.require_parent_terminal(f'28948425|{state}|{code}\n')


def test_exact_original_parent_and_one_completed_row_required():
    control.require_parent_terminal('28948425|COMPLETED|0:0\n28948425.batch|COMPLETED|0:0\n')
    for raw in ('28947723|COMPLETED|0:0', '28948425|COMPLETED|0:0\n28948425|COMPLETED|0:0'):
        with pytest.raises(ValueError): control.require_parent_terminal(raw)


def test_exact_held_resource_metadata_required():
    raw = ('JobId=123 JobState=PENDING Reason=JobHeldUser Partition=hfacnormal01 '
           'NumCPUs=32 NumTasks=1 CPUs/Task=32 NumNodes=1 TimeLimit=06:00:00 JobName=hfo2-E273-connection')
    control.validate_held(raw, '123')
    for changed in (raw.replace('NumCPUs=32','NumCPUs=64'), raw.replace('NumNodes=1','NumNodes=2'),
                    raw.replace('JobHeldUser','Dependency'), raw.replace('06:00:00','24:00:00')):
        with pytest.raises(ValueError): control.validate_held(changed, '123')


def test_inert_package_contains_exact_pinned_seed_and_no_submission(tmp_path):
    package = load('package'); output = tmp_path/'package'
    result = package.package(CASE, output)
    assert result['new_DFT_calls'] == result['scheduler_mutations'] == 0
    assert package.sha(output/'bundle/bundle_files.json') == result['bundle_pin']
    manifest = json.loads((output/'bundle/bundle_files.json').read_bytes())
    assert set(manifest) == {'control.py','run.slurm','protocol.md',
        'prepared/manifest.json','prepared/factory_parameters.json','prepared/connection_replay.json',
        'prepared/seed.traj','prepared/initial.vasp','prepared/final.vasp','prepared/substrate.vasp'}
    for name, digest in manifest.items(): assert package.sha(output/'bundle'/name) == digest
    script = (output/'bundle/run.slurm').read_text()
    assert 'sbatch' not in script and 'scontrol' not in script
    assert '--steps 24' in script and '--no-climb' in script and '--fmax .10' in script
    assert 'mode-subspace' not in script and 'mpirun -np 32' in script
    assert 'make_clamped_resume_cached_factory' in script


def test_run_requires_accepted_identity_and_two_allocation_cap(tmp_path, monkeypatch):
    monkeypatch.setattr(control, 'ROOT', tmp_path)
    (tmp_path/'accepted.json').write_text(json.dumps(dict(job_id='123',bundle_pin='p')))
    (tmp_path/'preflight.json').write_text(json.dumps(dict(bundle_pin='p',new_DFT_calls=0,decks=[{}]*9)))
    active = ['123','456']
    monkeypatch.setattr(control, 'validate', lambda pin: SimpleNamespace(queue=lambda: active))
    env = dict(RUN_DFT='1', SLURM_JOB_ID='123', SLURM_JOB_NAME='hfo2-E273-connection',
        SLURM_JOB_PARTITION='hfacnormal01', SLURM_JOB_NUM_NODES='1', SLURM_NTASKS='1', SLURM_CPUS_PER_TASK='32',
        OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    for key, value in env.items(): monkeypatch.setenv(key,value)
    control.run_check('p')
    active.append('789')
    with pytest.raises(ValueError, match='exact accepted32CPU'): control.run_check('p')
    active.pop(); monkeypatch.setenv('SLURM_CPUS_PER_TASK','64')
    with pytest.raises(ValueError): control.run_check('p')


def test_terminal_export_is_not_a_live_or_failed_batch_promotion():
    path = CASE.with_name('preserving_pattern_branches_E272_20261018')/'collect_terminal.py'
    spec = importlib.util.spec_from_file_location('E272_export_admission_test',path)
    collector = importlib.util.module_from_spec(spec); spec.loader.exec_module(collector)
    raw = '28948425|COMPLETED|0:0|32|hfacnormal01|iai806|1|5000|start|end'
    assert collector.terminal_row(raw)[7] == '5000'
    for wrong in (raw.replace('COMPLETED','RUNNING'), raw.replace('COMPLETED','TIMEOUT'),
                  raw.replace('|32|','|64|'), raw.replace('0:0','1:0')):
        with pytest.raises(ValueError): collector.terminal_row(wrong)
