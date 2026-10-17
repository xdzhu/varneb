"""Pre-label source/admission checks; no scheduler writes or DFT."""
from pathlib import Path
from copy import deepcopy
import importlib.util
import json
from types import SimpleNamespace
import numpy as np
import pytest

BASE=Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/G3_P8_execution_E266_20261018'


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


@pytest.fixture(scope='module')
def worker():return load('E266_source_admission',BASE/'worker.py')


def test_original_frozen_targets_policy_and_forecasts_replay(worker):
    actual=worker.prepare(BASE);saved=json.loads((BASE/'config.json').read_bytes())
    assert actual==saved
    assert [p['id'] for p in actual['points']]==list(range(1,9))
    assert actual['policy']['shards']==[[1,2,3,4],[5,6,7,8]]
    assert actual['policy']['maximum_new_attempts_cohort']==96
    assert actual['policy']['gradient_tolerance_eV_A']==.03
    assert actual['policy']['zeta_gradient_tolerance_eV_A']==1e-4
    assert actual['physical_input_sha256']['INPUT']=='dc6684ffa709bf3d953c647272589b05d844c1293477bd144226d495207cb9df'
    assert not actual['physical_barrier_panel_frozen'] and not actual['stable_conditional_model_admitted']


@pytest.mark.parametrize('field,value',[('maximum_new_attempts_cohort',192),('zeta_gradient_tolerance_eV_A',.03),('MPI_ranks',1)])
def test_changed_budget_gradient_or_MPI_contract_refuses(worker,field,value):
    config=deepcopy(worker.prepare(BASE));config['policy'][field]=value
    with pytest.raises(ValueError,match='immutable local'):worker.validate(config,BASE)


def test_unallocated_run_refuses_before_any_calculator_or_directory(worker,monkeypatch):
    monkeypatch.delenv('RUN_DFT',raising=False)
    monkeypatch.setattr(worker,'validate',lambda *a,**k:pytest.fail('must refuse before native/DFT access'))
    with pytest.raises(ValueError,match='authorized one-node32MPI'):worker.run({},BASE,0)


def test_actual_eight_native_inert_decks_keep_all_six_inputs_without_new_DFT(worker):
    report=json.loads((BASE/'native_preflight.json').read_bytes())
    old,_=worker.old_worker()
    assert report['new_DFT']==0 and report['original_runtime_and_assets_rechecked']
    assert [r['id'] for r in report['decks']]==list(range(1,9))
    assert report['config_sha256']==worker.sha(BASE/'config.json')
    assert all({k:r['input_sha256'][k] for k in old.study.CONTRACT}==old.study.CONTRACT for r in report['decks'])


@pytest.mark.parametrize('nodes',['1','1-1'])
def test_actual_slurm_one_node_range_schema_accepted(nodes):
    dispatcher=load('E266_admit_dispatch_'+nodes,BASE/'dispatch.py')
    raw=f'JobId=123 JobState=PENDING Reason=JobHeldUser Partition=hfacnormal01 NumCPUs=32 NumTasks=1 CPUs/Task=32 NumNodes={nodes} TimeLimit=04:00:00 JobName=hfo2-G3-E266-P0'
    assert dispatcher.validate_held(raw,0,'123')['NumCPUs']=='32'


@pytest.mark.parametrize('replacement',[('NumCPUs=32','NumCPUs=64'),('NumNodes=1','NumNodes=2'),('JobHeldUser','Dependency'),('04:00:00','24:00:00')])
def test_wrong_held_resources_or_state_refuse(replacement):
    dispatcher=load('E266_bad_held_dispatch',BASE/'dispatch.py')
    raw='JobId=123 JobState=PENDING Reason=JobHeldUser Partition=hfacnormal01 NumCPUs=32 NumTasks=1 CPUs/Task=32 NumNodes=1 TimeLimit=04:00:00 JobName=hfo2-G3-E266-P0'
    with pytest.raises(ValueError,match='held resources'):dispatcher.validate_held(raw.replace(*replacement),0,'123')


def test_four_target_controller_and_native_reader_routing_with_analytic_oracle(worker,tmp_path,monkeypatch):
    """Synthetic routing only: no physical chart, MPI or material validation."""
    from vcneb.branch_relaxation import relax_branch_candidate
    real,_=worker.old_worker();config=json.loads((BASE/'config.json').read_bytes())
    bundle=tmp_path/'bundle';bundle.mkdir();(bundle/'config.json').write_text(json.dumps(config))
    monkeypatch.setattr(worker,'EXEC',tmp_path)
    (tmp_path/'accepted_0.json').write_text(json.dumps(dict(job_id='123',config_sha256=worker.sha(bundle/'config.json'))))
    for k,v in dict(RUN_DFT='1',SLURM_JOB_ID='123',SLURM_JOB_PARTITION='hfacnormal01',SLURM_JOB_NUM_NODES='1',
        SLURM_NTASKS='1',SLURM_CPUS_PER_TASK='32',SLURM_JOB_NAME='hfo2-G3-E266-P0',
        OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1').items():monkeypatch.setenv(k,v)
    sentinel=object()
    class Oracle:
        def __init__(self,*a,workdir,**kw):
            self.directory=workdir;workdir.mkdir();self.n_new_evaluations=0;self.n_evaluation_attempts=0
        def __call__(self,x):
            self.n_new_evaluations+=1;self.n_evaluation_attempts+=1
            path=self.directory/'eval-000000';path.mkdir()
            (path/'result.json').write_text(json.dumps(dict(coordinates=x.tolist(),gradient=[0.]*36,enthalpy_eV=-1.)))
            return -1.,np.zeros(36)
    class Reader:
        def native_record(self,w,data,x,path,directory):
            assert w is sentinel
            return json.loads(path.read_bytes()),dict(result_path=str(path))
    data=dict(chart=object(),basis=np.eye(36),center_energy=0.,center_gradient=np.zeros(36))
    old=SimpleNamespace(study=SimpleNamespace(make_factory=lambda **kw:object(),ASSETS=tmp_path,COMMAND='not-called',
        result_validator=None,exclusive_json=real.study.exclusive_json),geometry=lambda *a:None,
        faces=SimpleNamespace(canary=SimpleNamespace(CalculatorModeEvaluator=Oracle,InternalEvaluator=lambda e,*a:e)),
        solver=SimpleNamespace(relax_branch_candidate=relax_branch_candidate))
    monkeypatch.setattr(worker,'validate',lambda *a,**kw:(old,(Reader(),sentinel,data,None,None,np.eye(36),np.eye(36),None)))
    monkeypatch.setattr(worker,'queue',lambda *a,**kw:dict(ready=True))
    worker.run(config,bundle,0)
    summary=json.loads((tmp_path/'shard_0/completed.json').read_bytes())
    assert [r['id'] for r in summary['rows']]==[1,2,3,4]
    assert sum(r['reserved_attempts'] for r in summary['rows'])==4
    assert all(r['status']=='gradient_candidate_stability_unchecked' for r in summary['rows'])
    assert summary['new_independent_targets_processed']==4 and not summary['physical_plus005_generated_or_read']
    with pytest.raises(ValueError,match='fresh exact'):worker.run(config,bundle,0)
