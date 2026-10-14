"""Bounded frozen controls; injected evaluators below are explicitly not DFT."""
from copy import deepcopy
import importlib.util
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'benchmarks/hfo2_channels/20261008/G3_frozen_R_controls_E166_20261015'
spec=importlib.util.spec_from_file_location('hfo2_G3_frozen_R_case',CASE/'worker.py')
worker=importlib.util.module_from_spec(spec);spec.loader.exec_module(worker)


@pytest.fixture(scope='module')
def config():return worker.prepare()


def test_actual_eight_prepared_geometries_keep_Q_freeze_R_and_preserve_control(config):
    data=worker.validate(config)
    assert config['policy']['training_slots_committed_after_reservation_per_face']==29
    assert config['policy']['training_slots_remaining_per_face']==52
    assert config['policy']['independent_slots_remaining_per_face']==8
    assert config['zero_reference_new_SCFs']==0
    assert not config['physical_reserved_geometry_generated_or_read']
    for channel,points in config['points'].items():
        assert len(points)==4 and [p['id'] for p in points]==[1,2,3,4]
        d=data[channel];b,a=d['basis'],d['axes']
        for p in points:
            x=np.asarray(p['coordinates_full_A']);u=b.T@x
            np.testing.assert_allclose(a.T@u,p['q_A'],atol=2e-12,rtol=0)
            np.testing.assert_allclose(u-a@(a.T@u),0.,atol=2e-12,rtol=0)
            assert x[-1]==0 and worker.faces.canary.geometry_guard(d['chart'],b,u)
            atoms=d['chart'].to_atoms(x)
            np.testing.assert_array_equal(atoms.cell.array[:2],d['chart'].reference_atoms.cell.array[:2])


def test_prepared_policy_and_physical_hash_maps_do_not_alias_runtime_globals():
    c=worker.prepare();c['policy']['channel_order'].reverse();c['original_six_input_sha256']['INPUT']='0'*64
    assert worker.POLICY['channel_order']==['preserving','T']
    assert worker.study.CONTRACT['INPUT']=='dc6684ffa709bf3d953c647272589b05d844c1293477bd144226d495207cb9df'


@pytest.mark.parametrize('fault',['extra_point','nonzero_R','control','wrong_Q','cutoff_hash','budget','source','counter_type'])
def test_changed_geometry_physics_or_budget_is_refused(config,fault):
    c=deepcopy(config)
    if fault=='extra_point':c['points']['T'].append(deepcopy(c['points']['T'][0]))
    elif fault=='nonzero_R':c['points']['T'][0]['coordinates_full_A'][0]+=.001
    elif fault=='control':c['points']['T'][0]['coordinates_full_A'][-1]=1e-15
    elif fault=='wrong_Q':c['points']['T'][0]['q_A']=[.2,0.]
    elif fault=='cutoff_hash':c['original_six_input_sha256']['INPUT']='0'*64
    elif fault=='budget':c['policy']['attempts_per_point']=2
    elif fault=='source':c['normalized_source_sha256'][worker.REL+'/worker.py']='0'*64
    else:c['policy']['attempts_per_point']=True
    with pytest.raises(ValueError):worker.validate(c)


def parents(t=('COMPLETED','0:0'),p=('RUNNING','0:0')):
    return {worker.PARENTS['T']:t,worker.PARENTS['preserving']:p}


def row(job='28862441',name='hfo2-G3-E161-preserving'):
    return dict(job_id=job,name=name,state='RUNNING',CPUs=32,nodes=1)


def test_admission_waits_for_real_parent_completion_and_counts_pending_allocations():
    assert not worker.admission_ready(parents(t=('RUNNING','0:0')),[row(),row('28862319','hfo2-G3-E161-T')])
    assert worker.admission_ready(parents(),[row()])
    assert not worker.admission_ready(parents(),[row(),row('new','hfo2-G3-E166-frozen')])
    assert worker.admission_ready(json.loads(json.dumps(parents())),[row()])


@pytest.mark.parametrize('terminal',['FAILED','TIMEOUT','CANCELLED','COMPLETING'])
def test_non_normal_parent_requires_review_not_automatic_submission(terminal):
    assert not worker.admission_ready(parents(t=(terminal,'0:0')),[row()])


def test_missing_parent_or_wrong_resource_is_not_assumed_safe():
    with pytest.raises(ValueError,match='handles'):worker.admission_ready({},[])
    wrong=row();wrong['CPUs']=1
    with pytest.raises(ValueError,match='resource'):worker.admission_ready(parents(),[wrong])
    assert not worker.admission_ready(parents(t=('COMPLETED','1:0')),[row()])


def test_running_control_must_be_visible_and_not_a_third_study_allocation():
    own=row('new',worker.POLICY['job_name'])
    assert worker.admission_ready(parents(),[row(),own],'new')
    assert not worker.admission_ready(parents(),[row()],'new')
    assert not worker.admission_ready(parents(),[row(),own,row('third','hfo2-G3-unexpected')],'new')


def test_toy_compute_runs_exactly_eight_single_calls_without_optimizer_or_sign_clipping(config,tmp_path,monkeypatch):
    calls=[]
    class Toy:
        def __init__(self,chart,factory,*,workdir,max_evaluation_attempts,**kw):
            assert max_evaluation_attempts==1
            workdir.mkdir(parents=True);self.n_new_evaluations=self.n_evaluation_attempts=0
        def __call__(self,x):
            calls.append(x.copy());self.n_new_evaluations=self.n_evaluation_attempts=1
    monkeypatch.setattr(worker.faces.canary,'CalculatorModeEvaluator',Toy)
    monkeypatch.setattr(worker.study,'make_factory',lambda **kw:None)
    monkeypatch.setattr(worker.faces.canary,'relax_orthogonal_at_q',lambda *a,**kw:pytest.fail('frozen controls must not relax'))
    monkeypatch.setattr(worker,'native_point',lambda data,p,d:dict(point_id=p['id'],q_A=p['q_A'],
        frozen_minus_released_meV_cell=-1.,new_SCFs=1,attempts=1,explicit_toy_not_native=True))
    result=worker.compute(config,{c:dict(chart=None) for c in worker.POLICY['channel_order']},tmp_path/'frozen')
    assert len(calls)==result['new_SCFs']==result['total_attempts']==8
    assert all(r['frozen_minus_released_meV_cell']==-1. for r in result['rows'])
    assert [r['channel'] for r in result['rows']]==['preserving']*4+['T']*4
    assert result['cached_original_zero_new_SCFs']==0
    assert not result['independent_material_prediction_or_released_stability_or_TS_certified']


def test_extra_controls_refuse_before_creating_a_workdir(config,tmp_path):
    c=deepcopy(config);c['points']['T'].append(c['points']['T'][0])
    with pytest.raises(ValueError,match='eight'):worker.compute(c,{},tmp_path/'frozen')
    assert not (tmp_path/'frozen').exists()


def test_login_refuses_before_scheduler_or_calculator(config,tmp_path,monkeypatch):
    monkeypatch.delenv('RUN_DFT',raising=False)
    monkeypatch.setattr(worker,'scheduler_snapshot',lambda *a:pytest.fail('must reject login first'))
    with pytest.raises(ValueError,match='allocation'):worker.run(config,tmp_path/'unused')


def test_native_preflight_refuses_existing_preview_before_source_reads(tmp_path,monkeypatch):
    path=tmp_path/'old';path.mkdir()
    monkeypatch.setattr(worker,'validate',lambda *a:pytest.fail('must reject output first'))
    with pytest.raises(ValueError,match='protected'):worker.native_preflight({},path)


def test_existing_configuration_refuses_before_preparation(tmp_path,monkeypatch):
    path=tmp_path/'old.json';path.write_text('preserve')
    monkeypatch.setattr('sys.argv',['frozen','prepare','--config',str(path)])
    monkeypatch.setattr(worker,'prepare',lambda:pytest.fail('must refuse first'))
    with pytest.raises(FileExistsError):worker.main()
    assert path.read_text()=='preserve'


def test_run_CLI_invokes_worker_exactly_once(tmp_path,monkeypatch):
    monkeypatch.setattr('sys.argv',['frozen','run','--config','unused','--config-sha256','unused',
        '--workdir',str(tmp_path/'frozen')])
    monkeypatch.setattr(worker.faces.canary,'pinned',lambda *a:{})
    calls=[]
    def run(*a):calls.append(True);return dict(status='explicit_toy',new_SCFs=8)
    monkeypatch.setattr(worker,'run',run);worker.main()
    assert calls==[True]


def test_fixed_wrapper_is_inert_32MPI_one_hour_and_never_submits():
    text=(CASE/'run.slurm').read_text()
    for line in ('#SBATCH -p hfacnormal01','#SBATCH -N 1','#SBATCH --ntasks=1',
                 '#SBATCH --cpus-per-task=32','#SBATCH --time=01:00:00'):
        assert line in text
    assert 'sbatch ' not in text and 'srun ' not in text
    assert worker.study.COMMAND.startswith('mpirun -np 32 ')


def test_frozen_actual_configuration_has_strict_bytes_and_one_unique_budget_reservation():
    raw=(CASE/'config.json').read_bytes()
    assert hashlib.sha256(raw).hexdigest()=='5ffb1fe0d707aec5d0a6546a44e0ad0002cd0f0280f1ba17ca2d955ac2d4c930'
    c=json.loads(raw);worker.validate(c)
    assert c['policy']['reservation_id']=='E166-four-frozen-R-axial-training-controls-per-face'
    assert c['policy']['total_new_SCF_attempts']==8
    assert c['policy']['training_slots_committed_after_reservation_per_face']+c['policy']['training_slots_remaining_per_face']==81


def test_actual_eight_inert_native_decks_pass_physical_bytes_and_do_not_claim_DFT():
    raw=(CASE/'native_preflight.json').read_bytes()
    assert hashlib.sha256(raw).hexdigest()=='3ec6390084497dc3f3a23a34b8e02baaef72b80ffc773cbb2437ea2128329b12'
    report=json.loads(raw)
    assert report['status']=='frozen_R_native_preflight_no_DFT' and report['new_DFT_calls']==0
    assert len(report['preview_decks'])==len({r['input_sha256']['STRU'] for r in report['preview_decks']})==8
    for row in report['preview_decks']:
        assert row['minimum_distance_A']>=1.6 and row['original_asset_symlinks_not_copies']
        assert all(row['input_sha256'][n]==h for n,h in worker.study.CONTRACT.items())
    assert not report['admission']['ready'] and report['admission']['scheduler_writes']==0
    assert {r['job_id'] for r in report['admission']['committed_G3_queue']}==set(worker.PARENTS.values())
    assert all(r['state']=='RUNNING' and r['CPUs']==32 and r['nodes']==1 for r in report['admission']['committed_G3_queue'])
    assert not report['physical_reserved_geometry_generated_or_read']
