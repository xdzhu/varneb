"""Bounded5x5 orchestration; toy evaluators below are explicitly NOT native DFT."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from scripts import hfo2_G3_conditional_faces as faces


@pytest.mark.parametrize('channel',['T','preserving'])
def test_actual_native_candidate_preparation_covers_25_safe_training_starts(channel):
    config=faces.prepare(channel);data=faces.load_inputs(channel)
    assert len(config['points'])==25 and config['points'][0]['q_A']==[0.,0.]
    assert len({tuple(p['q_A']) for p in config['points']})==25
    assert config['policy']['maximum_new_attempts_per_face']==24*12
    assert config['policy']['remaining_training_node_allowance']==81-25
    assert config['policy']['remaining_independent_point_allowance']==8
    assert not config['physical_reserved_geometry_generated_or_read']
    assert not config['final_predictive_representative_or_domain_locked']
    assert data['origin_point']['orthogonal_gradient_norm_eV_A']<=.03
    faces.validate(config,channel)


@pytest.mark.parametrize('fault',['budget','coordinates','source','pressure','counter_type'])
def test_changed_policy_or_source_refuses_before_native_run(fault):
    config=faces.prepare('T')
    if fault=='budget':config['policy']['maximum_new_attempts_per_face']+=1
    elif fault=='coordinates':config['points'][1]['q_A'][0]=.3
    elif fault=='source':config['source_text_sha256']['scripts/hfo2_G3_conditional_faces.py']='0'*64
    elif fault=='pressure':config['policy']['pressure_GPa']=1.
    elif fault=='counter_type':config['policy']['independent_test_points_used']=False
    with pytest.raises(ValueError,match='grid/source/physics/budget'):faces.validate(config,'T')


def test_prepared_policy_does_not_alias_global_settings():
    config=faces.prepare('T');config['policy']['axis_A'][0]=99
    assert faces.POLICY['axis_A'][0]==-.2


@pytest.fixture
def toy(monkeypatch):
    data=faces.load_inputs('T');b=data['basis'];a=data['axes']
    r=np.linalg.svd(a.T,full_matrices=True)[2][2:].T[:,0]
    h=np.eye(36)+.6*(np.outer(a[:,0],r)+np.outer(r,a[:,0]))
    data.update(hessian=h,center_energy=0.,center_gradient=np.zeros(40),
        origin_point=dict(coordinates_internal_A=np.zeros(36).tolist(),coordinates_full_A=np.zeros(40).tolist(),
                          energy_eV_cell=0.,orthogonal_gradient_norm_eV_A=0.,callback_evaluations=0),
        origin_record=dict(gradient=np.zeros(40).tolist(),validation=dict(explicit_toy_not_native=True)))
    created=[]
    class ToyEvaluator:
        fixed_coordinate_covectors=np.eye(40)[:,-1:]
        def __init__(self,chart,factory,*,workdir,max_evaluation_attempts,**kwargs):
            self.path=workdir;self.path.mkdir();self.cap=max_evaluation_attempts
            self.n_evaluation_attempts=0;self.n_new_evaluations=0;self.cache={};created.append(self)
        def __call__(self,x):
            key=hashlib.sha256(x.tobytes()).hexdigest()
            if key in self.cache:return self.cache[key]
            if self.n_evaluation_attempts==self.cap:
                raise RuntimeError('mode evaluation attempt budget exhausted; calculator was not called')
            index=self.n_evaluation_attempts;self.n_evaluation_attempts+=1
            u=b.T@x;e=float(.5*u@h@u);g=b@(h@u)
            directory=self.path/f'eval-{index:06d}-{key[:12]}';directory.mkdir()
            row=dict(status='explicit_toy_not_native',coordinates=x.tolist(),gradient=g.tolist(),enthalpy_eV=e,
                validation=dict(elapsed_seconds=1.,explicit_toy_not_native=True))
            faces.study.exclusive_json(directory/'result.json',row)
            self.n_new_evaluations+=1;self.cache[key]=(e,g)
            return e,g.copy()
    monkeypatch.setattr(faces.canary,'CalculatorModeEvaluator',ToyEvaluator)
    monkeypatch.setattr(faces.study,'make_factory',lambda **kw:None)
    return data,created,ToyEvaluator


def test_full_toy_grid_reuses_origin_and_preserves_nonzero_Q_and_attempt_caps(tmp_path,toy):
    data,created,_=toy
    config=dict(channel='T',points=faces.ordered_points())
    result=faces.compute_face(config,data,tmp_path/'face')
    assert result['gradient_candidates']==25 and len(created)==24
    assert result['cached_origin_new_SCFs']==0
    assert 0<result['new_reserved_attempts']<=288
    assert all(e.n_evaluation_attempts<=12 for e in created)
    assert len(list((tmp_path/'face').glob('point_*.json')))==25
    for row in result['rows']:
        assert row['point']['orthogonal_gradient_norm_eV_A']<=.03
        np.testing.assert_allclose(data['axes'].T@row['point']['coordinates_internal_A'],row['q_A'],atol=2e-12,rtol=0)
        assert row['point']['coordinates_full_A'][-1]==0
    assert not result['released_stability_global_branch_TS_contour_or_prediction_certified']


def test_unsafe_continuation_skips_without_constructing_an_evaluator(tmp_path,toy,monkeypatch):
    data,created,_=toy
    monkeypatch.setattr(faces.canary,'geometry_guard',lambda *a:False)
    result=faces.compute_face(dict(channel='T',points=faces.ordered_points()),data,tmp_path/'face')
    assert result['gradient_candidates']==1 and result['new_reserved_attempts']==0 and not created
    assert all(r['point'] is None for r in result['rows'][1:])


def test_real_backend_error_category_stops_the_rest_of_the_face(tmp_path,toy,monkeypatch):
    data,created,cls=toy
    def fail(self,x):
        directory=self.path/'eval-failure';directory.mkdir()
        faces.study.exclusive_json(directory/'failure.json',dict(explicit_injected_backend_failure=True))
        raise RuntimeError('explicit injected backend failure')
    monkeypatch.setattr(cls,'__call__',fail)
    with pytest.raises(RuntimeError,match='injected backend'):
        faces.compute_face(dict(channel='T',points=faces.ordered_points()),data,tmp_path/'face')
    assert len(created)==1 and not (tmp_path/'face/face.json').exists()
    assert len(list((tmp_path/'face').glob('node_*/eval-*/failure.json')))==1


def test_full_per_node_budgets_leave_holes_not_fake_converged_points(tmp_path,toy,monkeypatch):
    data,created,_=toy
    transverse=np.linalg.svd(data['axes'].T,full_matrices=True)[2][2:].T[:,0]
    def exhausted(plane,q,evaluate,**kwargs):
        for i in range(13):evaluate(data['axes']@q+transverse*(i+1)*.001)
        raise AssertionError('thirteenth fresh call must have refused')
    monkeypatch.setattr(faces.canary,'relax_orthogonal_at_q',exhausted)
    result=faces.compute_face(dict(channel='T',points=faces.ordered_points()),data,tmp_path/'face')
    assert result['new_reserved_attempts']==result['new_completed_SCFs']==288
    assert len(created)==24 and all(e.n_evaluation_attempts==12 for e in created)
    assert result['gradient_candidates']==1
    assert all(r['point'] is None and r['status']=='optimizer_not_converged' for r in result['rows'][1:])


def test_continuation_uses_only_successes_and_fixed_order_for_distance_ties(toy):
    data,_,_=toy
    origin=dict(id=0,q_A=[0.,0.],point=data['origin_point'])
    positive=deepcopy(origin);positive.update(id=1,q_A=[.1,0.])
    negative=deepcopy(origin);negative.update(id=2,q_A=[-.1,0.])
    selected=faces.choose_start(data,[negative,positive],np.array([0.,0.]))
    assert selected[0]==1
    np.testing.assert_allclose(data['axes'].T@selected[1],[0.,0.],atol=2e-12,rtol=0)


def test_login_cannot_reach_native_preflight(tmp_path,monkeypatch):
    monkeypatch.delenv('RUN_DFT',raising=False);calls=[]
    monkeypatch.setattr(faces,'native_preflight',lambda *a:calls.append(True))
    with pytest.raises(ValueError,match='allocation'):faces.run({},'T',tmp_path/'face')
    assert not calls and not (tmp_path/'face').exists()


def test_existing_configuration_refuses_before_any_source_read(tmp_path,monkeypatch):
    path=tmp_path/'original.json';path.write_text('preserve')
    monkeypatch.setattr('sys.argv',['faces','prepare','--channel','T','--config',str(path)])
    calls=[];monkeypatch.setattr(faces,'prepare',lambda *a:calls.append(True))
    with pytest.raises(FileExistsError):faces.main()
    assert not calls and path.read_text()=='preserve'


def test_fixed_queue_wrapper_never_submits_or_uses_an_alternate_launcher():
    text=(faces.study.ROOT/'cluster/hf_hfo2_G3_conditional_faces.slurm').read_text()
    for line in ('#SBATCH -p hfacnormal01','#SBATCH -N 1','#SBATCH --ntasks=1',
                 '#SBATCH --cpus-per-task=32','#SBATCH --time=12:00:00'):
        assert line in text
    assert 'sbatch ' not in text and 'srun ' not in text and faces.study.COMMAND.startswith('mpirun -np 32 ')
