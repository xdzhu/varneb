"""Two existing released bases, exact stencils, no physical backend in tests."""
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import importlib.util,json
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'benchmarks/hfo2_channels/20261008/G3_released_curvature_E181_20261015'


def load(name):
    s=importlib.util.spec_from_file_location('E181_'+name,CASE/(name+'.py'))
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


w=load('worker');submission=load('submit_once')


@pytest.fixture(scope='module')
def prepared():
    c=w.prepare(CASE);return c,w.validate(c,CASE)


def test_actual_frozen_configuration_not_only_regenerated_fixture_passes():
    frozen=json.loads((CASE/'config.json').read_bytes());w.validate(frozen,CASE)


def test_separately_amended_budget_does_not_claim_original_or_new_surface_points(prepared):
    p=prepared[0]['policy']
    assert p['new_attempts_per_center']==136 and p['new_attempts_total']==272
    assert p['additional_to_original_raw_peak_attempts']==296
    assert p['additional_verification_base_geometries']==2
    assert p['new_Q_training_nodes_per_face']==0
    assert p['training_committed']==dict(T=34,preserving=35)
    assert p['independent_face_points_untouched_each']==8 and p['automatic_retries']==0


@pytest.mark.parametrize('channel',['T','preserving'])
def test_complete_both_steps_and_exact_native_anchors_all_original_guards(prepared,channel):
    config,data=prepared;c=config['channels'][channel];d=data[channel]
    base=np.asarray(c['literal_native_anchor_40_A']);r=np.asarray(c['literal_directions_40_by_34'])
    np.testing.assert_array_equal(base,d['origin_record']['coordinates'])
    assert c['base_released_gradient_norm_eV_A']<=.03
    assert len(c['probes'])==136
    assert len({np.asarray(p['coordinates_full_A']).tobytes() for p in c['probes']})==136
    for h in (.01,.02):
        points=[p for p in c['probes'] if p['step_A']==h]
        assert [(p['direction'],p['sign']) for p in points]==[(i,s) for i in range(34) for s in (1,-1)]
        for p in points:
            full=base+r[:,p['direction']]*(p['sign']*h)
            np.testing.assert_array_equal(full,p['coordinates_full_A'])
            w.safe_probe(d,full)


@pytest.mark.parametrize('fault',['cutoff','budget','holdout','new_center','extra_point','coordinate','rank','source'])
def test_changes_refused_before_DFT(prepared,fault):
    c=deepcopy(prepared[0])
    if fault=='cutoff':c['physical_input_sha256']['INPUT']='0'*64
    elif fault=='budget':c['policy']['new_attempts_total']=273
    elif fault=='holdout':c['physical_reserved_geometry_generated_or_read']=True
    elif fault=='new_center':c['channels']['T']['literal_native_anchor_40_A'][0]+=.01
    elif fault=='extra_point':c['channels']['T']['probes'].append(c['channels']['T']['probes'][0])
    elif fault=='coordinate':c['channels']['T']['probes'][0]['coordinates_full_A'][0]=np.nextafter(c['channels']['T']['probes'][0]['coordinates_full_A'][0],1.)
    elif fault=='rank':c['channels']['T']['released_basis_36_by_34'][0][0]+=.1
    else:c['worker_text_sha256']='0'*64
    with pytest.raises(ValueError):w.validate(c,CASE)


def test_literal_full_directions_allow_roundoff_reconstruction_not_changed_probe_bytes(prepared):
    config,data=prepared;c=config['channels']['T'];d=deepcopy(data['T'])
    d['basis'][0,0]=np.nextafter(d['basis'][0,0],1.)
    actual=w.make_channel(d,c['released_basis_36_by_34'],declared=c)
    assert actual==c


@pytest.mark.parametrize('fault',['radius','control','Q','nonfinite'])
def test_unsafe_probe_is_rejected_not_moved_or_repaired(prepared,fault):
    d=prepared[1]['T'];x=np.asarray(d['origin_record']['coordinates']).copy()
    if fault=='radius':x+=d['basis'][:,5]*2
    elif fault=='control':x[-1]=.005*w.study.SCALE
    elif fault=='Q':x+=d['basis']@d['axes'][:,0]*.1
    else:x[0]=np.nan
    before=x.copy()
    with pytest.raises(ValueError):w.safe_probe(d,x)
    np.testing.assert_array_equal(x,before)


def test_full_gradient_central_differences_preserve_negative_sign_and_both_steps(prepared):
    spec=prepared[0]['channels']['T'];directions=np.asarray(spec['literal_directions_40_by_34'])
    base=np.asarray(spec['literal_native_anchor_40_A'])
    eigen=np.arange(1,35,dtype=float);eigen[0]=-2
    h=directions@np.diag(eigen)@directions.T;g=np.linspace(-.2,.2,40)
    gradients=[g+h@(np.asarray(p['coordinates_full_A'])-base) for p in spec['probes']]
    blocks=w.blocks_from_gradients(directions,gradients)
    for b in blocks:
        assert b['complete_pairs']==34 and b['symmetric_projected_eigenvalues_eV_A2'][0]<0
        np.testing.assert_allclose(b['raw_projected_eV_A2'],np.diag(eigen),atol=2e-12,rtol=0)
    with pytest.raises(ValueError,match='complete136'):w.blocks_from_gradients(directions,gradients[:-1])
    gradients[3][2]=np.nan
    with pytest.raises(ValueError,match='complete136'):w.blocks_from_gradients(directions,gradients)


def test_pending_and_held_allocations_count_and_unknown_layout_refuses(monkeypatch):
    def queue(s):monkeypatch.setattr(w.subprocess,'run',lambda *a,**k:SimpleNamespace(stdout=s))
    queue('28867177|hfo2-G3-E175-atoms|PENDING|32|1\n');assert w.scheduler_snapshot()['ready']
    queue('28867177|hfo2-G3-E175-atoms|PENDING|32|1\n2|hfo2-G3-E181-release-Hrr|PENDING|32|1\n')
    assert not w.scheduler_snapshot()['ready'] and w.scheduler_snapshot('2')['ready']
    queue('2|hfo2-G3-E181-release-Hrr|RUNNING|64|1\n')
    with pytest.raises(ValueError):w.scheduler_snapshot('2')
    queue('garbage')
    with pytest.raises(ValueError):w.scheduler_snapshot()


def test_held_job_must_have_exact_requested_resources_command_and_time():
    exact='JobState=PENDING Reason=JobHeldUser Partition=hfacnormal01 NumNodes=1 NumCPUs=32 NumTasks=1 CPUs/Task=32 TimeLimit=12:00:00 Command='+str(submission.BUNDLE/'run.slurm')
    submission.validate_allocation(exact)
    for old,new in [('NumCPUs=32','NumCPUs=64'),('12:00:00','24:00:00'),('JobHeldUser','Resources'),('hfacnormal01','other')]:
        with pytest.raises(ValueError):submission.validate_allocation(exact.replace(old,new))


def test_run_cannot_start_from_login_environment_or_emit_completions(prepared,tmp_path,monkeypatch):
    monkeypatch.delenv('RUN_DFT',raising=False)
    with pytest.raises(ValueError,match='authorized1node32MPI'):
        w.run(prepared[0],tmp_path/'bundle',tmp_path/'released_curvature')
    assert not (tmp_path/'released_curvature').exists()


@pytest.mark.parametrize('fail',[False,True])
def test_finite_orchestration_and_failure_charge_no_retry_or_partial_promotion(prepared,tmp_path,monkeypatch,fail):
    config,data=prepared;bundle=tmp_path/'bundle';bundle.mkdir();calls=[]
    for key,value in dict(RUN_DFT='1',SLURM_JOB_PARTITION='hfacnormal01',SLURM_JOB_NUM_NODES='1',
            SLURM_NTASKS='1',SLURM_CPUS_PER_TASK='32',SLURM_JOB_ID='2',
            SLURM_JOB_NAME=w.POLICY['job_name'],OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1').items():
        monkeypatch.setenv(key,value)
    monkeypatch.setattr(w.study,'RUNS',tmp_path)
    monkeypatch.setattr(w,'preflight',lambda *a:(data,{'new_DFT':0}))
    monkeypatch.setattr(w,'scheduler_snapshot',lambda *a:dict(ready=True))
    class Fake:
        def __init__(self,*args,workdir,max_evaluation_attempts,**kwargs):
            assert max_evaluation_attempts==136
            self.workdir=workdir;workdir.mkdir(parents=True)
            self.n_evaluation_attempts=0;self.n_new_evaluations=0
        def __call__(self,x):
            calls.append(x.copy());i=self.n_evaluation_attempts;self.n_evaluation_attempts+=1
            (self.workdir/f'attempt-{i:06d}.json').write_text('{}')
            p=self.workdir/f'eval-{i:06d}-toy';p.mkdir()
            if fail:
                (p/'failure.json').write_text('{}')
                raise RuntimeError('injected native failure, charged once')
            self.n_new_evaluations+=1
            (p/'result.json').write_text(json.dumps({'validation':{'elapsed_seconds':1.}}))
            return 0.,x.copy()  # explicit synthetic full-space quadratic, not DFT
    monkeypatch.setattr(w.faces.canary,'CalculatorModeEvaluator',Fake)
    if fail:
        with pytest.raises(RuntimeError,match='charged once'):
            w.run(config,bundle,tmp_path/'released_curvature')
        assert len(calls)==1
        assert len(list((tmp_path/'released_curvature').glob('*/evaluations/attempt-*.json')))==1
        assert not list((tmp_path/'released_curvature').glob('*/curvature.json'))
        assert not (tmp_path/'released_curvature/completed.json').exists()
    else:
        w.run(config,bundle,tmp_path/'released_curvature')
        assert len(calls)==272
        summary=json.loads((tmp_path/'released_curvature/completed.json').read_text())
        assert summary['new_native_SCFs']==272 and summary['new_center_SCFs']==0
        for channel in ('T','preserving'):
            result=json.loads((tmp_path/f'released_curvature/{channel}/curvature.json').read_text())
            assert result['new_native_SCFs']==136 and not result['physical_stability_TS_or_independent_prediction_certified']
            for block in result['blocks']:
                np.testing.assert_allclose(block['raw_projected_eV_A2'],np.eye(34),atol=2e-12,rtol=0)
